pacotes_cran <- c("sf", "geobr", "sidrar", "dplyr", "tidyr", "ggplot2",
                  "spdep", "spatialreg", "remotes")
faltando <- setdiff(pacotes_cran, rownames(installed.packages()))
if (length(faltando)) install.packages(faltando)

if (!requireNamespace("microdatasus", quietly = TRUE)) {
  remotes::install_github("rfsaldanha/microdatasus", upgrade = "never")
}

library(sf)
library(dplyr)
library(ggplot2)
library(spdep)
library(spatialreg)

UF_SIGLA  <- "BA"
UF_CODIGO <- 29
ANOS      <- 2020:2022

carregar_malha <- function(uf) {
  geobr::read_municipality(code_muni = uf, year = 2022, showProgress = FALSE) |>
    mutate(cod6 = substr(as.character(code_muni), 1, 6)) |>
    select(code_muni, cod6, name_muni)
}

# Homicídios: causas básicas X85 a Y09 da CID-10 (agressões)
carregar_homicidios <- function(uf, anos) {
  sim <- microdatasus::fetch_datasus(
    year_start = min(anos), year_end = max(anos),
    uf = uf, information_system = "SIM-DO"
  )
  sim |>
    mutate(cid3 = substr(CAUSABAS, 1, 3)) |>
    filter(cid3 >= "X85", cid3 <= "Y09") |>
    count(cod6 = as.character(CODMUNRES), name = "obitos")
}

carregar_populacao <- function(uf_codigo) {
  sidrar::get_sidra(
    x = 4714, variable = c(93, 614), period = "2022",
    geo = "City", geo.filter = list(State = uf_codigo)
  ) |>
    select(cod7 = `Município (Código)`, variavel = `Variável`, valor = Valor) |>
    tidyr::pivot_wider(names_from = variavel, values_from = valor) |>
    transmute(cod6 = substr(cod7, 1, 6),
              populacao = `População residente`,
              densidade = `Densidade demográfica`)
}

# Suavização bayesiana empírica: estabiliza as taxas de municípios pequenos
preparar_base <- function(malha, homicidios, populacao, n_anos) {
  base <- malha |>
    left_join(populacao,  by = "cod6") |>
    left_join(homicidios, by = "cod6") |>
    mutate(obitos     = coalesce(obitos, 0L),
           pessoa_ano = populacao * n_anos,
           taxa_bruta = obitos / pessoa_ano * 1e5,
           log_dens   = log(densidade))

  eb <- EBest(base$obitos, base$pessoa_ano)
  base$taxa_eb <- eb$estmm * 1e5
  base
}

base <- preparar_base(carregar_malha(UF_SIGLA),
                      carregar_homicidios(UF_SIGLA, ANOS),
                      carregar_populacao(UF_CODIGO),
                      length(ANOS))
stopifnot(!anyNA(base$populacao))
summary(base[, c("taxa_bruta", "taxa_eb")])

vizinhos <- poly2nb(base, queen = TRUE)
pesos    <- nb2listw(vizinhos, style = "W", zero.policy = TRUE)
summary(vizinhos)

moran_global <- moran.mc(base$taxa_eb, pesos, nsim = 999, zero.policy = TRUE)
moran_global

classificar_lisa <- function(base, pesos, alfa = 0.05) {
  lisa <- localmoran(base$taxa_eb, pesos, zero.policy = TRUE)
  quadrante <- attr(lisa, "quadr")$mean
  p_valor   <- lisa[, "Pr(z != E(Ii))"]

  base$lisa <- ifelse(p_valor < alfa, as.character(quadrante), "Não significativo")
  base$lisa <- factor(base$lisa, levels = c("High-High", "Low-Low", "High-Low",
                                            "Low-High", "Não significativo"))
  base
}

base <- classificar_lisa(base, pesos)
table(base$lisa)

mapa_taxa <- ggplot(base) +
  geom_sf(aes(fill = taxa_eb), color = NA) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "Homicídios\npor 100 mil") +
  labs(title = "Taxa suavizada de homicídios", subtitle = "Bahia, 2020-2022") +
  theme_void()

cores_lisa <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
                "High-Low" = "#fdae61", "Low-High" = "#abd9e9",
                "Não significativo" = "grey90")

mapa_lisa <- ggplot(base) +
  geom_sf(aes(fill = lisa), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = cores_lisa, name = "Cluster LISA", drop = FALSE) +
  labs(title = "Clusters locais de homicídios") +
  theme_void()

ggsave("homicidios_taxa.png", mapa_taxa, width = 7, height = 7, dpi = 150)
ggsave("homicidios_lisa.png", mapa_lisa, width = 7, height = 7, dpi = 150)

ols <- lm(taxa_eb ~ log_dens, data = base)
summary(ols)

lm.morantest(ols, pesos, zero.policy = TRUE)
lm.RStests(ols, pesos, test = "all", zero.policy = TRUE)

sar <- lagsarlm(taxa_eb ~ log_dens, data = base, listw = pesos, zero.policy = TRUE)
sem <- errorsarlm(taxa_eb ~ log_dens, data = base, listw = pesos, zero.policy = TRUE)

data.frame(modelo = c("OLS", "SAR", "SEM"),
           AIC    = c(AIC(ols), AIC(sar), AIC(sem)))

impacts(sar, listw = pesos)
