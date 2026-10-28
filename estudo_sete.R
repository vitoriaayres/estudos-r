pacotes <- c("sf", "geobr", "abjData", "dplyr", "tidyr", "ggplot2", "spdep", "cluster")
faltando <- setdiff(pacotes, rownames(installed.packages()))
if (length(faltando)) install.packages(faltando)

library(sf)
library(dplyr)
library(ggplot2)
library(spdep)

UF_SIGLA  <- "BA"
UF_CODIGO <- 29
SEMENTE   <- 42

# Indicadores do Atlas do Desenvolvimento Humano (Censo 2010)
VARIAVEIS <- c("t_analf15m", "pind", "gini", "mort1", "agua_esgoto", "rdpc")

carregar_malha <- function(uf) {
  geobr::read_municipality(code_muni = uf, year = 2010, showProgress = FALSE) |>
    st_transform(5880) |>
    mutate(code_muni = as.character(code_muni)) |>
    select(code_muni, name_muni)
}

carregar_indicadores <- function(uf_codigo, variaveis) {
  atlas <- abjData::pnud_muni
  ausentes <- setdiff(variaveis, names(atlas))
  if (length(ausentes)) {
    stop("Colunas ausentes em pnud_muni: ", paste(ausentes, collapse = ", "),
         ". Confira names(abjData::pnud_muni).")
  }
  atlas |>
    filter(ano == 2010, substr(as.character(codmun7), 1, 2) == as.character(uf_codigo)) |>
    transmute(code_muni = as.character(codmun7), across(all_of(variaveis)))
}

base <- carregar_malha(UF_SIGLA) |>
  left_join(carregar_indicadores(UF_CODIGO, VARIAVEIS), by = "code_muni")
stopifnot(!anyNA(st_drop_geometry(base)[VARIAVEIS]))

# Renda invertida para que valores maiores signifiquem sempre mais vulnerabilidade
montar_matriz <- function(base) {
  dados <- st_drop_geometry(base) |>
    mutate(renda_inv = -log(rdpc)) |>
    select(-code_muni, -name_muni, -rdpc)
  scale(as.matrix(dados))
}

X <- montar_matriz(base)
round(cor(X), 2)

avaliar_k <- function(X, ks = 2:8) {
  distancias <- dist(X)
  resultados <- lapply(ks, function(k) {
    km <- kmeans(X, centers = k, nstart = 25)
    sil <- cluster::silhouette(km$cluster, distancias)
    data.frame(k = k, wss = km$tot.withinss, silhueta = mean(sil[, "sil_width"]))
  })
  do.call(rbind, resultados)
}

set.seed(SEMENTE)
avaliar_k(X)

K <- 4

set.seed(SEMENTE)
km <- kmeans(X, centers = K, nstart = 50)
base$kmeans <- km$cluster

# A árvore geradora do SKATER exige grafo conexo; ilhas são ligadas ao vizinho mais próximo
garantir_conectividade <- function(vizinhos, geometria) {
  coords <- st_coordinates(st_point_on_surface(st_geometry(geometria)))
  k <- 1
  while (n.comp.nb(vizinhos)$nc > 1) {
    knn <- knn2nb(knearneigh(coords, k = k), row.names = attr(vizinhos, "region.id"))
    vizinhos <- union.nb(vizinhos, make.sym.nb(knn))
    k <- k + 1
  }
  vizinhos
}

vizinhos <- poly2nb(base, queen = TRUE) |> garantir_conectividade(base)

custos <- nbcosts(vizinhos, X)
pesos  <- nb2listw(vizinhos, custos, style = "B")
arvore <- mstree(pesos)

regioes <- skater(edges = arvore[, 1:2], data = X, ncuts = K - 1, crit = 10)
base$skater <- regioes$groups

soma_quadrados_intra <- function(X, grupos) {
  sum(sapply(split(seq_len(nrow(X)), grupos), function(i) {
    sum(scale(X[i, , drop = FALSE], scale = FALSE)^2)
  }))
}

silhueta_media <- function(X, grupos) {
  mean(cluster::silhouette(grupos, dist(X))[, "sil_width"])
}

contar_fragmentos <- function(vizinhos, grupos) {
  sum(sapply(unique(grupos), function(g) n.comp.nb(subset(vizinhos, grupos == g))$nc))
}

comparar_metodos <- function(X, vizinhos, base, metodos = c("kmeans", "skater")) {
  do.call(rbind, lapply(metodos, function(m) {
    grupos <- base[[m]]
    data.frame(metodo     = m,
               ssw        = soma_quadrados_intra(X, grupos),
               silhueta   = silhueta_media(X, grupos),
               fragmentos = contar_fragmentos(vizinhos, grupos))
  }))
}

comparar_metodos(X, vizinhos, base)

# Renumera os clusters do menos (1) ao mais vulnerável (K) para comparar os mapas
ordenar_por_vulnerabilidade <- function(grupos, indice) {
  posicao <- rank(tapply(indice, grupos, mean))
  factor(unname(posicao[as.character(grupos)]), levels = seq_along(posicao))
}

indice_vulnerabilidade <- rowMeans(X)
base$kmeans <- ordenar_por_vulnerabilidade(base$kmeans, indice_vulnerabilidade)
base$skater <- ordenar_por_vulnerabilidade(base$skater, indice_vulnerabilidade)

base |>
  st_drop_geometry() |>
  group_by(skater) |>
  summarise(municipios = n(), across(all_of(VARIAVEIS), \(v) round(mean(v), 2)))

mapa_comparativo <- base |>
  select(name_muni, kmeans, skater) |>
  tidyr::pivot_longer(c(kmeans, skater), names_to = "metodo", values_to = "cluster") |>
  ggplot() +
  geom_sf(aes(fill = cluster), color = "white", linewidth = 0.05) +
  facet_wrap(~ metodo, labeller = as_labeller(c(kmeans = "k-means", skater = "SKATER"))) +
  scale_fill_brewer(palette = "YlOrRd", name = "Vulnerabilidade\n(1 = menor)") +
  labs(title = "Agrupamentos de vulnerabilidade social", subtitle = "Bahia, Censo 2010") +
  theme_void()

ggsave("vulnerabilidade_clusters.png", mapa_comparativo, width = 11, height = 6, dpi = 150)
