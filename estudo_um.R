# Inferência causal com os dados de LaLonde (NSW + CPS)
# Requer o pacote causaldata: install.packages('causaldata')
library(causaldata)
set.seed(2024)
#=================================================================================================

adicionar_desemprego <- function(tabela) {
  tabela <- as.data.frame(tabela)
  tabela$u74 <- as.integer(tabela$re74 == 0)
  tabela$u75 <- as.integer(tabela$re75 == 0)
  tabela
}
experimento <- adicionar_desemprego(nsw_mixtape)
cps <- adicionar_desemprego(cps_mixtape)
observacional <- rbind(experimento[experimento$treat == 1, ], cps)

covariaveis <- c('age', 'educ', 'black', 'hisp', 'marr', 'nodegree', 're74', 're75', 'u74', 'u75')

stopifnot(nrow(experimento) == 445, all(cps$treat == 0), sum(observacional$treat == 1) == 185)
cat(sprintf("Experimento: %d pessoas | Observacional: %d tratados + %d controles do censo\n",
            nrow(experimento), sum(observacional$treat == 1), sum(observacional$treat == 0)))
#--------------------------------------------------------------------------------------------------

diferenca_medias <- function(tabela) {
  mean(tabela$re78[tabela$treat == 1]) - mean(tabela$re78[tabela$treat == 0])
}
efeito_experimental <- diferenca_medias(experimento)
efeito_ingenuo <- diferenca_medias(observacional)
stopifnot(abs(efeito_experimental - 1794) < 1)

cat(sprintf("Efeito real (experimento com sorteio):  US$ %8.0f\n", efeito_experimental))
cat(sprintf("Efeito ingênuo (tratados vs. censo):    US$ %8.0f\n", efeito_ingenuo))

media_por_grupo <- function(tabela) colMeans(tabela[c(covariaveis, 're78')])
comparacao <- data.frame(
  participantes = media_por_grupo(experimento[experimento$treat == 1, ]),
  controle_sorteado = media_por_grupo(experimento[experimento$treat == 0, ]),
  censo = media_por_grupo(cps)
)
print(round(comparacao, 2))

options(repr.plot.width= 12, repr.plot.height= 6)
barplot(t(as.matrix(comparacao[c('re74', 're75', 're78'), ])), beside= TRUE,
        names.arg= c('RENDA 1974', 'RENDA 1975', 'RENDA 1978'),
        col= c('TOMATO', 'GREY60', 'STEELBLUE'), ylab= 'US$',
        main= 'QUEM ENTROU NO PROGRAMA JÁ GANHAVA MUITO MENOS ANTES DELE',
        legend.text= c('PARTICIPANTES', 'CONTROLE SORTEADO', 'CENSO (CPS)'),
        args.legend= list(x= 'topleft', bty= 'n'))

#===========================================================================================

diferenca_padronizada <- function(tabela, pesos= rep(1, nrow(tabela))) {
  tratado <- tabela$treat == 1
  sapply(covariaveis, function(variavel) {
    x <- tabela[[variavel]]
    media_tratados <- weighted.mean(x[tratado], pesos[tratado])
    media_controles <- weighted.mean(x[!tratado], pesos[!tratado])
    (media_tratados - media_controles) / sd(x[tratado])
  })
}
print(round(diferenca_padronizada(observacional), 2))
#===============================================================================

formula_resultado <- re78 ~ treat + age + I(age^2) + educ + black + hisp + marr +
  nodegree + re74 + re75 + u74 + u75
regressao <- lm(formula_resultado, data = observacional)
cat(sprintf("Efeito pela regressão: US$ %.0f (IC 95%%: %.0f a %.0f)\n",
            coef(regressao)[["treat"]],
            confint(regressao)["treat", 1], confint(regressao)["treat", 2]))
#====================================================================================

formula_propensao <- treat ~ age + I(age^2) + I(age^3) + educ + I(educ^2) + marr +
  nodegree + black + hisp + re74 + re75 + u74 + u75 + educ:re74

estimar_propensao <- function(tabela) {
  modelo <- suppressWarnings(glm(formula_propensao, data = tabela, family = binomial))
  as.numeric(fitted(modelo))
}

observacional$propensao <- estimar_propensao(observacional)

logito <- qlogis(pmin(pmax(observacional$propensao, 1e-6), 1 - 1e-6))
quebras <- seq(min(logito), max(logito), length.out = 50)
histograma_controles <- hist(logito[observacional$treat == 0], breaks= quebras, plot= FALSE)
histograma_tratados <- hist(logito[observacional$treat == 1], breaks= quebras, plot= FALSE)

plot(histograma_controles, freq = FALSE, col = adjustcolor("grey40", 0.5), border = NA,
     ylim = c(0, max(histograma_controles$density, histograma_tratados$density)),
     xlab = "Escore de propensão (escala logito)", ylab = "Densidade",
     main = "Sobreposição: poucos controles se parecem com os participantes")
plot(histograma_tratados, freq = FALSE, col = adjustcolor("tomato", 0.6), border = NA, add = TRUE)
legend("topright", c("Censo (controles)", "Participantes"),
       fill = c(adjustcolor("grey40", 0.5), adjustcolor("tomato", 0.6)), bty = "n")

minimo_tratados <- min(observacional$propensao[observacional$treat == 1])
cat(sum(observacional$treat == 0 & observacional$propensao >= minimo_tratados), "de",
    sum(observacional$treat == 0), "controles têm escore dentro da faixa dos participantes\n")
#========================================================================================================

estimar_pareamento <- function(tabela, caliper= 0.2) {
  logito <- qlogis(pmin(pmax(tabela$propensao, 1e-8), 1 - 1e-8))
  tratados <- which(tabela$treat == 1)
  controles <- which(tabela$treat == 0)
  limite <- caliper * sd(logito)
  par_de <- sapply(tratados, function(i) {
    distancias <- abs(logito[controles] - logito[i])
    j <- which.min(distancias)
    if (distancias[j] <= limite) controles[j] else NA
  })
  com_par <- !is.na(par_de)
  list(efeito = mean(tabela$re78[tratados[com_par]] - tabela$re78[par_de[com_par]]),
       tratados = tratados[com_par], controles = par_de[com_par])
}
calcular_pesos_pareamento <- function(tabela, pareamento) {
  pesos <- numeric(nrow(tabela))
  pesos[pareamento$tratados] <- 1
  usos <- table(pareamento$controles)
  pesos[as.integer(names(usos))] <- as.numeric(usos)
  pesos
}

pareamento <- estimar_pareamento(observacional)
pesos_pareamento <- calcular_pesos_pareamento(observacional, pareamento)
cat(sprintf('EFEITO PELO PAREAMENTO: US$ %.0f\n', pareamento$efeito))
cat(sprintf('PARTICIPANTES PAREADOS: %d de 185 | CONTROLES DISTINTOS USADOS: %d\n',
            length(pareamento$tratados), length(unique(pareamento$controles))))
#================================================================================================

estimar_ipw <- function(tabela) {
  p <- pmin(tabela$propensao, 0.99)
  tratado <- tabela$treat == 1
  pesos <- ifelse(tratado, 1, p / (1 - p))
  list(efeito = mean(tabela$re78[tratado]) - weighted.mean(tabela$re78[!tratado], pesos[!tratado]),
       pesos= pesos)
}
estimar_duplamente_robusto <- function(tabela, pesos) {
  tabela$pesos_ipw <- pesos
  coef(lm(formula_resultado, data = tabela, weights = pesos_ipw))[['treat']]
}
ipw <- estimar_ipw(observacional)
efeito_duplamente_robusto <- estimar_duplamente_robusto(observacional, ipw$pesos)
pesos_controles <- ipw$pesos[observacional$treat == 0]

cat(sprintf("Efeito por IPW:                 US$ %.0f\n", ipw$efeito))
cat(sprintf("Efeito duplamente robusto:      US$ %.0f\n", efeito_duplamente_robusto))
cat(sprintf("Tamanho efetivo dos controles:  %.0f (de %d)\n",
            sum(pesos_controles)^2 / sum(pesos_controles^2), length(pesos_controles)))
#================================================================================================

balanco <- data.frame(
  antes = abs(diferenca_padronizada(observacional)),
  pareamento = abs(diferenca_padronizada(observacional, pesos_pareamento)),
  ipw = abs(diferenca_padronizada(observacional, ipw$pesos))
)
balanco <- balanco[order(balanco$antes), ]
posicoes <- seq_len(nrow(balanco))

par(mar= c(5, 8, 4, 2))
plot(balanco$antes, posicoes, xlim= c(0, max(balanco)), yaxt= 'n', pch= 19, cex= 1.4,
     col= 'grey50', xlab= 'DIFERENÇA MÉDIA PADRONIZADA (ABSOLUTA)', ylab= ' ',
     main= 'BALANÇO DAS COVARIÁVEIS ANTES E DEPOIS DO AJUSTE')
points(balanco$pareamento, posicoes, pch= 17, cex= 1.4, col= 'steelblue')
points(balanco$ipw, posicoes, pch= 15, cex= 1.4, col= 'darkorange')
axis(2, at= posicoes, labels= rownames(balanco), las= 1)
abline(v= 0.1, lty= 2)
legend('bottomright', c('ANTES', 'PAREAMENTO', 'IPW'), pch= c(19, 17, 15),
       col= c('grey50', 'steelblue', 'darkorange'), bty= 'n')
par(mar= c(5, 4, 4, 2))
print(round(balanco, 3))
#================================================================================================

estimar_tudo <- function(tabela) {
  tabela$propensao <- estimar_propensao(tabela)
  ipw <- estimar_ipw(tabela)
  c(ingenuo = diferenca_medias(tabela),
    regressao = coef(lm(formula_resultado, data = tabela))[["treat"]],
    pareamento = estimar_pareamento(tabela)$efeito,
    ipw = ipw$efeito,
    duplamente_robusto = estimar_duplamente_robusto(tabela, ipw$pesos))
}

reamostrar <- function(tabela) {
  tratados <- which(tabela$treat == 1)
  controles <- which(tabela$treat == 0)
  tabela[c(sample(tratados, replace = TRUE), sample(controles, replace = TRUE)), ]
}

inicio <- Sys.time()
estimativas <- estimar_tudo(observacional)
amostras_bootstrap <- t(replicate(200, estimar_tudo(reamostrar(observacional))))
intervalos <- apply(amostras_bootstrap, 2, quantile, probs = c(0.025, 0.975))
cat("Tempo:", round(as.numeric(difftime(Sys.time(), inicio, units = "secs"))), "segundos\n\n")

tabela_final <- data.frame(
  estimativa = round(estimativas),
  ic_inferior = round(intervalos[1, ]),
  ic_superior = round(intervalos[2, ]),
  distancia_do_real = round(estimativas - efeito_experimental)
)
print(tabela_final)
#====================================================================================================

rotulos <- c(ingenuo = 'COMPARAÇÃO INGÊNUA', regressao= 'REGRESSÃO LINEAR',
             pareamento= 'PAREAMENTO POR ESCORE', ipw= 'PONDERAÇÃO (IPW)',
             duplamente_robusto= 'DUPLAMENTE ROBUSTO')

desenhar_resultado <- function() {
  ordem <- rev(names(rotulos))
  posicoes <- seq_along(ordem)
  par(mar = c(5, 15, 5, 2))
  plot(estimativas[ordem], posicoes, xlim= range(c(intervalos, efeito_experimental, 0)),
       ylim= c(0.5, length(ordem) + 0.5), yaxt= 'n', pch= 19, cex= 1.8,
       xlab= 'EFEITO ESTIMADO DO PROGRAMA SOBRE A RENDA DE 1978 (US$)', ylab= ' ',
       main= 'MESMOS DADOS, CONCLUSÕES OPOSTAS:\nO MÉTODO DE INFERÊNCIA CAUSAL FAZ A DIFERENÇA')
  segments(intervalos[1, ordem], posicoes, intervalos[2, ordem], posicoes, lwd= 3)
  abline(v= efeito_experimental, col= 'forestgreen', lwd= 3, lty= 2)
  abline(v= 0, col= 'grey60')
  axis(2, at= posicoes, labels= rotulos[ordem], las= 1)
  legend('bottomleft', legend= sprintf('EFEITO REAL (EXPERIMENTO): US$ %.0f', efeito_experimental),
         col= 'forestgreen', lty= 2, lwd= 3, bty= 'n')
  par(mar= c(5, 4, 4, 2))
}
png("inferencia_causal_resultado.png", width = 1800, height = 1000, res = 160)
desenhar_resultado()
invisible(dev.off())
cat("Gráfico salvo em inferencia_causal_resultado.png\n")
