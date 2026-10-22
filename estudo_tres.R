set.seed(123)
taxa_real <- c(A = 0.040, B= 0.046)
dias <- 30
visitantes_por_dia <- 1000

simular_experimento <- function(taxas, dias, visitantes_por_dia) {
  data.frame(
    dia= seq_len(dias),
    visitantes_a= visitantes_por_dia, 
    conversoes_a= rbinom(dias, visitantes_por_dia, taxas[1]), 
    visitantes_b= visitantes_por_dia,
    conversoes_b= rbinom(dias, visitantes_por_dia, taxas[2])
  )
}
acumular <- function(experimento) {
  data.frame(
    dia= experimento$dia ,
    n_a= cumsum(experimento$visitantes_a), c_a= cumsum(experimento$conversoes_a),
    n_b= cumsum(experimento$visitantes_b), c_b= cumsum(experimento$conversoes_b)
  )
}
experimento <- simular_experimento(taxa_real, dias, visitantes_por_dia)
acumulado <- acumular(experimento)
final <- acumulado[dias, ]
stopifnot(is.data.frame(acumulado), nrow(acumulado) == dias)

cat(sprintf("Conversão final: A = %.2f%% | B = %.2f%% (%d visitantes por variante)\n",
            100 * final$c_a / final$n_a, 100 * final$c_b / final$n_b, final$n_a))

options(repr.plot.width = 10, repr.plot.height = 5)
plot(acumulado$dia, 100 * acumulado$c_a / acumulado$n_a, type = "l", lwd = 3, col = "grey40",
     ylim = range(100 * c(acumulado$c_a / acumulado$n_a, acumulado$c_b / acumulado$n_b)),
     xlab = "Dia", ylab = "Conversão acumulada (%)", main = "Conversão acumulada ao longo do teste")
lines(acumulado$dia, 100 * acumulado$c_b / acumulado$n_b, lwd = 3, col = "#377EB8")
legend("topright", c("A (atual)", "B (novo botão)"), col = c("grey40", "#377EB8"), lwd = 3, bty = "n")


#=======================================================================================================

p_valor_z <- function(c_a, n_a, c_b, n_b) {
  p_combinada <- (c_a + c_b) / (n_a + n_b)
  z <- (c_b / n_b - c_a / n_a) / sqrt(p_combinada * (1 - p_combinada) * (1 / n_a + 1 / n_b))
  2 * pnorm(-abs(z))
}

teste_classico <- prop.test(c(final$c_b, final$c_a), c(final$n_b, final$n_a), correct = FALSE)
stopifnot(abs(teste_classico$p.value - p_valor_z(final$c_a, final$n_a, final$c_b, final$n_b)) < 1e-10)

cat(sprintf("p-valor: %.4f\n", teste_classico$p.value))
cat(sprintf("Diferença B - A: %.2f p.p. (IC 95%%: %.2f a %.2f)\n",
            100 * (final$c_b / final$n_b - final$c_a / final$n_a),
            100 * teste_classico$conf.int[1], 100 * teste_classico$conf.int[2]))

#========================================================================================================
prior <- c(alfa = 1, beta = 1)

posterior <- function(conversoes, visitantes, prior) {
  c(alfa = prior[["alfa"]] + conversoes, beta = prior[["beta"]] + visitantes - conversoes)
}

resumo_bayesiano <- function(c_a, n_a, c_b, n_b, prior) {
  alfa_a <- prior[["alfa"]] + c_a
  beta_a <- prior[["beta"]] + n_a - c_a
  alfa_b <- prior[["alfa"]] + c_b
  beta_b <- prior[["beta"]] + n_b - c_b
  diferenca <- alfa_b / (alfa_b + beta_b) - alfa_a / (alfa_a + beta_a)
  desvio <- sqrt(alfa_a * beta_a / ((alfa_a + beta_a)^2 * (alfa_a + beta_a + 1)) +
                   alfa_b * beta_b / ((alfa_b + beta_b)^2 * (alfa_b + beta_b + 1)))
  z <- diferenca / desvio
  list(prob_b_melhor = pnorm(z),
       perda_escolher_b = desvio * dnorm(z) - diferenca * pnorm(-z),
       perda_escolher_a = desvio * dnorm(z) + diferenca * pnorm(z))
}

post_a <- posterior(final$c_a, final$n_a, prior)
post_b <- posterior(final$c_b, final$n_b, prior)

limites <- range(qbeta(c(0.0005, 0.9995), post_a[["alfa"]], post_a[["beta"]]),
                 qbeta(c(0.0005, 0.9995), post_b[["alfa"]], post_b[["beta"]]))
curve(dbeta(x, post_a[["alfa"]], post_a[["beta"]]), from = limites[1], to = limites[2], n = 500,
      lwd = 3, col = "grey40", xlab = "Taxa de conversão", ylab = "Densidade",
      main = "O que sabemos sobre cada taxa depois de 30 dias")
curve(dbeta(x, post_b[["alfa"]], post_b[["beta"]]), add = TRUE, n = 500, lwd = 3, col = "#377EB8")
abline(v = taxa_real, lty = 2, col = c("grey40", "#377EB8"))
legend("topright", c("A", "B", "Taxas reais (tracejado)"), col = c("grey40", "#377EB8", "black"),
       lwd = c(3, 3, 1), lty = c(1, 1, 2), bty = "n")
#===========================================================================================================
amostras_a <- rbeta(1e5, post_a[['alfa']], post_a[['beta']])
amostras_b <- rbeta(1e5, post_b[['alfa']], post_b[['beta']])
lift <- amostras_b / amostras_a - 1

prob_monte_carlo <- mean(amostras_b > amostras_a)
perda_monte_carlo <- mean(pmax(amostras_a - amostras_b, 0))
aproximacao <- resumo_bayesiano(final$c_a, final$n_a, final$c_b, final$n_b, prior)
stopifnot(abs(prob_monte_carlo - aproximacao$prob_b_melhor) < 0.01,
          abs(perda_monte_carlo - aproximacao$perda_escolher_b) < 5e-5)
cat(sprintf("P(B > A): %.1f%%\n", 100 * prob_monte_carlo))
cat(sprintf("Lift de B: %.1f%% (intervalo de 95%%: %.1f%% a %.1f%%)\n",
            100 * median(lift), 100 * quantile(lift, 0.025), 100 * quantile(lift, 0.975)))
cat(sprintf("Perda esperada se escolher B: %.4f p.p. | se escolher A: %.4f p.p.\n",
            100 * aproximacao$perda_escolher_b, 100 * aproximacao$perda_escolher_a))

hist(100 * lift, breaks = 80, col = "#377EB8", border = "white", freq = FALSE,
     xlab = "Lift de B sobre A (%)", ylab = "", yaxt = "n",
     main = sprintf("Distribuição do lift: P(B > A) = %.1f%%", 100 * prob_monte_carlo))
abline(v = 0, lwd = 2)
abline(v = 15, lwd = 2, lty = 2, col = "#E41A1C")
legend("topright", c("Sem efeito", "Lift real (15%)"), lty = c(1, 2), col = c("black", "#E41A1C"),
       lwd = 2, bty = "n")
#====================================================================================================
limiar_perda <- 0.0002

evolucao <- with(acumulado, data.frame( 
  dia, 
  as.data.frame(resumo_bayesiano(c_a, n_a, c_b, n_b, prior)), 
  p_valor= p_valor_z(c_a, n_a, c_b, n_b)
  ))
menor_perda <- pmin(evolucao$perda_escolher_b, evolucao$perda_escolher_a)
dia_decisao <- which(menor_perda < limiar_perda)[1]

if (is.na(dia_decisao)) {
  cat("A regra de perda esperada não atingiu o limiar em", dias, "dias\n")
} else {
  cat(sprintf("Regra de perda esperada: decide no dia %d, escolhendo %s\n", dia_decisao,
              ifelse(evolucao$perda_escolher_b[dia_decisao] < evolucao$perda_escolher_a[dia_decisao], "B", "A")))
}

options(repr.plot.width = 16, repr.plot.height = 5)
par(mfrow = c(1, 3))
plot(evolucao$dia, 100 * evolucao$prob_b_melhor, type = "b", pch = 19, ylim = c(0, 100),
     xlab = "Dia", ylab = "P(B > A) (%)", main = "Probabilidade de B ser melhor")
abline(h = 95, lty = 2, col = "#E41A1C")
plot(evolucao$dia, 100 * menor_perda, type = "b", pch = 19, log = "y",
     xlab = "Dia", ylab = "Perda esperada (p.p., escala log)", main = "Custo de decidir agora")
abline(h = 100 * limiar_perda, lty = 2, col = "#E41A1C")
plot(evolucao$dia, evolucao$p_valor, type = "b", pch = 19, ylim = c(0, 1),
     xlab = "Dia", ylab = "p-valor", main = "p-valor ao longo do teste")
abline(h = 0.05, lty = 2, col = "#E41A1C")
par(mfrow = c(1, 1))
#==============================================================================================
simular_regras <- function(taxas, n_simulacoes = 2000) {
  resultados <- replicate(n_simulacoes, {
    acum <- acumular(simular_experimento(taxas, dias, visitantes_por_dia))
    p <- with(acum, p_valor_z(c_a, n_a, c_b, n_b))
    bayes <- with(acum, resumo_bayesiano(c_a, n_a, c_b, n_b, prior))
    parada <- which(pmin(bayes$perda_escolher_b, bayes$perda_escolher_a) < limiar_perda)[1]
    if (is.na(parada)) parada <- dias
    c(p_valor_fixo = p[dias] < 0.05,
      p_valor_espiando = any(p < 0.05),
      prob_95_espiando = any(bayes$prob_b_melhor > 0.95 | bayes$prob_b_melhor < 0.05),
      dia_parada = parada,
      escolheu_b = bayes$perda_escolher_b[parada] < bayes$perda_escolher_a[parada])
  })
  as.data.frame(t(resultados))
}
set.seed(2024)
inicio <- Sys.time()
sem_diferenca <- simular_regras(c(0.040, 0.040))
com_diferenca <- simular_regras(taxa_real)
cat('Tempo: ', round(as.numeric(difftime(Sys.time(), inicio, units= 'secs'))), 'segundos\n\n')

tabela_regras <- data.frame(
  regra = c("p-valor só no fim (horizonte fixo)", "p-valor espiando todo dia", "P(B>A) > 95% espiando todo dia"),
  falso_positivo_aa = round(100 * colMeans(sem_diferenca[1:3]), 1),
  poder_ab = round(100 * colMeans(com_diferenca[1:3]), 1),
  row.names = NULL
)
print(tabela_regras)

cat(sprintf("\nRegra de perda esperada:\n  Teste A/B: decide em média no dia %.1f e escolhe B em %.1f%% dos casos\n",
            mean(com_diferenca$dia_parada), 100 * mean(com_diferenca$escolheu_b)))
cat(sprintf("  Teste A/A: decide em média no dia %.1f (escolher qualquer uma custa zero)\n",
            mean(sem_diferenca$dia_parada)))
#========================================================================================================================
priors <- list(
  uniforme = c(alfa= 1, beta= 1), 
  historica = c(alfa= 40, beta= 960),
  cetica_forte = c(alfa = 400, beta= 9600)
)
cores_prior <- c("black", "#377EB8", "#E41A1C")

options(repr.plot.width = 10, repr.plot.height = 5)
plot(NULL, xlim= c(1, dias), ylim= c(0, 100), xlab= 'Dia', ylab= 'P(B > A) (%)',
     main= 'Com dados suficientes, a escolha da prior deixa de importar')
for (i in seq_along(priors)) {
  prob <- with(acumulado, resumo_bayesiano(c_a, n_a, c_b, n_b, priors[[i]])$prob_b_melhor)
  lines(acumulado$dia, 100 * prob, lwd= 3, col= cores_prior[i])
}
abline(h= 50, col= 'grey70')
legend('bottomright', names(priors), col= cores_prior, lwd= 3, bty= 'n')
#===============================================================================================
salvar_png <- function(arquivo, desenhar, largura= 2000, altura= 900) {
  png(arquivo, width= largura, height = altura, res= 160)
  desenhar()
  invisible(dev.off())
  cat('Salvo: ', arquivo, '\n')
}
desenhar_resumo <- function() {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  hist(100 * lift, breaks = 80, col = "#377EB8", border = "white", freq = FALSE, yaxt = "n",
       xlab = "Lift de B sobre A (%)", ylab = "",
       main = sprintf("Teste A/B bayesiano: P(B > A) = %.1f%%", 100 * prob_monte_carlo))
  abline(v = 0, lwd = 2)
  barplot(tabela_regras$falso_positivo_aa,
          names.arg = c("p-valor\nsó no fim", "p-valor\nespiando", "P(B>A) > 95%\nespiando"),
          col = c("#4DAF4A", "#E41A1C", "#E41A1C"), border = NA, ylab = "Falsos positivos em testes A/A (%)",
          ylim = c(0, max(tabela_regras$falso_positivo_aa) * 1.2),
          main = "Espiar o teste todo dia multiplica os falsos positivos")
  abline(h = 5, lty = 2)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
}
options(repr.plot.width = 16, repr.plot.height = 7)
desenhar_resumo()
salvar_png('teste_ab_bayesiano.png', desenhar_resumo)
#=============================================================================================================================
