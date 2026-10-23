set.seed(42)
tamanho_vocabulario <- 10
comprimento <- 8

gerar_entradas <- function(n) {
  matrix(sample(tamanho_vocabulario, n * comprimento, replace = TRUE), n, comprimento)
}

gerar_lote_inversao <- function(n) {
  entradas <- gerar_entradas(n)
  list(entradas = entradas, saidas = entradas[, comprimento:1, drop = FALSE])
}

gerar_lote_maximo <- function(n) {
  entradas <- gerar_entradas(n)
  list(entradas = entradas, saidas = matrix(apply(entradas, 1, max), n, comprimento))
}

exemplo <- gerar_lote_inversao(2)
stopifnot(all(exemplo$saidas[1, ] == rev(exemplo$entradas[1, ])))
cat("Inversão  | entrada:", exemplo$entradas[1, ] - 1, "| saída:", exemplo$saidas[1, ] - 1, "\n")
exemplo <- gerar_lote_maximo(2)
cat("Máximo    | entrada:", exemplo$entradas[1, ] - 1, "| saída:", exemplo$saidas[1, ] - 1, "\n")

#=============================================================================================================================
softmax_linhas <- function(m) {
  m <- m - apply(m, 1, max)
  exponenciais <- exp(m)
  exponenciais / rowSums(exponenciais)
}

iniciar_modelo <- function(tamanho_vocabulario, comprimento, dimensao) {
  aleatoria <- function(linhas, colunas) {
    matrix(rnorm(linhas * colunas, sd = 1 / sqrt(dimensao)), linhas, colunas)
  }
  list(E = aleatoria(tamanho_vocabulario, dimensao), P = aleatoria(comprimento, dimensao),
       Wq = aleatoria(dimensao, dimensao), Wk = aleatoria(dimensao, dimensao),
       Wv = aleatoria(dimensao, dimensao), Wo = aleatoria(dimensao, tamanho_vocabulario),
       bo = rep(0, tamanho_vocabulario))
}

propagar_sequencia <- function(modelo, tokens) {
  escala <- 1 / sqrt(ncol(modelo$Wq))
  X <- modelo$E[tokens, , drop = FALSE] + modelo$P
  Q <- X %*% modelo$Wq
  K <- X %*% modelo$Wk
  V <- X %*% modelo$Wv
  A <- softmax_linhas(Q %*% t(K) * escala)
  H <- A %*% V
  prob <- softmax_linhas(sweep(H %*% modelo$Wo, 2, modelo$bo, "+"))
  list(X = X, Q = Q, K = K, V = V, A = A, H = H, prob = prob, escala = escala)
}

calcular_perda <- function(modelo, entradas, saidas) {
  perdas <- sapply(seq_len(nrow(entradas)), function(i) {
    prob <- propagar_sequencia(modelo, entradas[i, ])$prob
    -mean(log(prob[cbind(seq_len(ncol(saidas)), saidas[i, ])] + 1e-12))
  })
  mean(perdas)
}

prever <- function(modelo, tokens) {
  max.col(propagar_sequencia(modelo, tokens)$prob, ties.method = "first")
}

avaliar <- function(modelo, lote) {
  acertos <- t(sapply(seq_len(nrow(lote$entradas)), function(i) {
    prever(modelo, lote$entradas[i, ]) == lote$saidas[i, ]
  }))
  c(por_token = mean(acertos), por_sequencia = mean(apply(acertos, 1, all)))
}

modelo_teste <- iniciar_modelo(tamanho_vocabulario, comprimento, 16)
saida_teste <- propagar_sequencia(modelo_teste, 1:8)
stopifnot(all(abs(rowSums(saida_teste$A) - 1) < 1e-12), all(dim(saida_teste$prob) == c(8, 10)))
cat("Camada de atenção verificada\n")
#================================================================================================================
calcular_gradientes <- function(modelo, entradas, saidas) {
  gradientes <- lapply(modelo, function(parametro) parametro * 0)
  n_total <- length(saidas)
  for (i in seq_len(nrow(entradas))) {
    tokens <- entradas[i, ]
    f <- propagar_sequencia(modelo, tokens)
    
    dZ <- f$prob
    corretos <- cbind(seq_along(tokens), saidas[i, ])
    dZ[corretos] <- dZ[corretos] - 1
    dZ <- dZ / n_total
    gradientes$Wo <- gradientes$Wo + t(f$H) %*% dZ
    gradientes$bo <- gradientes$bo + colSums(dZ)
    
    dH <- dZ %*% t(modelo$Wo)
    dA <- dH %*% t(f$V)
    dV <- t(f$A) %*% dH
    dS <- f$A * (dA - rowSums(dA * f$A))
    dQ <- dS %*% f$K * f$escala
    dK <- t(dS) %*% f$Q * f$escala
    
    gradientes$Wq <- gradientes$Wq + t(f$X) %*% dQ
    gradientes$Wk <- gradientes$Wk + t(f$X) %*% dK
    gradientes$Wv <- gradientes$Wv + t(f$X) %*% dV
    
    dX <- dQ %*% t(modelo$Wq) + dK %*% t(modelo$Wk) + dV %*% t(modelo$Wv)
    gradientes$P <- gradientes$P + dX
    por_token <- rowsum(dX, tokens)
    linhas <- as.integer(rownames(por_token))
    gradientes$E[linhas, ] <- gradientes$E[linhas, ] + por_token
  }
  gradientes
}

verificar_gradientes <- function() {
  modelo <- iniciar_modelo(5, 4, 3)
  lote <- list(entradas = matrix(sample(5, 12, replace= TRUE), 3, 4), 
               saidas= matrix(sample(5, 12, replace= TRUE), 3, 4))
  analitico <- calcular_gradientes(modelo, lote$entradas, lote$saidas)
  erros <- c()
  for (nome in names(modelo)) {
    for (j in seq_along(modelo[[nome]])) {
      mais <- modelo
      mais[[nome]][[j]] <- mais[[nome]][[j]] + 1e-5 
      menos <- modelo
      menos[[nome]][[j]] <- menos[[nome]][[j]] - 1e-5
      numerico <- (calcular_perda(mais, lote$entradas, lote$saidas) - 
                     calcular_perda(menos, lote$entradas, lote$saidas)) / 2e-5
      derivada <- analitico[[nome]][[j]]
      erros <- c(erros, abs(numerico - derivada) / max(1e-8, abs(numerico) + abs(derivada)))
    }
  }
  erros
}
erros_relativos <- verificar_gradientes()
cat(sprintf('PARAMETROS VERIFICADOS: %d | MAIOR ERRO RELATIVO: %.2e\n', 
            length(erros_relativos), max(erros_relativos)))
stopifnot(max(erros_relativos) < 1e-4)
#================================================================================================================
treinar <- function(modelo, gerar_lote, passos = 2000, tamanho_lote = 32, taxa = 0.01,
                    congelados = character(0), registrar_a_cada = 50) {
  momento1 <- lapply(modelo, function(parametro) parametro * 0)
  momento2 <- momento1
  historico <- data.frame()
  for (passo in seq_len(passos)) {
    lote <- gerar_lote(tamanho_lote)
    gradientes <- calcular_gradientes(modelo, lote$entradas, lote$saidas)
    for (nome in setdiff(names(modelo), congelados)) {
      momento1[[nome]] <- 0.9 * momento1[[nome]] + 0.1 * gradientes[[nome]]
      momento2[[nome]] <- 0.999 * momento2[[nome]] + 0.001 * gradientes[[nome]]^2
      m_corrigido <- momento1[[nome]] / (1 - 0.9^passo)
      v_corrigido <- momento2[[nome]] / (1 - 0.999^passo)
      modelo[[nome]] <- modelo[[nome]] - taxa * m_corrigido / (sqrt(v_corrigido) + 1e-8)
    }
    if (passo %% registrar_a_cada == 0) {
      validacao <- gerar_lote(200)
      historico <- rbind(historico, data.frame(
        passo = passo, perda = calcular_perda(modelo, validacao$entradas, validacao$saidas),
        t(avaliar(modelo, validacao))))
    }
  }
  list(modelo = modelo, historico = historico)
}

dimensao <- 32
set.seed(1)
inicio <- Sys.time()
inversao <- treinar(iniciar_modelo(tamanho_vocabulario, comprimento, dimensao), gerar_lote_inversao)
cat(sprintf("Treino em %.0f segundos\n", as.numeric(difftime(Sys.time(), inicio, units = "secs"))))
stopifnot(is.list(inversao), nrow(inversao$historico) > 0)

resultado_inversao <- avaliar(inversao$modelo, gerar_lote_inversao(1000))
cat(sprintf("Inversão | acerto por dígito: %.1f%% | sequências inteiras corretas: %.1f%%\n",
            100 * resultado_inversao[["por_token"]], 100 * resultado_inversao[["por_sequencia"]]))

options(repr.plot.width = 10, repr.plot.height = 5)
plot(inversao$historico$passo, 100 * inversao$historico$por_sequencia, type = "l", lwd = 3,
     col = "#377EB8", ylim = c(0, 100), xlab = "Passo de treino",
     ylab = "Sequências 100% corretas (%)", main = "Aprendendo a inverter sequências")
#=============================================================================================================
desenhar_atencao <- function(modelo, tokens, titulo) {
  A <- propagar_sequencia(modelo, tokens)$A
  L <- length(tokens)
  image(1:L, 1:L, t(A[L:1, ]), col = hcl.colors(50, "Blues 3", rev = TRUE), zlim = c(0, 1),
        axes = FALSE, xlab = "Posição consultada", ylab = "Posição que consulta", main = titulo)
  axis(1, at = 1:L, labels = tokens - 1, tick = FALSE)
  axis(2, at = L:1, labels = tokens - 1, las = 1, tick = FALSE)
}
exemplos <- gerar_lote_inversao(3)$entradas
options(repr.plot.width = 15, repr.plot.height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (i in 1:3) {
  previsto <- prever(inversao$modelo, exemplos[i, ]) - 1
  desenhar_atencao(inversao$modelo, exemplos[i, ], paste("Saída:", paste(previsto, collapse = " ")))
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
#========================================================================================================
set.seed(2)
maximo <- treinar(iniciar_modelo(tamanho_vocabulario, comprimento, dimensao), gerar_lote_maximo)

resultado_maximo <- avaliar(maximo$modelo, gerar_lote_maximo(1000))
lote_referencia <- gerar_lote_maximo(1000)
sempre_nove <- mean(lote_referencia$saidas[, 1] == tamanho_vocabulario)
cat(sprintf("Máximo | acerto: %.1f%% | referência (chutar sempre 9): %.1f%%\n",
            100 * resultado_maximo[["por_sequencia"]], 100 * sempre_nove))

exemplos_maximo <- rbind(c(3, 1, 4, 1, 5, 9, 2, 6), c(7, 2, 8, 1, 8, 2, 8, 4), c(1, 4, 1, 4, 2, 1, 3, 5)) + 1
options(repr.plot.width= 15, repr.plot.height= 5)
par(mfrow= c(1, 3), mar= c(4.5, 4.5, 3, 1))
for (i in 1:3) {
  previsto <- prever(maximo$modelo, exemplos_maximo[i, ])[1] - 1
  desenhar_atencao(maximo$modelo, exemplos_maximo[i, ], paste('MAXIMO PREVISTO: ', previsto))
}
par(mfrow= c(1, 1), mar= c(5, 4 , 4, 2))
#=============================================================================================================
set.seed(1)
modelo_sem_posicao <- iniciar_modelo(tamanho_vocabulario, comprimento, dimensao)
modelo_sem_posicao$P[] <- 0 
sem_posicao <- treinar(modelo_sem_posicao, gerar_lote_inversao, congelados ='P')

resultado_sem_posicao <- avaliar(sem_posicao$modelo, gerar_lote_inversao(1000))
cat(sprintf('SEM POSICAO | ACERTO POR DIGITO: %.1f%% | SEQUENCIA CORRETAS: %.1f%%\n', 
            100 * resultado_sem_posicao[['por_token']], 100 * resultado_sem_posicao[['por_sequencia']]))

options(repr.plot.width = 10, repr.plot.height = 5)
plot(inversao$historico$passo, 100 * inversao$historico$por_token, type = "l", lwd = 3,
     col = "#377EB8", ylim = c(0, 100), xlab = "Passo de treino", ylab = "Acerto por dígito (%)",
     main = "Inversão com e sem embeddings de posição")
lines(sem_posicao$historico$passo, 100 * sem_posicao$historico$por_token, lwd = 3, col = "#E41A1C")
abline(h = 10, lty = 2, col = "grey50")
legend("right", c("Com posição", "Sem posição", "Chute aleatório (10%)"),
       col = c("#377EB8", "#E41A1C", "grey50"), lwd = c(3, 3, 1), lty = c(1, 1, 2), bty = "n")
#=============================================================================================================
salvar_png <- function(arquivo, desenhar, largura= 2400, altura= 800 ) { 
  png(arquivo, width= largura, height= altura, res= 160)
  desenhar()
  invisible(dev.off())
  cat('Salvo:', arquivo, '\n')
}
desenhar_resumo <- function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 4, 1))
  desenhar_atencao(inversao$modelo, exemplos[1, ], "Inverter: atenção por POSIÇÃO")
  desenhar_atencao(maximo$modelo, exemplos_maximo[1, ], "Encontrar o máximo: atenção por CONTEÚDO")
  plot(inversao$historico$passo, 100 * inversao$historico$por_token, type = "l", lwd = 3,
       col = "#377EB8", ylim = c(0, 100), xlab = "Passo de treino", ylab = "Acerto por dígito (%)",
       main = "Sem embeddings de posição, inverter é impossível")
  lines(sem_posicao$historico$passo, 100 * sem_posicao$historico$por_token, lwd = 3, col = "#E41A1C")
  legend("right", c("Com posição", "Sem posição"), col = c("#377EB8", "#E41A1C"), lwd = 3, bty = "n")
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
}
options(repr.plot.width= 18, repr.plot.height= 6)
desenhar_resumo()
salvar_png('atencao_do_zero.png', desenhar_resumo)
