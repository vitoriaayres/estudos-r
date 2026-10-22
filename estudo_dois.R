# Rede neural do zero em R base, só com álgebra de matrizes
set.seed(7)
cores <- c("#E41A1C", "#377EB8", "#4DAF4A")

gerar_espiral <- function(n_por_classe, n_classes, ruido) {
  x <- matrix(0, n_por_classe * n_classes, 2)
  y <- integer(n_por_classe * n_classes)
  for (k in seq_len(n_classes)) {
    linhas <- ((k - 1) * n_por_classe + 1):(k * n_por_classe)
    raio <- seq(0.05, 1, length.out = n_por_classe)
    angulo <- seq((k - 1) * 4, k * 4, length.out = n_por_classe) + rnorm(n_por_classe, sd = ruido)
    x[linhas, ] <- cbind(raio * sin(angulo), raio * cos(angulo))
    y[linhas] <- k
  }
  list(x = x, y = y)
}

treino <- gerar_espiral(200, 3, 0.25)
teste <- gerar_espiral(200, 3, 0.25)
stopifnot(nrow(treino$x) == 600, all(table(treino$y) == 200))

options(repr.plot.width = 7, repr.plot.height = 7)
plot(treino$x, col = cores[treino$y], pch = 19, cex = 0.7, asp = 1,
     xlab = "", ylab = "", main = "Três espirais entrelaçadas")
#======================================================================================================

relu <- function(z) pmax(z, 0)
one_hot <- function(y, n_classes) {
  matriz <- matrix(0, length(y), n_classes)
  matriz[cbind(seq_along(y), y)] <- 1
  matriz
}
softmax <- function(z) {
  z <- z - apply(z, 1, max)
  exponenciais <- exp(z)
  exponenciais / rowSums(exponenciais)
}
entropia_cruzada <- function(prob, y_one_hot) {
  -mean(rowSums(y_one_hot * log(prob + 1e-12)))
}

stopifnot(all(one_hot(c(2, 1), 3) == rbind(c(0, 1, 0), c(1, 0, 0))))
stopifnot(all(abs(rowSums(softmax(matrix(c(1000, 0, 0, 1000), 2))) - 1) < 1e-12))
stopifnot(all(relu(matrix(c(-1, 2), 1)) == c(0, 2)))
cat("Blocos básicos verificados\n")
#=======================================================================================================

iniciar_rede <- function(tamanhos) {
  lapply(seq_len(length(tamanhos) - 1), function(i) {
    list(W = matrix(rnorm(tamanhos[i] * tamanhos[i + 1], sd = sqrt(2 / tamanhos[i])),
                    tamanhos[i], tamanhos[i + 1]),
         b = rep(0, tamanhos[i + 1]))
  })
}

propagar <- function(rede, x) {
  ativacoes <- list(x)
  for (i in seq_along(rede)) {
    z <- sweep(ativacoes[[i]] %*% rede[[i]]$W, 2, rede[[i]]$b, "+")
    ativacoes[[i + 1]] <- if (i < length(rede)) relu(z) else softmax(z)
  }
  ativacoes
}

retropropagar <- function(rede, ativacoes, y_one_hot, lambda) {
  n <- nrow(y_one_hot)
  gradientes <- vector("list", length(rede))
  delta <- (ativacoes[[length(ativacoes)]] - y_one_hot) / n
  for (i in rev(seq_along(rede))) {
    gradientes[[i]] <- list(W = t(ativacoes[[i]]) %*% delta + lambda * rede[[i]]$W,
                            b = colSums(delta))
    if (i > 1) delta <- (delta %*% t(rede[[i]]$W)) * (ativacoes[[i]] > 0)
  }
  gradientes
}

calcular_perda <- function(rede, x, y_one_hot, lambda) {
  prob <- propagar(rede, x)[[length(rede) + 1]]
  penalidade <- sum(sapply(rede, function(camada) sum(camada$W^2)))
  entropia_cruzada(prob, y_one_hot) + lambda / 2 * penalidade
}

prever_classe <- function(rede, x) {
  max.col(propagar(rede, x)[[length(rede) + 1]], ties.method = "first")
}

acuracia <- function(rede, dados) mean(prever_classe(rede, dados$x) == dados$y)

rede_teste <- iniciar_rede(c(2, 8, 3))
stopifnot(length(rede_teste) == 2, all(dim(rede_teste[[1]]$W) == c(2, 8)))
stopifnot(all(dim(propagar(rede_teste, treino$x)[[3]]) == c(600, 3)))
cat("Estrutura da rede verificada\n")
#====================================================================================================

verificar_gradientes <- function() {
  rede <- iniciar_rede(c(2, 5, 3))
  x <- matrix(rnorm(20), 10, 2)
  y <- one_hot(sample(3, 10, replace= TRUE), 3)
  lambda <- 0.01
  analitico <- retropropagar(rede, propagar(rede, x), y, lambda)
  erros <- c()
  for (i in seq_along(rede)) {
    for (parametro in c('W', 'b')) {
      for (j in seq_along(rede[[i]][[parametro]])) {
        mais <- rede
        mais[[i]][[parametro]][j] <- mais[[i]][[parametro]][j] + 1e-5
        menos <- rede
        menos[[i]][[parametro]][j] <- menos[[i]][[parametro]][j] - 1e-5
        numerico <- (calcular_perda(mais, x, y, lambda) - calcular_perda(menos, x, y, lambda)) / 2e-5
        derivada <- analitico[[i]][[parametro]][j]
        erros <- c(erros, abs(numerico - derivada) / max(1e-8, abs(numerico) + abs(derivada)))
      }
    }
  }
  erros
}
erros_relativos <- verificar_gradientes()
cat(sprintf('PARÂMETROS VERIFICADOS: %d | MAIOR ERRO RELATIVO: %.2e\n',
            length(erros_relativos), max(erros_relativos)))
stopifnot(max(erros_relativos) < 1e-4)
#============================================================================================================

treinar <- function(tamanhos, dados, epocas= 400, taxa= 0.01, lambda= 1e-4,
                    tamanho_lote= 64, otimizador= 'adam', salvar_em= integer(0)) {
  rede <- iniciar_rede(tamanhos)
  y <- one_hot(dados$y, tamanhos[length(tamanhos)])
  momento1 <- lapply(rede, function(camada) lapply(camada, function(p) p * 0))
  momento2 <- momento1
  passo <- 0
  historico <- data.frame()
  instantaneos <- list()
  for (epoca in seq_len(epocas)) {
    ordem <- sample(nrow(dados$x))
    for (inicio in seq(1, length(ordem), by = tamanho_lote)) {
      lote <- ordem[inicio:min(inicio + tamanho_lote - 1, length(ordem))]
      ativacoes <- propagar(rede, dados$x[lote, , drop = FALSE])
      gradientes <- retropropagar(rede, ativacoes, y[lote, , drop = FALSE], lambda)
      passo <- passo + 1
      for (i in seq_along(rede)) {
        for (p in c("W", "b")) {
          g <- gradientes[[i]][[p]]
          if (otimizador == "adam") {
            momento1[[i]][[p]] <- 0.9 * momento1[[i]][[p]] + 0.1 * g
            momento2[[i]][[p]] <- 0.999 * momento2[[i]][[p]] + 0.001 * g^2
            m_corrigido <- momento1[[i]][[p]] / (1 - 0.9^passo)
            v_corrigido <- momento2[[i]][[p]] / (1 - 0.999^passo)
            rede[[i]][[p]] <- rede[[i]][[p]] - taxa * m_corrigido / (sqrt(v_corrigido) + 1e-8)
          } else {
            rede[[i]][[p]] <- rede[[i]][[p]] - taxa * g
          }
        }
      }
    }
    if (epoca == 1 || epoca %% 10 == 0) {
      historico <- rbind(historico, data.frame(
        epoca = epoca, perda = calcular_perda(rede, dados$x, y, lambda),
        acuracia = acuracia(rede, dados)))
    }
    if (epoca %in% salvar_em) instantaneos[[as.character(epoca)]] <- rede
  }
  list(rede = rede, historico = historico, instantaneos = instantaneos)
}

set.seed(1)
inicio_treino <- Sys.time()
profunda <- treinar(c(2, 64, 64, 3), treino, salvar_em = c(1, 20, 80, 400))
cat(sprintf("Treino em %.0f segundos | Acurácia no treino: %.1f%% | no teste: %.1f%%\n",
            as.numeric(difftime(Sys.time(), inicio_treino, units = "secs")),
            100 * acuracia(profunda$rede, treino), 100 * acuracia(profunda$rede, teste)))
stopifnot(is.list(profunda), nrow(profunda$historico) > 0)

#=======================================================================================================
desenhar_fronteira <- function(rede, dados, titulo) {
  eixo <- seq(-1.2, 1.2, length.out = 150)
  grade <- as.matrix(expand.grid(eixo, eixo))
  classe <- prever_classe(rede, grade)
  image(eixo, eixo, matrix(classe, 150, 150), col = adjustcolor(cores, 0.35),
        breaks = c(0.5, 1.5, 2.5, 3.5), xlab = "", ylab = "", main = titulo, axes = FALSE, asp = 1)
  points(dados$x, col = cores[dados$y], pch = 19, cex = 0.5)
}

set.seed(1)
linear <- treinar(c(2, 3), treino, epocas = 200)

options(repr.plot.width = 14, repr.plot.height = 7)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
desenhar_fronteira(linear$rede, teste,
                   sprintf("Modelo linear: %.0f%% no teste", 100 * acuracia(linear$rede, teste)))
desenhar_fronteira(profunda$rede, teste,
                   sprintf("Rede com 2 camadas ocultas: %.0f%% no teste", 100 * acuracia(profunda$rede, teste)))
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

#========================================================================================================================

salvar_png <- function(arquivo, desenhar, largura= 2000, altura= 1100) {
  png(arquivo, width= largura, height = altura, res= 160)
  desenhar()
  invisible(dev.off())
  cat('Gráfico salvo em', arquivo, '\n')
}
desenhar_evolucao <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5, 5, 5), 2, byrow = TRUE), heights = c(1.3, 1))
  par(mar = c(1, 1, 3, 1))
  for (epoca in names(profunda$instantaneos)) {
    rede <- profunda$instantaneos[[epoca]]
    desenhar_fronteira(rede, treino, sprintf("Época %s: %.0f%% de acerto", epoca, 100 * acuracia(rede, treino)))
  }
  par(mar = c(4.5, 4.5, 2, 1))
  plot(profunda$historico$epoca, profunda$historico$perda, type = "l", lwd = 3, log = "y",
       col = "#377EB8", xlab = "Época", ylab = "Perda (escala log)",
       main = "Rede neural em R base: backpropagation e Adam implementados do zero")
  layout(1)
  par(mar = c(5, 4, 4, 2))
}

options(repr.plot.width = 16, repr.plot.height = 9)
desenhar_evolucao()
salvar_png("rede_neural_evolucao.png", desenhar_evolucao)
#==================================================================================================================
set.seed(1)
com_adam <- treinar(c(2, 64, 64, 3), treino, taxa= 0.01, otimizador = 'adam')
set.seed(1)
com_sgd <- treinar(c(2, 64, 64, 3), treino, taxa= 0.3, otimizador= 'sgd')

options(repr.plot.width = 10, repr.plot.height = 6)
plot(com_adam$historico$epoca, com_adam$historico$perda, type = "l", lwd = 3, log = "y",
     col = "#377EB8", ylim = range(c(com_adam$historico$perda, com_sgd$historico$perda)),
     xlab = "Época", ylab = "Perda (escala log)", main = "Mesma rede, mesma inicialização: Adam contra SGD")
lines(com_sgd$historico$epoca, com_sgd$historico$perda, lwd = 3, col = "#E41A1C")
legend("topright", sprintf("%s: %.1f%% no teste", c("Adam (taxa 0,01)", "SGD (taxa 0,3)"),
                           100 * c(acuracia(com_adam$rede, teste), acuracia(com_sgd$rede, teste))),
       col = c("#377EB8", "#E41A1C"), lwd = 3, bty = "n")
#====================================================================================================================
set.seed(3)
poucos_dados <- gerar_espiral(25, 3, 0.35)

set.seed(1)
sem_regularizacao <- treinar(c(2, 128, 128, 3), poucos_dados, epocas= 1500, lambda = 0)
set.seed(1)
com_regularizacao <- treinar(c(2, 128, 128, 3), poucos_dados, epocas= 1500, lambda = 3e-3)

descrever <- function(modelo) {
  sprintf('treino %.0f%%  teste %.0f%%',
          100 * acuracia(modelo$rede, poucos_dados), 100 * acuracia(modelo$rede, teste))
}

options(repr.plot.width = 14, repr.plot.height = 7)
par(mfrow= c(1, 2), mar= c(1, 1, 3, 1))
desenhar_fronteira(sem_regularizacao$rede, poucos_dados, paste('SEM REGULARIZAÇÃO:', descrever(sem_regularizacao)))
desenhar_fronteira(com_regularizacao$rede, poucos_dados, paste('COM L2:', descrever(com_regularizacao)))
par(mfrow= c(1, 1), mar= c(5, 4, 4, 2))
#=====================================================================================================================

desenhar_neuronios <- function(rede, camada, neuronios) {
  eixo <- seq(-1.2, 1.2, length.out= 120)
  grade <- as.matrix(expand.grid(eixo, eixo))
  ativacao <- propagar(rede, grade)[[camada + 1]]
  for (neuronio in neuronios) {
    image(eixo, eixo, matrix(ativacao[, neuronio], 120, 120), col= hcl.colors(50, 'Inferno'),
          axes= FALSE, xlab= ' ', ylab= ' ', asp= 1,
          main= sprintf('Camada %d, neurônio %d', camada, neuronio))
  }
}
ativos_camada2 <- order(-colMeans(propagar(profunda$rede, treino$x)[[3]]))[1:4]

options(repr.plot.width= 16, repr.plot.height= 8)
par(mfrow= c(2, 4), mar= c(1, 1, 3, 1))
desenhar_neuronios(profunda$rede, 1, 1:4)
desenhar_neuronios(profunda$rede, 2, ativos_camada2)
par(mfrow= c(1, 1), mar= c(5, 4, 4, 2))
