library(janeaustenr)

livros <- austen_books()
texto <- tolower(paste(livros$text, collapse = ' '))
texto <- gsub('_', ' ', texto)

tokens <- unlist(strsplit(texto, "[^a-z']+"))
tokens <- gsub("'s$", "", tokens)
tokens <- gsub("^'+|'+$", "", tokens)
tokens <- tokens[tokens != ""]

stopifnot(length(tokens) > 600000)
cat(sprintf("Livros: %d | Palavras no corpus: %d | Palavras distintas: %d\n",
            length(unique(livros$book)), length(tokens), length(unique(tokens))))
head(tokens, 30)
#=========================================================================================
tamanho_vocabulario <- 3000 
frequencias <- sort(table(tokens), decreasing= TRUE)
vocabulario <- names(frequencias)[seq_len(tamanho_vocabulario)]
ids <- match(tokens, vocabulario)


contar_coocorrencias <- function(ids, tamanho_vocabulario, janela) {
  n_celulas <- tamanho_vocabulario^2
  contagem <- numeric(n_celulas)
  n <- length(ids)
  for (distancia in seq_len(janela)) {
    centro <- ids[1:(n - distancia)]
    contexto <- ids[(1 + distancia):n]
    validos <- !is.na(centro) & !is.na(contexto)
    chave <- (centro[validos] - 1) * tamanho_vocabulario + contexto[validos]
    contagem <- contagem + tabulate(chave, nbins = n_celulas) / distancia
  }
  matriz <- matrix(contagem, tamanho_vocabulario, tamanho_vocabulario)
  matriz + t(matriz)
}
coocorrencia <- contar_coocorrencias(ids, tamanho_vocabulario, janela= 5)
dimnames(coocorrencia) <- list(vocabulario, vocabulario)

stopifnot(isSymmetric(coocorrencia), all(rowSums(coocorrencia) > 0))
cat(sprintf("Cobertura do vocabulário: %.1f%% dos tokens\n", 100 * mean(!is.na(ids))))
coocorrencia[c("elizabeth", "darcy", "happy", "money"), c("mr", "miss", "very", "pounds")]
#===============================================================================================
calcular_ppmi <- function(contagens, suavizacao= 0.75) { 
  total <- sum(contagens)
  prob_palavras <- rowSums(contagens) / total
  contexto_suavisado <- colSums(contagens)^suavizacao
  prob_contexto <- contexto_suavisado / sum(contexto_suavisado)
  pmi <- log((contagens / total) / outer(prob_palavras, prob_contexto))
  pmax(pmi, 0)
  }
ppmi <- calcular_ppmi(coocorrencia)
stopifnot(all(is.finite(ppmi)), all(ppmi >= 0))

cat(sprintf('Proporção de pares com associação positiva: %.1f%%\n', 100 * mean(ppmi > 0)))
round(ppmi[c('elizabeth', 'darcy', 'happy', 'money'), c('mr', 'miss', 'very', 'pounds')], 2)
#===============================================================================================
dimensoes <- 100 
inicio <- Sys.time()
decomposicao <- svd(ppmi, nu= dimensoes, nv= 0)
cat('SVD em', round(as.numeric(difftime(Sys.time(), inicio, units= 'secs'))), 'segundos\n')

embeddings <- sweep(decomposicao$u, 2, sqrt(decomposicao$d[seq_len(dimensoes)]), '*')
embeddings <- embeddings / sqrt(rowSums(embeddings^2))
rownames(embeddings) <- vocabulario
stopifnot(all(is.finite(embeddings)), all(abs(rowSums(embeddings^2) - 1) < 1e-10))

options(repr.plot.width = 10, repr.plot.height = 5)
plot(decomposicao$d[1:300], type = "l", lwd = 3, log = "y", col = "#377EB8",
     xlab = "Componente", ylab = "Valor singular (escala log)",
     main = "A informação se concentra nas primeiras dimensões")
abline(v = dimensoes, lty = 2)
#==================================================================================================
vizinhos <- function(palavra, vetores, n = 7) {
  if (!palavra %in% rownames(vetores)) stop("Palavra fora do vocabulário: ", palavra)
  similaridade <- as.numeric(vetores %*% vetores[palavra, ])
  names(similaridade) <- rownames(vetores)
  similaridade <- similaridade[names(similaridade) != palavra]
  head(sort(similaridade, decreasing = TRUE), n)
}
vetores_contagem <- coocorrencia / sqrt(rowSums(coocorrencia^2))

palavras_teste <- c('happy', 'letter', 'sister', 'morning', 'money', 'walk', 'handsome')
palavras_teste <- palavras_teste[palavras_teste %in% vocabulario]

for (palavra in palavras_teste) {
  cat(sprintf("\n%s\n", toupper(palavra)))
  cat("  contagem bruta: ", paste(names(vizinhos(palavra, vetores_contagem)), collapse = ", "), "\n")
  cat("  embeddings:     ", paste(names(vizinhos(palavra, embeddings)), collapse = ", "), "\n")
}
#============================================================================================================
analogia <- function(a, b, c, n = 3) {
  palavras <- c(a, b, c)
  ausentes <- palavras[!palavras %in% rownames(embeddings)]
  if (length(ausentes) > 0) return(paste("fora do vocabulário:", paste(ausentes, collapse = ", ")))
  alvo <- embeddings[b, ] - embeddings[a, ] + embeddings[c, ]
  similaridade <- as.numeric(embeddings %*% alvo)
  names(similaridade) <- rownames(embeddings)
  similaridade <- similaridade[!names(similaridade) %in% palavras]
  paste(names(head(sort(similaridade, decreasing = TRUE), n)), collapse = ", ")
}
analogias <- data.frame(
  a= c('he', 'brother', 'father', 'his', 'day', 'mr', 'good'), 
  b= c('she', 'sister', 'mother', 'her', 'days', 'mrs', 'better'), 
  c= c('him', 'son', 'uncle', 'himself', 'year', 'sir', 'bad'), 
  esperado= c('her', 'daughter', 'aunt', 'herself', 'years', 'lady', 'worse')
)
analogias$top3 <- mapply(analogia, analogias$a, analogias$b, analogias$c)
analogias$acertou <- mapply(function(top3, esperado) esperado %in% strsplit(top3, ', ')[[1]], 
                            analogias$top3, analogias$esperado)
print(analogias, row.names= FALSE)
cat(sprintf("\nAcertos no top 3: %d de %d\n", sum(analogias$acertou), nrow(analogias)))
#=================================================================================================
grupos <- list( 
  familia= c('father', 'mother', 'brother', 'sister', 'son', 'daughter', 'uncle', 'aunt', 'husband', 'wife'), 
  sentimentos= c('happy', 'glad', 'pleased', 'sorry', 'unhappy', 'miserable', 'angry', 'delighted'), 
  tempo= c('morning', 'evening', 'night', 'day', 'week', 'month', 'year', 'hour'), 
  lugares= c('house', 'room', 'park', 'town', 'village', 'church', 'garden', 'street'), 
  dinheiro= c('money', 'fortune', 'income', 'pounds', 'rich', 'poor', 'thousand', 'hundred'))

cores_grupo <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

palavras_mapa <- unlist(grupos)
grupo_de <- rep(seq_along(grupos), lengths(grupos))
no_vocabulario <- palavras_mapa %in% vocabulario
palavras_mapa <- palavras_mapa[no_vocabulario]
grupo_de <- grupo_de[no_vocabulario]

desenhar_mapa <- function() {
  coordenadas <- prcomp(embeddings[palavras_mapa, ])$x[, 1:2]
  plot(coordenadas, type = "n", xlab = "Componente 1", ylab = "Componente 2",
       main = "Mapa semântico aprendido de Jane Austen, sem nenhum rótulo")
  text(coordenadas, labels = palavras_mapa, col = cores_grupo[grupo_de], cex = 0.9, font = 2)
  legend("bottomright", names(grupos), text.col = cores_grupo, bty = "n")
}
options(repr.plot.width =  11, repr.plot.height = 8)
desenhar_mapa()
#=========================================================================================================
pares_genero <- rbind(c("he", "she"), c("him", "her"), c("man", "woman"),
                      c("mr", "mrs"), c("brother", "sister"), c("father", "mother"))
pares_genero <- pares_genero[pares_genero[, 1] %in% vocabulario & pares_genero[, 2] %in% vocabulario, ]

eixo_genero <- colMeans(embeddings[pares_genero[, 2], ] - embeddings[pares_genero[, 1], ])
eixo_genero <- eixo_genero / sqrt(sum(eixo_genero^2))

adjetivos <- c("handsome", "pretty", "beautiful", "elegant", "clever", "sensible", "agreeable",
               "amiable", "proud", "rich", "foolish", "silly", "charming", "lively", "sweet",
               "gentlemanlike", "respectable", "accomplished", "good", "fine")
adjetivos <- adjetivos[adjetivos %in% vocabulario]
projecao <- sort(as.numeric(embeddings[adjetivos, ] %*% eixo_genero))
names(projecao) <- adjetivos[order(as.numeric(embeddings[adjetivos, ] %*% eixo_genero))]

desenhar_vies <- function() {
  par(mar = c(5, 9, 4, 2))
  barplot(projecao, horiz = TRUE, las = 1, border = NA,
          col = ifelse(projecao > 0, "#E41A1C", "#377EB8"),
          xlab = "<- contexto masculino        contexto feminino ->",
          main = "Viés de gênero que o modelo aprendeu do corpus")
  abline(v = 0)
  par(mar = c(5, 4, 4, 2))
}

desenhar_vies()
#=======================================================================================
salvar_png <- function(arquivo, desenhar, largura= 2200, altura= 1000) { 
  png(arquivo, width= largura, height= altura, res= 160 )
  desenhar()
  invisible(dev.off())
  cat('Salvo: ', arquivo, '\n')
}
desenhar_resumo <- function() {
  par(mfrow= c(1, 2))
  desenhar_mapa()
  desenhar_vies()
  par(mfrow= c(1, 1))
}
options(repr.plot.width= 18, repr.plot.height = 8)
desenhar_resumo()
salvar_png('embeddings_jane_austen.png', desenhar_resumo)
#=======================================================================================
