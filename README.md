# Estudos em R

Os cinco primeiros estudos são escritos em R base, sem bibliotecas de modelagem prontas. A partir do sexto, a série entra em estatística espacial e machine learning geoespacial com dados públicos brasileiros, usando o ecossistema `sf` e `spdep`.

## 1. Inferência causal com os dados de LaLonde (`estudo_um.R`)

Um programa de treinamento profissional (NSW) foi avaliado por sorteio, o que dá o efeito real sobre a renda: **US$ 1.794**. O estudo troca o grupo de controle sorteado por pessoas do censo (CPS) e mostra como diferentes métodos se saem ao tentar recuperar esse efeito a partir de dados observacionais:

- comparação ingênua de médias;
- regressão linear com covariáveis;
- pareamento por escore de propensão (vizinho mais próximo com caliper);
- ponderação pelo inverso da probabilidade (IPW);
- estimador duplamente robusto (regressão ponderada pelo IPW).

Os intervalos de confiança vêm de 200 reamostragens bootstrap.

![Comparação dos métodos de inferência causal](inferencia_causal_resultado.png)

| Método | Estimativa (US$) | Distância do efeito real |
|---|---:|---:|
| Comparação ingênua | -8.498 | -10.292 |
| Regressão linear | 1.137 | -657 |
| Pareamento por escore | 952 | -842 |
| IPW | 1.782 | -12 |
| Duplamente robusto | 1.702 | -92 |

## 2. Rede neural do zero (`estudo_dois.R`)

Uma rede neural totalmente conectada implementada só com álgebra de matrizes, classificando três espirais entrelaçadas:

- propagação com ReLU e softmax, perda de entropia cruzada com regularização L2;
- backpropagation verificado contra gradientes numéricos (erro relativo máximo na ordem de 1e-9);
- otimizadores Adam e SGD com mini-lotes;
- comparação entre modelo linear e rede profunda, Adam e SGD, e com e sem regularização;
- visualização das fronteiras de decisão e das ativações dos neurônios ocultos.

![Evolução da fronteira de decisão durante o treino](rede_neural_evolucao.png)

## 3. Teste A/B bayesiano (`estudo_tres.R`)

Um teste A/B simulado com efeito conhecido (conversão de 4,0% contra 4,6%), analisado pelas abordagens clássica e bayesiana:

- teste de duas proporções e p-valor;
- modelo Beta-Binomial conjugado, probabilidade de B ser melhor e distribuição do lift por Monte Carlo;
- perda esperada em forma fechada como regra de decisão, conferida contra a simulação;
- 2.000 testes A/A e A/B simulados para medir quanto "espiar" o resultado todo dia inflaciona os falsos positivos;
- sensibilidade à prior (uniforme, histórica e cética).

## 4. Word embeddings do zero (`estudo_quatro.R`)

Vetores de palavras construídos a partir dos seis romances de Jane Austen, sem rede neural:

- matriz de coocorrência com janela de 5 palavras, ponderada pela distância;
- PPMI com suavização de contexto e redução para 100 dimensões com SVD;
- vizinhos semânticos (contagem bruta contra embeddings) e analogias do tipo he : she :: him : ?;
- mapa semântico em 2D por PCA;
- medição do viés de gênero que o modelo aprende do corpus.

## 5. Atenção do zero (`estudo_cinco.R`)

O mecanismo de self-attention do Transformer implementado só com multiplicação de matrizes, treinado em duas tarefas escolhidas para tornar a atenção visível:

- embeddings de dígitos e de posição, projeções Q, K e V e softmax escalado por √d;
- retropropagação escrita à mão, verificada contra gradientes numéricos;
- treino com Adam em lotes gerados a cada passo;
- **inverter sequências:** a atenção aprende uma diagonal invertida, ou seja, atenção por posição;
- **encontrar o máximo:** todas as posições olham para o maior dígito, ou seja, atenção por conteúdo;
- ablação sem embeddings de posição: a atenção é invariante a permutações e não consegue inverter a sequência.

## 6. Autocorrelação espacial de homicídios (`estudo_seis.R`)

Homicídios por município na Bahia entre 2020 e 2022 (SIM/DATASUS, CID-10 X85 a Y09), com população e densidade do Censo 2022:

- taxas por 100 mil habitantes com suavização bayesiana empírica, para não confundir ruído de municípios pequenos com padrão espacial;
- vizinhança por contiguidade (queen) e pesos padronizados por linha;
- Moran global com teste de permutação (999 simulações);
- LISA para localizar hotspots (alto-alto), coldspots (baixo-baixo) e outliers espaciais;
- regressão OLS diagnosticada com testes de Rao (LM) e comparada com os modelos SAR (defasagem espacial) e SEM (erro espacial), incluindo os impactos diretos e indiretos do SAR.

Gera `homicidios_taxa.png` e `homicidios_lisa.png`.

## 7. Regiões de vulnerabilidade social: k-means contra SKATER (`estudo_sete.R`)

Seis indicadores do Atlas do Desenvolvimento Humano (Censo 2010) para os municípios da Bahia: analfabetismo, pobreza extrema, Gini, mortalidade infantil, água e esgoto inadequados e renda per capita:

- indicadores orientados na mesma direção e padronizados por z-score;
- escolha do número de grupos pelo cotovelo e pela silhueta;
- k-means, que agrupa só por semelhança, contra SKATER, que corta a árvore geradora mínima do grafo de vizinhança e produz regiões contíguas;
- tratamento de municípios-ilha (Itaparica e Vera Cruz), que deixariam o grafo desconexo;
- comparação do trade-off: homogeneidade interna e silhueta contra número de fragmentos espaciais.

Gera `vulnerabilidade_clusters.png`, com os dois mapas lado a lado.

## Como executar

Requer R 4.x. O primeiro estudo usa os dados do pacote `causaldata`, e o quarto, os do `janeaustenr`:

```r
install.packages(c("causaldata", "janeaustenr"))
```

Os estudos 6 e 7 instalam os próprios pacotes quando necessário (`sf`, `geobr`, `spdep`, `spatialreg`, `sidrar`, `abjData` e `microdatasus`, este último pelo GitHub) e precisam de internet para baixar malhas e dados do IBGE, DATASUS e Atlas Brasil.

```bash
Rscript estudo_um.R     # cerca de 1 minuto por causa do bootstrap
Rscript estudo_dois.R   # cerca de 1 minuto
Rscript estudo_tres.R   # menos de 1 minuto
Rscript estudo_quatro.R # depende do SVD de uma matriz 3000 x 3000, de 1 a alguns minutos
Rscript estudo_cinco.R  # três treinos de 2.000 passos, alguns minutos
Rscript estudo_seis.R   # o download do SIM pelo FTP do DATASUS domina o tempo
Rscript estudo_sete.R   # menos de 1 minuto depois de instalados os pacotes
```

Cada script gera o seu gráfico principal em PNG na pasta atual.
