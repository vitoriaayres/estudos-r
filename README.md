# Estudos em R

Dois estudos escritos em R base, sem bibliotecas de modelagem prontas.

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

## Como executar

Requer R 4.x. O primeiro estudo usa os dados do pacote `causaldata`:

```r
install.packages("causaldata")
```

```bash
Rscript estudo_um.R    # cerca de 1 minuto por causa do bootstrap
Rscript estudo_dois.R  # cerca de 1 minuto
```

Cada script gera o seu gráfico principal em PNG na pasta atual.
