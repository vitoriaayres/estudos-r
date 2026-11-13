# Estudo 9: detecção de fraude com classe rara (SMOTE, PR-AUC e ajuste de limiar)

# 1. Pacotes
install.packages(c("tidymodels", "themis", "xgboost", "probably"))

library(tidymodels)
library(themis)
library(probably)
tidymodels_prefer()
set.seed(2024)

# 2. Simulação de transações com ~2% de fraude
n <- 30000
dados <- tibble(
  valor                  = rlnorm(n, meanlog = 4, sdlog = 1),
  hora                   = sample(0:23, n, replace = TRUE),
  distancia_casa_km      = rexp(n, rate = 1 / 10),
  transacoes_ultima_hora = rpois(n, lambda = 1),
  idade_conta_dias       = rexp(n, rate = 1 / 500),
  exterior               = rbinom(n, 1, 0.05)
)

risco <- with(
  dados,
  0.5 * log(valor) + 0.02 * distancia_casa_km + 0.6 * transacoes_ultima_hora -
    0.002 * idade_conta_dias + 1.2 * (hora < 5) +
    1.0 * exterior * (log(valor) > 5)
)

# Ajusta o intercepto para a prevalência-alvo de 2%
intercepto <- uniroot(
  function(b) mean(plogis(b + risco)) - 0.02,
  interval = c(-30, 5)
)$root

dados$classe <- factor(
  ifelse(rbinom(n, 1, plogis(intercepto + risco)) == 1, "fraude", "legitima"),
  levels = c("fraude", "legitima")   # o primeiro nível é o evento de interesse
)

dados |> count(classe) |> mutate(prop = n / sum(n))

# 3. Divisão e validação cruzada
divisao <- initial_split(dados, prop = 0.8, strata = classe)
treino  <- training(divisao)
teste   <- testing(divisao)
folds   <- vfold_cv(treino, v = 5, strata = classe)

# 4. Duas receitas: sem e com SMOTE (SMOTE só dentro dos folds de treino)
receita_base <- recipe(classe ~ ., data = treino) |>
  step_log(valor) |>
  step_normalize(all_numeric_predictors())

receita_smote <- receita_base |>
  step_smote(classe, over_ratio = 0.5)

# 5. Modelos
modelo_logistico <- logistic_reg() |>
  set_engine("glm")

modelo_xgb <- boost_tree(
  trees = tune(), tree_depth = tune(), learn_rate = tune()
) |>
  set_engine("xgboost") |>
  set_mode("classification")

# 6. Comparação: 2 receitas x 2 modelos
conjunto <- workflow_set(
  preproc = list(base = receita_base, smote = receita_smote),
  models  = list(logistico = modelo_logistico, xgb = modelo_xgb)
)

resultados <- workflow_map(
  conjunto,
  fn        = "tune_grid",
  resamples = folds,
  grid      = 10,
  metrics   = metric_set(pr_auc, roc_auc),
  seed      = 1,
  verbose   = TRUE
)

autoplot(resultados, metric = "pr_auc")
rank_results(resultados, rank_metric = "pr_auc", select_best = TRUE)

# 7. Finalizar o melhor e avaliar no teste
id_melhor <- rank_results(resultados, rank_metric = "pr_auc", select_best = TRUE) |>
  slice(1) |>
  pull(wflow_id)

melhores_params <- resultados |>
  extract_workflow_set_result(id_melhor) |>
  select_best(metric = "pr_auc")

ajuste_final <- resultados |>
  extract_workflow(id_melhor) |>
  finalize_workflow(melhores_params) |>
  last_fit(divisao, metrics = metric_set(pr_auc, roc_auc))

collect_metrics(ajuste_final)

previsoes <- collect_predictions(ajuste_final)

grafico_pr <- previsoes |>
  pr_curve(truth = classe, .pred_fraude) |>
  autoplot()

ggsave("fraude_curva_pr.png", grafico_pr, width = 6, height = 5, dpi = 150)

# 8. Ajuste do limiar de decisão (o padrão de 0.5 raramente é o ideal)
limiares <- previsoes |>
  threshold_perf(
    truth      = classe,
    estimate   = .pred_fraude,
    thresholds = seq(0.05, 0.95, by = 0.05),
    metrics    = metric_set(precision, recall, f_meas)
  )

grafico_limiares <- limiares |>
  ggplot(aes(x = .threshold, y = .estimate, color = .metric)) +
  geom_line() +
  labs(x = "Limiar", y = "Valor", color = "Métrica")

ggsave("fraude_limiares.png", grafico_limiares, width = 6, height = 4, dpi = 150)

# 9. Matriz de confusão no limiar com melhor F1
melhor_limiar <- limiares |>
  filter(.metric == "f_meas") |>
  slice_max(.estimate, n = 1, with_ties = FALSE) |>
  pull(.threshold)

previsoes |>
  mutate(
    predita = factor(
      ifelse(.pred_fraude >= melhor_limiar, "fraude", "legitima"),
      levels = c("fraude", "legitima")
    )
  ) |>
  conf_mat(truth = classe, estimate = predita)
