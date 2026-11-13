# Estudo 8: preço de imóveis com tidymodels (Lasso contra random forest)

# 1. Pacotes
install.packages(c("tidymodels", "modeldata", "glmnet", "ranger", "vip"))

library(tidymodels)
library(vip)
tidymodels_prefer()
set.seed(123)

# 2. Dados (o alvo em log10 estabiliza a variância)
data(ames, package = "modeldata")
ames <- ames |> mutate(Sale_Price = log10(Sale_Price))

# 3. Divisão estratificada e validação cruzada (só no treino)
divisao <- initial_split(ames, prop = 0.8, strata = Sale_Price)
treino  <- training(divisao)
teste   <- testing(divisao)
folds   <- vfold_cv(treino, v = 10, strata = Sale_Price)

# 4. Receita de pré-processamento
receita <- recipe(Sale_Price ~ ., data = treino) |>
  step_other(all_nominal_predictors(), threshold = 0.01) |>
  step_dummy(all_nominal_predictors()) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# 5. Modelos
modelo_lasso <- linear_reg(penalty = tune(), mixture = 1) |>
  set_engine("glmnet")

modelo_floresta <- rand_forest(trees = 500, mtry = tune(), min_n = tune()) |>
  set_engine("ranger", importance = "permutation") |>
  set_mode("regression")

# 6. Conjunto de workflows e ajuste de hiperparâmetros
conjunto <- workflow_set(
  preproc = list(base = receita),
  models  = list(lasso = modelo_lasso, floresta = modelo_floresta)
)

resultados <- workflow_map(
  conjunto,
  fn        = "tune_grid",
  resamples = folds,
  grid      = 15,
  metrics   = metric_set(rmse, rsq),
  seed      = 1,
  verbose   = TRUE
)

autoplot(resultados)
rank_results(resultados, rank_metric = "rmse", select_best = TRUE)

# 7. Finalizar o melhor workflow e avaliar no teste (uma única vez)
id_melhor <- rank_results(resultados, rank_metric = "rmse", select_best = TRUE) |>
  slice(1) |>
  pull(wflow_id)

melhores_params <- resultados |>
  extract_workflow_set_result(id_melhor) |>
  select_best(metric = "rmse")

ajuste_final <- resultados |>
  extract_workflow(id_melhor) |>
  finalize_workflow(melhores_params) |>
  last_fit(divisao)

collect_metrics(ajuste_final)

# 8. Predito vs. observado
collect_predictions(ajuste_final) |>
  ggplot(aes(x = Sale_Price, y = .pred)) +
  geom_abline(lty = 2, color = "gray50") +
  geom_point(alpha = 0.4) +
  coord_obs_pred() +
  labs(x = "Observado (log10)", y = "Predito (log10)")

# 9. Importância das variáveis
ajuste_final |>
  extract_fit_parsnip() |>
  vip(num_features = 15)
