# Estudo 10: interpretabilidade com SHAP e modelo empacotado com vetiver

# 1. Pacotes
install.packages(c("tidymodels", "modeldata", "xgboost", "shapviz", "vetiver", "plumber"))

library(tidymodels)
library(shapviz)
library(vetiver)
library(plumber)
tidymodels_prefer()
set.seed(42)

# 2. Dados (níveis de Attrition: "No", "Yes")
data(attrition, package = "modeldata")

divisao <- initial_split(attrition, prop = 0.8, strata = Attrition)
treino  <- training(divisao)
teste   <- testing(divisao)
folds   <- vfold_cv(treino, v = 5, strata = Attrition)

# 3. Receita: ordinais viram números e nominais viram one-hot (SHAP por categoria)
receita <- recipe(Attrition ~ ., data = treino) |>
  step_ordinalscore(all_ordered_predictors()) |>
  step_dummy(all_nominal_predictors(), one_hot = TRUE) |>
  step_zv(all_predictors())

# 4. Modelo e ajuste de hiperparâmetros
modelo <- boost_tree(
  trees = tune(), tree_depth = tune(), learn_rate = tune(), min_n = tune()
) |>
  set_engine("xgboost") |>
  set_mode("classification")

fluxo <- workflow() |>
  add_recipe(receita) |>
  add_model(modelo)

ajuste_grade <- tune_grid(
  fluxo,
  resamples = folds,
  grid      = 15,
  metrics   = metric_set(roc_auc)
)

show_best(ajuste_grade, metric = "roc_auc", n = 5)

# 5. Modelo final e avaliação no teste
fluxo_final <- fluxo |>
  finalize_workflow(select_best(ajuste_grade, metric = "roc_auc"))

ajuste_final <- last_fit(fluxo_final, divisao)
collect_metrics(ajuste_final)

# 6. SHAP: extrai o booster nativo e a matriz de preditores já transformada
fluxo_ajustado <- extract_workflow(ajuste_final)
booster        <- extract_fit_engine(fluxo_ajustado)
receita_prep   <- extract_recipe(fluxo_ajustado)

X_teste <- bake(receita_prep, new_data = teste, all_predictors(), composition = "matrix")

shap <- shapviz(booster, X_pred = X_teste, X = X_teste)

# Importância global (média do |SHAP|) e distribuição por variável
sv_importance(shap, kind = "bar", max_display = 15)
grafico_beeswarm <- sv_importance(shap, kind = "beeswarm", max_display = 15)
ggsave("shap_importancia.png", grafico_beeswarm, width = 7, height = 5, dpi = 150)

# Explicação de um caso individual
grafico_waterfall <- sv_waterfall(shap, row_id = 1)
ggsave("shap_waterfall.png", grafico_waterfall, width = 7, height = 5, dpi = 150)

# Efeito de uma variável numérica
sv_dependence(shap, v = "MonthlyIncome", color_var = "auto")

# 7. Empacotamento com vetiver (sem gravar em disco)
v <- vetiver_model(fluxo_ajustado, model_name = "modelo_attrition")
v

# Predição local, simulando o que a API faria
predict(bundle::unbundle(v$model), new_data = head(teste, 3), type = "prob")

# Monta a API sem iniciá-la
api <- pr() |> vetiver_api(v)
api
