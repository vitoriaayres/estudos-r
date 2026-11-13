# Estudo 11: previsão de demanda com séries temporais (naive, ARIMA, Prophet e XGBoost)

# 1. Pacotes
install.packages(c("tidymodels", "tidyverse", "modeltime", "timetk", "prophet", "xgboost"))

library(tidymodels)
library(tidyverse)
library(modeltime)
library(timetk)
tidymodels_prefer()
set.seed(7)

# 2. Dados: aluguéis diários
dados <- bike_sharing_daily |>
  select(data = dteday, alugueis = cnt)

dados |> plot_time_series(data, alugueis, .interactive = FALSE)

# 3. Separação temporal: últimos 90 dias para teste
divisao <- time_series_split(dados, assess = "90 days", cumulative = TRUE)

divisao |>
  tk_time_series_cv_plan() |>
  plot_time_series_cv_plan(data, alugueis, .interactive = FALSE)

treino <- training(divisao)
teste  <- testing(divisao)

# 4. Modelos
mod_naive <- naive_reg() |>
  set_engine("naive") |>
  fit(alugueis ~ data, data = treino)

mod_arima <- arima_reg() |>
  set_engine("auto_arima") |>
  fit(alugueis ~ data, data = treino)

mod_prophet <- prophet_reg(seasonality_yearly = TRUE) |>
  set_engine("prophet") |>
  fit(alugueis ~ data, data = treino)

# XGBoost com features de calendário
receita_ml <- recipe(alugueis ~ data, data = treino) |>
  step_timeseries_signature(data) |>
  step_rm(matches("(.iso$)|(.xts$)")) |>
  step_zv(all_predictors()) |>
  step_normalize(matches("(index.num$)|(_year$)")) |>
  step_dummy(all_nominal(), one_hot = TRUE) |>
  step_fourier(data, period = 365, K = 5)

mod_xgb <- workflow() |>
  add_model(
    boost_tree(mode = "regression", trees = 300, learn_rate = 0.05) |>
      set_engine("xgboost")
  ) |>
  add_recipe(receita_ml |> update_role(data, new_role = "ID")) |>
  fit(treino)

# 5. Tabela de modelos e calibração no teste
tabela <- modeltime_table(mod_naive, mod_arima, mod_prophet, mod_xgb)

calibrado <- tabela |> modeltime_calibrate(new_data = teste)

calibrado |> modeltime_accuracy()

calibrado |>
  modeltime_forecast(new_data = teste, actual_data = dados) |>
  plot_modeltime_forecast(.interactive = FALSE)

# 6. Validação por janela deslizante (5 janelas de 60 dias)
janelas <- time_series_cv(
  dados,
  initial     = "1 year",
  assess      = "60 days",
  skip        = "60 days",
  cumulative  = FALSE,
  slice_limit = 5
)

janelas |>
  plot_time_series_cv_plan(data, alugueis, .interactive = FALSE)

resultados_janelas <- tabela |>
  modeltime_fit_resamples(resamples = janelas, control = control_resamples(verbose = TRUE))

resultados_janelas |> modeltime_resample_accuracy(summary_fns = list(mean = mean, sd = sd))

# 7. Reajuste com todos os dados e previsão dos próximos 60 dias
calibrado |>
  modeltime_refit(data = dados) |>
  modeltime_forecast(h = "60 days", actual_data = dados) |>
  plot_modeltime_forecast(.interactive = FALSE)
