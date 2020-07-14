## SKG
## June 18, 2020
## Analyzing the loglike and best model

devtools::load_all()
library(knitr)
library(tidyverse)
library(kableExtra)

files <- list.files()
fn_rds <- grep("RDS", files, value = TRUE)
fn <- fn_rds[4]

data_list <- readRDS(fn)

## loglike table

loglike_df <- data_list$loglike_df

loglike_df$vars <- c("Smear",
                     "Smear/HIV",
                     "Smear/HIV/rel. time",
                     "Smear/HIV/rel. time/race")

loglike_df %>% select(Model, vars, n_params, loglike, aic) %>%
    kable(format = "latex", digits = 2,
          caption = "Reported log likelihood for the Maryland TB data for the four different models described in Eq. \\eqref{eq:data-models}.",
          label = "data-models",
          col.names = c("Model", "Covariates", "# parameters", "Log like", "AIC"),
          booktabs = TRUE) %>%
    kable_styling(latex_options = "striped", position = "center") %>%
    row_spec(3, background = "yellow")

## Likelihood ratio tests
## 1 vs 3
my_chisq <- qchisq(.95, df = 3)
lrt_stat <- 2 * (loglike_df$loglike[3] -
                 loglike_df$loglike[1])
c(lrt_stat, my_chisq)


## Results for best model
## which is model 3

best_mod <- data_list$output_list[[5]]
est_pars <- best_mod$est_pars
rownames(est_pars) <- c("Intercept", "Smear",
                        "HIV neg.: HIV pos.",
                        "HIV unknown: HIV pos.",
                        "Relative Time")
est_pars %>%
    kable(format = "latex", booktabs = TRUE,
          col.names = c("Mean Est.", "Lower", "Upper"),
          row.names = TRUE,
          caption = "Coefficient estimates for the best selected model.  The columns `Lower' and `Upper' comprise the lower and upper bound estimates, respectively,  of the 95\\% CI.  The highlighted rows are significant coefficients at the $\\alpha = .05$-level.",
          label = "best-model-ests",
          digits = 2) %>%
    kable_styling(latex_options = "striped", position = "center") %>%
    row_spec(c(3,5), background = "yellow")


## odds
odds <- as.data.frame(est_pars) %>%
    mutate(odds_mean = exp(mean),
           odds_lower = exp(lower),
           odds_upper = exp(upper))
           
