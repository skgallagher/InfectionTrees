## SKG
## August 14, 2020
## Combining the tables together from three model types

devtools::load_all()
library(knitr)
library(tidyverse)
library(kableExtra)



## loglike table BASE MODEL
files <- list.files()
fn_rds <- grep("RDS", files, value = TRUE)
fn <- grep("base", fn_rds, value = TRUE)
data_list_base <- readRDS(fn)

loglike_df_base <- data_list_base$loglike_df

loglike_df_base$vars <- c("Null",
                          "Smear",
                     "Smear/HIV",
                     "Smear/HIV/rel. time",
                     "Smear/HIV/rel. time/race")
loglike_df_base$type <- "Base"

############################################
## OUTSIDER MODEL

files_o <- list.files()
fn_rds_o <- grep("RDS", files_o, value = TRUE)
fn_rds2_o <- fn_rds_o[grep("outsider", fn_rds_o)]
fn_o <- fn_rds2_o[5] ## 8-11-22 is date this should be/  I used last

data_list_o <- readRDS(fn_o)


loglike_df_o <- data_list_o$loglike_df

loglike_df_o$vars <- c("Null",
                       "Smear",
                       "Smear/HIV",
                       "Smear/HIV/rel. time",
                       "Smear/HIV/rel. time/race")
loglike_df_o$type <- "Multiple outside"


combined_loglike_df <- dplyr::bind_rows(loglike_df_base,
                                        loglike_df_o) %>%
    select(type, everything())

combined_loglike_df %>%
    select(type, Model, vars, n_params, loglike, aic) %>%
    kable(format = "latex", digits = 2,
          caption = "Reported log likelihood for the Maryland TB data for the five different models described in Eq. \\eqref{eq:data-models} for both the base and multiple outside transmissions model types.",
          label = "data-models",
          col.names = c("Type", "Model #", "Covariates", "# parameters", "Log like", "AIC"),
          booktabs = TRUE,
          linesep = "") %>%
    kable_styling(latex_options = "striped", position = "center") %>%
    row_spec(c(4,9), background = "yellow") %>%
    row_spec(c(5), extra_latex_after = "\\midrule") 


###############################################################
## BEST MODELS
## smear/HIV/rel time


## BASE
best_mod_base2 <- (data_list_base$output_list[[4]])$est_pars
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_base2) <- NULL
best_mod_base <- data.frame(Type = "Base",
                            Coeff = coeff_names,
                            best_mod_base2)

## OUTSIDER
best_mod_o2 <- (data_list_o$output_list[[4]])$est_pars
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_o2) <- NULL
best_mod_o <- data.frame(Type = "Multiple outside",
                            Coeff = coeff_names,
                         best_mod_o2)


##  NAIVE MODEL
best_mod_naive2 <- data.frame(readRDS("naive-tab.RDS"))
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_naive2) <- NULL
best_mod_naive <- data.frame(Type = "Naive",
                            Coeff = coeff_names,
                         best_mod_naive2)

########
## Putting it all together
combined_best_mod <- rbind(best_mod_base,
                           best_mod_o,
                           best_mod_naive)








combined_best_mod %>%
    kable(format = "latex", booktabs = TRUE,
          col.names = c("Type", "Coeff.",
                        "Mean Est.", "Lower", "Upper", "SE"),
          caption = "Coefficient estimates for the best selected model.  The columns `Lower' and `Upper' comprise the lower and upper bound estimates, respectively,  of the 95\\% CI, which come from likelihood profiling.  The highlighted rows are significant coefficients at the $\\alpha = .05$-level.  We also report the standard error of the estimate (SE).",
          label = "best-model-ests",
          digits = 2,
          linesep = "") %>%
    kable_styling(latex_options = "striped", position = "center") %>%
    row_spec(c(3, 5, 8, 9, 10, 13), background = "yellow") %>%
    row_spec(c(5, 10), extra_latex_after = "\\midrule") 
