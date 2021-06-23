## SKG
## August 14, 2020
## Combining the tables together from three model types
## Updated June 22, 2021

devtools::load_all()
library(knitr)
library(tidyverse)
library(kableExtra)



## loglike table BASE MODEL
files <- list.files()
fn_rds <- grep("RDS", files, value = TRUE)
fn <- fn_rds[1]

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
fn_rds2_o <- fn_rds_o[grep("mot", fn_rds_o)]
fn_o <- fn_rds2_o[1] 

data_list_o <- readRDS(fn_o)


loglike_df_o <- data_list_o$loglike_df

loglike_df_o$vars <- c("Null",
                       "Smear",
                       "Smear/HIV",
                       "Smear/HIV/rel. time",
                       "Smear/HIV/rel. time/race")
loglike_df_o$type <- "Multiple outside"
#########################################################################
## TABLE 5 - loglike for models
#####################################################
combined_loglike_df <- dplyr::bind_rows(loglike_df_base,
                                        loglike_df_o) %>%
    select(type, everything())

combined_loglike_df %>%
    select(type, model, vars, n_params, loglike, aic) %>%
    kable(format = "latex", digits = 2,
          caption = "Reported log likelihood for the Maryland TB data for the five different models described in Eq. \\eqref{eq:data-models} for both the base and multiple outside transmissions model types.",
          label = "data-models",
          col.names = c("Type", "Model #", "Covariates", "# parameters", "Log like", "AIC"),
          booktabs = TRUE,
          linesep = "") %>%
    kable_styling(latex_options = "striped", position = "center",
                  stripe_index = c(1, 3, 5, 7)) %>%
    row_spec(c(4), background = "yellow") %>%
    row_spec(c(9), background = "yellow") %>%
    row_spec(c(5), extra_latex_after = "\\midrule") 


#########################################################
###############################################################
## TABLE 6 - best model coefficients
## BEST MODELS
## smear/HIV/rel time


## BASE
best_mod_base2 <- data_list_base$beta_list[[4]][,c(1, 4)]
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_base2) <- NULL
best_mod_base <- data.frame(Type = "Base",
                            Coeff = coeff_names,
                            best_mod_base2)

## OUTSIDER
best_mod_o2 <- data_list_o$beta_list[[4]][,c(1,4)]
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_o2) <- NULL
best_mod_o <- data.frame(Type = "Multiple outside",
                            Coeff = coeff_names,
                         best_mod_o2)


## NAIVE MODEL
best_mod_naive2 <- data.frame(readRDS("naive-tab.RDS"))
coeff_names <- c("Intercept", "Smear pos.",
                 "HIV neg.: HIV pos.",
                 "HIV unk.: HIV pos.",
                 "Relative time")
rownames(best_mod_naive2) <- NULL
best_mod_naive <- data.frame(Type = "Naive",
                            Coeff = coeff_names,
                            best_mod_naive2[, c(1, 4)]) %>%
    rename("Est." = Mean)

########
## Putting it all together
combined_best_mod <- rbind(best_mod_base,
                           best_mod_o,
                           best_mod_naive) %>%
    mutate(or = exp(`Est.`),
           lower = exp(`Est.`- 1.96 * SE),
           upper = exp(`Est.` + 1.96 * SE)) %>%
    mutate(pretty_or = paste0(round(or, 2), " [",
                              round(lower,2), ", ",
                              round(upper, 2), "]"))
           








combined_best_mod %>%
    select(-c(lower, upper, or)) %>%
    kable(format = "latex", booktabs = TRUE,
          col.names = c("Type", "Coeff.",
                        "$\\hat{\\beta}$",  "$\\widehat{SE}(\\hat{\\beta})$",
                        "Odds Ratio and 95\\% CI"),
          caption = "Coefficient estimates for the best selected model.  We report the mean estimate and standard error of the $\\beta$ values in addition to the odds ratio and 95\\% CI.",
          label = "best-model-ests",
          digits = 2,
          linesep = "",
          escape = FALSE,
          align = c("l", "l", "c", "c", "c")) %>%
    kable_styling(latex_options = "striped", position = "center",
                  stripe_index = c(1,3, 7, 9, 11, 13)) %>%
    row_spec(c(5, 10, 15), background = "yellow") %>%
    row_spec(c(5, 10), extra_latex_after = "\\midrule") 
