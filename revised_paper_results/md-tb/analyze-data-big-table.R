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
fn_rds2_o <- grep("mot", fn_rds_o, value = TRUE)

data_list_o <- readRDS(fn_rds2_o)


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
## add boot strap SE
base_boot <- readRDS("biowulf-results/bootstrap_sims_base_2021-06-23.RDS")
best_mod_base$se_boot <- base_boot$se_vec

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
## add boot strap SE
mot_boot <- readRDS("biowulf-results/bootstrap_sims_mot_2021-06-24.RDS")
best_mod_o$se_boot <- mot_boot$se_vec


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
## Se boot
naive_boot <- readRDS("naive-tab-boot.RDS")
best_mod_naive$se_boot <- naive_boot$se_boot

########
## Putting it all together
combined_best_mod <- rbind(best_mod_base,
                           best_mod_o,
                           best_mod_naive) %>%
    mutate(or = exp(`Est.`),
           lower = exp(`Est.`- 1.96 * se_boot),
           upper = exp(`Est.` + 1.96 * se_boot)) %>%
    mutate(pretty_or = sprintf("%.2f [%.2f, %.2f]", or, lower, upper)) %>%
    mutate(combined_se = sprintf("(%.2f, %.2f)", SE, se_boot))









combined_best_mod %>%
    select(c("Type", "Coeff", "Est.", "combined_se", "pretty_or")) %>%
    kable(format = "latex", booktabs = TRUE,
          col.names = c("Type", "Coeff.",
                        "$\\hat{\\beta}$",  "$\\widehat{SE}_F(\\hat{\\beta}), \\widehat{SE}_B(\\hat{\\beta})$",
                        "Odds Ratio [95\\% CI]"),
          caption = paste("Coefficient estimates for the best selected model.",
                          "We report the mean estimate and two estimates of the standard error of the $\\beta$ values",
                          "in addition to the odds ratio and 95\\% CI.  $SE_F$ represents standard error",
                          "derived from the Fisher information, and $SE_B$ represents the bootstrap standard error.",
                          "The 95\\% CI is estimated using the bootstrap SE."),
          label = "best-model-ests",
          digits = 2,
          linesep = "",
          escape = FALSE,
          align = c("l", "l", "c", "c", "c")) %>%
    kable_styling(latex_options = "striped", position = "center",
                  stripe_index = c(1, 7,  11)) %>%
    row_spec(c(3:5, 9,10, 13:15), background = "yellow") %>%
    row_spec(c(5, 10), extra_latex_after = "\\midrule")
