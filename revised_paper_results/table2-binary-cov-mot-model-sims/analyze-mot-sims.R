## SKG
## July 31, 2020
## Analyze simulations
## JUNE 15, 2021 UPDATES FOR JCGS REVISION

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(knitr)
library(kableExtra)
library(viridis)


## Biowulf output from outsider-sims.R
folder_name <- "biowulf-results"

filenames <- list.files(file.path(folder_name))

out_list <- lapply(file.path(folder_name, filenames),
                      readRDS)[[1]]

## Parameter sets

## Simulations to try
M <- 100 # number of simulations for each parameter configuration
B <- 10000 # number of MC samples
K <- 1000 # number of clusters in a data set


par_df1 <- data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.5, -2.5, -2.5,
                         -1.5, -1.5, -1.5),
                     beta_1 = c(.5, 0, -.5,
                                1, 0, -1),
                     gamma = .5)

par_df2 <- data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.5, -2.5, -2.5,
                         -1.5, -1.5, -1.5),
                     beta_1 = c(.5, 0, -.5,
                                .5, 0, -.5),
                     gamma = .7)

par_df3 <-  data.frame(
                     M = M,
                     B = B,
                     beta_0 = c(-2.75, -2.75, -2.75),
                     beta_1 = c(2.25, 0, -.25),
                     gamma = .1)
###########################################################################

## data frame of simulation parameter configurations
par_df <- dplyr::bind_rows(par_df1, par_df2, par_df3)


sizes <- as.data.frame(do.call('rbind', lapply(out_list, "[[", 7)))
cov <- as.data.frame(do.call('rbind', lapply(out_list, "[[", 8))[, 1:2])
colnames(cov) <- c("fi_beta0", "fi_beta1")

mean_pars <- as.data.frame(do.call('rbind', lapply(
                                   lapply(out_list, "[[", 1),
                                   colMeans)))
mean_pars <- as.data.frame(do.call('rbind', lapply(
    lapply(out_list, "[[", 1),
    function(mat){
        apply(mat, 2, mean)
    }))) %>%
    rename(mean_beta0 = V1,
           mean_beta1 = V2)

med_pars <- as.data.frame(do.call('rbind', lapply(
    lapply(out_list, "[[", 1),
    function(mat){
        apply(mat, 2, mean)
    }))) %>%
    rename(med_beta0 = V1,
           med_beta1 = V2)

iqr_pars <- as.data.frame(do.call('rbind', lapply(
    lapply(out_list, "[[", 1),
    function(mat){
        apply(mat, 2, function(x){
            quantile(x, prob = .75) - quantile(x, prob = .25)
        })
    }))) %>%
    rename(IQR_beta0 = V1,
           IQR_beta1 = V2)

se_pars <- as.data.frame(do.call('rbind', lapply(
    lapply(out_list, "[[", 1),
    function(mat){
        apply(mat, 2, function(x){
            sd(x)
        })
    }))) %>%
    rename(se_beta0 = V1,
           se_beta1 = V2)

summary_df <- cbind(par_df,
                    mean_pars,
                    med_pars,
                    iqr_pars,
                    se_pars,
                    cov,
                    sizes) %>%
    arrange(max, q90)

pretty_df <- summary_df %>%
    mutate(sizes = paste0("(", `max`,
                          ", ", `q90`,
                          ", ", `med`,
                          ")"),
           prob_pos = gamma,
           betas = paste0("(", beta_0,
                          ", ", beta_1,
                          ")"),
           est_betas1 = paste0("(", round(mean_beta0, 2),
                               " (", round(se_beta0, 2),
                               "), ", round(mean_beta1, 2),
                               " (", round(se_beta1, 2), "))"),
           est_betas2 = paste0("(", round(med_beta0, 2),
                               " (", round(IQR_beta0, 2),
                               "), ", round(med_beta1, 2),
                               " (", round(IQR_beta1, 2), "))"),
           cov = paste0("[", fi_beta0 * 100, ", ",
                        fi_beta1 * 100, "]\\%")
           ) %>%
    select(sizes, prob_pos, betas,
           est_betas1,
           est_betas2,
           cov)


pretty_df %>%
    knitr::kable("latex", digits = 2, booktabs = TRUE, align = "c",
                 linesep = "",
                 col.names = c("$|\\mathcal{C}|$ (Max, 90Q, 50Q)",
                               "P(+)",
                               "$(\\beta_0$, $\\beta_1)$",
                               "(Mean $\\hat{\\beta}_0$, (SE), Mean $\\hat{\\beta}_1$ (SE))",
                               "(Median $\\hat{\\beta}_0$, (IQR), Median $\\hat{\\beta}_1$ (IQR))",
                               "Coverage [$\\beta_0$, $\\beta_1$]"
                               ), escape = FALSE,
                 caption = paste0("Results of simulations where we first generate data from Model in Eq. \\eqref{eq:sim-mod}",
                                  "and then exclude the seed individual to reflect the outsider generation process",
                 " For each row in the table, we generate 100 sets of outbreak data each with 1000 transmission trees for a given $P(+)$, $\\beta_0$ and $\\beta_1$. ",
                 "to estimate the MLE of $\\bm{\\beta}$.  However, fewer than 1000 transmission trees are observed once we exclude the first generation (the outside infection).",
                 "The mean and bootstrap standard error (SE) along with the median and inner quartile range (IQR) ",
                 "are reported over the 100 sets of outbreak data. ",
                 "Additionally we report the coverage where a 95\\% CI was derived ",
                 "using a Fisher information derived estimate of the SE."
                 ),
                 label = "sim-results-mot") %>%
    kableExtra::kable_styling(latex_options = c("scale_down",
                                                "striped"))
