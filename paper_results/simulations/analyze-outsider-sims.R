## SKG
## Aug 31, 2020
## Analyze outsider simulations from Biowulf
##



## Biowulf output from paper_results/simulations/outsider-sims.R


folder_name <- "single_output_outsider_all15sets"

filenames <- list.files(file.path(folder_name))

output_list <- lapply(file.path(folder_name, filenames),
                      readRDS)


df <- do.call('rbind', lapply(output_list, function(ll){
    betas_mean <- colMeans(ll$best_params_mat)
    betas_med <- apply(ll$best_params_mat, 2, median)
    lower <- apply(ll$best_params_mat, 2, quantile, prob = c(.025))
    upper <- apply(ll$best_params_mat, 2, quantile, prob = c(.975))
    data.frame(max = ll$cluster_sizes[2],
               med = ll$cluster_sizes[1],
               q90 = ll$cluster_sizes[3],
               p_pos = ll$gamma,
               beta0 = ll$par_set[1],
               beta1 = ll$par_set[2],
               beta0_mean = betas_mean[1],
               beta1_mean = betas_mean[2],
               beta0_med = betas_med[1],
               beta1_med = betas_med[2],
               beta0_lower =  lower[1],
               beta1_lower = upper[1],
               beta0_upper = lower[2],
               beta1_upper = upper[2]
               )}))
rownames(df) <- NULL

library(tidyverse)
library(knitr)
library(kableExtra)

df %>% arrange(max) %>%
  kable("latex", digits = 2, booktabs = TRUE, align = "l",
             linesep = "",
             col.names = c("Max $|C_m|$",
                           "Median $|C_m|$",
                           "90Q $|C_m|$",
                           "P(+)",
                           "$\\beta_0$",
                           "$\\beta_1$",
                           "Mean $\\hat{\\beta}_0$",
                           "Mean $\\hat{\\beta}_1$",
                           "Median $\\hat{\\beta}_0$",
                           "Median $\\hat{\\beta}_1$",
                           "2.5Q $\\hat{\\beta}_0$",
                           "97.5Q  $\\hat{\\beta}_0$",
                           "2.5Q $\\hat{\\beta}_1$",
                           "97.5Q $\\hat{\\beta}_1$"
             ), escape = FALSE,
             caption = "Results of simulations where we generate from the multiple outside transmissions model.  For each row in the table, we generate 100 sets of outbreak data each with 1000 transmission trees.  To estimate the MLE of $\\bm{\\beta}$, To obtain the confidence interval, we take the 2.5 and 97.5 quantiles over the 100 estimates of $\\bm{\\beta}$.",
             label = "sim-results-mot") %>%
  kableExtra::kable_styling(latex_options = c("scale_down",
                                              "striped"))
