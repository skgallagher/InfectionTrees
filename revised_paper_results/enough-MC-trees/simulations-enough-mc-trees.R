## ----setup, include=FALSE--------------------------------------------------------------------------------------

devtools::load_all()
library(InfectionTrees)
library(tidyr)
library(dplyr)
library(kableExtra)
library(ggplot2)
theme_set(theme_bw() + theme(axis.title = element_text()))


## --------------------------------------------------------------------------------------------------------------
set.seed(2020)
## Generate possible covariates
sample_covariates_df <- data.frame(x = c(1,0))

beta0 <- -2
beta1 <- 1.5
M <- 1000
data <- simulate_bp(K = M,
                                  inf_params = c(beta0, beta1),
                                  sample_covariates_df = sample_covariates_df,
                                  covariate_names = "x",
                                  max_size = 50)

 data %>% group_by(cluster_id) %>%
  summarize(cluster_size = n()) %>%
  ungroup() %>%
   group_by(cluster_size) %>%
   summarize(freq = n()) %>%
   ungroup() %>%
  kable(type = "html") %>%
   kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F,
                 position = "center")




## --------------------------------------------------------------------------------------------------------------
## Summarize data set into a table
binary_cluster_summaries <- summarize_binary_clusters(data)

B <- 5
K <- 5000

beta_mat <- matrix(0, nrow = B, ncol = 2)
var_mat <- matrix(0, nrow = B, ncol = 2)
seed <- 9152020

for(bb in 1:B){

    ## Get different seeds for MC sets
  seed <- seed + bb * 10
  set.seed(seed)

  ## Sample MC trees
  mc_trees <- sample_mc_binary_cov(K,
                                 observed_cluster_summaries = binary_cluster_summaries)

  ## Fit model
  inf_params_0 <- c(0,0)
  best_params <- optim(inf_params_0, fn = bp_loglike_binary_cov,
                       obs_data_summary = binary_cluster_summaries,
                       mc_samples_summary = mc_trees,
                       return_neg = TRUE,
                       lower = c(-5, -5),
                       upper = c(0, 5),
                       method = "L-BFGS-B",
                       hessian = TRUE)



  beta_mat[bb,] <- best_params$par
  var_mat[bb,] <- diag(solve(best_params$hessian))


}



## --------------------------------------------------------------------------------------------------------------
se_mc <- apply(beta_mat, 2, sd)
wald_scores_fisher <- beta_mat / sqrt(var_mat)
wald_se_fisher <- apply(wald_scores_fisher, 2, sd)




wald_se_mat <- rbind(wald_se_fisher, se_emp)
rownames(wald_se_mat) <- c("Fisher", "Emp.")

kable(wald_se_mat, digits = 2) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))


