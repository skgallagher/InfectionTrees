## SKG
## June 15, 2021
## Comparing our like function to the cluster generated in
## Fig. 2

rep_from_lib <- TRUE

if(rep_from_lib){
    library(InfectionTrees)
} else {
    devtools::load_all()
}


## Calculate the true likelihood for BP in fig. 2
true_like_31 <- function(beta0, beta1){
    p_pos <- 1 / (1 + exp(-beta0 - beta1))
    p_neg <- 1 / (1 + exp(-beta0))
    (
        (1-p_pos) * (1-p_neg)^2 *
        (p_pos^2 + p_neg^2 +
         p_neg^2 + p_pos * p_neg +
         p_neg * p_pos + p_neg^2)
    ) / 6
}


## "observed data of one cluster - size 3 with 1 pos
obs_data_summary <- data.frame(freq = 1, cluster_size = 3,
                                        x_pos = 1, x_neg = 2)

B <- 5000
mc_trees <- sample_mc_binary_cov(B = B,
                                 observed_cluster_summaries =
                                     obs_data_summary) %>%
    dplyr::select(-freq)

mc_loglike <-function(x){
    sapply(x, function(x){
        bp_loglike_binary_cov(inf_params = c(-.406, x),
                              obs_data_summary = obs_data_summary,
                                    mc_samples_summary =
                                        mc_trees)
    })
}

df <- data.frame(x = seq(-2, 2, by = .1))
df$true_loglike <- log(true_like_31(-.406, df$x))
df$mc_loglike <- mc_loglike(df$x)


library(ggplot2)
ggplot(data = df %>%
           tidyr::pivot_longer(-x),
       aes(x = x, y = value, col = name)) +
    geom_point()


mc_loglike_beta0 <-function(x){
    sapply(x, function(x){
        bp_loglike_binary_cov(inf_params = c(x,0),
                              obs_data_summary = obs_data_summary,
                                    mc_samples_summary =
                                        mc_trees)
    })
}

df <- data.frame(x = seq(-2, 2, by = .1))
df$true_loglike <- log(true_like_31(df$x, 0))
df$mc_loglike <- mc_loglike_beta0(df$x)


library(ggplot2)
ggplot(data = df %>%
           tidyr::pivot_longer(-x),
       aes(x = x, y = value, col = name)) +
    geom_point()


## go for a heat map oooh
df <- expand.grid(beta0 = seq(-2, 2, by = .1),
                  beta1 = seq(-2, 2, by = .1))
df$true_loglike <- NA
df$mc_loglike <- NA

for(row in 1:nrow(df)){
    df$true_loglike[row] <- log(true_like_31(df$beta0[row], df$beta1[row]))
    df$mc_loglike[row] <-  bp_loglike_binary_cov(inf_params =
                                                     c(df$beta0[row],
                                                            df$beta1[row]),
                                                       obs_data_summary =
                                                           obs_data_summary,
                                                       mc_samples_summary =
                                                           mc_trees)
}
                  

ggplot(data = df,
       aes(x = beta0,
           y = beta1,
           fill = exp(true_loglike) - exp(mc_loglike))) +
geom_tile() +
    scale_fill_distiller(palette = "RdBu",
                         limits = c(-0.0004, .0004)) +
    labs(x = latex2exp::TeX("$\\beta_0$"),
         y = latex2exp::TeX("$\\beta_1$"),
         title = "Difference between True and Sampled Likelihood",
         subtitle = "For a cluster of size 3, with 1 smear positive",
         fill = "Difference") +
    geom_vline(xintercept = -.406, linetype = "dashed") +
    geom_hline(yintercept = -0.00005,
               linetype = "dashed") +
    theme_bw(base_size = 14) +
    geom_contour(aes(z = true_loglike),
                 bins = 30)
       
           
