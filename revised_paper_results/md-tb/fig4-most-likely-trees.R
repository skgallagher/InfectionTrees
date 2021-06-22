## SKG
## June 22, 2021
## JCGS revisiosn to reproduce code for most likely trees


#################################################
## A couple helper functions

#' Assign x-coordinates based on generation and index for a single cluster
#'
#' @param gen generation number
#' @return generation number
assign_coords_x <- function(gen){
    return(gen)
}


#' Assign x-coordinates based on generation and index for a single cluster
#'
#' @param index index ONLY for those in a single generation
#' @param ymin default is -1
#' @param ymax default is 1
#' @return index space, evenly spaced by number in generation
assign_coords_y <- function(index,
                            ymin = -1, ymax = 1){
    gen_ranks <- rank(index)
    n_in_gen <- length(gen_ranks)
    if(n_in_gen == 1){
        return(0)
    }
    y <- -1 + (gen_ranks - 1) /(.5 * (length(gen_ranks)-1))
    return(y)

}

###########################################################


devtools::load_all()
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)


## Subset data to cluster 27, it's cool
tb_ex <- tb_clean %>%
    filter(group == 27) %>%
    mutate(rel_time = as.numeric(rel_time) / 365) %>%
    arrange(rel_time) %>%
    mutate(order = 1:n()) %>%
    dplyr::mutate(smear = ifelse(spsmear == "Positive",
                                 1, 0),
                  cluster_id = group,
                  hiv_f = ifelse(hivstatus == "Negative", "neg",
                          ifelse(hivstatus == "Positive", "pos",
                                 "unk"))) %>%
    dplyr::mutate(hiv_neg_pos = ifelse(hiv_f == "neg", 1, 0),
                  hiv_unk_pos = ifelse(hiv_f == "unk", 1, 0)) %>%
    ungroup()

## Load in best results from base model
fns <- list.files()
fns_rds <- grep(".RDS", fns, value = TRUE)
base_model_results <- readRDS(grep("base", fns_rds, value = TRUE))
inf_params <- base_model_results$beta_list[[4]][,1]

## Sample a bunch of MC trees
B <- 100000
set.seed(622021)
covariate_names <- c("smear",
                         "hiv_neg_pos", "hiv_unk_pos",
                         "rel_time")
mc_trees <- sample_mc_trees(tb_ex,
                            B = B,
                            multiple_outside_transmissions = FALSE,
                            covariate_names = c(covariate_names,
                                                "order"))
## get the covariate matrix
cov_mat <- covariate_df_to_mat(mc_trees,
                                cov_names = covariate_names)
mc_trees$prob_inf <- 1 / (1 + exp(-(cov_mat %*% inf_params)))
cluster_id <- n_inf <- orig_id <- prob_inf <- NULL
my_prob_inf <-  1 / (1 + exp(-(cov_mat %*% inf_params)))

## Turning mc_trees to a data table
mc_trees.dt <- data.table::as.data.table(mc_trees)
mc_trees.dt <- mc_trees.dt[, prob_inf := my_prob_inf]
cluster_id <- NULL

## Getting the likelihood for each sampled cluster
like_df <- mc_trees.dt[,
                    .(like = general_cluster_like.dt(prob_inf, n_inf)),
                    by = .(orig_id, cluster_id)]

mc_trees_like <- left_join(mc_trees, like_df,
                           by = c(orig_id, cluster_id)) %>%
    mutate(prob = like / sum(like))



df <- mc_trees_like %>% group_by(cluster_id) %>%
    mutate(x = assign_coords_x(gen)) %>%
    group_by(cluster_id, gen) %>%
    mutate(y = assign_coords_y(n_in_gen)) %>%
    ungroup() %>%
    group_by(cluster_id) %>%
    mutate(gen_size = max(gen)) %>%
    arrange(desc(prob))

top_groups_by_gen <- df %>%
    ungroup() %>%
    group_by(gen_size, cluster_id) %>%
    summarize(prob = min(prob),
              .groups = "drop_last") %>%
    slice_max(order_by = prob, n = 3, with_ties = FALSE)

top_clusts <- top_groups_by_gen$cluster_id


df_sub <- df %>% filter(cluster_id %in% top_clusts)
df_sub <-  df_sub %>%
    mutate(facet = paste(gen_size, cluster_id, sep = "-"))


inf_df <- df_sub %>%
    select(id, cluster_id, x, y) %>%
    rename(x_to = x, y_to = y,
           inf_id = id)

jsub <- left_join(df_sub, inf_df, by = c("cluster_id",  "inf_id"))
fctr_sub <- df_sub %>% group_by(cluster_id) %>%
    summarize(like = prob[1])
jsub$factor_id <- factor(jsub$cluster_id,
                         labels = paste0("Sample ",
                                         fctr_sub$cluster_id,
                                         ": P(C) = ",
                                         formatC(fctr_sub$like,
                                                 format = "e",
                                                 digits = 2)))
jsub$factor_gen <- factor(jsub$gen_size,
                          levels = 1:max(jsub$gen_size),
                          labels = paste0("Gen. ", 1:max(jsub$gen_size)))



fctr_sub2 <- df_sub %>% group_by(facet, gen_size, cluster_id) %>%
    summarize(like = prob[1]) %>% 
              mutate(id = stringr::str_sub(cluster_id, start = -3))
jsub$facet_lab <- factor(jsub$facet,
                         levels = unique(fctr_sub2$facet),
                         labels = paste0("# Gens: ",
                                         fctr_sub2$gen_size,
                                         "; ID: ",
                                         fctr_sub2$id,
                                         "; P(C) = ",
                                         formatC(fctr_sub2$like,
                                                 format = "e",
                                                 digits = 2)))

my_orig_id <- 27
                        
ggplot(data = jsub,
       aes(x = x, y = y,
           group = cluster_id,
           col = factor(order),
           ##  shape = factor(hiv_neg_pos)
           )) +
    geom_curve(aes(xend = x_to, yend = y_to),
               size = 2, curvature = -0,
               col = "black") +
    geom_point(size = 6, stroke = 4, col = "darkgray") +
    geom_point(size = 5, stroke = 4) +
    facet_wrap(~facet_lab, ncol = 3) +
    scale_color_brewer(palette = "Set1", name = "Detection Order")  +
    ##  scale_shape_manual(values = c(3, 16),
    ##                 labels = c("Pos./Unk.", "Neg."),
    ##                name = "HIV Status") +
    xlim(.8, 7.2) +
    ylim(-1.2, 1.2) +
    theme_bw(base_size = 18) +
    theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Generation in transmission tree",
         y = latex2exp::TeX("Order of infection in gen. among 'siblings' $\\rightarrow$"),
         title = "Most likely sampled trees by number of generations",
         subtitle = paste0("Cluster ", my_orig_id)) +
    theme(legend.position = "bottom")
   
ggsave("trees-7.pdf", height = 15, width = 13)
