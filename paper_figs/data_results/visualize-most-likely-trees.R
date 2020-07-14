## SKG
## Basic coding to visualize clusters
## June 19, 2020

library(tidyverse)
library(data.table)
library(stringr)
devtools::load_all()


files <- list.files()
fn_rds <- grep('*outsider.*RDS', files,
               value = TRUE)
fn <- max(fn_rds)

## Let's look at a cluster of size x
data_list <- readRDS(fn)

sampled_data <- data_list$sampled_data
covariate_names <- c("smear", "hiv_neg_pos",
                     "hiv_unk_pos", "rel_time")
best_model_number <- 3 

my_orig_id <- "27" #121 is for 5 # 27 is for a cool one with 8
## 27 is 8
## 12 is for 4 people
my_orig_id <- "27"

#dog <- sampled_data %>% filter(cluster_size == 4)

sub_data <- sampled_data %>%
    filter(orig_id == my_orig_id)



joined_data <- sub_data


cov_mat <- covariate_df_to_mat(joined_data, cov_names = covariate_names)
## hard coded best params
best_pars_list <- data_list$output_list[[best_model_number]]
best_pars <- best_pars_list$est_pars[,1]
rownames(best_pars) <- NULL
joined_data$prob_trans <- 1 / (1 + exp(-cov_mat %*% best_pars))

           
joined_dt <- data.table::as.data.table(joined_data)
like_df <- joined_dt[, .(like = general_cluster_like.dt(prob_trans,
                                                     n_inf)),
                  by = .(cluster_id)]
## join likelihood back
joined_data2 <- left_join(joined_data, like_df, by = "cluster_id")
joined_data2$like2 <- joined_data2$like / sum(joined_data2$like)
joined_data3 <- joined_data2 %>%
    filter(gen!= 1)

df <- joined_data3 %>% group_by(cluster_id) %>%
    mutate(x = assign_coords_x(gen)) %>%
    group_by(cluster_id, gen) %>%
    mutate(y = assign_coords_y(n_in_gen)) %>%
    ungroup() %>%
    group_by(cluster_id) %>%
    mutate(gen_size = max(gen)) %>%
    arrange(desc(like2))

top_groups_by_gen <- df %>%
    ungroup() %>%
    group_by(gen_size, cluster_id) %>%
    summarize(prob = min(like2),
              .groups = "drop_last") %>%
    top_n(n = 3, wt = prob)

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
    summarize(like = like2[1])
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
    summarize(like = like2[1]) %>% 
              mutate(id = str_sub(cluster_id, start = -3))
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


                        
ggplot(data = jsub,
       aes(x = x, y = y,
           group = cluster_id,
           col = feature_id,
           shape = factor(hiv_neg_pos))) +
    geom_curve(aes(xend = x_to, yend = y_to),
               size = 2, curvature = -0,
               col = "black") +
    geom_point(size = 6, stroke = 4, col = "darkgray") +
    geom_point(size = 5, stroke = 4) +
    facet_wrap(~facet_lab, ncol = 3) +
    scale_color_brewer(palette = "Set1", name = "Person ID")  +
    scale_shape_manual(values = c(3, 16),
                       labels = c("Pos./Unk.", "Neg."),
                       name = "HIV Status") +
    xlim(1.8, 8.2) +
    ylim(-1.2, 1.2) +
    theme_bw(base_size = 18) +
    theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Generation in transmission tree",
         y = latex2exp::TeX("Order of infection in gen. among 'siblings' $\\rightarrow$"),
         title = "Most likely sampled trees by number of generations",
         subtitle = paste0("Cluster ", my_orig_id)) 
   
ggsave("trees-7.pdf", height = 15, width = 13)
