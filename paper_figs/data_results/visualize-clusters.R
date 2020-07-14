## SKG
## Basic coding to visualize clusters
## June 19, 2020

library(tidyverse)
library(data.table)
devtools::load_all()


files <- list.files()
fn_rds <- grep('^(?!.*outsider).*RDS', files,
               value = TRUE, perl = TRUE)
fn <- max(fn_rds)

## Let's look at a cluster of size x
data_list <- readRDS(fn)

sampled_data <- data_list$sampled_data
covariate_names <- c("scaled_rel_time")
best_model_number <- 8

my_orig_id <- 27 #121 is for 5 # 27 is for a cool one with 7

sub_data <- sampled_data %>%
    filter(cluster_size == 5)

sub_data <- sampled_data %>% filter(orig_id == my_orig_id)

unique_covs <- sub_data %>%
    distinct(smear, hiv_pos, rel_time,
             scaled_rel_time, .keep_all = TRUE) %>%
    select(contains(c(covariate_names))) %>%
    mutate(person_id = row_number())

joined_data <- left_join(sub_data, unique_covs)


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



## DO next
## Get probability of transmission from MLE
## Get the likelihood for each cluster
## Make df with |person id| rel. time | cluster_like| person1_weighted_prob | ... | personK_weighted_prob|


## Ok but do it right htis time

mat <- matrix(0, ncol = length(unique(joined_data$person_id)),
              nrow = nrow(joined_data))
cluster_ids <- unique(joined_data2$cluster_id)
for(ii in 1:length(cluster_ids)){
    if((ii %% 100) == 0) print(ii)
    my_cluster_id <- cluster_ids[ii]
    sub_data <- joined_data2 %>% filter(cluster_id == my_cluster_id)
    for(jj in 1:nrow(sub_data)){
        id <- sub_data$id[jj]
        for(kk in 1:nrow(sub_data)){
            if(!is.na(sub_data$inf_id[kk]) &
               sub_data$inf_id[kk] == id){
                infected_person_id <- sub_data$person_id[kk]
                mat[(ii-1) * nrow(sub_data)  + jj,
                    infected_person_id] <- sub_data$prob_trans[jj] *
                    sub_data$like2[jj]
            }
        }
    }
}
mat <- mat / sum(mat)
colnames(mat) <- paste0("person_id", 1:ncol(mat))

## Putting it all together
fig_df <- cbind(joined_data2, mat) %>%
    rename(index_person = person_id) %>%
    select(index_person, contains("person_id")) 

mean_weights <- fig_df %>% group_by(index_person) %>%
    summarize_all(sum)

feature_df <- joined_data2 %>% distinct(person_id, rel_time,
                                        hiv_pos,
                                        smear,
                                        scaled_rel_time) %>%
    rename(index_person = person_id)


full_df <- left_join(feature_df, mean_weights,
                     by = "index_person")

long_df <- full_df %>%
    pivot_longer(cols = contains("person_id")) %>%
    mutate(infected_person = parse_number(name))



## LINEAR
feature_df$x <- feature_df$rel_time
feature_df$y <- 0

ggplot(feature_df,
       aes(x = x, y = y)) + geom_point(size = 5) +
    geom_label(aes(label = round(rel_time, 2)),
               nudge_x = .1) +
    ylim(-1, 1)



long_df2 <- left_join(long_df, feature_df,
                      by = c("index_person", "smear", "hiv_pos", "rel_time",
                             "scaled_rel_time")) %>%
    rename(x_from = "x",
           y_from = "y")

long_df3 <- left_join(long_df2, feature_df %>% select(index_person, x, y),
                      by = c("infected_person" = "index_person")) %>%
    rename(x_to = "x",
           y_to = "y")

final_df <-  long_df3 %>% filter(index_person != infected_person)


ggplot(final_df,
       aes(x = x_from, y = y_from)) +  
    ## xlim(-1, 1.2) +
    ylim(-.5, .5) +
    geom_hline(yintercept = 0, linetype  = "dashed") + 
    geom_curve(aes(x = x_from, y = y_from,
                   xend = x_to, yend = y_to,
                   col = value * 100),
               arrow = arrow(length = unit(0.03, "npc"),
                             type = "closed"),
               curvature = .5,
               alpha = 1,
               size = 2) +
    scale_color_viridis_c(end = .9,
                          limits = c(0, 4),
                          name = latex2exp::TeX("Prob. A infects B (%)")) +
      geom_label(aes(label = paste(round(rel_time, 2), "years")),
                 nudge_y = -.02, size = 8) +
    geom_point(size = 5, shape = 22,
               aes(fill = factor(hiv_pos))) +
    scale_fill_manual(values = c("red", "black"),
                      labels = c("HIV +", "HIV -/unk"),
                      name = "HIV Status") +
    theme_void(base_size = 25) +
    labs(title = "Estimated transmission probability of TB from person A to B",
         subtitle = sprintf("Individuals in cluster %d", my_orig_id)) +
    annotate("text", x = .6, y = .02,
             label =  latex2exp::TeX("Rel. time of detection from min time $\\rightarrow$"),
             size = 6)
   # theme(legend.position = "bottom")

ggsave(paste0("graph_probs_network_cluster_", my_orig_id, ".pdf"),
       width = 14)





##################################3
## dumb way to get person_inf_id
## very slow
old_cluster_id <- 0
joined_data$inf_person_id <- NA
for(ii in 1:nrow(joined_data)){
    my_cluster_id <- joined_data$cluster_id[ii]
    inf_id <- joined_data$inf_id[ii]
    if(old_cluster_id != my_cluster_id){
        sub_data <- joined_data %>% filter(cluster_id == my_cluster_id)
    }
    ind <- which(sub_data$id == inf_id)
    if(length(ind) > 0){
        joined_data$inf_person_id[ii] <- sub_data$person_id[ind]
    }

}
