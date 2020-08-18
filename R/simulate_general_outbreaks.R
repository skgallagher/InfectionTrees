## SKG
## May 29, 2020
##
## Generate clusters based on an outsider
## Basically we will relabel the no-outsider simulations
## All generations go down by one number
## Will take out the index person
## update i accordingly




#' Simulate the branching process of flipping until failure for K clusters
#'
#' @param K number of total clusters to simulate
#' @param inf_params vector with beta coefficients to use in logistic function for probability of transmission
#' @param sample_covariates_df Data frame of covariates to sample from
#' @param covariate_names names of the covariates.  Must match size of inf_params - 1.
#' @param covariate_weights default is NULL which draws uniformly at random with replacement from the sample_covariates_df.  Otherwise, the weights are used.
#' @param max_size maximum size a cluster can be
#' @return data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{person_id}{order of infection in the cluster}
#' \item{gen}{generation number (>=0)}
#' \item{inf_id}{ID of the infector}
#' \item{n_inf}{number of people infected by person}
#' \item{censored}{whether the cluster end was censored or not}
#' \item{cluster_size}{size of the cluster}
#' \item{covariates}{covariates of the individuals}
#' }
#' @details Generate a branching process according to the following process.
#' First a root infector is drawn covariates \eqn{X}
#'  from some distribution $F$ (given by the set of covariates in \code{sample_covariates_df}) and has probability of
#'   transmission according to a logit function.
#'   The number of infections produced by the root node $N_{(1,1)}
#'    is a geometric random variable with probability $p_{(1,1)}$
#'   where the indexing represents $(g=$, generation, $i=$ index).
#'    If $N_{(1,1)} > 0$, then the $N_{(1,1)}$ infections are added to
#'     the cluster and assigned to generation $g=2$ with indices $i=1,
#'      \dots, N_{(1,1)}$ and covariats are drawn for these new infections.
#'      The infection process continues with individuals $(2, 1)$ through $(2, $N_{(1,1)})$
#'       where new infections are added, in order to the subsequent generation.
#'        The process terminates when either there are no new infections or the
#'         maximum number of infections specified in \code{max_size} is reached.
#' \deqn{X_{(g,i)} \sim F}
#' \deqn{p_{(g,i)} = logit^{-1}\left ( X_{(g,i)} \beta\right )}
#' \deqn{N_{(g,i)} \sim Geometric(p_{(g,i)})}
#' @export
#' @examples
#' set.seed(2020)
#' inf_params <- c("beta_0" = -2, "beta_1" = 2)
#' df <- data.frame(x= c(0, 1))
#' branching_processes <- simulate_bp(K = 10,
#' inf_params = inf_params,
#' covariate_names = "x",
#' sample_covariates_df = df)
#' head(branching_processes)
#' table(branching_processes$cluster_size) /
#' sort(unique(branching_processes$cluster_size))
simulate_bp <- function(K,
                              inf_params,
                              sample_covariates_df,
                              covariate_names,
                              covariate_weights = NULL,
                              max_size = 50){

    sample_covariates_df <- sample_covariates_df %>%
        dplyr::select(tidyselect::all_of(covariate_names))

    stopifnot((length(inf_params) - 1) == ncol(sample_covariates_df))

    cov_mat <- as.matrix(sample_covariates_df)
    cov_mat <- cbind(rep(1, nrow(cov_mat)), cov_mat)

    sample_covariates_df$prob_inf <- 1 / (1 + exp(-cov_mat %*% inf_params))

    cluster_list <- vector(mode = "list", length = K)
    for(ii in 1:K){
        cluster_list[[ii]] <- simulate_general_outbreak_inner(cluster_id = ii,
                                                              sample_covariates_df,
                                                              covariate_weights,
                                                              max_size = max_size)
    }

    clusters <- dplyr::bind_rows(cluster_list)

    return(clusters)

}

#' Generate a single cluster of infections
#'
#' @param cluster_id ID of cluster
#' @param sample_covariates_df data frame of covariates with pre-calculated probability of infection \code{prob_inf}
#' @param covariate_weights default is NULL which draws uniformly at random with replacement from the sample_covariates_df.  Otherwise, the weights are used.
#' @param max_size maximum size of cluster
#' @return data frame with the following columns
#' \describe{
#' \item{cluster_id}{unique cluster ID}
#' \item{person_id}{order of infection in the cluster}
#' \item{gen}{generation number (>=0)}
#' \item{inf_id}{ID of the infector}
#' \item{n_inf}{number of people infected by person}
#' \item{cluster_pos}{number of positive smears in the cluster}
#' \item{cluster_size}{number in cluster}
#' \item{censored}{whether the cluster end was censored or not}
#' \item{covariates}{other covariates in sample_covariates_df}
#' }
#' @details randomly assigns covariates from sample_df.   breadth not depth.  Generate generation by generation as opposed to going up the branch til termination.
simulate_general_outbreak_inner <- function(cluster_id,
                                            sample_covariates_df,
                                            covariate_weights = NULL,
                                            max_size = 50){

    df <- data.frame(cluster_id = rep(cluster_id, max_size),
                     person_id = NA,
                     gen = NA,
                     inf_id = NA,
                     n_inf = NA,
                     cluster_size = 1,
                     censored = FALSE)

    censored <- FALSE


    cluster_size <- 1
    gen_list <- vector(mode = "list", length = max_size)
    cov_1 <- draw_covariate_rows(1, sample_covariates_df, covariate_weights)
    gen_list[[1]] <- cbind(data.frame(cluster_id = cluster_id,
                                person_id = paste0("C", cluster_id, "-G1-N1"),
                                gen = 1,
                                inf_id = NA,
                                n_inf = NA,
                                cluster_size = 1,
                                stringsAsFactors = FALSE),
                           cov_1)
    for(gg in 2:(max_size - 1)){
        next_gen_output <- general_generation_infection(cluster_id = cluster_id,
                                                        sample_covariates_df,
                                                        covariate_weights,
                                                        gen = gg,
                                                        prev_gen = gen_list[[gg-1]],
                                                        cluster_size = cluster_size,
                                                        max_size = max_size)
        if(cluster_size == next_gen_output$new_cluster_size){ # no new generations
            break
        }
        cluster_size <- next_gen_output$new_cluster_size
        gen_list[[gg]] <- next_gen_output$cur_gen
        gen_list[[gg-1]]$n_inf <- next_gen_output$n_inf
        if(max_size <= cluster_size){ # Hit the max size
            censored <- TRUE
            break
        }



    }
    sim_df <- dplyr::bind_rows(gen_list)
    sim_df$cluster_size <- nrow(sim_df)
    sim_df$censored <- censored
    sim_df$cluster_size <- cluster_size
    return(sim_df %>% dplyr::select(-.data$prob_inf))


}





#' Infect the current generation given the previous generation
#'
#' @param cluster_id unique cluster ID
#' @param df covariate data frame.  One is prob_inf
#' @param covariate_weights either NULL or weight of covariates
#' @param gen new generation number
#' @param prev_gen data frame with at least smear status (0/1) and person_id
#' @param cluster_size size of previous cluster
#' @param max_size maximum size of cluster
#' @return list with the following entries
#' \describe{
#' \item{new_cluster_size}{new cluster size with the current generation}
#' \item{cur_gen}{data frame of current generation (new infections)}
#' \item{n_inf}{number of infections for each person in the previous generation}
#' }
#' @importFrom rlang .data
general_generation_infection <- function(cluster_id,
                                         df,
                                         covariate_weights = NULL,
                                         gen,
                                         prev_gen,
                                         cluster_size,
                                         max_size){


    n_new_inf <- rep(0, nrow(prev_gen))
    p <- prev_gen$prob_inf
    p <- ifelse(p == 1, .99999999, p) # pretty hacky
    for(ii in 1:nrow(prev_gen)){
        n_new_inf[ii] <- stats::rgeom(n = 1, prob = 1 - p[ii])
        new_cluster_size <- cluster_size +
            sum(n_new_inf[1:ii])
        if(new_cluster_size >= max_size){
                                        #  browser()
            if(ii > 1){
                n_new_inf[ii] <- (max_size -  cluster_size - sum(n_new_inf[1:(ii-1)]))
            } else {
                n_new_inf[ii] <- max_size - cluster_size
            }
            break  # Have found enough
        }
    }
                                        #  browser()
    n <- sum(n_new_inf)
    if(n == 0){ # No new infections in next generations
        return(list(new_cluster_size = cluster_size,
                    cur_gen = NULL,
                    n_inf = n_new_inf))
    }
    new_cluster_size <- n + cluster_size
   # browser()

    cov_df <- draw_covariate_rows(n, df, covariate_weights)
    cur_gen <- cbind(data.frame(cluster_id = rep(cluster_id, n),
                          person_id = paste0("C", cluster_id, "-G", gen, "-N", 1:n),
                          gen = gen,
                          inf_id = rep(prev_gen$person_id, times = n_new_inf),
                          n_inf = NA,
                          cluster_size = NA,
                          stringsAsFactors = FALSE),
                     cov_df)

    return(list(new_cluster_size = n + cluster_size,
                cur_gen = cur_gen,
                n_inf = n_new_inf))

}


#' Randomly sample covariate rows
#'
#' @param n number rows to draw
#' @param df data frame to draw from
#' @param covariate_weights default is NULL which draws uniformly at random with replacement from the sample_covariates_df.  Otherwise, the weights are used.
#' @return new df
draw_covariate_rows <- function(n, df,
                                covariate_weights = NULL){

    if(is.null(covariate_weights)){
        row_inds <- sample(1:nrow(df), size = n, replace = TRUE)
        new_df <- df[row_inds,]
        rownames(new_df) <- NULL
    } else {

        weights <- covariate_weights / sum(covariate_weights)
        row_inds <- sample(1:nrow(df), size = n, replace = TRUE,
                           prob = weights)
        new_df <- df[row_inds,]

        rownames(new_df) <- NULL

    }
    return(new_df)
}
