## May 27, 2020

test_that("general_loglike", {

    ## with data.table
     observed_data <- data.frame(cluster_id = c(1,
                                               2, 2,
                                               3, 3),
                                x = c(1,
                                      0, 0,
                                      1, 1))
    B <- 1
    covariate_names <- "x"

    sampled_data <- sample_permuted_cond_trees(observed_data,
                                      B = B,
                                      covariate_names = "x")
    inf_params <- c(-2, 1)
    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = "x")
    out <- general_loglike(inf_params,
                           observed_data,
                           data.table::as.data.table(sampled_data),
                           return_neg = FALSE,
                           cov_mat = cov_mat,
                           cov_names = "x")

    p_neg <- 1/ (1 + exp(2))
    p_pos <- 1 / (1 + exp(1))
    exp_out <- log((1 - p_pos) *
                   (1-p_neg)^2 * p_neg *
                   (1-p_pos)^2 * p_pos
                   )
    expect_equal(exp_out, out)



    #######################################################3
    ## Uses original data
    ## Permutes the clusters so better representation of the average loglike for a  particular cluster as opposed to sampling from all clusters
    observed_data <- data.frame(cluster_id = c(1,
                                               2, 2,
                                               3, 3),
                                x = c(1,
                                      0, 0,
                                      1, 1))
    B <- 1
    covariate_names <- "x"

    sampled_data <- sample_permuted_cond_trees(observed_data,
                                      B = B,
                                      covariate_names = "x")
    inf_params <- c(-2, 1)
    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = "x")
    out <- general_loglike(inf_params,
                           observed_data,
                           sampled_data,
                           return_neg = FALSE,
                           cov_mat = cov_mat,
                           cov_names = "x")

    p_neg <- 1/ (1 + exp(2))
    p_pos <- 1 / (1 + exp(1))
    exp_out <- log((1 - p_pos) *
                   (1-p_neg)^2 * p_neg *
                   (1-p_pos)^2 * p_pos
                   )
    expect_equal(exp_out, out)

    ## two variables
    observed_data <- data.frame(cluster_id = c(1, 1,
                                               2, 2),
                                x = c(1, 1, 0, 0),
                                y = c(1, 1, 1, 1))
    B <- 10
    covariate_names <- c("x", "y")

    sampled_data <- sample_permuted_cond_trees(observed_data,
                                      B = B,
                                      covariate_names = c("x", "y"))
    inf_params <- c(-2, .5, .5)
    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = c("x", "y"))
    out <- general_loglike(inf_params,
                           observed_data,
                           sampled_data,
                           return_neg = FALSE,
                           cov_mat = cov_mat,
                           cov_names = c("x", "y"))

    p_11 <- 1 / ( 1 + exp(1))
    p_01 <- 1 / ( 1 + exp(1.5))
    exp_out <- log((1 - p_11)^2 * p_11)  +
        log((1 - p_01)^2 * p_01)

    expect_equal(exp_out, out)


})

test_that("general_loglike_summarized", {
    inf_params <- c(1, 0)
    obs_summarized_data <- data.frame(cluster_size = 1,
                                      freq = 3)

    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                         cluster_size = 4,
                         prob_inf = c(.8, .1, .1, .1),
                         n_inf = c(2, 0, 1, 0))
    clust_3 <- data.frame(cluster_id = "c",
                          cluster_size = 1,
                          prob_inf = .2,
                          n_inf = 0)
    sampled_data <- rbind(clust_1, clust_2, clust_3)

    sampled_data$x <- rnorm(nrow(sampled_data),
                            mean = 0,
                            sd = .5)

    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = c("x"))
    expect_equal(cov_mat[,2], sampled_data$x)

    out <- general_loglike_summarized(inf_params,
                           obs_summarized_data,
                           sampled_data,
                           return_neg = FALSE,
                           cov_mat = cov_mat)
    expect_equal(length(out), 1)  ### Ok this tells like nothing lol
    p <- 1 / (1 + exp(- cov_mat %*% inf_params))
    expect_equal(3 * log(1 - p[1]), out)
    ## ###############################
    ##  # now x is 1/0

    inf_params <- c(0, 1)
    obs_summarized_data <- data.frame(cluster_size = 1,
                                      freq = 3)

    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                          cluster_size = 4,
                          prob_inf = c(.8, .1, .1, .1),
                          n_inf = c(2, 0, 1, 0))
    clust_3 <- data.frame(cluster_id = "c",
                          cluster_size = 1,
                          prob_inf = .2,
                          n_inf = 0)
    sampled_data <- rbind(clust_1, clust_2, clust_3)

    sampled_data$x <- c(1,
                        1, 0, 1, 0,
                        0)

    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = c("x"))
    expect_equal(cov_mat[,2], sampled_data$x)

    out <- general_loglike_summarized(inf_params,
                           obs_summarized_data,
                           sampled_data,
                           return_neg = FALSE,
                           cov_mat = cov_mat)
    expect_equal(length(out), 1)
    p_mat <- 1 / (1 + exp(- cov_mat %*% inf_params))
    exp_out <- 3 * log(1/2 * ((1 - p_mat[1]) +(1 - p_mat[6])))
    expect_equal(exp_out, out)

    ## ################################3
    ## # And more than one freq
      ##  # now x is 1/0

    inf_params <- c(0, 1)
    obs_summarized_data <- data.frame(cluster_size = c(1, 4),
                                      freq = c(3, 1))

    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                          cluster_size = 4,
                          prob_inf = c(.8, .1, .1, .1),
                          n_inf = c(2, 0, 1, 0))
    clust_3 <- data.frame(cluster_id = "c",
                          cluster_size = 1,
                          prob_inf = .2,
                          n_inf = 0)
    sampled_data <- rbind(clust_1, clust_2, clust_3)

    sampled_data$x <- c(1,
                        1, 0, 1, 0,
                        0)

    cov_mat <- covariate_df_to_mat(sampled_data,
                                   cov_names = c("x"))
    expect_equal(cov_mat[,2], sampled_data$x)

    out <- general_loglike_summarized(inf_params,
                           obs_summarized_data,
                           sampled_data,
                           return_neg = FALSE,
                           cov_mat = cov_mat)
    expect_equal(length(out), 1)
    p_mat <- 1 / (1 + exp(- cov_mat %*% inf_params))
    exp_out <- 3 * log(1/2 * ((1 - p_mat[1]) +(1 - p_mat[6]))) +
        log((1-p_mat[1]) * p_mat[1]^2 *
            (1-p_mat[6]) *
            (1-p_mat[1]) * p_mat[1] *
            (1-p_mat[6]))

    expect_equal(exp_out, out)

})



test_that("general_loglike_sampled_data", {
    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                         cluster_size = 4,
                         prob_inf = c(.8, .1, .1, .1),
                         n_inf = c(2, 0, 1, 0))
    sampled_data <- rbind(clust_1, clust_2)

    out <- general_loglike_sampled_data(sampled_data)
    exp_out1 <- log(.6)
    exp_out2 <- log(.2 * .8^2 *
                    .9 *
                    .9 * .1 *
                    .9)
    expect_equal(out$avg_loglike[1], exp_out1)
    expect_equal(out$avg_loglike[2], exp_out2)

    ## # averages yay
    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                         cluster_size = 4,
                         prob_inf = c(.8, .1, .1, .1),
                         n_inf = c(2, 0, 1, 0))
    clust_3 <- data.frame(cluster_id = "c",
                          cluster_size = 1,
                          prob_inf = .2,
                          n_inf = 0)
    sampled_data <- rbind(clust_1, clust_2, clust_3)

    out <- general_loglike_sampled_data(sampled_data)
    exp_out1 <- log((.6 + .8) / 2)
    exp_out2 <- log(.2 * .8^2 *
                    .9 *
                    .9 * .1 *
                    .9)
    expect_equal(out$avg_loglike[1], exp_out1)
    expect_equal(out$avg_loglike[2], exp_out2)

    ## #One more for good measure
    ## # averages yay
    clust_1 <- data.frame(cluster_id = "a",
                          cluster_size = 1,
                          prob_inf = .4,
                          n_inf = 0)
    clust_2 <- data.frame(cluster_id = "b",
                         cluster_size = 4,
                         prob_inf = c(.8, .1, .1, .1),
                         n_inf = c(2, 0, 1, 0))
    clust_3 <- data.frame(cluster_id = "c",
                          cluster_size = 4,
                          prob_inf = c(.2, .3, .2, .2),
                          n_inf = c(1, 2, 0, 0))
    sampled_data <- rbind(clust_1, clust_2, clust_3)

    out <- general_loglike_sampled_data(sampled_data)
    exp_out1 <- log(.6)
    exp_out2 <- log((.2 * .8^2 *
                    .9 *
                    .9 * .1 *
                    .9 +
                    .8 * .2^1 *
                    .7 * .3^2 *
                    .8 *
                    .8
                    )/2)
    expect_equal(out$avg_loglike[1], exp_out1)
    expect_equal(out$avg_loglike[2], exp_out2)


})



test_that("general_cluster_like", {

    data <- data.frame(prob_inf = .3,
                       n_inf = 0)

    out <- general_cluster_like(data)
    expect_equal(.7, out)
    ##
    data <- data.frame(prob_inf = c(.3, .2),
                       n_inf = c(1, 0))
    out <- general_cluster_like(data)
    expect_equal(.7 * .3 * .8 , out)

})
