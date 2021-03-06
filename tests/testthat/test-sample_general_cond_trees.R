test_that("sample_mc_tree", {
## the outsider version

  observed_data <- data.frame(cluster_id = c(1,
                                             2, 2,
                                             3, 3),
                              x = c(1,
                                    0, 1,
                                    1, 1))
  B <- 3
  covariate_names <- "x"
  multiple_outside_transmissions <- TRUE
  out <- sample_mc_trees(observed_data,
                         B = B,
                         covariate_names = covariate_names,
                         multiple_outside_transmissions = multiple_outside_transmissions)

  expect_equal(nrow(out), 8 * B)


})

test_that("sample_mc_trees", {
    observed_data <- data.frame(cluster_id = c(1,
                                               2, 2,
                                               3, 3),
                                x = c(1,
                                      0, 1,
                                      1, 1))
    B <- 3
    covariate_names <- "x"

    out <- sample_mc_trees(observed_data,
                                      B = B,
                                      covariate_names = "x")
    expect_equal(nrow(out), B * nrow(observed_data))
    expect_true(all(1:3 %in% out$orig_id))
    df3 <- out %>% dplyr::filter(orig_id == 3) %>%
                 dplyr::select(x)
    expect_equal(unique(df3$x),
                 1)


    ## Check that each row is used once
    observed_data <- data.frame(cluster_id = rep(1, 8),
                                x = 1:8)
    sampled_data <- sample_mc_trees(observed_data,
                                      B = 10,
                                      covariate_names = "x")
    expect_true(all(table(sampled_data$x) == 10))

})


test_that("draw_features", {
    tree_df <- data.frame(id = 1:5)
    feature_type <- "empirical"
    params_list <- list(weights = NULL,
                        covariate_df = data.frame(x = c(1:3)))
    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_equal(nrow(out), 5)

    ##########################################3
    ## draw features
    ## "binary_cov"
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov"
    params_list <- list(x_pos = 3,
                        x_neg = 0,
                        root_node = NULL)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_equal(out$x, rep(1, nrow(tree_df)))
    expect_equal(nrow(out), nrow(tree_df))
    ##
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov"
    params_list <- list(x_pos = 2,
                        x_neg = 1,
                        root_node = NULL)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_equal(sum(out$x), params_list$x_pos * 4)
    expect_equal(nrow(out), nrow(tree_df))
     ##
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov"
    params_list <- list(x_pos = 1,
                        x_neg = 2,
                        root_node = NULL)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_equal(sum(out$x), params_list$x_pos * 4)
    expect_equal(nrow(out), nrow(tree_df))
                                        #
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov"
    params_list <- list(x_pos = 0,
                        x_neg = 3,
                        root_node = NULL)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_equal(sum(out$x), params_list$x_pos * 4)
    expect_equal(nrow(out), nrow(tree_df))
    ##
    ## binary_cov_out
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov_out"
    params_list <- list(x_pos = 2,
                        x_neg = 0,
                        root_node = 1)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_true(all(out$x[c(1, 4, 7, 10)] == 1))
    expect_equal(nrow(out), nrow(tree_df))
    ##
    tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov_out"
    params_list <- list(x_pos = 1,
                        x_neg = 1,
                        root_node = 1)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_true(all(out$x[c(1, 4, 7, 10)] == 1))
    expect_equal(nrow(out), nrow(tree_df))
    ##
       tree_df <- data.frame(cluster_id = c(1, 1, 1,
                                         2, 2, 2,
                                         3, 3, 3,
                                         4, 4, 4))
    feature_type <- "binary_cov_out"
    params_list <- list(x_pos = 0,
                        x_neg = 2,
                        root_node = 1)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_true(all(out$x[c(1, 4, 7, 10)] == 1))
    expect_equal(nrow(out), nrow(tree_df))
    ##
    tree_df <- data.frame(cluster_id = c(1, 1,
                                         2, 2,
                                         3, 3))
    feature_type <- "binary_cov_out"
    params_list <- list(x_pos = 1,
                        x_neg = 0,
                        root_node = 1)

    out <- draw_features(tree_df,
                         feature_type,
                         params_list)
    expect_true(all(out$x[c(1, 3, 5)] == 1))
    expect_equal(nrow(out), nrow(tree_df))




})



test_that("sample_general_cond_trees", {
    n_vec <- 1:50
    B <- 10
    out <- sample_general_cond_trees(n_vec,
                                     B)
    expect_equal(nrow(out), sum(n_vec) * B)
    expect_equal(length(unique(out$cluster_id)),
                 length(n_vec) * B)

    ## empirical
    n_vec <- 5
    B <- 10
    feature_type <- "empirical"
    params_list <- list(weights = NULL,
                        covariate_df = data.frame(x = c(1:3)))
    out <- sample_general_cond_trees(n_vec,
                                     B,
                                     feature_type = "empirical",
                                     params_list = params_list)
    expect_equal(nrow(out), sum(n_vec) * B)
    expect_equal(length(unique(out$cluster_id)),
                 length(n_vec) * B)

    ## empirical 2 variables whoa
    n_vec <- 5
    B <- 10
    feature_type <- "empirical"
    params_list <- list(weights = NULL,
                        covariate_df = data.frame(x = c(0,1, 0, 1),
                                                  y = c(0, 0, 1, 1)))
    out <- sample_general_cond_trees(n_vec,
                                     B,
                                     feature_type = "empirical",
                                     params_list = params_list)
    expect_equal(nrow(out), sum(n_vec) * B)
    expect_equal(length(unique(out$cluster_id)),
                 length(n_vec) * B)




})





test_that("sample_unif_trees_no_features", {
    n_vec <- 1:5
    B <- 10
    out <- sample_unif_trees_no_features(n_vec, B)
    expect_equal(nrow(out), sum(n_vec) * B)
    expect_equal(length(unique(out$cluster_id)),
                 length(n_vec) * B)
})

test_that("general_tree_sampler", {
    n <- 1
    B <- 10
    out <- general_tree_sampler(n, B)
    expect_equal(nrow(out), B)
    ## ######################3
    n <- 3
    B <- 10
    out <- general_tree_sampler(n, B)
    expect_equal(nrow(out), B * n)
    expect_equal(length(unique(out$cluster_id)), B)
    expect_equal(sum(is.na(out$inf_id)), B)
    ## ############################
        ## ######################3
    n <- 5
    B <- 10
    out <- general_tree_sampler(n, B)
    expect_equal(nrow(out), B * n)
    expect_equal(length(unique(out$cluster_id)), B)
    expect_equal(sum(is.na(out$inf_id)), B)
})







test_that("sample_tree_perm", {
    gen_sizes <- c(1, 1, 1)
    out <- sample_general_tree_perm(gen_sizes)
    expect_equal(nrow(out), sum(gen_sizes))
    expect_equal(out$n_inf[1], 1)
    expect_equal(out$n_inf[3], 0)
    ## ##########################
    gen_sizes <- c(1, 2)
    out <- sample_general_tree_perm(gen_sizes)
    expect_equal(nrow(out), sum(gen_sizes))
    expect_equal(out$n_inf[1], 2)
    expect_equal(out$n_inf[3], 0)
    ## ##########################
    gen_sizes <- c(1, 2, 1)
    out <- sample_general_tree_perm(gen_sizes)
    expect_equal(nrow(out), sum(gen_sizes))
    expect_equal(out$n_inf[1], 2)
    expect_equal(out$n_inf[sum(gen_sizes)], 0)

})

test_that("sample_trees_fixed_g", {

    perm_mat <- matrix(c(1,1,
                         1, 1,
                         1, 1), ncol = 3)
    out <- sample_trees_fixed_g(perm_mat)
    expect_equal(nrow(out), 6)
    expect_equal(length(unique(out$cluster_id)), 3)
    ## #######################3

    perm_mat <- matrix(c(1, 1, 2,
                         1, 2, 1,
                         1, 2, 1), nrow = 3)
    out <- sample_trees_fixed_g(perm_mat)
    expect_equal(nrow(out), 12)
    expect_equal(length(unique(out$cluster_id)), 3)
    expect_equal(sum(is.na(out$inf_id)), 3)

})

