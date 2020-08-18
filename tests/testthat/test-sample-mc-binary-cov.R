
test_that("summarize_binary_trees", {
    mc_trees <- data.frame(cluster_id = c(1, 1, 1,
                                          2, 2, 2,
                                          3, 3, 3),
                           n_inf = c(2, 0, 0,
                                     1, 1, 0,
                                     1, 0, 1),
                           x = c(1, 0, 0,
                                 0, 1, 0,
                                 0, 0, 1),
                           gen = c(1, 2, 2,
                                   1, 2, 3,
                                   1, 3, 2))
    multiple_outside_transmissions <- FALSE

    out <- summarize_binary_cov_trees(mc_trees,
                                  multiple_outside_transmissions =
                                      multiple_outside_transmissions)
    expect_true(all(out$x_pos == 1))
    expect_true(all(out$x_neg == 2))
    expect_equal(as.numeric(out[1, 3:5]), c(1, 1, 2))
    expect_equal(as.numeric(out[2, 3:5]), c(2, 0, 1))

    ## tests for outside
    ## TODO

})


test_that("sample_mc_binary_cov_inner", {
    x_pos <- 2
    x_neg <- 1
    B <- 5
    root_node <- NULL
    out <- sample_mc_binary_cov_inner(x_pos,
                                        x_neg,
                                        B,
                                      root_node)
    expect_equal(nrow(out), B * (x_pos + x_neg))
    expect_equal(x_pos * B, sum(out$x))


    ## tests for outside
    ## TODO

})


test_that("sample_mc_binary_cov", {
    B <- 5
    observed_cluster_summaries <- data.frame(freq = c(4, 2),
                                             cluster_size = c(1, 1),
                                             x_pos = c(1,0),
                                             x_neg = c(0, 1))

    out <- sample_mc_binary_cov(B,
                                observed_cluster_summaries)
    expect_equal(out$mc_freq, c(5, 5))
    expect_equal(out$freq, c(4, 2))
    ##
    B <- 10
    observed_cluster_summaries <- data.frame(freq = c(4),
                                             cluster_size = c(3),
                                             x_pos = c(2),
                                             x_neg = c(1))

    out <- sample_mc_binary_cov(B,
                                observed_cluster_summaries)
    expect_equal(sum(out$mc_freq), B)
})
