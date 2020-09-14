test_that("mc_trees_A_infects_B",{
    cluster_id <- "A"
    par_ests <- c("intercept" = -2,
                  "x" = 1)
    set.seed(42)
    observed_data <- data.frame(cluster_id = c("A", "A", "A", "A"),
                                x = c(1, 0, 0, 1))
    mc_trees <- sample_mc_trees(observed_data = observed_data,
                                B = 10,
                                covariate_names = "x")

    out <- mc_trees_to_A_infects_B(my_id = "A",
                                   par_ests = par_ests,
                                   mc_trees = mc_trees)
    expect_true(nrow(out) == 4)


})


test_that("compute_A_infects_B", {

    my_df <- data.frame(cluster_id = c("A", "A", "A",
                                       "B", "B", "B"),
                        prob_trans = c(.2, .8, .4,
                                       .2, .8, .4),
                        like = c(.3, .3, .3,
                                 .2, .2, .2),
                        id = c("n1", "n2", "n3",
                               "n1", "n2", "n3"),
                        inf_id = c(NA, "n1", "n2",
                                   "n3", "n1", NA))
    my_df$feature_id <- factor(my_df$prob_trans)
    my_df$like2 <- c(.6, .6, .6,
                     .4, .4, .4)
    like2 <- my_df$like2
    out <- compute_A_infects_B(df = my_df)
    exp_out <- matrix(0, nrow = 3, ncol = 3)
    exp_out[1, 3] <- .2 * like2[2] + .2 * like2[5]
    exp_out[3, 2] <- .8 * like2[3]
    exp_out[2, 1] <- .4 * like2[4]
    expect_equal(as.numeric(out),
                 as.numeric(exp_out))
    
    

})


test_that("calculate_transmission_prob", {
  df <- data.frame(x = c(1, 0, 1, 0),
                     y = c(-3, -11, 2, 4))
  inf_params <- c(-2, 1)
  cov_names <- c("x")

  out <- calculate_transmission_prob(data = df,
                                     inf_params = inf_params,
                                     covariate_names = cov_names)

  exp_out <- 1 / (1 + exp(2 - c(1, 0, 1, 0)))
  expect_equal(out, exp_out)
###########################################3
  ### using named inf_params
  df <- data.frame(x = c(1, 0, 1, 0),
                   y = c(-3, -11, 2, 4))
  inf_params <- c("intercept" = -2, "x" = 1)

  out <- calculate_transmission_prob(data = df,
                                     inf_params = inf_params)

  exp_out <- 1 / (1 + exp(2 - c(1, 0, 1, 0)))
  expect_equal(out, exp_out)
  ## ############################################################
  ## the dumb  case where you only have the intercept
  df <- data.frame(x = c(1, 0, 1, 0),
                   y = c(-3, -11, 2, 4))
  inf_params <- c("intercept" = -2)

  out <- calculate_transmission_prob(data = df,
                                     inf_params = inf_params)

  exp_out <- 1 / (1 + exp(2 - c(0, 0, 0, 0)))
  expect_equal(out, exp_out)

  ## ###############################################3
  ## x and y
  df <- data.frame(x = c(1, 0, 1, 0),
                   y = c(-3, -11, 2, 4))
  inf_params <- c(-2, 1, -.1)
  cov_names <- c("x", "y")

  out <- calculate_transmission_prob(data = df,
                                     inf_params = inf_params,
                                     covariate_names = cov_names)

  exp_out <- 1 / (1 + exp(2 - c(1, 0, 1, 0) +
                            .1 * c(-3, -11, 2, 4)))
  expect_equal(out, exp_out)
})
