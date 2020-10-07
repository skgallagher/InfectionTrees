
# test_that("real data with mc_trees_A_infects_B", {
#
#     data(tb_clean)
#     clusters <- tb_clean %>%
#         dplyr::mutate(smear = ifelse(spsmear == "Positive",
#                                      1, 0),
#                       cluster_id = group,
#                       hiv_f = hivstatus) %>%
#         dplyr::mutate(hiv_neg_pos = ifelse(hiv_f == "neg", 1, 0),
#                hiv_unk_pos = ifelse(hiv_f == "unk", 1, 0)) %>%
#         dplyr::group_by(cluster_id) %>%
#         dplyr::mutate(rel_time = rel_time / 365) %>%
#         dplyr::mutate(cluster_size = dplyr::n()) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(race_f = forcats::fct_collapse(race,
#                                      white = "White",
#                                      black = "Black or African American",
#                                      asian = "Asian")) %>%
#         dplyr::mutate(race_asian_white = ifelse(race_f == "asian", 1, 0),
#                race_black_white = ifelse(race_f == "black", 1, 0)) %>%
#         dplyr::select(cluster_id, smear,
#                hiv_neg_pos,
#                hiv_unk_pos,
#                rel_time,
#                race_asian_white,
#                race_black_white,
#                cluster_size)
#
#     par_ests <- c("Intercept" = -.72690792,
#               "smear" = -.09322607,
#               "hiv_neg_pos" = -0.36415841,
#               "hiv_unk_pos" = -0.56817810,
#               "rel_time" = .345552628)
#
#
#     set.seed(42)
#     mc_trees <-  sample_mc_trees(clusters %>%
#                                  dplyr::filter(cluster_id == 27),
#                                  B = 100,
#                                  multiple_outside_transmissions = FALSE,
#                                  covariate_names = names(par_ests))
#
#
#     plotting_df <- mc_trees_to_A_infects_B(my_id = 27,
#                                    par_ests = par_ests,
#                                    mc_trees = mc_trees,
#                                    output = "tidy")
#     expect_true(sum(plotting_df$prob == 0) <
#                 sum(plotting_df$prob != 0))
#
#
#     adj_mat <- mc_trees_to_A_infects_B(my_id = 27,
#                                    par_ests = par_ests,
#                                    mc_trees = mc_trees,
#                                    output = "matrix")
#     expect_true(all(diag(adj_mat) == 0))
#
# })


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
