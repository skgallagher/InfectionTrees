test_that("bp_loglike_binary_cov", {

  inf_params <- c(-1, 1)
  mc_samples_summary <- data.frame(freq = c(10, 5, 6, 1, 1),
                                   mc_freq = c(10, 10, 10, 3, 7),
                                   cluster_size = c(1, 1, 2, 2, 2),
                                   x_pos = c(1, 0, 2, 1, 1),
                                   x_neg = c(0, 1, 0, 1, 1),
                                   x_pos_trans = c(0, 0, 1, 1, 0),
                                   x_neg_trans = c(0, 0, 0, 0, 1))
  p_neg <- 1 / (1 + exp(1))
  p_pos <- .5
  prob_inf <- with(mc_samples_summary,
                   (1-p_pos)^x_pos * (1-p_neg)^x_neg *
                     p_pos^x_pos_trans * p_neg^x_neg_trans)

  expect_equal(prob_inf, c(.5,
                 (1-p_neg),
                 (1-p_pos)^2 * p_pos^1,
                 (1-p_pos) * (1-p_neg) * p_pos,
                 (1-p_pos) * (1-p_neg) * p_neg))


  avg_like <- c(  prob_inf[1],
                 prob_inf[2],
                 prob_inf[3],
                 (3 * prob_inf[4] + 7 * prob_inf[5]) /10 )
  out <- bp_loglike_binary_cov(inf_params = inf_params,
                               mc_samples_summary = mc_samples_summary,
                               return_neg = FALSE,
                               multiple_outsider_transmission = FALSE)
  expect_equal(sum(c(10, 5, 6, 2) * log(avg_like)), out)
})
