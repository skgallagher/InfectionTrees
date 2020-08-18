test_that("summarize_binary_clusters", {
  example_cluster <- data.frame(cluster_id = c(1, 1, 1,
                                                2, 2,
                                                3, 3, 3, 3,
                                               4),
                                               x = c(0, 1, 1,
                                               0, 0,
                                               1, 0, 1, 1,
                                               0))
  out <- summarize_binary_clusters(example_cluster)
  expect_equal(out$freq, c(1, 1, 1, 1))
  expect_equal(out$cluster_size, 1:4)
  expect_equal(out$x_pos, c(0, 0, 2, 3))
  expect_equal(out$x_neg, out$cluster_size - out$x_pos)

  ########################################################################
  example_cluster <- data.frame(cluster_id = c(1, 1, 1,
                                               2, 2,
                                               3, 3, 3,
                                               4),
                                x = c(0, 1, 1,
                                      0, 0,
                                      1, 0, 1,
                                      0))
  out <- summarize_binary_clusters(example_cluster)
  expect_equal(out$freq, c(1, 1, 2))
  expect_equal(out$cluster_size, 1:3)
  expect_equal(out$x_pos, c(0, 0, 2))
  expect_equal(out$x_neg, out$cluster_size - out$x_pos)
  expect_true("binary_cov" %in% class(out))

})
