test_that("bootstrap_clusters", {
  data <- data.frame(cluster_id = c("a", "a", "a",
                                    "b", "b",
                                    "c", "c",
                                    "d"),
                     x = c(1, 0, 1,
                           1, 1,
                           0, 1,
                           0))

  boot_clusters <- bootstrap_clusters(clusters = data)
  expect_equal(length(unique(boot_clusters$cluster_id)),
               length(unique(data$cluster_id)))
  expect_true(all(boot_clusters$original_id %in%
                    data$cluster_id))
})
