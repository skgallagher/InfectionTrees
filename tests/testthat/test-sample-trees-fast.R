

test_that("sample_connections", {
  gen_sizes <- c(1)
  out <- sample_connections(gen_sizes)
  expect_true(is.na(out$inf_id))
  #############
  ######
  gen_sizes <- c(1, 2, 2)
  out <- sample_connections(gen_sizes)
  expect_equal(sum(is.na(out$inf_id)), 1)
  inf_id_g <- as.integer(substr(out$inf_id, 1, 1))
  expect_true(all(is.na(inf_id_g) | inf_id_g == (out$gen - 1)))
  ######
  gen_sizes <- c(1, 2, 1)
  out <- sample_connections(gen_sizes)
  expect_equal(sum(is.na(out$inf_id)), 1)
  inf_id_g <- as.integer(substr(out$inf_id, 1, 1))
  expect_true(all(is.na(inf_id_g) | inf_id_g == (out$gen - 1)))
})




test_that("sample_unique_perms", {
  g <- 1
  n <- 1
  B <- 2
  out <- sample_unique_perms(g, n, B)
  expect_equal(as.numeric(out), rep(1, B))
  ################
  g <- 2
  n <- 2
  B <- 2
  out <- sample_unique_perms(g, n, B)
  expect_equal(dim(out), c(g, B))
  expect_true(all(out[1,] == 1))
  expect_true(all(colSums(out) == n))
  ################
  g <- 2
  n <- 2
  B <- 2
  out <- sample_unique_perms(g, n, B)
  expect_equal(dim(out), c(g, B))
  expect_true(all(out[1,] == 1))
  expect_true(all(colSums(out) == n))
  ################
  g <- 4
  n <- 10
  B <- 100
  out <- sample_unique_perms(g, n, B)
  expect_equal(dim(out), c(g, B))
  expect_true(all(out[1,] == 1))
  expect_true(all(colSums(out) == n))
})



