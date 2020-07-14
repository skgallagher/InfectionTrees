## Tests for the outsider general permutations
## Very different from original file
## Don't think I need recursion??

test_that("sample_gen_out_cond", {
    df1 <- data.frame(cluster_id = "a",
                      x = 1:5,
                      y = c(1, 1, 1, 0, 0))
    df2 <- data.frame(cluster_id = "b",
                      x = c(2,3),
                      y = 0)
    observed_data <- rbind(df1, df2)
    observed_data <- observed_data %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::mutate(person_id = paste(cluster_id,
                                        dplyr::row_number(),
                                        sep = "-")) %>%
        dplyr::ungroup()




    B <- 3
    out <- sample_gen_out_cond(observed_data,
                               B = B)
    expect_true(all(as.numeric(table(out$orig_id)) == c((nrow(df1) + 1) * B,
                                                        (nrow(df2) + 1) * B)))
    expect_equal(nrow(out), B * (nrow(df1) + 1 + nrow(df2) + 1))
})

test_that("sample_single_gen_out_cond", {

    observed_data <- data.frame(cluster_id = "a",
                                x = 1:5,
                                y = c(1, 1, 1, 0, 0))
    observed_data <- observed_data %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::mutate(person_id = paste(cluster_id,
                                        dplyr::row_number(),
                      sep = "-")) %>%
        dplyr::ungroup()

    B <- 5
    person_ids <- 1:12

    out <- sample_single_gen_out_cond(observed_data$person_id,
                                      B = B,
                                      person_ids,
                                      orig_id = "a")

    expect_true(all(out$x[1:(nrow(observed_data)) + 1] %in%
                    observed_data$x))
    expect_true(all(table(out$x) >= B))


})
