
test_that("Widest path", {
  # Using the example from
  # https://en.wikipedia.org/wiki/Widest_path_problem
  direct <- matrix(0, 7, 7)
  node_names <- c("B","C","D","F","H","M","T")
  rownames(direct) <- node_names
  colnames(direct) <- node_names
  direct["B","D"] <- 15
  direct["B","H"] <- 40
  direct["B","F"] <- 46
  direct["C","H"] <- 17
  direct["C","T"] <- 29
  direct["C","M"] <- 40
  direct["D","H"] <- 53
  direct["F","T"] <- 3
  direct["F","M"] <- 11
  direct["H","T"] <- 31
  direct["M","T"] <- 8

  #Make the matrix symmetrical
  direct[lower.tri(direct)] <-  t(direct)[lower.tri(direct)]

  expected_result <- matrix(c(46, 29, 40, 46, 40, 29, 31,
                              29, 40, 29, 29, 29, 40, 29,
                              40, 29, 53, 40, 53, 29, 31,
                              46, 29, 40, 46, 40, 29, 31,
                              40, 29, 53, 40, 53, 29, 31,
                              29, 40, 29, 29, 29, 40, 29,
                              31, 29, 31, 31, 31, 29, 31), 7, 7)
  rownames(expected_result) <- node_names
  colnames(expected_result) <- node_names

  expect_equal( expected_result, pairwise_widest_path_distances_from_matrix( direct ))
})

test_that("connectivity for single group", {
  # Each point has one sample close and two samples far
  #                         point1        point2       point3
  base_points <-  matrix( c(0,     1,     0,     0,    1,     0    ), nrow = 3, ncol = 2, byrow = TRUE)
  samples <- list(
                  matrix( c(0.1,   1.2,  -1.8,   0.3,  0.5,   0.9  ), nrow = 3, ncol = 2, byrow = TRUE),
                  matrix( c(0.9,   1.1,   0.2,   0.18, 1.4,   -1.1 ), nrow = 3, ncol = 2, byrow = TRUE),
                  matrix( c(-2.9,  0.17,  0.9,   1.04, 1.25,  -0.1 ), nrow = 3, ncol = 2, byrow = TRUE)
  )
  rownames(base_points) <- as.character(1:3)

  samples_further_asymmetric <- matrix(c(3, 1, 1,    2, 3, 2,   0, 2, 3), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equivalent( connectivity_matrix_direct_assymmetric(base_points, samples), samples_further_asymmetric / length(samples))

  samples_further <- matrix(c(3, 2, 1,    2, 3, 2,   1, 2, 3), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equivalent( connectivity_matrix_direct(base_points, samples), samples_further / length(samples))

  stats_group <- connectivity_stats_group(base_points, samples)
  expect_equal( stats_group$connectivity_min, 2/3)
})


test_that("Connectivity works with jackknife samples", {
  set.seed(20191011)
  base_points <-  matrix( rnorm(500), nrow = 50, ncol = 10)
  rownames(base_points) <- as.character(1:50)
  samples <- sample_posterior_jackknife_observations(base_points)
  res <- connectivity_stats_all_groups(base_points, samples, c(1:3, sample(1:3, size = 7, replace = TRUE)))
  expect_true(nrow(res) == 3)
  expect_true(all(!is.na(res$group)))
  expect_true(all(!is.na(res$connectivity_average)))
  expect_true(all(!is.na(res$connectivity_min)))
  expect_true(all(!is.na(res$bottleneck_a)))
  expect_true(all(!is.na(res$bottleneck_b)))
})
