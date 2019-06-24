
test_that("connectivity for single group", {
  # Each point has one sample close and two samples far
  #                         point1        point2       point3                           
  base_points <-  matrix( c(0,     1,     0,     0,    1,     0    ), nrow = 3, ncol = 2, byrow = TRUE)
  samples <- list(
                  matrix( c(0.1,   1.2,  -1.8,   0.3,  0.5,   0.9  ), nrow = 3, ncol = 2, byrow = TRUE),
                  matrix( c(0.9,   1.1,   0.2,   0.18, 1.4,   -1.1 ), nrow = 3, ncol = 2, byrow = TRUE),
                  matrix( c(-2.9,  0.17,  0.9,   1.04, 1.25,  -0.1 ), nrow = 3, ncol = 2, byrow = TRUE)
  )
  
  samples_further_asymmetric <- matrix(c(3, 1, 1,    2, 3, 2,   0, 2, 3), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equal( connectivity_direct_assymmetric(base_points, samples), samples_further_asymmetric / length(samples))
  
  samples_further <- matrix(c(3, 2, 1,    2, 3, 2,   1, 2, 3), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equal( connectivity_direct(base_points, samples), samples_further / length(samples))
  
  expect_equal( connectivity_group(base_points, samples), 2/3)
})
