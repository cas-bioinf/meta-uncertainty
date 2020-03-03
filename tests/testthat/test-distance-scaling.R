test_that("distnace_scaling", {
  points_x <- c(0, 1, 0, -1)
  points_y <- c(1, 0, -1, 0)
  points <- cbind(points_x,points_y)
  expect_equivalent(distance_scaling(points), 1)

  expect_equivalent(distance_scaling(points * 3), 3)
  expect_equivalent(distance_scaling(points * -2), 2)
  expect_equivalent(distance_scaling(points + 5), 1)

  rand_points <- matrix(rnorm(16), ncol = 2)
  expect_equivalent(distance_scaling(rand_points), distance_scaling(rand_points + 1))
  expect_equivalent(distance_scaling(rand_points), 0.5 * distance_scaling(rand_points * 2))
})
