test_that("consistency_location", {
  tol = 1e-14
  base_points <- cbind(c(0, 1, 0),  c(0, 0, 1))
  expect_equivalent(consistency_location_per_point(base_points, list(base_points, base_points)),
                    c(1,1,1), tolerance = tol)

  expect_equivalent(consistency_location_per_point(base_points, list(base_points * 0.5, base_points * 0.3)),
                    1 - sqrt(c(0, 0.5 ^ 2 + 0.7^2,0.5^2 + 0.7^2) / (2 * 2)), tolerance = tol)

})
