test_that("consistency_location", {
  tol = 1e-14
  base_points <- cbind(c(0, 1, 0),  c(0, 0, 1))
  expect_equivalent(consistency_location_per_point(base_points, list(base_points, base_points)),
                    c(1,1,1), tolerance = tol)

  # I previously computed the test case for max distance, but now use distance scaling
  scale_factor <- distance_scaling(base_points)
  max_distance_squared <- max(dist(base_points, method = "euclidean")) ^ 2
  expect_equivalent(consistency_location_per_point(base_points, list(base_points * 0.5, base_points * 0.3)),
                    1 - sqrt(c(0, 0.5 ^ 2 + 0.7^2,0.5^2 + 0.7^2) * max_distance_squared / (2 * 2 * scale_factor)), tolerance = tol)

})


test_that("consistency_location_missing", {
  tol = 1e-14
  base_points <- cbind(c(0, 1, 0),  c(0, 0, 1))
  rownames(base_points) <- as.character(1:3)

  expect_equivalent(consistency_location_per_point(base_points, list(base_points[c(1,2),], base_points[c(1,3),], base_points[c(2,3),])),
                    c(1,1,1), tolerance = tol)


})
