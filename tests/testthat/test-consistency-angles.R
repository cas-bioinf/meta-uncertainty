
test_that("compute_angles", {
  points_x <- c(0, 1, 1, 0)
  points_y <- c(0, 0, 1, 1)
  points <- cbind(points_x,points_y)
  expect_error(compute_angles(points[1:2,], " < 3 points"))

  tol <- 1e-14
  expect_equivalent(compute_angles(points[1:3,]), c(pi / 4, pi / 2, pi / 4), tolerance = tol)
  angles_4 <- compute_angles(points)
  expect_true(all(!is.na(angles_4)))
  for(i in 1:4) {
    slice <- ((i - 1) * 3 + 1) : (i * 3)
    angles_subset <- angles_4[slice]
    expect_equivalent(sort(angles_subset), c(pi / 4, pi / 4, pi / 2), tolerance = tol)
  }
  angles_4_2 <- compute_angles(points + 5)
  expect_equivalent(angles_4, angles_4_2, tolerance = tol)

  points_hexagon_x <- c( -1, -0.5, 0.5, 1, 0.5, -0.5)
  points_hexagon_y <- c( 0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2)
  points_hexagon <- cbind(points_hexagon_x, points_hexagon_y)
  angles_hexagon <- compute_angles(points_hexagon)
  expect_true(all(!is.na(angles_hexagon)))
  for(i in 1:6) {
    slice <- ((i - 1) * 10 + 1) : (i * 10)
    angles_subset <- angles_hexagon[slice]
    expect_equivalent(sort(angles_subset), pi * c(1/6, 1/6, 1/6, 1/6, 1/3, 1/3, 1/3, 1/2, 1/2, 2/3 ), tolerance = tol)
  }

  #Degenerate angles should return 0
  expect_equivalent(compute_angles(cbind( c(1,0,0),c(1,1,1))), c(0,0,0),  tolerance = tol)
})

test_that("consistency_angles", {
  tol = 1e-14
  base_points <- cbind(c(0, 1, 0),  c(0, 0, 1))
  expect_equivalent(consistency_angles_per_point(base_points, list(base_points, base_points + 3, base_points * 5 + 1)), c(1,1,1), tolerance = tol)

  shifted_points <- cbind(c(0, 1, 0.5), c(0, 0, sqrt(3) / 2))
  expect_equivalent(consistency_angles_per_point(base_points, list(shifted_points)), 1 - sqrt(c(pi/2 - pi/3, pi/3 - pi/4,pi/3 - pi/4)^2 / pi^2), tolerance = tol)

  points_5 <- cbind(rep(0,5), c(0,1,0,1,0))
  points_shifted_5 <- points_5
  points_shifted_5[1,] <- c(1,0)
  change_1 <- 4 * (1/4)^2
  change_even <- 2* (1/4)^2
  change_odd <- 2 * (1/2)^2
  expect_equivalent(consistency_angles_per_point(points_5, list(points_shifted_5)),
                     1 - sqrt((c(change_1, change_even, change_odd, change_even, change_odd)) / 6), tolerance = tol)

  # points_hexagon_x <- c( -1, -0.5, 0.5, 1, 0.5, -0.5)
  # points_hexagon_y <- c( 0, sqrt(3)/2, sqrt(3)/2, 0, -sqrt(3)/2, -sqrt(3)/2)
  # points_hexagon <- cbind(points_hexagon_x, points_hexagon_y)
  #
  # shifted_points_6 <- points_hexagon
  # shifted_points_6[1,] <- 0
  # expect_equivalent(consistency_angles_per_point(points_hexagon, list(shifted_points_6)),
  #                   c(1 - sqrt(sum( (c(0, 1/3,1/3,1/3,1/6,1/6,1/2,1/2,1/2,1/2) * pi)^2) / (10 * pi^2)) , rep(1, 5)), tolerance = tol)
})
