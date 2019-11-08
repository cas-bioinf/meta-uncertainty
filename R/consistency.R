#' @export
consistency_location_per_point <- function(base_points, aligned_samples_points) {
#distance / max distance within original mds
  check_point_matrix(base_points)
  n_points <- nrow(base_points)
  n_samples <- length(aligned_samples_points)

  max_distance_squared <-  max(dist(base_points, method = "euclidean"))^2

  sum_squared_distances <- array(0, n_points)
  n_squared_distances <- array(0, n_points)
  names(sum_squared_distances) <- rownames(base_points)

  for(i in 1:n_samples) {
    check_point_matrix(aligned_samples_points[[i]])
    if(n_points == nrow(aligned_samples_points[[i]])) {
      if(!is.null(rownames(base_points)) && !is.null(rownames(aligned_samples_points[[i]])) &&
         !all(rownames(base_points) == rownames(aligned_samples_points[[i]]))) {
        stop("Row names are inconsistent between base_points and aligned_samples_points")
      }
      points_present <- 1:n_points
    } else {
      points_present = rownames(aligned_samples_points[[i]])
      if(!identical(intersect(points_present, rownames(base_points)), points_present)) {
        stop("Some points in aligned_samples_points not found in base_points")
      }
    }
    squared_distances <- rowSums((base_points[points_present,] - aligned_samples_points[[i]]) ^ 2)
    sum_squared_distances[points_present] <- sum_squared_distances[points_present] + squared_distances
    n_squared_distances[points_present] <- n_squared_distances[points_present] + 1
  }
  1 - sqrt(sum_squared_distances / (n_squared_distances * max_distance_squared))
}

#' @export
consistency_location <- function(base_points, aligned_samples_points) {
  mean(consistency_location_per_point(base_points, aligned_samples_points))
}

# All angles between threes of points, in radians
#' @export
compute_angles <- function(points) {
  check_point_matrix(points)
  n_points <- nrow(points)
  if(n_points < 3) {
    stop("Can't compute angles with < 3 points")
  }
  n_angles <- (n_points * (n_points - 1) * (n_points - 2)) / 2
  centers <- array(NA_integer_, n_angles)
  edges1 <- array(NA_integer_, n_angles)
  edges2 <- array(NA_integer_, n_angles)

  next_angle <- 1
  for(center in 1:n_points) {
    for(edge1_raw in 2:(n_points - 1)) {
      if(edge1_raw >= center) {
        edge1 <- edge1_raw + 1
      }
      else {
        edge1 <- edge1_raw
      }

      segment_end <- next_angle + edge1_raw - 2
      centers[next_angle:segment_end] <- center
      edges1[next_angle:segment_end] <- edge1

      if(center == 1) {
        edges2[next_angle:segment_end] <- 2:edge1_raw
      } else if(center > edge1_raw - 1) {
        edges2[next_angle:segment_end] <- 1:(edge1_raw - 1)
      } else {
        edges2[next_angle:(next_angle + center - 2)] <- 1:(center-1)
        edges2[(next_angle + center - 1):segment_end] <- (center + 1):(edge1 - 1)
      }

      next_angle <- segment_end + 1
    }
  }

  v1 <- points[edges1,] - points[centers,]
  v2 <- points[edges2,] - points[centers,]
  norm <- sqrt(rowSums(v1 * v1) * rowSums(v2 * v2))
  normalized_dot_product <- ifelse(norm == 0, 1, rowSums(v1 * v2) / norm)


  acos(normalized_dot_product)
}

#' @export
consistency_angles_per_point <- function(base_points, aligned_samples_points) {
  check_point_matrix(base_points)
  n_points <- nrow(base_points)

  n_samples <- length(aligned_samples_points)

  base_angles <- compute_angles(base_points)

  sum_squared_distances <- array(0, n_points)
  n_squared_distances <- array(0, n_points)
  names(sum_squared_distances) <- rownames(base_points)

  for(i in 1:n_samples) {
    check_point_matrix(aligned_samples_points[[i]])
    if(n_points == nrow(aligned_samples_points[[i]])) {
      if(!is.null(rownames(base_points)) && !is.null(rownames(aligned_samples_points[[i]])) &&
         !all(rownames(base_points) == rownames(aligned_samples_points[[i]]))) {
        stop("Row names are inconsistent between base_points and aligned_samples_points")
      }
      points_present <- 1:n_points
      points_present_base_id <- 1:n_points
    } else {
      points_present = rownames(aligned_samples_points[[i]])
      if(!identical(intersect(points_present, rownames(base_points)), points_present)) {
        stop("Some points in aligned_samples_points not found in base_points")
      }
      points_present_base_id <- integer(length(points_present))
      for(p in 1:length(points_present)) {
        points_present_base_id[p] <- which(rownames(base_points) == points_present[p])
      }
    }
    sample_angles <- compute_angles(aligned_samples_points[[i]])
    n_sample_points <- length(points_present)
    n_angles_per_point <- (n_sample_points - 1) * (n_sample_points - 2) / 2

    for(p in 1:n_sample_points) {
      sample_point_angles <- sample_angles[ ((p - 1) * n_angles_per_point + 1) : (p * n_angles_per_point)]
      base_point_angles <- base_angles[ ((points_present_base_id[p] - 1) * n_angles_per_point + 1) : (points_present_base_id[p] * n_angles_per_point)]
      # print(sample_point_angles)
      # print(base_point_angles)
      sum_squared_distances[points_present_base_id[p]] <- sum_squared_distances[points_present_base_id[p]] + sum((sample_point_angles - base_point_angles) ^ 2)
      n_squared_distances[points_present_base_id[p]] <- n_squared_distances[points_present_base_id[p]] + n_angles_per_point
    }
  }
  1 - sqrt(sum_squared_distances / (n_squared_distances * pi^2))
}

#' @export
consistency_angles <- function(base_points, aligned_samples_points) {
  mean(consistency_angles_per_point(base_points, aligned_samples_points))
}
