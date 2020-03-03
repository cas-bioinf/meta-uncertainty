#' @export
distance_scaling <- function(base_points) {
  check_point_matrix(base_points)
  centroid <- base_points %>% colMeans()
  centroid_matrix <- matrix(rep(centroid, each = nrow(base_points)), ncol = 2)
  # The denominator akin sd / DoF
  sum(sqrt(rowSums(((base_points - centroid_matrix) ^ 2)))) / (nrow(base_points) - 1)
}

aggregate_consistencies <- function(consistencies) {
  #TODO choose
  #1 - sqrt(mean((1 - consistencies) ^ 2))
  1 - sqrt(sum((1 - consistencies) ^ 2) / (length(consistencies) - 1 ))
}

transform_mse_to_consistency <- function(mean_squared_error, scale_factor) {
  #TODO choose
  pmax(0, 1 - sqrt(mean_squared_error) / scale_factor)
  #exp( - exp(1) * (mean_squared_error / scale_factor^2))
  #sqrt(mean_squared_error) / scale_factor
}

#' @export
consistency_location_per_point <- function(base_points, aligned_bootstrap_points) {
#distance / max distance within original mds
  check_point_matrix(base_points)
  n_points <- nrow(base_points)
  n_draws <- length(aligned_bootstrap_points)

  #max_distance_squared <-  max(dist(base_points, method = "euclidean"))^2
  #mean_distance_squared <-  mean(dist(base_points, method = "euclidean"))^2
  scale_factor <- distance_scaling(base_points)

  sum_squared_distances <- array(0, n_points)
  n_squared_distances <- array(0, n_points)
  names(sum_squared_distances) <- rownames(base_points)
  for(i in 1:n_draws) {
    check_point_matrix(aligned_bootstrap_points[[i]])
    if(n_points == nrow(aligned_bootstrap_points[[i]])) {
      if(!is.null(rownames(base_points)) && !is.null(rownames(aligned_bootstrap_points[[i]])) &&
         !all(rownames(base_points) == rownames(aligned_bootstrap_points[[i]]))) {
        stop("Row names are inconsistent between base_points and aligned_bootstrap_points")
      }
      points_present <- 1:n_points
      points_present_base_id <- 1:n_points
    } else {
      points_present = rownames(aligned_bootstrap_points[[i]])
      if(!identical(intersect(points_present, rownames(base_points)), points_present)) {
        stop("Some points in aligned_bootstrap_points not found in base_points")
      }
      points_present_base_id <- integer(length(points_present))
      for(p in 1:length(points_present)) {
        points_present_base_id[p] <- which(rownames(base_points) == points_present[p])
      }

    }
    squared_distances <- rowSums((base_points[points_present_base_id,] - aligned_bootstrap_points[[i]]) ^ 2)
    sum_squared_distances[points_present_base_id] <- sum_squared_distances[points_present_base_id] + squared_distances
    n_squared_distances[points_present_base_id] <- n_squared_distances[points_present_base_id] + 1
  }
  transform_mse_to_consistency(sum_squared_distances / n_squared_distances , scale_factor)
}

#' @export
consistency_location <- function(base_points, aligned_bootstrap_points) {
  aggregate_consistencies(consistency_location_per_point(base_points, aligned_bootstrap_points))
}

#' @export
consistency_distances_per_point <- function(base_points, aligned_bootstrap_points) {
  check_point_matrix(base_points)
  n_points <- nrow(base_points)

  n_draws <- length(aligned_bootstrap_points)

  scale_factor <- distance_scaling(base_points)
  base_dist <- dist(base_points, method = "euclidean")
  base_distances <- as.matrix(base_dist)

  sum_squared_differences <- array(0, n_points)
  n_squared_differences <- array(0, n_points)
  names(sum_squared_differences) <- rownames(base_points)

  for(i in 1:n_draws) {
    check_point_matrix(aligned_bootstrap_points[[i]])
    if(n_points == nrow(aligned_bootstrap_points[[i]])) {
      if(!is.null(rownames(base_points)) && !is.null(rownames(aligned_bootstrap_points[[i]])) &&
         !all(rownames(base_points) == rownames(aligned_bootstrap_points[[i]]))) {
        stop("Row names are inconsistent between base_points and aligned_bootstrap_points")
      }
      points_present <- 1:n_points
      points_present_base_id <- 1:n_points
    } else {
      points_present = rownames(aligned_bootstrap_points[[i]])
      if(!identical(intersect(points_present, rownames(base_points)), points_present)) {
        stop("Some points in aligned_bootstrap_points not found in base_points")
      }
      points_present_base_id <- integer(length(points_present))
      for(p in 1:length(points_present)) {
        points_present_base_id[p] <- which(rownames(base_points) == points_present[p])
      }
    }
    draw_distances <- as.matrix(dist(aligned_bootstrap_points[[i]], method = "euclidean"))
    n_draw_points <- length(points_present)
    n_angles_per_point <- (n_draw_points - 1) * (n_draw_points - 2) / 2

    diff_matrix <- draw_distances - base_distances[points_present_base_id, points_present_base_id]
    sum_squared_differences[points_present_base_id] <-
      sum_squared_differences[points_present_base_id] + rowSums(diff_matrix ^ 2)

    n_squared_differences[points_present_base_id] <-
      #n_squared_differences[points_present_base_id] + n_points - 1
      #TODO Attempting something like degrees of freedom
      n_squared_differences[points_present_base_id] + n_points - 2
  }
  transform_mse_to_consistency(sum_squared_differences / n_squared_differences, scale_factor)
}

consistency_distances <- function(base_points, aligned_bootstrap_points) {
  aggregate_consistencies(consistency_distances_per_point(base_points, aligned_bootstrap_points))
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

  # Numerical inaccuracies can get me out of acos bounds
  normalized_dot_product[normalized_dot_product < 0 & normalized_dot_product > -1e-10] <- 0
  normalized_dot_product[normalized_dot_product > 1 & normalized_dot_product < 1 + 1e-10] <- 1

  ret <- acos(normalized_dot_product)

  if(any(is.na(ret))) {
    print(v1[is.na(ret)])
    print(v2[is.na(ret)])
    print(norm[is.na(ret)])
    print(normalized_dot_product[is.na(ret)])
    stop("NA angle")
  }

  ret
}

get_n_angles_per_point <- function(n_points) {
  (n_points - 1) * (n_points - 2) / 2
}

get_angles_per_point <- function(all_angles, n_points, point_id) {
  n_angles_per_point <- get_n_angles_per_point(n_points)
  if(length(all_angles) != n_points * n_angles_per_point) {
    stop("Inconsistent arguments for get_angles_per_point")
  }
  all_angles[ ((point_id - 1) * n_angles_per_point + 1) : (point_id * n_angles_per_point)]
}

#' @export
consistency_angles_per_point <- function(base_points, aligned_bootstrap_points) {
  check_point_matrix(base_points)
  n_points <- nrow(base_points)

  n_draws <- length(aligned_bootstrap_points)


  sum_squared_distances <- array(0, n_points)
  n_squared_distances <- array(0, n_points)
  names(sum_squared_distances) <- rownames(base_points)

  base_angles_all <- compute_angles(base_points)


  for(i in 1:n_draws) {
    check_point_matrix(aligned_bootstrap_points[[i]])
    if(n_points == nrow(aligned_bootstrap_points[[i]])) {
      if(!is.null(rownames(base_points)) && !is.null(rownames(aligned_bootstrap_points[[i]])) &&
         !all(rownames(base_points) == rownames(aligned_bootstrap_points[[i]]))) {
        stop("Row names are inconsistent between base_points and aligned_bootstrap_points")
      }
      points_present <- 1:n_points
      points_present_base_id <- 1:n_points
      matched_base_angles <- base_angles_all
    } else {
      points_present = rownames(aligned_bootstrap_points[[i]])
      if(!identical(intersect(points_present, rownames(base_points)), points_present)) {
        stop("Some points in aligned_bootstrap_points not found in base_points")
      }
      points_present_base_id <- integer(length(points_present))
      for(p in 1:length(points_present)) {
        points_present_base_id[p] <- which(rownames(base_points) == points_present[p])
      }

      matched_base_angles <- compute_angles(base_points[points_present_base_id,])
    }

    draw_angles <- compute_angles(aligned_bootstrap_points[[i]])
    n_draw_points <- length(points_present)
    n_angles_per_point <- get_n_angles_per_point(n_draw_points)

    if(any(is.na(draw_angles))){
      stop(paste0("NA angle in draw ", i))
    }

    for(p in 1:n_draw_points) {
      draw_point_angles <- get_angles_per_point(draw_angles, n_draw_points, p)
      base_point_angles <- get_angles_per_point(matched_base_angles, n_draw_points, p)
      # Note that since the angles are computed in [0, pi] I don't need to worry about periodicity 2pi ~ 0
      sum_squared_distances[points_present_base_id[p]] <- sum_squared_distances[points_present_base_id[p]] + sum((draw_point_angles - base_point_angles) ^ 2)
      #n_squared_distances[points_present_base_id[p]] <- n_squared_distances[points_present_base_id[p]] + n_angles_per_point
      #TODO Attempting something like "degrees of freedom"
      n_squared_distances[points_present_base_id[p]] <- n_squared_distances[points_present_base_id[p]] + n_angles_per_point - n_points + 2
    }
  }
  transform_mse_to_consistency(sum_squared_distances / n_squared_distances, 0.5 * pi)
}

#' @export
consistency_angles <- function(base_points, aligned_bootstrap_points) {
  aggregate_consistencies(consistency_angles_per_point(base_points, aligned_bootstrap_points))
}




