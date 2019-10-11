#' Asymetric direct connectivity between pair of observations:
#'     fraction of resamples of source further from source than the target
#' Value: array with connectivy from `x` to `y` as `result[x,y]`
#' @export
connectivity_matrix_direct_assymmetric <- function(base_points, aligned_samples_points) {
  base_distances <-  dist(base_points, method = "euclidean") %>% as.matrix()
  n_points <- nrow(base_points)
  n_samples <- length(aligned_samples_points)
  n_samples_far <- matrix(0, n_points, n_points)
  rownames(n_samples_far) <- rownames(base_points)
  colnames(n_samples_far) <- rownames(base_points)
  for(i in 1:n_samples) {
    if(nrow(base_points) == nrow(aligned_samples_points[[i]])) {
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
    dist_to_samples <- sqrt(rowSums((base_points[points_present,] - aligned_samples_points[[i]]) ^ 2))
    for(p in points_present) {
      is_further <- dist_to_samples >= base_distances[points_present, p]
      n_samples_far[points_present,p] <- n_samples_far[points_present, p] + is_further
    }
  }
  n_samples_far / n_samples
}

#' Direct connectivity between pair: maximum of asymetric direct connectivity
#' @export
connectivity_matrix_direct <- function(base_points, aligned_samples_points) {
  assymetric <- connectivity_matrix_direct_assymmetric(base_points, aligned_samples_points)
  pmax(assymetric, t(assymetric))
}

pairwise_widest_path_distances_from_matrix <- function(direct) {
  n_points <- nrow(direct)

  # Modified Floyd-Warshall for widest path
  # Probably suboptimal, since I could use maximal spanning tree for this, but easiest to implement

  # Init with edges
  widest_matrix <- direct

  # 2) Main loop
  for(k in 1:n_points) {
    for(i in 1:n_points) {
      for(j in 1:n_points) {
        candidate_width <- min(widest_matrix[i,k], widest_matrix[k, j])
        widest_matrix[i, j] <- max(widest_matrix[i, j], candidate_width)
      }
    }
  }

  widest_matrix
}

#' Connectivity betweeen pair: width of the widest path in the graph of direct connectivity
#' Connectivity for a group: minimum of connectivity between pairs across the group
#'
#' Since widest path between any two poitns has to be completely included in the maximum spanning tree,
#' connectivity for a group can be found as
#' the minimum edge on the maximum spanning tree.
#'
#' @export
connectivity_stats_group <- function(base_points, aligned_samples_points) {
  direct <- connectivity_matrix_direct(base_points, aligned_samples_points)

  #I am looking for maximum spanning tree, but spantree searches for minimum, hence the `1 - x` transform
  spanning_tree <- vegan::spantree(1 - direct)

  inverted_dist <- 1 - spanning_tree$dist
  min_index <- which.min(inverted_dist)
  connectivity_min = inverted_dist[min_index]

  min_parent <- spanning_tree$labels[min_index + 1]
  min_kid <- spanning_tree$labels[spanning_tree$kid[min_index]]

  all_distances <- pairwise_widest_path_distances_from_matrix(direct)

  connectivity_average <- mean(sqrt(all_distances[lower.tri(all_distances)])) ^ 2

  adjacency_matrix <- direct > 0
  n_components <- igraph::components(igraph::graph.adjacency(adjacency_matrix))$no



  data.frame(n_components,
             connectivity_average,
             connectivity_min, bottleneck_a = min_parent, bottleneck_b = min_kid)
}

#' @importFrom tidyr %>%
#' @export
connectivity_stats_all_groups <- function(base_points, aligned_samples_points, groups) {
  get_group_indices <- function(aligned_sample, group) {
    observation_names <- intersect(rownames(base_points)[groups == group], rownames(aligned_sample))
    aligned_sample[observation_names, ]
  }

  groups %>% unique() %>%
    purrr::map(function(g) {
      filtered_base <- base_points[groups == g, ]
      filtered_samples <- aligned_samples_points %>% purrr::map(get_group_indices, group = g)
      result <- connectivity_stats_group(filtered_base, filtered_samples) %>%
        dplyr::mutate(group = g) %>% dplyr::select(group, dplyr::everything())
      result$group <- g
      result
    }) %>%
  do.call(rbind, .)
}
