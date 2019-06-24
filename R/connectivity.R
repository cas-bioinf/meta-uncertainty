#' Asymetric direct connectivity between pair of observations: 
#'     fraction of resamples of source further from source than the target
#' Value: array with connectivy from `x` to `y` as `result[x,y]`
#' @export
connectivity_direct_assymmetric <- function(base_points, aligned_samples_points) {
  base_distances <-  dist(base_points, method = "euclidean") %>% as.matrix()
  n_points <- nrow(base_points)
  n_samples <- length(aligned_samples_points)
  n_samples_far <- matrix(0, n_points, n_points)
  for(i in 1:n_samples) {
    dist_to_samples <- sqrt(rowSums((base_points - aligned_samples_points[[i]]) ^ 2))
    for(p in 1:n_points) {
      is_further <- dist_to_samples >= base_distances[,p]
      n_samples_far[,p] <- n_samples_far[,p] + is_further
    }
  }
  n_samples_far / n_samples
}

#' Direct connectivity between pair: maximum of asymetric direct connectivity
#' @export
connectivity_direct <- function(base_points, aligned_samples_points) {
  assymetric <- connectivity_direct_assymmetric(base_points, aligned_samples_points)
  pmax(assymetric, t(assymetric))
}


#' Connectivity betweeen pair: width of the widest path in the graph of direct connectivity
#' Connectivity for a group: minimum of connectivity between pairs across the group
#' 
#' Since widest path between any two poitns has to be completely included in the maximum spanning tree, 
#' connectivity for a group can be found as
#' the minimum edge on the maximum spanning tree.
#' 
#' @export
connectivity_group <- function(base_points, aligned_samples_points) {
  direct <- connectivity_direct(base_points, aligned_samples_points)
  
  #I am looking for maximum spanning tree, but spantree searches for minimum, hence the `1 - x` transform
  spanning_tree <- vegan::spantree(1 - direct) 
 
  min(1 - spanning_tree$dist) 
}

connectivity_all_groups <- function(base_points, aligned_samples_points, groups) {
  groups %>% unique() %>%
    purrr::map(function(g) {
      filtered_base <- base_points[groups == g, ]
      filtered_samples <- aligned_samples_points %>% map( ~ .x[groups == g, ])
      data.frame(group = g, connectivity = connectivity_group(filtered_base, filtered_samples))
    }) %>%
  do.call(rbind, .)
}