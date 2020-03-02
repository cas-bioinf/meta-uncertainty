metaMDS_per_draw <- function(bootstrap, cores = parallel::detectCores(), ...) {
  cl <- parallel::makePSOCKcluster(cores);
  currentLibPaths <- .libPaths()
  tryCatch({
    parallel::clusterExport(cl, "currentLibPaths", environment())
    parallel::clusterEvalQ(cl, .libPaths(currentLibPaths))
    parallel::clusterExport(cl, "bootstrap", environment())
    parallel::parLapplyLB(cl, 1:get_n_bootstrap_draws(bootstrap),
                          fun = function(i) {vegan::metaMDS(metagenboot::get_bootstrap_draw(bootstrap,i), parallel = 1, ...)},
                          chunk.size = 1)
  }, finally = {
    parallel::stopCluster(cl)
  })
}

align_single_MDS <- function(mds_to_align, base_mds, allow_unaligned = FALSE) {
  if(!identical(rownames(vegan::scores(base_mds)), rownames(vegan::scores(mds_to_align)))) {
    if(!allow_unaligned) {
      stop("Observations are not aligned")
    }
    all_obs_base <- rownames(vegan::scores(base_mds))
    all_obs_to_align <- rownames(vegan::scores(mds_to_align))
    selected_observations <- all_obs_base[all_obs_base %in% all_obs_to_align]
    vegan::procrustes(X = vegan::scores(base_mds)[selected_observations,], Y = vegan::scores(mds_to_align)[selected_observations,])
  } else {
    vegan::procrustes(X = base_mds, Y = mds_to_align)
  }
}


bootstrap_func_dm <- function(N_draws, N_reads = "original", prior = default_dm_prior) {
  function(observed_matrix) {
    bootstrap_reads_dm(N = N_draws, observed = observed_matrix, N_reads = N_reads, prior = prior)
  }
}

mds_sensitivity_check_compute <- function(observed_matrix, bootstrap_func = bootstrap_func_dm(20),
                                          ...) {
  if(!is.matrix(observed_matrix)) {
    stop("observed_matrix has to be a matrix")
  }

  base_mds <- observed_matrix %>% vegan::metaMDS(...)

  bootstrap_func <- rlang::as_function(bootstrap_func)
  draws <- bootstrap_func(observed_matrix)

  aligned_bootstrap_mds <- draws %>%
    metaMDS_per_draw(...) %>% purrr::map(align_single_MDS, base_mds = base_mds,
                                           allow_unaligned = is_observation_subset(draws))

  list(observed_matrix = observed_matrix,
       base_mds = base_mds,
       aligned_bootstrap_mds = aligned_bootstrap_mds)
}


mds_sensitivity_check_eval <- function(compute_result, mapping, group_column = Group,
                                       id_column = Sample) {
  if(!is.data.frame(mapping) && !is.tibble(mapping)) {
    stop("mapping has to be a data.frame or a tibble")
  }

  group_column <- enquo(group_column)
  id_column <- enquo(id_column)

  if(! (rlang::as_name(group_column) %in% names(mapping))) {
    stop("`group_column` does not appear to be a member of `mapping`")
  }
  if(! (rlang::as_name(id_column) %in% names(mapping))) {
    stop("`id_column` does not appear to be a member of `mapping`")
  }

  group_values_tmp <- mapping %>% pull(!!group_column)
  names(group_values_tmp) <- mapping %>% pull(!!id_column)
  group_values <- group_values_tmp[rownames(compute_result$observed_matrix)]


  base_mds <- compute_result$base_mds

  aligned_bootstrap_points <- compute_result$aligned_bootstrap_mds %>% purrr::map(~ .x$Yrot)
  res_connectivity_stats <- connectivity_stats_all_groups(base_mds$points,
                                                          aligned_bootstrap_points,
                                                          group_values)

  res_consistency_location <- consistency_location(base_mds$points, aligned_bootstrap_points)
  res_consistency_angles <- consistency_angles(base_mds$points, aligned_bootstrap_points)
  res_consistency_distances <- consistency_distances(base_mds$points, aligned_bootstrap_points)

  #Per-point stats

  res_consistency_location_per_point <- consistency_location_per_point(base_mds$points, aligned_bootstrap_points)
  res_consistency_angles_per_point <- consistency_angles_per_point(base_mds$points, aligned_bootstrap_points)
  res_consistency_distances_per_point <- consistency_distances_per_point(base_mds$points, aligned_bootstrap_points)

  structure(
    list(base_mds = base_mds,
         aligned_bootstrap_mds = compute_result$aligned_bootstrap_mds,
         mapping = mapping,
         group_column = group_column,
         group_values = group_values,
         connectivity_stats = res_connectivity_stats,
         consistency_stats = data.frame(consistency_location = res_consistency_location,
                                        consistency_angles = res_consistency_angles,
                                        consistency_distances = res_consistency_distances),
         per_point_consistency = data.frame(location = res_consistency_location_per_point,
                                            angles = res_consistency_angles_per_point,
                                            distances = res_consistency_distances_per_point)
    ),
    class = "mds_sensitivity"
  )
}

mds_sensitivity_check <- function(observed_matrix, mapping,
                                  group_column = Group,
                                  id_column = Sample,
                                  bootstrap_func = bootstrap_func_dm(20),
                                  ...) {

  compute_result <- mds_sensitivity_check_compute(observed_matrix, bootstrap_func, ...)

  mds_sensitivity_check_eval(compute_result, mapping, group_column, id_column)
}

