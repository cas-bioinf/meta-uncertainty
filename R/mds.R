metaMDS_per_sample <- function(posterior, cores = parallel::detectCores(), ...) {
  cl <- parallel::makePSOCKcluster(cores);
  currentLibPaths <- .libPaths()
  tryCatch({
    parallel::clusterExport(cl, "currentLibPaths", environment())
    parallel::clusterEvalQ(cl, .libPaths(currentLibPaths))
    parallel::clusterExport(cl, "posterior", environment())
    parallel::parLapplyLB(cl, 1:get_n_posterior_samples(posterior),
                          fun = function(i) {vegan::metaMDS(metagenboot::get_posterior_sample(posterior,i), parallel = 1, ...)},
                          chunk.size = 1)
  }, finally = {
    parallel::stopCluster(cl)
  })
}

align_single_MDS <- function(mds_sample, base_mds, allow_unmatched_observations = FALSE) {
  if(!identical(rownames(vegan::scores(base_mds)), rownames(vegan::scores(mds_sample)))) {
    if(!allow_unmatched_observations) {
      stop("Observations are not aligned")
    }
    all_obs_base <- rownames(vegan::scores(base_mds))
    all_obs_sample <- rownames(vegan::scores(mds_sample))
    selected_observations <- all_obs_base[all_obs_base %in% all_obs_sample]
    vegan::procrustes(X = vegan::scores(base_mds)[selected_observations,], Y = vegan::scores(mds_sample)[selected_observations,])
  } else {
    vegan::procrustes(X = base_mds, Y = mds_sample)
  }
}


mds_sensitivity_dm_check <- function(N_samples, observed_matrix, mapping, group_column = Group, id_column = Sample, N_draws = "original",
                                     prior = default_prior, ...) {
  if(!is.matrix(observed_matrix)) {
    stop("observed_matrix has to be a matrix")
  }
  if(!is.data.frame(mapping) || !is.tibble(mapping)) {
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
  group_values <- group_values_tmp[rownames(observed_matrix)]

  base_mds <- observed_matrix %>% vegan::metaMDS(...)

  samples_dm <- sample_posterior_dm(N_samples, observed_matrix, prior = prior, N_draws = N_draws)

  resampled_aligned_mds <- samples_dm %>%
    metaMDS_per_sample(...) %>% purrr::map(vegan::procrustes, X = base_mds)

  connectivity_stats <- connectivity_stats_all_groups(base_mds$points,
                                                      resampled_aligned_mds %>% purrr::map(~ .x$Yrot),
                                                      group_values)

  structure(
    list(base_mds = base_mds,
         samples_dm = samples_dm,
         resampled_aligned_mds = resampled_aligned_mds,
         mapping = mapping,
         group_column = group_column,
         group_values = group_values,
         connectivity_stats = connectivity_stats
         ),
    class = "mds_sensitivity"
  )
}

