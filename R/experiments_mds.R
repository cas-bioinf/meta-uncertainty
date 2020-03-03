run_mds_experiment <- function(otu_tab, bootstrap_func, n_observations_to_use_min, n_observations_to_use_max,
                               n_steps, mapping,
                               accept_subset_func,
                               trymax = 20, maxit = 50,
                               ...) {

  step_results <- list()
  full_mds <- otu_tab %>% metaMDS(trymax = trymax, maxit = maxit, trace = 0)


  for(n in 1:n_steps) {
    n_observations_to_use <- rdunif(1, n_observations_to_use_min, n_observations_to_use_max)
    repeat {
      observations_to_use <- sample(rownames(otu_tab), size = n_observations_to_use)
      otu_tab_filtered <- otu_tab[observations_to_use,]
      mapping_filtered <- mapping %>% filter(Sample %in% observations_to_use)

      if(accept_subset_func(mapping_filtered)){
        break;
      }
    }

    compute_result <- mds_sensitivity_check_compute(
                        observed_matrix = otu_tab_filtered, bootstrap_func = bootstrap_func,
                        trace = 0, trymax = trymax, maxit = maxit, ...)


    step_results[[n]] <- list(compute_result = compute_result,
                              observations_to_use = observations_to_use,
                              otu_tab_filtered = otu_tab_filtered,
                              mapping_filtered = mapping_filtered)

    if(n %% 10 == 0) {
      cat("Step", n, " completed\n")
    }

  }

  list(full_mds = full_mds, step_results = step_results)
}


#' @export
eval_mds_experiment <- function(run_result, group_column) {
  stats_list <- list()
  per_point_stats_list <- list()
  sens_check_list <- list()

  full_mds <- run_result$full_mds

  step_results <- run_result$step_results
  n_steps <- length(step_results)

  for(n in 1:n_steps) {
    sens_check <- mds_sensitivity_check_eval(step_results[[n]]$compute_result,
                                             step_results[[n]]$mapping_filtered,
                                            group_column = {{ group_column }})

    sens_check_list[[n]] <- sens_check

    observations_to_use <- step_results[[n]]$observations_to_use
    n_observations_to_use <- length(observations_to_use)

    full_mds_subset <- full_mds$points[observations_to_use,]
    procrustes_to_full <- vegan::procrustes(full_mds_subset, sens_check$base_mds)


    location_error_per_point <- sqrt(rowSums((full_mds_subset - procrustes_to_full$Yrot) ^ 2))

    rmse_distances_per_point <- sqrt(rowMeans((
      as.matrix(dist(full_mds_subset, method = "euclidean")) -
        as.matrix(dist(procrustes_to_full$Yrot, method = "euclidean"))) ^ 2))

    rmse_distances <- sqrt( mean(
      (dist(full_mds_subset, method = "euclidean") - dist(procrustes_to_full$Yrot, method = "euclidean")) ^ 2 ))


    angles_full_subset <- compute_angles(full_mds_subset)
    angles_step <- compute_angles(sens_check$base_mds$points)

    rmse_angles_per_point <- numeric(n_observations_to_use)
    for(p in 1:n_observations_to_use){
      rmse_angles_per_point[p] <- sqrt(mean(
        (get_angles_per_point(angles_full_subset, n_observations_to_use, p) -
          get_angles_per_point(angles_step, n_observations_to_use, p)) ^ 2))
    }

    rmse_angles <- sqrt( mean(
      (angles_full_subset - angles_step) ^ 2 ))



    stats_list[[n]] <- sens_check$connectivity_stats %>%
      summarise(min_connectivity_min = min(connectivity_min),
                min_connectivity_average = min(connectivity_average),
                avg_connectivity_min = mean(connectivity_min),
                avg_connectivity_average = mean(connectivity_average)) %>%
      mutate( step = n,
              n_observations = n_observations_to_use,
              rmse_location = sqrt(procrustes_to_full$ss / n_observations_to_use),
              rmse_distances = rmse_distances,
              rmse_angles = rmse_angles
              ) %>%
      crossing(sens_check$consistency_stats)

    per_point_stats_list[[n]] <- cbind(sens_check$per_point_consistency, n_observations = n_observations_to_use,
                                       data.frame(step = n, observation = observations_to_use,
                                                  location_error = location_error_per_point,
                                                  rmse_distances = rmse_distances_per_point,
                                                  rmse_angles = rmse_angles_per_point))

    if(n %% 10 == 0) {
      cat("Step", n, " completed\n")
    }

  }

  list(global = stats_list %>% do.call(rbind, .),
       per_point = per_point_stats_list %>% do.call(rbind, .),
       sens_checks = sens_check_list,
       step_results = step_results)
}
