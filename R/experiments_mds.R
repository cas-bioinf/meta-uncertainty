#' @export
run_mds_experiment <- function(otu_tab, n_samples_to_use, n_steps, mapping, group_column,
                               accept_samples_func, n_sensitivity_samples = 20,
                               trymax = 20, maxit = 50,
                               ...) {
  stats_list <- list()

  full_mds <- otu_tab %>% metaMDS(trymax = trymax, maxit = maxit, trace = 0)


  for(n in 1:n_steps) {
    n_samples_to_use <- 15
    repeat {
      samples_to_use <- sample(rownames(otu_tab), size = n_samples_to_use)
      otu_tab_filtered <- otu_tab[samples_to_use,]
      mapping_filtered <- mapping %>% filter(Sample %in% samples_to_use)

      if(accept_samples_func(mapping_filtered)){
        break;
      }
    }

    sens_check <- mds_sensitivity_check(n_sensitivity_samples, otu_tab_filtered, mapping_filtered,
                                        group_column = {{ group_column }}, trace = 0,
                                        trymax = trymax, maxit = maxit,
                                        ...)

    full_mds_subset <- full_mds$points[samples_to_use,]
    procrustes_to_full <- vegan::procrustes(full_mds_subset, sens_check$base_mds)

    rmse_angles <- sqrt( mean(
      (compute_angles(full_mds_subset) - compute_angles(sens_check$base_mds$points)) ^ 2 ))

    rmse_distances <- sqrt( mean(
      (dist(full_mds_subset, method = "euclidean") - dist(procrustes_to_full$Yrot, method = "euclidean")) ^ 2 ))

    stats_list[[n]] <- sens_check$connectivity_stats %>%
      summarise(min_connectivity_min = min(connectivity_min),
                min_connectivity_average = min(connectivity_average),
                avg_connectivity_min = mean(connectivity_min),
                avg_connectivity_average = mean(connectivity_average)) %>%
      mutate( n_samples = n_samples_to_use,
              rmse_location = sqrt(procrustes_to_full$ss / n_samples_to_use),
              rmse_angles = rmse_angles,
              rmse_distances = rmse_distances
              ) %>%
      crossing(sens_check$consistency_stats)

    if(n %% 10 == 0) {
      cat("Step", n, " completed\n")
    }

  }

  stats_list %>% do.call(rbind, .)
}
