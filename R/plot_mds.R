plot_mds <- function(base_mds, mapping, aligned_samples = NULL,  color_aes = NULL, shape_aes = NULL, show_paths = "all", sample_point_alpha = 0.3) {
  color_aes <- enquo(color_aes)
  shape_aes <- enquo(shape_aes)
  
  
  n_samples <- length(aligned_samples)
  if(is.null(aligned_samples) || show_paths == "none" || show_paths == 0) {
    path_geom <- NULL
    path_frac <- 0
    path_group_indices = rep(1, n_samples)
  } else if(show_paths == "all" || show_paths == 1) {
    path_geom <- geom_path(aes(group = observation), size = 1, alpha = 0.1)
    path_frac <- 1
    path_group_indices <- rep(1, n_samples)
    print(path_group_indices)
  } else if(is.numeric(show_paths) && show_paths > 0 && show_paths < 1) {
    path_geom <- geom_path(aes(group = paste0(observation,"___",path_group)), size = 1, alpha = 0.1)
    path_frac <- show_paths
    n_path_groups <- round(n_samples * (1 - path_frac))
    if(n_path_groups == 0) {
      path_group_indices <- 2:(n_samples + 1)
    } else {
      path_group_indices <- rep(1:(n_path_groups), length.out = n_samples)
    }
  } else {
    stop("show_paths has to be a number between 0 and 1 or 'none' or 'all'")
  }
  
  to_data <- function(x) { 
    x %>% as.data.frame() %>% set_names("MDS1","MDS2") %>% rownames_to_column("observation") 
  }
  base_points <- base_mds$points %>% to_data() %>% mutate(sample = "Orig")
  
  
  sampled_points <- aligned_samples %>% 
    purrr::imap(~ .x$Yrot %>% to_data() %>% 
                  mutate(sample = paste0("S_",.y))
               ) %>%
    do.call(rbind, .)
  
  rbind(base_points, sampled_points) %>%
    inner_join(mapping %>% mutate(Sample = as.character(Sample)), by =c("observation" = "Sample")) %>%
    mutate(is_original = sample == "Orig") %>%
    group_by(observation, is_original) %>%
    mutate(path_group = ifelse(is_original, 1, sample(path_group_indices))) %>% #deliberately not using if_else as I need to hack around a bit
    ungroup() %>%
    sample_frac() %>% #resamples everything - just reorders, to have random paths
    ggplot(aes(MDS1, MDS2, shape = !!shape_aes, color = !!color_aes, size = is_original, alpha = is_original)) + 
    path_geom +
    geom_point() + 
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) + 
    scale_alpha_manual(values = c("FALSE" = sample_point_alpha, "TRUE" = 1))
}