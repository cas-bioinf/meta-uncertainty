plot_mds <- function(base_mds, mapping, aligned_samples = NULL,  color_aes = NULL, shape_aes = NULL) {
  color_aes <- enquo(color_aes)
  shape_aes <- enquo(shape_aes)
  to_data <- function(x) { 
    x %>% as.data.frame() %>% set_names("MDS1","MDS2") %>% rownames_to_column("observation") 
  }
  base_points <- base_mds$points %>% to_data() %>% mutate(sample = "Orig")
  sampled_points <- aligned_samples %>% purrr::imap(~ .x$Yrot %>% to_data() %>% mutate(sample = paste0("S_",.y))) %>%
    do.call(rbind, .)
  
  #TODO randomize to have different sample order for each drawn path.
  if(is.null(aligned_samples)) {
    path_geom <- NULL
  } else {
    path_geom <- geom_path(aes(group = observation), size = 1, alpha = 0.1)
  }
  
  rbind(base_points, sampled_points) %>%
    inner_join(mapping %>% mutate(Sample = as.character(Sample)), by =c("observation" = "Sample")) %>%
    mutate(is_original = sample == "Orig") %>%
    ggplot(aes(MDS1, MDS2, shape = !!shape_aes, color = !!color_aes, size = is_original, alpha = is_original)) + 
    path_geom +
    geom_point() + 
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 3)) + 
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1))
}