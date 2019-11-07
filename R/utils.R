#' @export
is_point_matrix <- function(points) {
  dims <- dim(points)
  is.numeric(points) && !is.null(dims) && length(dims == 2) && dims[2] == 2
}

#' @importFrom rlang !!
#' @export
check_point_matrix <- function(points) {
  name <- rlang::as_label(dplyr::enquo(points))
  if(!is_point_matrix(points)) {
    stop(paste0("The parameter `", name, "` is not a matrix (numeric matrix with dimensions [X, 2])"))
  }
}
