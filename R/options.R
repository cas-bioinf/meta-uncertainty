# Using the settings package to manage options
# Variable, global to package's namespace.
# This function is not exported to user space and does not need to be documented.
MYPKGOPTIONS <- settings::options_manager(cores = parallel::detectCores())

#' Set or get options for the metagenboot package
#'
#' @param ... Option names to retrieve option values or \code{[key]=[value]} pairs to set options.
#'
#' @section Supported options:
#' The following options are supported
#' \itemize{
#'  \item{\code{cores}}{(\code{numeric};\code{parallel::detectCores()})
#'  number of cores to use for metaMDS and other multithreading operations }
#' }
#'
#' @export
metagenboot_options <- function(...){
  # protect against the use of reserved words.
  settings::stop_if_reserved(...)
  MYPKGOPTIONS(...)
}
