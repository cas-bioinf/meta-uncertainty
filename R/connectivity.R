#' Asymetric direct connectivity between pair of observations: 
#'     fraction of resamples of source further from source than the target
#'     
#' Direct connectivity between pair: maximum of asymetric direct connectivity
#' 
#' Connectivity betweeen pair: width of the widest path in the graph of direct connectivity
#' 
#' Connectivity for a group: minimum of connectivity between pairs across the group
#' 
#' Since widest path has to be completely included in the maximum spanning tree, connectivity for a group can be found as
#' the minimum edge on the maximum spanning tree.