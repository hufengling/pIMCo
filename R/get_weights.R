#' Assign weights to each voxel in neighborhood
#'
#' Assign weights to each voxel in neighborhood for weighted PCA/regression.
#' Weights are assigned according to sigma (generated from FWHM)
#'
#' @param offsets numerical vector of voxel distances from central voxel
#' @param voxel_dims numerical vector of voxel dimensions
#' @param sigma numeric defined from fwhm
#'
#' @return numerical vector of weights defined based on sigma
#' @noRd
#' @examples
#' \dontrun{
#' get_weights(offsets = c(3, 1, 3), voxelDims = c(2, 2, 2), sigma = sigma)
#' }
get_weights <- function(offsets, voxel_dims, sigma) {
  dists <- t(apply(offsets, 1, function(x, y) x * y, y = voxel_dims))
  sqNorm <- apply(dists * dists, 1, function(x) sum(x))
  # Constant doesn't matter because we would rescale by
  # max weight so that center voxel would be be 1
  return(exp(-1 * sqNorm / (2 * sigma^2)))
}
