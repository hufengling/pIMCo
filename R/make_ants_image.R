#' Makes ANTs Image in mask out of vector of intensities
#'
#' @param vec numerical vector of intensities with 0s removed
#' @param mask_indices vector of indices at which mask is 1
#' @param reference antsImage to supply metadata
#'
#' @return antsImage
#' @export
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @importFrom ANTsRCore as.antsImage
make_ants_image <- function(vec, mask_indices, reference) {
  arr <- array(0, dim = dim(reference))
  arr[mask_indices] <- vec
  x <- ANTsRCore::as.antsImage(arr, reference = reference)
  return(x)
}
