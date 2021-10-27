#' Calculate size of neighborhood
#'
#' @param fwhm numeric value of full width at half maximum
#' @param vDims numeric vector of voxel dimensions
#' @param brainMask antsImage
#' @param verbose TRUE or FALSE
#'
#' @return list of voxel properties
#'
#' @examples
#' \dontrun{
#' get_nhood_size(fwhm = 3, vDims = c(2, 2, 2),
#' brainMask = "inst/extdata/gm10_pcals_rest.nii.gz", verbose = TRUE)
#' }
#' @importFrom fslr fslsmooth
#' @importFrom extrantsr check_ants ants2oro
get_nhood_size <- function(fwhm = NULL, vDims = NULL, brainMask, verbose = TRUE) {
  if (is.null(vDims)) {
    stop("Must specify a vector of voxel size in mm")
  }
  if (is.null(fwhm)) {
    stop("Must specify fwhm")
  }
  if (!is.numeric(fwhm)) {
    stop("fwhm must be numeric")
  }
  if (fwhm <= 0) {
    stop("fwhm must be positive")
  }

  # Load brain mask
  bMask <- extrantsr::check_ants(brainMask)

  # FWHM => sigma
  # Note: We specify FWHM in mm.
  # Using 2.35482 which is what FSL uses
  # Note 2.35482 \approx 2*sqrt(2*log(2))
  sigma <- fwhm / 2.35482
  bMaskOro <- extrantsr::ants2oro(brainMask)
  bMaskOro <- bMaskOro * 0
  midVoxel <- floor(dim(bMaskOro) / 2)
  bMaskOro[midVoxel[1], midVoxel[2], midVoxel[3]] <- 1
  sm <- fslr::fslsmooth(bMaskOro, sigma = sigma)
  # AMV implementation
  # Set diam to match vDims x,y,z
  # x direction length
  diamx <- sum(sm[, midVoxel[2], midVoxel[3]] != 0)
  # y direction length
  diamy <- sum(sm[midVoxel[1], , midVoxel[3]] != 0)
  # z direction length
  diamz <- sum(sm[midVoxel[1], midVoxel[2], ] != 0)
  # create radius vector
  radius <- c((diamx - 1) / 2, (diamy - 1) / 2, (diamz - 1) / 2)
  if (any(radius < 1)) {
    stop("FWHM is too small for at least one side of the neighborhood")
  }

  # Message the mm and voxel dimensions
  if (verbose == TRUE) {
    # Message the x,y,z diameter length in voxels
    message(paste0(
      "# Neighborhood diameters in voxels are as follows: \n",
      "# \t x=", diamx, "\n",
      "# \t y=", diamy, "\n",
      "# \t z=", diamz, "\n"
    ))
    # Message the dimension x,y,z diameter length in voxels
    message(paste0("# Neighborhood is x=", diamx, " by y=", diamy, " by z=", diamz, " voxels"))

    # Message the x,y,z diameter length in mm
    message(paste0(
      "# Neighborhood diameters in mm are as follows: \n",
      "# \t x=", diamx * vDims[1], "\n",
      "# \t y=", diamy * vDims[2], "\n",
      "# \t z=", diamz * vDims[3], "\n"
    ))

    # Message the dimension x,y,z diameter length in mm
    message(paste0("# Neighborhood is x=", diamx * vDims[1], " by y=", diamy * vDims[2], " by z=", diamz * vDims[3], " mm"))
  }

  res <- list(
    v_radius = radius,
    v_diamx = diamx,
    v_diamy = diamy,
    v_diamz = diamz,
    mm_radius = radius * vDims,
    mm_diamx = diamx * vDims[1],
    mm_diamy = diamy * vDims[2],
    mm_diamz = diamz * vDims[3],
    sigma = sigma
  )

  return(res)
}
