#' Performs PCA-based Intermodal Coupling
#'
#' @param files list of paths to niftis (string) or list of antsImages
#' @param brain_mask path to mask nifti or antsImage mask to apply before coupling (optional)
#' @param out_dir string path to output directory
#' @param out_name string prefix for file name. Final name is "'out_name'_coupling.nii.gz"
#' @param fwhm numerical value of full width at half maximum for calculating neighborhood size; default is 3
#' @param verbose TRUE or FALSE; default TRUE
#' @param prop_miss proportion of voxels in a neighborhood allowed to be missing; default is 0.9
#' @param use_ratio default FALSE for (logit of 1st eigenvalue scaled to 0-1), TRUE for (1st eigenvalue/total variance)
#' @param return_coupling default FALSE to write image to out_dir, TRUE to return coupling image without writing
#'
#' @return list of coupled antsImages and file paths
#' @export
#'
#' @examples
#' \dontrun{
#' mod1_filepath <- "/path/to/proj/dir/mod1_subj1.nii.gz"
#' mod2_filepath <- "/path/to/proj/dir/mod2_subj1.nii.gz"
#' mod3_filepath <- "/path/to/proj/dir/mod3_subj1.nii.gz"
#' imco(files = list(mod1_files, mod2_files, mod3_files),
#'     brain_mask = "/path/to/proj/dir/grey_matter_mask.nii.gz",
#'     out_dir = "/path/to/proj/dir/coupled_images",
#'     out_name = "subj-1",
#'     fwhm = 3)
#' }
#' @importFrom ANTsRCore antsGetSpacing getNeighborhoodInMask
#' @importFrom neurobase zscore_img
imco <- function(files, brain_mask,
                 out_dir = NULL, out_name = NULL,
                 fwhm = 3,
                 verbose = TRUE,
                 prop_miss = 0.9,
                 use_ratio = FALSE,
                 return_coupling = FALSE) {
  fileList <- extrantsr::check_ants(files)
  for (i in 2:length(files)) {
    if (!all(dim(fileList[[i - 1]]) == dim(fileList[[i]]))) {
      stop("Image dimensions do not match")
    }
    if (!all(ANTsRCore::antsGetSpacing(fileList[[i - 1]]) == ANTsRCore::antsGetSpacing(fileList[[i]]))) {
      stop("Voxel dimensions do not match")
    }
    if (!all(ANTsRCore::antsGetDirection(fileList[[i - 1]]) == ANTsRCore::antsGetDirection(fileList[[i]]))) {
      stop("Image directions do not match")
    }
    if (!all(ANTsRCore::antsGetOrigin(fileList[[i - 1]]) == ANTsRCore::antsGetOrigin(fileList[[i]]))) {
      stop("Image origins/locations do not match")
    }
  }
  # Read in brain mask
  mask <- extrantsr::check_ants(brain_mask)
  if (!all(dim(mask) == dim(fileList[[1]]))) {
    stop("Image dimensions do not match the brain mask")
  }
  if (!all(ANTsRCore::antsGetSpacing(mask) == ANTsRCore::antsGetSpacing(fileList[[1]]))) {
    stop("Voxel dimensions do not match the brain mask")
  }

  if (!is.numeric(fwhm)) {
    fwhm <- as.numeric(fwhm)
    if (!is.numeric(fwhm)) {
      stop("fwhm must be a numeric value")
    }
  }

  if (!is.numeric(prop_miss)) {
    prop_miss <- as.numeric(prop_miss)
    if (!is.numeric(prop_miss)) {
      stop("prop_miss must be a numeric value")
    }
  }

  if (!return_coupling) {
    if (is.null(out_dir) | is.null(out_name)) {
      stop("To write coupling image, out_dir and out_name must be specified.")
    }
  }

  fileList <- lapply(fileList, neurobase::zscore_img)
  fileList <- lapply(fileList, check_ants)

  # Dimension of each voxel in mm
  vDims <- ANTsRCore::antsGetSpacing(fileList[[1]])

  nhood_dims <- get_nhood_size(
    fwhm = fwhm,
    vDims = vDims,
    brainMask = mask,
    verbose = verbose
  )

  # Assign anything outside brain mask to NA
  fileList <- lapply(fileList, function(x) {
    newx <- x
    newx[mask == 0] <- NA
    return(newx)
  })

  if (verbose) {
    cat("# Extracting neighborhood data \n")
  }
  # Neighborhood data from each modality
  mask_indices <- which(as.array(mask) > 0)
  nhoods <- lapply(fileList, function(x) ANTsRCore::getNeighborhoodInMask(image = x, mask = mask, radius = nhood_dims$v_radius, spatial.info = TRUE))

  if (is.null(fwhm) == FALSE) {
    # Check that the dimension from getting neighborhoods is the same as
    # calculated diamx, diamy, and diamz (i.e. check volumes match)
    for (i in 1:length(fileList)) {
      volume <- dim(nhoods[[i]]$values)[1]
      if (volume != nhood_dims$v_diamx * nhood_dims$v_diamy * nhood_dims$v_diamz) {
        stop("Calculated diameters from get_nhood_size (v_diamx, v_diamy, v_diamz) and ANTsRCore::getNeighborhoodInMask do not match. Report bug.")
      }
    }
  }

  # Will use to map back to niftis
  inds <- nhoods[[1]]$indices
  offs <- nhoods[[1]]$offsets

  nhood_weights <- get_weights(offs, vDims, sigma = nhood_dims$sigma)
  coupling <- imco_pca(
    files = fileList,
    nhoods = nhoods,
    nhood_weights = nhood_weights,
    mask_indices = mask_indices,
    verbose = verbose,
    prop_miss = prop_miss,
    pca_type = "global_wcov",
    use_ratio = use_ratio
  )

  if (!return_coupling) {
    antsImageWrite(coupling, file.path(out_dir, paste0(out_name, "_coupling.nii.gz")))
    return(NULL)
  }
  else {
    return(coupling)
  }
}
