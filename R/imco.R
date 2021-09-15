#' Makes ANTs Image in mask out of vector of intensities
#'
#' @param vec numerical vector of intensities
#' @param mask_indices binary array - 1 in mask and 0 out of mask
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

#' Assign weights to each voxel in neighborhood
#'
#' Assign weights to each voxel in neighborhood for weighted PCA/regression.
#' Weights are assigned according to sigma (generated from FWHM)
#'
#' @param offsets numerical vector of voxel distances from central voxel
#' @param voxelDims numerical vector of voxel dimensions
#' @param sigma numeric defined from fwhm
#'
#' @return numerical vector of weights defined based on sigma
#' @export
#'
#' @examples
#' \dontrun{
#' get_weights(offsets = c(3, 1, 3), voxelDims = c(2, 2, 2), sigma = sigma)
#' }
get_weights <- function(offsets, voxelDims, sigma) {
  dists <- t(apply(offsets, 1, function(x, y) x * y, y = voxelDims))
  sqNorm <- apply(dists * dists, 1, function(x) sum(x))
  # Constant doesn't matter because we would rescale by
  # max weight so that center voxel would be be 1
  return(exp(-1 * sqNorm / (2 * sigma^2)))
}

#' Calculate size of neighborhood
#'
#' @param fwhm numeric value of full width at half maximum
#' @param vDims numeric vector of voxel dimensions
#' @param brainMask antsImage
#' @param verbose TRUE or FALSE
#'
#' @return list of voxel properties
#' @export
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

#' Performs Intermodal Coupling
#'
#' @param files list of files to couple
#' @param brainMask antsImage mask to apply before coupling (optional)
#' @param out_dir string for output directory
#' @param out_name TODO
#' @param fwhm numerical value of full width at half maximum for calculating neighborhood size
#' @param verbose TRUE or FALSE
#' @param full_pca_dir string directory if the full PCA extraction is desired. Otherwise, NULL by default and only coupled image is returned.
#' @param propMiss TODO
#' @param pcaType "global" or "unscaled"
#' @param matrixType "wcov" (weighted covariance) or "wcor" (weighted correlation)
#' @param returnNeighborhoods TODO
#'
#' @return list of coupled antsImages and file paths
#' @export
#'
#' @examples
#' \dontrun{
#' mod1_files <- list("mod1_subj1.nii.gz")
#' mod2_files <- list("mod2_subj1.nii.gz")
#' mod3_files <- list("mod3_subj1.nii.gz")
#' imco(files = list(mod1_files, mod2_files, mod3_files),
#'     brainMask = "inst/extdata/gm10_pcals_rest.nii.gz",
#'     out_dir = coupled_images, out_name = "", type = "pca", fwhm = 3)
#' }
#' @importFrom ANTsRCore antsGetSpacing getNeighborhoodInMask
#' @importFrom neurobase zscore_img
imco <- function(files, brainMask, out_dir, out_name,
                 fwhm = NULL,
                 verbose = TRUE,
                 full_pca_dir = NULL, propMiss = 1,
                 pcaType = NULL, matrixType = NULL,
                 returnNeighborhoods = FALSE) {
  nf <- length(files)
  fileList <- extrantsr::check_ants(files)
  for (i in 2:length(files)) {
    if (!all(dim(fileList[[i - 1]]) == dim(fileList[[i]]))) {
      stop("Image dimensions do not match")
    }
    if (!all(ANTsRCore::antsGetSpacing(fileList[[i - 1]]) == ANTsRCore::antsGetSpacing(fileList[[i]]))) {
      stop("Voxel dimensions do not match")
    }
  }
  # Read in brain mask
  bMask <- extrantsr::check_ants(brainMask)
  if (!all(dim(bMask) == dim(fileList[[1]]))) {
    stop("Image dimensions do not match the brain mask")
  }
  if (!all(ANTsRCore::antsGetSpacing(bMask) == ANTsRCore::antsGetSpacing(fileList[[1]]))) {
    stop("Voxel dimensions do not match the brain mask")
  }

  if (pcaType == "global") {
    fileList <- lapply(fileList, neurobase::zscore_img)
    fileList <- lapply(fileList, check_ants)
  }

  # Dimension of each voxel in mm
  vDims <- ANTsRCore::antsGetSpacing(fileList[[1]])

  nhood_dims <- get_nhood_size(
    fwhm = fwhm,
    vDims = vDims,
    brainMask = bMask,
    verbose = verbose
  )

  # Assign anything outside brain mask to NA
  fileList <- lapply(fileList, function(x) {
    newx <- x
    newx[bMask == 0] <- NA
    return(newx)
  })

  if (verbose) {
    cat("# Extracting neighborhood data \n")
  }
  # Neighborhood data from each modality
  mask_indices <- which(as.array(bMask) > 0)
  nhoods <- lapply(fileList, function(x) ANTsRCore::getNeighborhoodInMask(image = x, mask = bMask, radius = nhood_dims$v_radius, spatial.info = TRUE))

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
  # ANTsR version 0.3.2 used 0 first index instead of 1
  # This was changed here: https://github.com/stnava/ANTsR/commit/706148aa994d414c9efd76e9292ae99351e2c4be
  # Thus, have commented out the following line and require latest release of ANTsR version 0.3.3
  # inds = inds + 1
  # Will use to compute distances from center voxel
  nWts <- get_weights(offs, vDims, sigma = nhood_dims$sigma)
  pcaObj <- imco_pca(
    files = fileList,
    nhoods = nhoods,
    nWts = nWts,
    mask_indices = mask_indices,
    verbose = verbose,
    full_pca_dir = full_pca_dir,
    propMiss = propMiss,
    pcaType = pcaType,
    matrixType = matrixType,
    returnNeighborhoods = returnNeighborhoods
  )
  eigenvalues <- pcaObj$eigenValueImages
  eigenvalue_sum <- Reduce("+", eigenvalues)
  coupling <- eigenvalues[[1]] / eigenvalue_sum
  antsImageWrite(coupling, file.path(out_dir, paste0(out_name, "_coupling.nii.gz")))

  return(coupling)
}

#' Intermodal Coupling with PCA
#'
#' @param files TODO
#' @param nhoods TODO
#' @param nWts TODO
#' @param mask_indices TODO
#' @param verbose TRUE or FALSE
#' @param full_pca_dir TODO
#' @param propMiss TODO
#' @param pcaType "global" or "unscaled"
#' @param matrixType "wcov" or "wcor"
#' @param returnNeighborhoods TRUE or FALSE
#'
#' @return TODO
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @importFrom rlist list.rbind
#' @importFrom stats cov
#' @importFrom ANTsRCore antsImageWrite
imco_pca <- function(files, nhoods, nWts, mask_indices, verbose = TRUE, full_pca_dir = NULL, propMiss = NULL, pcaType = NULL, matrixType = NULL, returnNeighborhoods = FALSE) {
  # Restructure to get eigen decomp at each voxel
  imgVals <- lapply(nhoods, function(x) x$values)
  bigDf <- rlist::list.rbind(imgVals)
  matList <- lapply(split(bigDf, col(bigDf)), function(x) matrix(x, ncol = length(files)))
  rmnaList <- lapply(matList, function(x) {
    w <- nWts
    xRows <- apply(as.matrix(x), 1, function(z) {
      !any(is.na(z))
    })
    if (sum(xRows) > 2) {
      # If the proportion of missing voxels is greater than propMiss, return NA
      if (mean(!xRows) > propMiss) {
        return(NA)
      } else {
        return(cbind(w[xRows], as.matrix(x)[xRows, ]))
      }
    }
    return(NA)
  })
  if (returnNeighborhoods == TRUE) {
    return(rmnaList)
  }
  if (verbose) {
    cat("# Computing weighted covariances \n")
  }
  if (pcaType == "unscaled" | pcaType == "global") {
    rmnaListCenter <- lapply(rmnaList, function(x) {
      if (!is.na(x)[1]) {
        w <- x[, 1]
        newx <- scale(x[, -1], center = TRUE, scale = FALSE)
        return(cbind(w, newx))
      }
      return(NA)
    })
  }
  if (pcaType == "scaled") {
    rmnaListCenter <- lapply(rmnaList, function(x) {
      if (!is.na(x)[1]) {
        w <- x[, 1]
        newx <- scale(x[, -1], center = TRUE, scale = TRUE)
        return(cbind(w, newx))
      }
      return(NA)
    })
  }
  rm(rmnaList)
  # Weighted cov of each matrix in matList
  if (matrixType == "wcov") {
    wcovList <- lapply(rmnaListCenter, function(x) {
      if (!is.na(x)[1]) {
        w <- x[, 1]
        newx <- x[, -1]
        return(stats::cov.wt(newx, wt = w, center = FALSE)$cov)
      }
      return(NA)
    })
  }
  if (matrixType == "wcor") {
    wcovList <- lapply(rmnaListCenter, function(x) {
      if (!is.na(x)[1]) {
        w <- x[, 1]
        newx <- x[, -1]
        return(stats::cov.wt(newx, wt = w, center = FALSE, cor = T)$cor)
      }
      return(NA)
    })
  }
  if (matrixType == "unweighted") {
    wcovList <- lapply(rmnaListCenter, function(x) {
      if (!is.na(x)[1]) {
        newx <- x[, -1]
        return(cov(newx))
      }
      return(NA)
    })
  }
  wcovList_corrected <- lapply(wcovList, function(x) {
    diag <- diag(x)
    if (any(is.na(diag) | diag == 0)) {
      diag(x) <- rep(1, nrow(x))
      x[is.na(x)] <- 0
    }

    return(x)
  })
  if (verbose) {
    cat("# Computing weighted PCs \n")
  }
  eigenList <- lapply(wcovList_corrected, function(x) {
    if (!is.na(x)[1]) {
      return(eigen(x))
    } else {
      return(NA)
    }
  })
  rm(wcovList)
  if (verbose) {
    cat("# Extracting IMCo images \n")
  }
  evals <- list()
  components <- list()
  for (j in 1:length(files)) {
    tmp <- as.vector(rlist::list.rbind(lapply(
      eigenList,
      function(x) {
        if (!is.na(x)[1]) {
          return(x$values[j])
        } else {
          return(NA)
        }
      }
    )))
    evals[[j]] <- make_ants_image(vec = tmp, mask_indices = mask_indices, reference = files[[1]])
    if (!is.null(full_pca_dir)) {
      ANTsRCore::antsImageWrite(evals[[j]], file.path(full_pca_dir, paste0("eigenValue-", j, ".nii.gz")))
    }
    components[[j]] <- list()
    for (k in 1:length(files)) {
      tmp <- as.vector(rlist::list.rbind(lapply(eigenList, function(x) {
        if (!is.na(x)[1]) {
          return(x$vectors[k, j])
        } else {
          return(NA)
        } # Removed y = rmnaListCenter from function because it doesn't show up in the function.
      })))
      components[[j]][[k]] <- make_ants_image(vec = tmp, mask_indices = mask_indices, reference = files[[1]])
      if (!is.null(full_pca_dir)) {
        ANTsRCore::antsImageWrite(components[[j]][[k]], file.path(full_pca_dir, paste0("component", j, "-", k, ".nii.gz")))
      }
    }
  }
  rm(eigenList)

  evals <- lapply(evals, extrantsr::ants2oro)
  for (j in 1:length(files)) {
    temp <- lapply(components[[j]], extrantsr::ants2oro)
    components[[j]] <- lapply(components[[j]], extrantsr::ants2oro)
  }
  return(list("eigenValueImages" = evals, "eigenVectorImages" = components))
}


load_images <- function(dir, mask) {
  files <- list.files(dir)
  file_paths <- file.path(dir, files)
  image_vector_list <- lapply(file_paths,
                              function(x, mask) {
                                image <- extrantsr::check_ants(x)
                                image[mask == 0] <- NA
                                image_vector <- image %>% as.numeric()
                                image_vector_no_NAs <- image_vector[!is.na(image_vector)]
                                return(image_vector_no_NAs)
                              },
                              mask = mask
  )
  return(image_vector_list)
}

transpose_list <- function(list) {
  matrix <- list %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = length(list))
  transposed_list <- lapply(seq_len(ncol(matrix)), function(i) matrix[, i])
  return(transposed_list)
}

make_descriptive_images <- function(voxel_vector_list) {
  descriptive_vector <- lapply(voxel_vector_list, function(voxel_vector) {
    voxel_mean <- sum(voxel_vector)/length(voxel_vector)
    voxel_variance <- var(voxel_vector)

    return(c(voxel_mean, voxel_variance))
  }) %>% unlist()

  voxel_mean_vector <- descriptive_vector[c(T, F)]
  voxel_variance_vector <- descriptive_vector[c(F, T)]

  return(list(voxel_mean_vector, voxel_variance_vector))
}

#' Calculates pvals for each voxel
#'
#' @param voxel_vector vector of length 1xn, where n is number of subjects
#' @param predictors design matrix of predictors
#'
#' @return vector of p-values from various models
#' @export
#'
#' @examples
#' \dontrun{
#' get_pvals_by_voxel(1:5, matrix(1:10, nrow = 5))
#' }
#' @importFrom stats lm var
get_pvals_by_voxel <- function(voxel_vector, predictors) {
  if (length(voxel_vector) != nrow(predictors)) {
    stop("n doesn't match")
  }

  regression <- stats::lm(voxel_vector ~
                            sex + ageAtScan1 +
                            race2 + pcaslRelMeanRMSMotion + restRelMeanRMSMotion,
                          data = predictors) %>%
    summary()
  reg_pvals <- regression$coefficients[c(2, 3), 4]

  interaction_regression <- stats::lm(voxel_vector ~
                                        sex * ageAtScan1 +
                                        race2 + pcaslRelMeanRMSMotion + restRelMeanRMSMotion,
                                      data = predictors) %>%
    summary()
  int_pvals <- interaction_regression$coefficients[8, 4]

  pvals <- c(reg_pvals, int_pvals)
  return(pvals)
}

#' Writes p-values to antsImage
#'
#' @param image_list TODO
#' @param mask TODO
#' @param dir TODO
#' @param is_descriptive TODO
#' @param out_dir TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @importFrom stringr str_split
write_pvals <- function(image_list, mask, dir, is_descriptive = FALSE, out_dir) {
  mask_indices <- which(as.array(mask) > 0)
  reference <- extrantsr::check_ants(file.path(dir, list.files(dir)[[1]]))
  file_name <- (dir %>% stringr::str_split("/"))[[1]][2]

  if (is_descriptive) {
    dir.create("three_modality/descriptive_images", showWarnings = FALSE)
    names <- c("_mean", "_variance")
    for (i in 1:length(image_list)) {
      descriptive_image <- make_ants_image(image_list[[i]], mask_indices, reference)
      ANTsRCore::antsImageWrite(
        descriptive_image,
        file.path(
          out_dir,
          paste0(file_name, names[i], ".nii.gz")
        )
      )
    }

    return(NULL)
  }

  for (i in 1:length(image_list)) {
    pval_image <- make_ants_image(image_list[[i]], mask_indices, reference)
    ANTsRCore::antsImageWrite(
      pval_image,
      file.path(
        out_dir,
        paste0(file_name, "_pval_", i, ".nii.gz")
      )
    )
  }

  return(NULL)
}

hypothesis_test_voxels <- function(dir, mask, predictors, out_dir_descriptive, out_dir_pval) {
  image_vector_list <- load_images(dir, mask)
  voxel_vector_list <- transpose_list(image_vector_list)
  descriptive_list <- make_descriptive_images(voxel_vector_list)
  pvalbyvoxel_list <- parallel::mclapply(voxel_vector_list,
                                         get_pvals_by_voxel,
                                         predictors = predictors,
                                         mc.cores = as.numeric(Sys.getenv("LSB_DJOB_NUMPROC")))
  pvalbycoef_list <- transpose_list(pvalbyvoxel_list)

  write_pvals(descriptive_list, mask, dir, is_descriptive = TRUE, out_dir = out_dir_descriptive)
  write_pvals(pvalbycoef_list, mask, dir, out_dir = out_dir_pval)
}
