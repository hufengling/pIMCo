library(extrantsr)
library(ANTsRCore)

mask <- check_ants("../testdata/gm10_pcals_rest.nii.gz")
alff <- check_ants("../testdata/2647_alffStd.nii.gz")

mask_indices <- which(as.array(mask) == 1)
full_indices <- 1:length(as.array(mask) %>% as.vector())

mask_vec <- mask %>% as.array() %>% as.vector()
alff_vec <- alff %>% as.array() %>% as.vector()

mask_vec_in_mask <- mask_vec[mask_vec == 1]

test_that("make_ants_image works", {
  expect_equal(make_ants_image(mask_vec, full_indices, mask) %>% as.array(), mask %>% as.array())
  expect_equal(make_ants_image(alff_vec, full_indices, mask) %>% as.array(), alff %>% as.array())
  expect_equal(make_ants_image(mask_vec_in_mask, mask_indices, mask) %>% as.array(), mask %>% as.array())
  expect_equal(make_ants_image(mask_vec_in_mask, mask_indices, alff) %>% as.array(), mask %>% as.array())
})
