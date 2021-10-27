test_that("get_weights works", {
  expect_equal(get_weights(as.matrix(c(5, 3, 1), nrows = 1), voxel_dims = c(2, 2, 2), sigma = 3), get_weights(as.matrix(c(10, 6, 2), nrows = 1), voxel_dims = c(1, 1, 1), sigma = 3))
  expect_equal(sum(get_weights(as.matrix(c(5, 3)), voxel_dims = c(3, 1), sigma = 2)), sum(get_weights(as.matrix(c(3, 5)), voxel_dims = c(1, 3), sigma = 2)))
  expect_length(get_weights(as.matrix(1:6, nrows = 3), voxel_dims = c(1, 1), sigma = 2), 6)
  expect_gt(get_weights(as.matrix(1:6, nrows = 3), voxel_dims = c(1, 1), sigma = 2)[1], get_weights(as.matrix(1:6, nrows = 3), voxel_dims = c(1, 1), sigma = 2)[2])
})
