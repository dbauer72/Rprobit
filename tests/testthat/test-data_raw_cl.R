test_that("data_raw_cl object can be created", {
  data_raw <- data_raw_cl$new()
  expect_s3_class(data_raw, c("data_raw_cl", "R6"), exact = TRUE)
})
