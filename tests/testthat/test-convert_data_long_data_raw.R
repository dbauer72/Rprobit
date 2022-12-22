test_that("data in long format can be transformed", {
  data("TravelMode", package = "AER")
  TravelMode["avinc"] <- with(TravelMode, (mode == "air") * income)
  data_raw <- convert_data_long_data_raw(
    df = TravelMode,
    id_dec = "individual",
    id_choice ="individual",
    alternative = "mode",
    choice = "choice"
  )
  expect_s3_class(data_raw, c("data_raw_cl", "R6"), exact = TRUE)
})
