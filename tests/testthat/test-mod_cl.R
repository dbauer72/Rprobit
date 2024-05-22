test_that("misspecification of fL can be detected", {
  expect_warning(
    mod <- mod_cl$new(
      Hb       = diag(2),
      fb       = matrix(c(0, 0), ncol = 1),
      HO       = matrix(0, 0, 0),
      fO       = matrix(0, 0, 0),
      HL       = diag(3),
      ### fL is misspecified, should be vector of length 3
      fL       = as.matrix(c(0, 0), ncol = 1),
      alt      = 2,
      ordered  = FALSE,
      validate = TRUE
    ),
    "Dimensions of HL and fL do not match. Resizing fL."
  )
})
