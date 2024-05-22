test_that("model parameters can be built", {
  mod <- mod_cl$new(
    ### first beta coefficient fixed to 1
    Hb      = diag(2)[, -1, drop = FALSE],
    fb      = matrix(c(1, 0), ncol = 1),
    HO      = matrix(0, 0, 0),
    fO      = matrix(0, 0, 0),
    ### normalize utility level
    HL      = diag(6)[, -c(1, 2, 3)],
    fL      = as.matrix(c(0, 0, 0, 0, 0, 0), ncol = 1),
    alt     = 3,
    ordered = FALSE
  )
  theta <- rnorm(4)
  expect_type(build_par_from_mod(theta = theta, mod = mod), "list")
})
