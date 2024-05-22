test_that("parameter identifiability can be checked", {
  ### this model is not identified (beta_3 = beta_1 + beta_2, all estimated)
  mod <- mod_cl$new(
    Hb  = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 0), nrow = 3, ncol = 3),
    fb  = matrix(c(0, 0, 0), ncol = 1),
    HO  = matrix(0, 0, 0),
    fO  = matrix(0, 0, 0),
    HL  = matrix(0, nrow = 3, ncol = 0),
    fL  = matrix(c(0, 0, 1), ncol = 1),
    alt = 2
  )
  expect_warning(
    expect_false(check_identifiability(mod)),
    "Parameter set seems not identifed."
  )

  ### not estimating beta_3 makes the model identified
  mod2 <- mod$clone()
  mod2$Hb <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
  expect_true(check_identifiability(mod2))
})

test_that("identifiability for ordered probit can be checked", {
  ### this model is not identified (beta_3 = beta_1 + beta_2, all estimated)
  mod <- mod_cl$new(
    Hb  = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 0), nrow = 3, ncol = 3),
    fb  = matrix(c(0, 0, 0), ncol = 1),
    HO   = matrix(0, 0, 0),
    fO   = matrix(0, 0, 0),
    HL   = matrix(c(1, 0, 0, 1, 0, 1), nrow = 6, ncol = 1),
    fL   = matrix(0, 6, 1),
    alt  = 4,
    ordered = TRUE
  )
  expect_warning(
    expect_false(check_identifiability(mod)),
    "Parameter set seems not identifed."
  )

  ### not estimating beta_3 makes the model identified
  mod2 <- mod$clone()
  mod2$Hb <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
  expect_true(check_identifiability(mod2))
})

test_that("scale restriction can be checked", {
  ### this model is not identified (no scale restriction)
  mod <- mod_cl$new(
    Hb  = diag(2),
    fb  = matrix(c(0, 0), ncol = 1),
    HO  = matrix(0, 0, 0),
    fO  = matrix(0, 0, 0),
    HL  = diag(3),
    fL  = as.matrix(c(0, 0, 0), ncol = 1),
    alt = 2
  )
  expect_warning(
    expect_false(check_identifiability(mod)),
    "Scale does not change system."
  )

  ### this modification fixes one beta entry, the model becomes identified
  mod2 <- mod$clone()
  mod2$Hb <- diag(2)[, -1, drop = FALSE]
  mod2$fb <- matrix(c(1, 0), ncol = 1)
  expect_true(check_identifiability(mod2))

  ### this modification fixes one Sigma entry, the model becomes identified
  mod3 <- mod$clone()
  mod3$HL <- matrix(0, nrow = 3, ncol = 0)
  mod3$fL <- matrix(c(0, 0, 1), ncol = 1)
  expect_true(check_identifiability(mod3))
})
