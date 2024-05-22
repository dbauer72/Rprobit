test_that("test #1: latent class, CML", {
  form <- choice ~ V1 + V2 | 0
  re <- c("V1")
  mod <- mod_latclass_cl$new(
    Hb = diag(2),
    fb = as.matrix(c(0, 0), ncol = 1),
    HO = matrix(0, 1, 0),
    fO = matrix(0.15, 1, 1),
    HL = matrix(0, 6, 0),
    fL = matrix(0, 6, 1),
    ordered = FALSE,
  )
  mod$fL[c(1, 4, 6), 1] <- 1
  mod$num_class <- 3

  # dec per decider: 3.
  RNGkind(sample.kind = "Rejection")
  set.seed(1)
  control_simul <- list(Tp = rep(5, 500))

  # draw parameters
  theta_0 <- c(3, 0, 0, 3, 3, 3, 0, 0)
  seed <- 1
  # simulate data
  Rprobit_obj <- setup_Rprobit(
    form = form,
    mod = mod,
    re = re,
    seed = seed,
    theta_0 = theta_0,
    control = control_simul
  )

  Rprobit_obj$control$probit <- FALSE
  Rprobit_obj$control$control_weights$cml_pair_type <- 0 # full pairwise
  Rprobit_obj$theta <- theta_0 # start at true value

  # generate mod_cl object for use in ll_macml
  modc <- mod_cl$new(
    Hb = diag(2),
    fb = as.matrix(c(0, 0), ncol = 1),
    HO = matrix(0, 1, 0),
    fO = matrix(0.15, 1, 1),
    HL = matrix(0, 6, 0),
    fL = matrix(0, 6, 1),
    ordered = FALSE,
  )
  modc$fL[c(1, 4, 6), 1] <- 1


  # evaluate at true parameter vector
  Rprobit_obj$theta <- rnorm(2)
  Rprobit_obj$theta_0 <- c(3, 0, .1)
  Rprobit_obj$mod <- modc
  
  Rp <- fit_Rprobit(Rprobit_obj, init_method = "theta")
  expect_snapshot(round(Rp$theta, 2))

  Rprobit_obj$theta <- rnorm(length(theta_0))
  Rprobit_obj$theta_0 <- theta_0
  Rprobit_obj$control$el <- 1
  Rprobit_obj$mod <- mod

  Rp2 <- fit_LC_Rprobit(Rprobit_obj, init_method = "kmeans")
  expect_snapshot(round(Rp2$theta, 2))
  
})

test_that("test #2: latent class, EM", {
  form <- choice ~ V1 + V2 | 0
  re <- c("V1")
  mod <- mod_latclass_cl$new(
    Hb = diag(2),
    fb = as.matrix(c(0, 0), ncol = 1),
    HO = matrix(0, 1, 0),
    fO = matrix(0.15, 1, 1),
    HL = matrix(0, 6, 0),
    fL = matrix(0, 6, 1),
    ordered = FALSE,
  )
  mod$fL[c(1, 4, 6), 1] <- 1
  mod$num_class <- 3

  # dec per decider: 3.
  RNGkind(sample.kind = "Rejection")
  set.seed(1)
  control_simul <- list(Tp = rep(5, 100))

  # draw parameters
  theta_0 <- c(3, 0, 0, 3, 3, 3, 0, 0)
  seed <- 1
  # simulate data
  Rprobit_obj <- setup_Rprobit(
    form = form,
    mod = mod,
    re = re,
    seed = seed,
    theta_0 = theta_0,
    control = control_simul
  )

  Rprobit_obj$control$probit <- FALSE
  Rprobit_obj$control$control_weights$cml_pair_type <- 0 # full pairwise
  Rprobit_obj$theta <- theta_0 # start at true value

  # generate mod_cl object for use in ll_macml
  modc <- mod_cl$new(
    Hb = diag(2),
    fb = as.matrix(c(0, 0), ncol = 1),
    HO = matrix(0, 1, 0),
    fO = matrix(0.15, 1, 1),
    HL = matrix(0, 6, 0),
    fL = matrix(0, 6, 1),
    ordered = FALSE,
  )
  modc$fL[c(1, 4, 6), 1] <- 1


  # evaluate at true parameter vector
  Rprobit_obj$theta <- rnorm(2)
  Rprobit_obj$theta_0 <- c(3, 0, .1)
  Rprobit_obj$mod <- modc

  Rp <- fit_Rprobit(Rprobit_obj, init_method = "theta")
  expect_snapshot(round(Rp$theta, 2))

  Rprobit_obj$theta <- rnorm(length(theta_0))
  Rprobit_obj$theta_0 <- theta_0
  Rprobit_obj$control$el <- 1
  Rprobit_obj$mod <- mod

  control_EM <- list()
  control_EM$iter_lim <- 100
  control_EM$iter_lim_one_step <- 5
  control_EM$tol <- 0.00001
  control_EM$optim_all <- FALSE

  Rp2 <- fit_LC_EM_Rprobit(Rprobit_obj, init_method = "kmeans", control_EM = control_EM)
  expect_snapshot(round(Rp2$theta, 2))
})
