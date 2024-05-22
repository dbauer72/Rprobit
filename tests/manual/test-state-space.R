test_that("test #1: state space, CML, zero initial distribution", {
  mod <- mod_AR_cl$new(
    alt = 3,
    Hb = diag(2)[, -2, drop = FALSE],
    fb = matrix(0, 2, 1),
    HO = matrix(1, 1, 1),
    fO = matrix(0, 1, 1),
    HL = matrix(0, 1, 0),
    fL = matrix(1, 1, 1),
    ordered = TRUE
  )

  mod$alt <- 3
  mod$fb[2] <- 1
  mod$lag_length <- 2
  mod$stationary <- FALSE

  RNGkind(sample.kind = "Rejection")
  set.seed(1)

  theta_0 <- c(1, 0, 0, 0.2, 0.5, 0.2) + rnorm(6) * 0.05
  time <- c(1:5)

  N <- 1000
  time <- c(1:5)
  quest <- 1
  Tp <- length(time) * quest

  # full data set contains
  quest_df <- rep(c(1:quest), length(time))
  iota_q <- rep(1, quest)
  time_df <- matrix(0, Tp, 1)
  for (j in 1:length(time)) {
    time_df[(j - 1) * quest + c(1:quest), 1] <- iota_q * time[j]
  }

  # system matrices
  syst <- build_system_from_model_AR_R(theta_0, mod, time)
  tauk <- syst$tauk

  ##### data preparation
  # Cholesky factors for random draws
  LG <- t(chol(syst$GammaT))
  LO <- t(chol(syst$Omega))


  # first decider to set up the data.frame
  X <- matrix(rnorm(2 * Tp), ncol = 2)
  lRE <- dim(LO)[1]
  gam <- LO %*% rnorm(lRE)

  U <- X %*% syst$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
  iU <- cbind((U > tauk[1]), (U > tauk[2]))
  y <- apply(iU, 1, sum) + 1

  df <- data.frame(ID = 1, X1 = X[, 1], X2 = X[, 2], y = y, time = time_df, quest = quest_df)

  # cycle over deciders
  for (j in 2:N) {
    X <- matrix(rnorm(2 * Tp), ncol = 2)
    gam <- LO %*% rnorm(lRE)

    U <- X %*% syst$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    iU <- cbind((U > tauk[1]), (U > tauk[2]))
    y <- apply(iU, 1, sum) + 1
    df <- rbind(df, data.frame(ID = j, X1 = X[, 1], X2 = X[, 2], y = y, time = time_df, quest = quest_df))
  }
  # convert to data_raw_StSp_cl.
  data_raw <- data_raw_StSp_cl$new(df = df, alt_names = c("1", "2", "3"), id = "ID", choice = "y", ordered = TRUE, varying = "", dec_char = c("X1", "X2"))
  data_raw$set_time_col("time")
  data_raw$set_quest_col("quest")

  # set up Rprobit_obj
  form <- y ~ 0 | X1 + X2 | 0

  control_simul <- list(Tp = Tp)
  control <- list(control_simulation = control_simul)
  re <- c("X1")

  Rprobit_obj <- setup_Rprobit(
    form = form, data_raw = data_raw,
    mod = mod,
    re = re,
    seed = 17,
    theta_0 = theta_0,
    control = control_simul
  )

  # estimate parameters using fit_Rprobit
  Rprobit_obj$theta <- Rprobit_obj$theta_0
  set.seed(1)

  Rp <- fit_Rprobit(Rprobit_obj, init_method = "theta")
  # summary(Rp)

  expect_snapshot(round(Rp$theta, 2))
})

test_that("test #2: state space, CML, stationary initialisation", {
  mod <- mod_AR_cl$new(
    alt = 3,
    Hb = diag(2)[, -2, drop = FALSE],
    fb = matrix(0, 2, 1),
    HO = matrix(1, 1, 1),
    fO = matrix(0, 1, 1),
    HL = matrix(0, 1, 0),
    fL = matrix(1, 1, 1),
    ordered = TRUE
  )

  mod$alt <- 3
  mod$fb[2] <- 1
  mod$lag_length <- 2
  mod$stationary <- TRUE

  set.seed(1)
  theta_0 <- c(1, 0, 0, 0.2, 0.5, 0.2) + rnorm(6) * 0.05
  time <- c(1:5)

  N <- 1000
  time <- c(1:5)
  quest <- 1
  Tp <- length(time) * quest

  # full data set contains
  quest_df <- rep(c(1:quest), length(time))
  iota_q <- rep(1, quest)
  time_df <- matrix(0, Tp, 1)
  for (j in 1:length(time)) {
    time_df[(j - 1) * quest + c(1:quest), 1] <- iota_q * time[j]
  }

  # system matrices
  syst <- build_system_from_model_AR_R(theta_0, mod, time)
  tauk <- syst$tauk

  ##### data preparation
  # Cholesky factors for random draws
  LG <- t(chol(syst$GammaT))
  LO <- t(chol(syst$Omega))

  RNGkind(sample.kind = "Rejection")
  set.seed(12)
  # first decider to set up the data.frame
  X <- matrix(rnorm(2 * Tp), ncol = 2)
  lRE <- dim(LO)[1]
  gam <- LO %*% rnorm(lRE)

  U <- X %*% syst$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
  iU <- cbind((U > tauk[1]), (U > tauk[2]))
  y <- apply(iU, 1, sum) + 1

  df <- data.frame(ID = 1, X1 = X[, 1], X2 = X[, 2], y = y, time = time_df, quest = quest_df)

  # cycle over deciders
  for (j in 2:N) {
    X <- matrix(rnorm(2 * Tp), ncol = 2)
    gam <- LO %*% rnorm(lRE)

    U <- X %*% syst$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    iU <- cbind((U > tauk[1]), (U > tauk[2]))
    y <- apply(iU, 1, sum) + 1
    df <- rbind(df, data.frame(ID = j, X1 = X[, 1], X2 = X[, 2], y = y, time = time_df, quest = quest_df))
  }
  # convert to data_raw_StSp_cl.
  data_raw <- data_raw_StSp_cl$new(df = df, alt_names = c("1", "2", "3"), id = "ID", choice = "y", ordered = TRUE, varying = "", dec_char = c("X1", "X2"))
  data_raw$set_time_col("time")
  data_raw$set_quest_col("quest")

  # set up Rprobit_obj
  form <- y ~ 0 | X1 + X2 | 0

  control_simul <- list(Tp = Tp)
  control <- list(control_simulation = control_simul)
  re <- c("X1")

  Rprobit_obj <- setup_Rprobit(
    form = form, data_raw = data_raw,
    mod = mod,
    re = re,
    seed = 17,
    theta_0 = theta_0,
    control = control_simul
  )

  # estimate parameters using fit_Rprobit
  Rprobit_obj$theta <- Rprobit_obj$theta_0
  set.seed(1)

  suppressWarnings(
    Rp <- fit_Rprobit(Rprobit_obj, init_method = "theta")
  )
  # summary(Rp)

  expect_snapshot(round(Rp$theta, 2))
})
