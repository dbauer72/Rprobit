test_that("test #1: simple cross sectional probit with three choices", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1), ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = as.matrix(c(0,0,0,0,0,0)),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #2: same as test 1, but probit estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1), ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = as.matrix(c(0,0,0,0,0,0)),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #3: same as test 1, but panel with three choices per person", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(3, 1000))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #4: same as test 3, probit estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(3, 1000))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #5: same as test 3, Tp 1-3 per decision, macml estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  RNGkind(sample.kind = "Rejection")
  set.seed(1)
  control_simul = list(Tp = sample(1:3, 1000, replace = TRUE))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #6: same as test 5, probit estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c()
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  RNGkind(sample.kind = "Rejection")
  set.seed(17)
  control_simul = list(Tp = sample(1:3, 1000, replace = TRUE))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = c("SJ", "TVBS", "TVBSv2", "ME"),
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #7: same as test 5, 2 random coefficients, macml estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c("V1", "V2")
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = diag(3),
    fO   = matrix(0,3,1),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  RNGkind(sample.kind = "Rejection")
  set.seed(17)
  control_simul = list(Tp = sample(1:3, 1000, replace = TRUE))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #8: same as test 7, probit  estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  re <- c("V1", "V2")
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = diag(3),
    fO   = matrix(0,3,1),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  RNGkind(sample.kind = "Rejection")
  set.seed(1)
  control_simul = list(Tp = sample(1:3, 2000, replace = TRUE))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #9: same as test 7, only 1 random coefficient", {
  form <- choice ~ V1 + V2 + V3 | 0
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(1,1,1),
    fO   = matrix(0,1,1),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  re <- c("V2")
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #10: same as test 9, probit estimation", {
  form <- choice ~ V1 + V2 + V3 | 0
  mod <- mod_cl$new(
    Hb   = diag(3)[,-3],
    fb   = as.matrix(c(0,0,1),ncol=1),
    HO   = matrix(1,1,1),
    fO   = matrix(0,1,1),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  re <- c("V2")
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #11: V1 is alternative specific regressor, macml estimation", {
  form <-  choice ~ V3 + V2  | 0 | V1
  mod <- mod_cl$new(
    Hb   = diag(5)[,-1],
    fb   = as.matrix(c(1,0,0,0,0),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  re <- c("V2")
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #12: same as test 11, probit estimation", {
  form <-  choice ~ V3 + V2  | 0 | V1
  mod <- mod_cl$new(
    Hb   = diag(5)[,-1],
    fb   = as.matrix(c(1,0,0,0,0),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = diag(6)[,-c(1,2,3)],
    fL   = matrix(0,6,1),
    ordered = FALSE
  )
  control_simul = list(Tp = rep(1, 1000))
  re <- c("V2")
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #13: five categories, two regressors, one random coefficient", {
  form <-  choice ~ V1 + V2 | 0
  mod <- mod_cl$new(
    Hb   = diag(2),
    fb   = as.matrix(diag(2)[,1],ncol=1),
    HO   = matrix(1,1,1),
    fO   = matrix(0,1,1),
    HL   = diag(15)[,c(10,13,15)],
    fL   = matrix(0,15,1),
    ordered = FALSE
  )
  re <- c("V2")
  mod$fL[6] <- 1
  RNGkind(sample.kind = "Rejection")
  set.seed(1)
  control_simul = list(Tp = sample(1:4, 500, replace = TRUE))
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #15: ordered probit, estimated using the CML approach", {
  form <-  choice ~ 0 |V1 + V2 | 0
  mod <- mod_cl$new(
    Hb   = diag(2)[,-1,drop=FALSE],
    fb   = as.matrix(c(1,0),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = matrix(1,1,1),
    fL   = as.matrix(c(1,0,1),ncol=1),
    ordered = TRUE
  )
  Tp <- 3
  HL <- diag(Tp*(Tp+1)/2)
  HL <- as.matrix(apply(HL[,c(1,4,6),drop=FALSE],1,sum))
  mod$HL <- HL
  mod$fL <- matrix(0,6,1)
  control_simul <- list(Tp = rep(3,1000))
  control <- list(control_simulation = control_simul)
  re <- c()
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = FALSE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #16: same as test 15, probit estimation", {
  form <-  choice ~ 0 |V1 + V2 | 0
  mod <- mod_cl$new(
    Hb   = diag(2)[,-1,drop=FALSE],
    fb   = as.matrix(c(1,0),ncol=1),
    HO   = matrix(0,0,0),
    fO   = matrix(0,0,0),
    HL   = matrix(1,1,1),
    fL   = as.matrix(c(1,0,1),ncol=1),
    ordered = TRUE
  )
  Tp <- 3
  HL <- diag(Tp*(Tp+1)/2)
  HL <- as.matrix(apply(HL[,c(1,4,6),drop=FALSE],1,sum))
  mod$HL <- HL
  mod$fL <- matrix(0,6,1)
  control_simul <- list(Tp = rep(3,1000))
  control <- list(control_simulation = control_simul)
  re <- c()
  est <- suppressMessages(
    run_mini_example(
      form = form,
      mod = mod,
      re = re,
      approx_method = "SJ",
      show_theta = FALSE,
      seed = 17,
      at_true = TRUE,
      runs = 1,
      probit = TRUE,
      nCores = 1,
      print.level = 0,
      control_simul = control_simul,
      verbose = FALSE
    )
  )
  expect_snapshot(round(est[[1]], 2))
})

test_that("test #17: comparison to mlogit via Train dataset works", {

  data("Train", package = "mlogit")
  Train$choiceid <- 1:nrow(Train)

  # estimate a MNL model using mlogit
  Train.ml <- mlogit::mlogit.data(
    Train, choice = "choice", shape = "wide", varying = 4:11, sep = "_",
    opposite = c("price", "time", "change", "comfort")
  )
  Train.mod <- mlogit::mlogit(
    choice ~ price + time + change + comfort | 0,
    Train.ml, reflevel = "A", R = 500
  )

  # set up the probit model
  train_mod <- setup_Rprobit(
    form = choice ~ price + time + change + comfort | 0,
    data_raw = Train,
    ids = c("id")
  )

  # change scale fixation to price parameter.
  train_mod$mod$HL = matrix(0,3,1)
  train_mod$mod$HL[3,1] = 1
  train_mod$mod$fL = matrix(0,3,1)
  train_mod$mod$Hb = diag(4)[,-1]
  train_mod$mod$fb = as.matrix(c(-0.001,0,0,0),ncol=1)

  # initialize with MNL results
  train_mod$theta_0 <- -c(as.numeric(Train.mod$coefficients[2:4]),1)
  train_fit <- fit_Rprobit(
    Rprobit_obj = train_mod, init_method = "theta",
    control_nlm = list(
      approx_method = "SJ",
      probit = FALSE,
      hess = 0,
      print.level=0
    )
  )

  ### compare logit and probit results
  pr_par = c(-0.001,train_fit$theta[1:train_fit$mod$lthb])
  comparison = rbind(
    "logit" = (Train.mod$coefficients[1:4]),
    "probit" = pr_par,
    "logit_norm" = as.numeric(Train.mod$coefficients[1:4]) /
      Train.mod$coefficients[1],
    "probit_norm" = pr_par/pr_par[1]
  )
  expect_snapshot(round(comparison, 2))
})
