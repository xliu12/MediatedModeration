library(origami)
library(mvtnorm)
library(parallel)
library(tidyverse)
library(glue)
library(SuperLearner)

devtools::load_all("multiMedMod")

GenData <- function(
    sim.seed = 123,
    n = 500,
    Mfamily = "binomial", # "gaussian",
    Yfamily = "binomial", # "gaussian",
    num_c = 20
) {
  set.seed(seed = sim.seed)
  # baseline covariates ----
  C_dat <- rmvnorm(n, sigma = diag(1, nrow = num_c))
  C_dat[, 1] <- 1*(C_dat[, 1] >0)
  if (num_c > 3) {
    C_dat[, which(1:num_c %% 3 == 1)] <- 1*(C_dat[, which(1:num_c %% 3 == 1)] > 0)
  }
  dat <- data.frame(ID = 1:n)
  dat$C <- C_dat

  # Subgroup ----
  gen_r <- list(r_c = sqrt(0.15/num_c))

  r_given <- gen_r[["r_c"]] * rowSums(dat$C)

  dat$R <- rbinom(n, 1, pr(1, r_given, gen_r))

  # Treatment ------
  gen_tt <- list(tt_c = sqrt(0.15/num_c), tt_r = 0.2)

  tt_given <- gen_tt[["tt_c"]] * rowSums(dat$C) + gen_tt[["tt_r"]] * dat$R

  dat$tt <- rbinom(n, 1, ptt(1, tt_given, gen_tt))


  # Mediators Z = M1, M = M2  -----------------------
  if (Mfamily == "binomial") {
    zintercept <- -1
    mintercept <- -1
  }
  if (Mfamily == "gaussian") {
    zintercept <- 0
    mintercept <- 0
  }

  gen_z <- list(zintercept = zintercept,
                z_r = -0.05, z_tt = 0.1, z_rtt = 0.1,
                z_c = sqrt(0.15/num_c) )
  gen_m <- list(mintercept = mintercept,
                m_r = -0.05, m_tt = 0.1, m_z = 0.2,
                m_rtt = 0.1,
                m_rz = 0, m_ttz = 0,
                m_zrtt = -0.1,
                m_c = sqrt(0.15/num_c) )
  if (Mfamily == "binomial" ) {
    gen_z[["z_rtt"]] <- 0.9
    gen_m[["m_rtt"]] <- 0.9
  }

  z_given <- gen_z[["zintercept"]] + gen_z[["z_c"]] * rowSums(dat$C)
  m_given <- gen_m[["mintercept"]] + gen_m[["m_c"]] * rowSums(dat$C)

  if (Mfamily == "binomial") {
    dat$Z <- rbinom(n, 1, prob = pz(tt = dat$tt, r = dat$R, z_given, gen_z, if_binary = TRUE) )
    dat$M <- rbinom(n, 1, prob = pm(dat$Z, dat$tt, dat$R, m_given, gen_m, if_binary = TRUE) )
  }
  if (Mfamily == "gaussian") {
    dat$Z <- pz(tt = dat$tt, r = dat$R, z_given, gen_z, if_binary = FALSE) + rnorm(n)
    dat$M <- pm(dat$Z, dat$tt, dat$R, m_given, gen_m, if_binary = FALSE) + rnorm(n)
  }


  # Outcome -----
  if (Yfamily=="gaussian") { yintercept <- 1 }
  if (Yfamily!="gaussian") { yintercept <- -1 }

  gen_y <- list(yintercept = yintercept,
                y_r = 0.2, y_tt = 0.2, y_m = 1, y_z = 1,
                y_rtt = 1,
                y_rm = 0, y_ttm = 0, y_rz = 0, y_ttz = 0, y_zm = 1,
                y_rttm = 0.4, y_rttz = 0.4, y_rmz = 0, y_ttmz = 0,
                y_ttrzm = 1,
                y_c = sqrt(0.15 / num_c) )

  y_given <- gen_y[["yintercept"]] + gen_y[["y_c"]] * rowSums(dat$C)

  if_binary <- (Yfamily == "binomial")

  if (Yfamily == "gaussian") {
    condmy <- my(m = dat$M, dat$Z, dat$tt, dat$R, y_given, gen_y, if_binary = FALSE)
    dat$Y <- condmy + rnorm(n)
  }
  if (Yfamily == "binomial") {
    condmy <- my(m = dat$M, dat$Z, dat$tt, dat$R, y_given, gen_y, if_binary = TRUE)
    dat$Y <- rbinom(n, 1, condmy)
  }

  datobs <- do.call(data.frame, dat)


  trueVals <- function() {
    vals <- rbind(
      expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(0,1), r2 = c(0), rjo = NA),
      expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(1), r2 = c(0,1), rjo = NA),
      expand.grid(tt = c(0, 1), r0 = c(1), r1 = NA, r2 = NA, rjo = c(0,1)),
      expand.grid(tt = c(0, 1), r0 = c(0), r1 = NA, r2 = NA, rjo = c(0))
    )


    if_binaryY <- (Yfamily == "binomial")
    if_binaryZM <- (Mfamily == "binomial")
    truevals <- list()

    j <- 1
    for (j in 1:nrow(vals)) {
      tt <- vals$tt[j]
      r0 <- vals$r0[j]
      r1 <- vals$r1[j]
      r2 <- vals$r2[j]
      rjo <- vals$rjo[j]

      trueval <- true.val(tt, r0, r1, r2, rjo, y_given, gen_y, if_binaryZM, if_binaryY, m_given, z_given, gen_m, gen_z, linear = FALSE)
      if (!is.na(r1)) {
        truevals[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- trueval
      }
      if (!is.na(rjo)) {
        truevals[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]] <- trueval
        if (r0==rjo) {
          truevals[[ glue("delta(t{tt},r{r0})") ]] <- truevals[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]]
        }
      }
    }
    true_effs <- effect(truevals)

    return(true_effs)
  }

  truevals <- trueVals()

  out <- mget( ls(), envir = environment() )

  return(out)
}





# Simulation----------------------------

condition_all <- data.frame(expand.grid(
  n = c(200, 300, 500, 800, 1000),
  Mfamily = c("binomial", "gaussian"),
  Yfamily = c("gaussian", "binomial")
))

condition <- condition_all

datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 1



OneData <- function(iseed = 1, cond = 1){

  Yfamily <- as.character(condition$Yfamily[cond])
  Mfamily <- as.character(condition$Mfamily[cond])

  # true values
  gen_largeN <- GenData(sim.seed = datseeds[iseed],
                        n = 1e5,
                        Mfamily = Mfamily,
                        Yfamily = Yfamily)
  true_values <- gen_largeN$truevals

  gen_data <- GenData(
    sim.seed = datseeds[iseed],
    n = condition$n[cond],
    Mfamily = Mfamily,
    Yfamily = Yfamily
  )
  data_in <- gen_data$datobs

  Rname <- grep("^R", colnames(data_in), value = T)
  ttname <- grep("^tt", colnames(data_in), value = T)
  Yname <- grep("^Y", colnames(data_in), value = T)
  Cnames <- grep("^C", colnames(data_in), value = T)
  M1names <- grep("^Z", colnames(data_in), value = T)
  M2names <- grep("^M", colnames(data_in), value = T)

  # estimators ----

  learners <- c("SL.glm", "SL.nnet")
  num_folds <- ifelse(condition$n[cond]<=300, 5, 4)
  num_folds <- ifelse(condition$n[cond] >= 2000, 4, 4)

  out <- MedMod(data_in = data_in,
                Yname = Yname,
                M1names = M1names, M2names = M2names, ttname = ttname, Rname = Rname, Cnames = Cnames,
                learners = learners,
                num_folds = num_folds)

  estimates <- out
  res <- data.frame(condition[cond, ],
                    estimand = names(true_values),
                    true_val = unlist(true_values),
                    estimates,
                    row.names = NULL)


  return(res)

}


# Parallel setup -----
num_reps <- 20
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))
