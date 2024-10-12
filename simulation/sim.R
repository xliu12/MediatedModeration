library(hal9001)
library(origami)
library(mvtnorm)
library(parallel)
library(tidyverse)
library(glue)
library(SuperLearner)
library(lightgbm)
library(earth)

GenData <- function(
    sim.seed = 123,
    n = 500,
    Mfamily = "binomial", # "gaussian",
    Yfamily = "binomial", # "gaussian",
    quadratic.R = FALSE,
    quadratic.M = FALSE,
    quadratic.Y = FALSE,
    interact = FALSE,
    no_tt = FALSE,
    num_c = 3
) {
  set.seed(seed = sim.seed)
  # baseline covariates ----
  # num_c <- 5
  C_dat <- rmvnorm(n, sigma = diag(1, nrow = num_c))
  C_dat[, 1] <- 1*(C_dat[, 1] >0)
  if (num_c > 3) {
    C_dat[, which(1:num_c %% 3 == 1)] <- 1*(C_dat[, which(1:num_c %% 3 == 1)] > 0)
  }
  # C_dat <- 1*(C_dat > 0)
  dat <- data.frame(ID = 1:n)
  dat$C <- C_dat

  # Subgroup ----
  gen_r <- list(r_c = sqrt(0.15/num_c))
  if (quadratic.R == FALSE) {
    r_given <- gen_r[["r_c"]] * rowSums(dat$C)
  }
  if (quadratic.R == TRUE) {
    Cquad <- (dat$C ^ 2 - 1) / sqrt(2)
    # Cquad <- dat$C[, 1]*dat$C
    # gen_r[["r_c"]] <- sqrt(0.4/num_c)
    r_given <- gen_r[["r_c"]] * rowSums(Cquad)
  }
  dat$R <- rbinom(n, 1, pr(1, r_given, gen_r))

  # Treatment ------
  gen_tt <- list(tt_c = sqrt(0.15/num_c), tt_r = 0.2)

  if (quadratic.R == FALSE) {
    tt_given <- gen_tt[["tt_c"]] * rowSums(dat$C) + gen_tt[["tt_r"]] * dat$R
  }
  if (quadratic.R == TRUE) {
    Cquad <- (dat$C ^ 2 - 1) / sqrt(2)
    # Cquad <- dat$C[, 1]*dat$C
    # gen_tt[["tt_c"]] <- sqrt(0.4/num_c)
    tt_given <- gen_tt[["tt_c"]] * rowSums(Cquad) + gen_tt[["tt_r"]] * dat$R
  }

  dat$tt <- rbinom(n, 1, ptt(1, tt_given, gen_tt))



  # intermediate  -----------------------
  if (Mfamily == "binomial") {
    zintercept <- -1
    mintercept <- -1
  }
  if (Mfamily == "gaussian") {
    zintercept <- 0
    mintercept <- 0
  }
  if (interact == FALSE) {
    gen_z <- list(zintercept = zintercept,
                  z_r = 0.1, z_tt = 0.1, z_rtt =0.5,
                  z_c = sqrt(0.15/num_c) )
    gen_m <- list(mintercept = mintercept,
                  m_r = 0.1, m_tt = 0.1, m_z = 0.2,
                  m_rtt = 0.5,
                  m_rz = 0, m_ttz = 0,
                  m_zrtt = 0,
                  m_c = sqrt(0.15/num_c) )
  }
  if (interact == TRUE) {
    gen_z <- list(zintercept = zintercept,
                  z_r = 0.1, z_tt = 0.1, z_rtt =0.5,
                  z_c = sqrt(0.15/num_c) )
    gen_m <- list(mintercept = mintercept,
                  m_r = 0.1, m_tt = 0.1, m_z = 0.2,
                  m_rtt = 0.5,
                  m_rz = 0, m_ttz = 0,
                  m_zrtt = 0.2,
                  m_c = sqrt(0.15/num_c) )

    if (Mfamily == "gaussian" ) { #|Mfamily == "binomial"
      gen_z <- list(zintercept = zintercept,
                    z_r = -0.05, z_tt = 0.1, z_rtt = 0.1,
                    z_c = sqrt(0.15/num_c) )
      gen_m <- list(mintercept = mintercept,
                    m_r = -0.05, m_tt = 0.1, m_z = 0.2,
                    m_rtt = 0.1,
                    m_rz = 0, m_ttz = 0,
                    m_zrtt = -0.1,
                    m_c = sqrt(0.15/num_c) )
    }
    if (Mfamily == "binomial"& Yfamily== "binomial" ) { #|Mfamily == "binomial"
      gen_z[["z_rtt"]] <- 0.9
      gen_m[["m_rtt"]] <- 0.9
    }
  }
  if (quadratic.M == FALSE) {
    z_given <- gen_z[["zintercept"]] + gen_z[["z_c"]] * rowSums(dat$C)
    m_given <- gen_m[["mintercept"]] + gen_m[["m_c"]] * rowSums(dat$C)
  }
  if (quadratic.M == TRUE) {
    Cquad <- (dat$C ^ 2 - 1) / sqrt(2)
    # Cquad <- dat$C[, 1]*dat$C
    # gen_z[["z_c"]] <- sqrt(0.4/num_c)
    # gen_m[["m_c"]] <- sqrt(0.4/num_c)
    z_given <- gen_z[["zintercept"]] + gen_z[["z_c"]] * rowSums(Cquad)
    m_given <- gen_m[["mintercept"]] + gen_m[["m_c"]] * rowSums(Cquad)
  }

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

  if (interact == FALSE) {
    gen_y <- list(yintercept = yintercept,
                  y_r = 0.2, y_tt = 0.2, y_m = 1, y_z = 1,
                  y_rtt = 1,
                  y_rm = 0, y_ttm = 0, y_rz = 0, y_ttz = 0, y_zm = 0,
                  y_rttm = 0, y_rttz = 0, y_rmz = 0, y_ttmz = 0,
                  y_ttrzm = 0,
                  y_c = sqrt(0.15 / num_c) )
  }
  if (interact == TRUE) {
    # gen_y <- list(yintercept = yintercept,
    #               y_r = 0.2, y_tt = 0.2, y_m = 1, y_z = 1,
    #               y_rtt = 1,
    #               y_rm = 0, y_ttm = 0, y_rz = 0, y_ttz = 0, y_zm = 1,
    #               y_rttm = 0, y_rttz = 0, y_rmz = 0, y_ttmz = 0,
    #               y_ttrzm = 0.8,
    #               y_c = sqrt(0.15 / num_c) )
    #
    # if(Mfamily=="gaussian") {
    gen_y <- list(yintercept = yintercept,
                  y_r = 0.2, y_tt = 0.2, y_m = 1, y_z = 1,
                  y_rtt = 1,
                  y_rm = 0, y_ttm = 0, y_rz = 0, y_ttz = 0, y_zm = 1,
                  y_rttm = 0.4, y_rttz = 0.4, y_rmz = 0, y_ttmz = 0,
                  y_ttrzm = 1,
                  y_c = sqrt(0.15 / num_c) )
    # }

  }
  if (quadratic.Y == FALSE) {
    y_given <- gen_y[["yintercept"]] + gen_y[["y_c"]] * rowSums(dat$C)
  }
  if (quadratic.Y == TRUE) {
    Cquad <- (dat$C ^ 2 - 1) / sqrt(2)
    # Cquad <- dat$C[, 1]*dat$C
    # gen_y[["y_c"]] <- sqrt(0.4/num_c)
    y_given <- gen_y[["yintercept"]] + gen_y[["y_c"]] * rowSums(Cquad)
  }

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
    if (no_tt==FALSE) {
      true_effs <- effect(truevals)
    }
    if (no_tt==TRUE) {
      true_effs <- effect_nott(truevals)
    }
    return(true_effs)
  }

  truevals <- trueVals()

  out <- mget( ls(), envir = environment() )

  return(out)
}





# Simulation----------------------------

condition_all <- data.frame(expand.grid(
  n = c(200, 300, 500, 800, 1000, 2000, 3000, 5000),
  # generating data
  interact=c(F, T),
  quadratic.R = c(F, T),
  quadratic.M = c(F, T),
  quadratic.Y = c(F, T),
  Mfamily = c("binomial", "gaussian"), #"gaussian",
  Yfamily = c("gaussian", "binomial"),
  num_c = c(20)
))

condition <- condition_all %>%
  filter(as.numeric(quadratic.R + quadratic.M + quadratic.Y) %in% c(0),
         n==1000, interact==T,
         num_c == 20,
         Mfamily=="binomial",
         Yfamily=="gaussian")

methds_all <- data.frame(expand.grid(
  fity.interact = c(F),
  Fit = c("glm", "mlr")
))
methds <- methds_all %>%
  filter(Fit == "mlr")

set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 1



OneData <- function(iseed = 1, cond = 1){

  Yfamily <- as.character(condition$Yfamily[cond])
  Mfamily <- as.character(condition$Mfamily[cond])
  interact <- as.logical(condition$interact[cond])
  num_c <- condition$num_c[cond]

  # true values -----------------------------
  gen_largeN <- GenData(sim.seed = datseeds[iseed],
                        n = 1e5,
                        Mfamily = Mfamily,
                        Yfamily = Yfamily,
                        quadratic.R = condition$quadratic.R[cond],
                        quadratic.M = condition$quadratic.M[cond],
                        quadratic.Y = condition$quadratic.Y[cond],
                        interact = interact,
                        no_tt = FALSE,
                        num_c = num_c )
  true_values <- gen_largeN$truevals

  unlist(true_values)

  gen_data <- GenData(
    sim.seed = datseeds[iseed],
    n = condition$n[cond],
    Mfamily = Mfamily,
    Yfamily = Yfamily, # "binomial",
    quadratic.R = condition$quadratic.R[cond],
    quadratic.M = condition$quadratic.M[cond],
    quadratic.Y = condition$quadratic.Y[cond],
    interact = interact,
    no_tt = FALSE,
    num_c = num_c
  )
  data_in <- gen_data$datobs

  Rname <- grep("^R", colnames(data_in), value = T)
  ttname <- grep("^tt", colnames(data_in), value = T)
  Yname <- grep("^Y", colnames(data_in), value = T)
  Cnames <- grep("^C", colnames(data_in), value = T)
  Znames <- grep("^Z", colnames(data_in), value = T)
  Mnames <- grep("^M", colnames(data_in), value = T)
  Mediators <- list(M1=Znames, M2=Mnames)

  # estimators -----

  jj <- 1
  one.jj <- function (jj = 1) {

    # par_out <- par.MedMod(data_in,
    #                       Yname,
    #                       Rname,
    #                       ttname,
    #                       Mediators,
    #                       Cnames,
    #                       learners,
    #                       Yfamily = Yfamily, linear.Y = F,
    #                       TotMod_dr = FALSE)
    # respar_out <- data.frame(methds[jj, , drop=F],
    #                          estimand = names(true_values),
    #                          true_val = unlist(true_values),
    #                          par_out,
    #                          row.names = NULL) %>%
    #   mutate( rbias = across(contains("est_"), ~1*(abs(true_val)>1e-3)*abs(. / true_val -1)+1*(abs(true_val)<=1e-3)*abs(. - true_val) ),
    #           bias = across(contains("est_"), ~abs(. - true_val) ) ) %>%
    #   do.call(data.frame, .) %>%
    #   mutate(across(where(~is.numeric(.)), ~round(., 3)))


    Fit <- as.character(methds$Fit[jj])

    if (Fit == "glm") {
      learners <- learners_a <- learners_m <- learners_y <- c("SL.glm")
      num_folds <- 1
    }
    if (Fit == "mlr") {
      # SL.hal , "SL.lightgbm"
      # learners <- c("SL.glm", "SL.mean", "SL.nnet")
      # learners <- c("SL.glm", "SL.glm.interaction",  "SL.nnet")
      learners <- c("SL.glm", "SL.nnet")
      # learners_h <- learners # "SL.xgboost"
      # learners_mu <- learners
      num_folds <- ifelse(condition$n[cond]<=300, 5, 4)
      num_folds <- ifelse(condition$n[cond] >= 2000, 4, 4)
    }

    if (Mfamily != "binomial") {
      learners <- c("SL.glm",  "SL.nnet")
      learners_h <- learners # "SL.xgboost"
      learners_mu <- learners
      # if (condition$n[cond] >= 2000) {
      #   learners_h <- c("SL.glm", "SL.ranger", "SL.nnet")
      #  # learners_mu <- c("SL.glm", "SL.glm.interaction", "SL.nnet")
      # }

      hout <- h.MedMod(
        data_in,
        Yname,
        Rname,
        ttname,
        Mediators,
        Cnames,
        learners,
        learners_h,
        learners_mu,
        Yfamily,
        Mfamily,
        fity.interact = FALSE,
        fitm.interact = FALSE,
        num_folds,
        TotMod_dr = FALSE,
        full.sample = TRUE,
        no_tt = FALSE,
        out_eif = FALSE
      )
      estimates <- hout

      # reshout1 <- data.frame(methds[jj, , drop=F],
      #                        estimand = names(true_values),
      #                        true_val = unlist(true_values),
      #                        hout,
      #                        row.names = NULL) %>%
      #   mutate( rbias = across(contains("est_"), ~abs(. / true_val -1)),
      #           bias = across(contains("est_"), ~abs(. - true_val) ) ) %>%
      #   do.call(data.frame, .) %>%
      #   mutate(across(where(~is.numeric(.)), ~round(., 3))) %>%
      #   select(contains("bias"), everything())

    }

    # save(hout, file = "h_RData/hout.RData")


    # save(reshout, file = "h_RData/reshout.RData")

    #
    if (Mfamily == "binomial") {
      # learners <- c("SL.glm", "SL.ranger", "SL.nnet")

      binaryout <- MedMod(
        data_in,
        Yname,
        Rname,
        ttname,
        Mediators,
        Cnames,
        learners,
        Yfamily,
        fity.interact = TRUE, #FALSE,
        fitm.interact = TRUE, #FALSE,
        num_folds,
        TotMod_dr = FALSE,
        out_eif = FALSE
      )
      estimates <- binaryout
    }



    # resbinary_out <- data.frame(methds[jj, , drop=F],
    #                             estimand = names(true_values),
    #                             true_val = unlist(true_values),
    #                             binaryout,
    #                             row.names = NULL) %>%
    #   mutate( rbias = across(contains("est_"), ~1*(abs(true_val)>1e-3)*abs(. / true_val -1)+1*(abs(true_val)<=1e-3)*abs(. - true_val) ),
    #           bias = across(contains("est_"), ~abs(. - true_val) ) ) %>%
    #   do.call(data.frame, .) %>%
    #   mutate(across(where(~is.numeric(.)), ~round(., 3)))



    res <- data.frame(methds[jj, , drop=F],
                      estimand = names(true_values),
                      true_val = unlist(true_values),
                      estimates,
                      row.names = NULL)

    return(res)
  }
  # one_boot <- res
  oneboot <- lapply(1:nrow(methds), one.jj)
  one_boot <- do.call(rbind, oneboot)

  res <- data.frame(condition[cond, ],
                    one_boot,
                    row.names = NULL)

  # res1 <- res %>% mutate( rbias = across(contains("est_"), ~abs(. / true_val -1) ) )

  return(res)

}


# Parallel setup -----
num_reps <- 20
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))
