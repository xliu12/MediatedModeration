

mu.mzc <- function(r, tt, data_in, varnames, Yfamily = "gaussian", interact, folds, learners, bounded = FALSE, full.sample = FALSE) {

  if (full.sample == TRUE) {
    tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  }
  if (full.sample == FALSE) {
    tr_vals <- expand.grid(tt=tt, r=r)
  }

  mu <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(mu) <- (glue("mu(t{tr_vals$tt},r{tr_vals$r},m,z,c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(tr_vals),
                         function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- tr_vals$r[jj]
      valid_v[, c(varnames$tt)] <- tr_vals$tt[jj]
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })
    if (full.sample == TRUE) {
      if (interact == TRUE) {
        alist <- crossfit(train, valid_list,
                          varnames$Y,
                          c(varnames$tt, varnames$R, varnames$M, varnames$Z,
                            # two-way
                            varnames$Rtt,
                            varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
                            # three-way and four-way
                            varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ,
                            varnames$C), varnames,
                          type = Yfamily, learners, "ID", bounded)
      }
      if (interact == FALSE) {
        alist <- crossfit(train, valid_list,
                          varnames$Y,
                          c(varnames$tt, varnames$R, varnames$M, varnames$Z,
                            # # two-way
                            # varnames$Rtt,
                            # varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
                            # # three-way and four-way
                            # varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ,
                            varnames$C), varnames,
                          type = Yfamily, learners, "ID", bounded)
      }
    }
    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
        (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        varnames$Y,
                        c(varnames$M, varnames$Z,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:nrow(tr_vals)) {
      mu[folds[[v]]$validation_set, glue("mu(t{tr_vals$tt[jj]},r{tr_vals$r[jj]},m,z,c)")] <-  preds[, jj]
    }
  }
  mu_mzc <- mu
  mu_mzc
}



# partially marginalized outcome

muM.zc <- function(r, tt, mu_mzc, h_mstarstar, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  # h_mstarstar <- h.m.starstar(r, rstarstar, tt, h_zm, r_tc, r_mtc)
  # h_mstarstar <- h.zstar.mstarstar(r, r, rstarstar, tt, h_zm, r_tc, r_mtc, r_ztc)
  data_in[["mu*hm"]] <- mu_mzc[, glue("mu(t{tt},r{r},m,z,c)")] * h_mstarstar
  muM_zc <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muM_zc) <- glue("muM(z,c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- r
      valid_v[, c(varnames$tt)] <- tt
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu*hm",
                        c(varnames$tt, varnames$R, varnames$Z,
                          # varnames$Rtt, varnames$RZ, varnames$ttZ, varnames$ttRZ,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
              (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu*hm",
                        c(varnames$Z, varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:1) {
      muM_zc[folds[[v]]$validation_set, glue("muM(z,c)")] <- 1 - preds[, jj]
      muM_zc[folds[[v]]$validation_set, glue("muM(z,c)")] <- preds[, jj]
    }
  }

  muM_zc
}



muZ.mc <- function(r, tt, mu_mzc, h_zstar, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  # h_zstar <- h.z.star(r, rstar, tt, h_zm, r_tc, r_ztc)
  # h_zstar <- h.zstar.mstarstar(r, rstar, r, tt, h_zm, r_tc, r_mtc, r_ztc)
  data_in[["mu*hz"]] <- mu_mzc[, glue("mu(t{tt},r{r},m,z,c)")] * h_zstar

  # tvals <- expand.grid(tt=c(0,1))
  # muZ_mc <- matrix(nrow = nrow(data_in), ncol = nrow(tvals))
  # colnames(muZ_mc) <- glue("muM(t{tvals$tt},r{r},z,c)")
  muZ_mc <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muZ_mc) <- glue("muZ(m,c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- r
      valid_v[, c(varnames$tt)] <- tt
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })
    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu*hz",
                        c(varnames$tt, varnames$R, varnames$M,
                          # varnames$Rtt, varnames$RM, varnames$ttM, varnames$ttRM,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
         (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu*hz",
                        c(varnames$M, varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:1) {
      muZ_mc[folds[[v]]$validation_set, glue("muZ(m,c)")] <- 1 - preds[, jj]
      muZ_mc[folds[[v]]$validation_set, glue("muZ(m,c)")] <- preds[, jj]
    }
  }

  muZ_mc
}

muZM.c <- function(rstarstar, tt, muZ_mc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  muZM_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muZM_c) <- glue("muZM(c)")

  data_in[["muZ(m,c)"]] <- muZ_mc[, "muZ(m,c)"]
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    valid_list <- lapply(1:1, function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstarstar
      valid_v[, c(varnames$tt)] <- tt
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "muZ(m,c)",
                        c(varnames$tt, varnames$R,
                          # varnames$Rtt,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
        (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        "muZ(m,c)",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    muZM_c[folds[[v]]$validation_set, glue("muZM(c)")] <-  preds[, 1]
  }
  muZM_c
}

muMZ.c <- function(rstar, tt, muM_zc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  muMZ_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muMZ_c) <- glue("muMZ(c)")

  data_in[["muM(z,c)"]] <- muM_zc[, "muM(z,c)"]
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    valid_list <- lapply(1:1, function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstar
      valid_v[, c(varnames$tt)] <- tt
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "muM(z,c)",
                        c(varnames$tt, varnames$R,
                          # varnames$Rtt,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
        (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      # predict_subset <- ( (train[[varnames$tt]] == tt) & (train[[varnames$R]] == rstar) )
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "muM(z,c)",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    muMZ_c[folds[[v]]$validation_set, glue("muMZ(c)")] <-  preds[, 1]
  }
  muMZ_c
}





mu.joint.c <- function(r, rstar, tt, mu_mzc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  data_in[["mu"]] <- mu_mzc[, glue("mu(t{tt},r{r},m,z,c)")]

  # tvals <- expand.grid(tt=c(0,1))
  mu_joint_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(mu_joint_c) <- glue("mu_joint(c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstar
      valid_v[, c(varnames$tt)] <- tt
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })
    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu",
                        c(varnames$tt, varnames$R,
                          # varnames$Rtt,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }


    if (full.sample == FALSE) {
      # predict_subset <- ( (train[[varnames$tt]] == tt) & (train[[varnames$R]] == rstar) )
      predict_subset <- (train[[varnames$tt]] == unique(valid_list[[1]][[varnames$tt]])) &
        (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }


    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    mu_joint_c[folds[[v]]$validation_set, glue("mu_joint(c)")] <- preds[, 1]

  }

  mu_joint_c
}


# ............................-----------------------------------------------------------
# no_tt = TRUE -----------------------------------------------
mu.mzc_t0 <- function(r, data_in, varnames, Yfamily = "gaussian", interact, folds, learners, bounded = FALSE, full.sample = FALSE) {

  if (full.sample == TRUE) {
    tr_vals <- expand.grid(tt=c(0), r=c(0,1))
  }
  if (full.sample == FALSE) {
    tr_vals <- expand.grid(tt=0, r=r)
  }

  mu <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(mu) <- (glue("mu(t{tr_vals$tt},r{tr_vals$r},m,z,c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(tr_vals),
                         function(jj=1) {
                           valid_v <- valid
                           valid_v[, c(varnames$R)] <- tr_vals$r[jj]
                           valid_v[, c(varnames$tt)] <- tr_vals$tt[jj]
                           # update the interactions
                           # two-way
                           valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
                           valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
                           valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
                           valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
                           valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
                           valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
                           # three-way and four-way
                           valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
                           valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
                           valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
                           valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
                           valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
                           valid_v
                         })
    if (full.sample == TRUE) {
      if (interact == TRUE) {
        alist <- crossfit(train, valid_list,
                          varnames$Y,
                          c(varnames$R, varnames$M, varnames$Z,
                            # two-way
                            varnames$RM, varnames$RZ, varnames$MZ,
                            # three-way and four-way
                            varnames$RMZ,
                            varnames$C), varnames,
                          type = Yfamily, learners, "ID", bounded)
      }
      if (interact == FALSE) {
        alist <- crossfit(train, valid_list,
                          varnames$Y,
                          c(varnames$R, varnames$M, varnames$Z,
                            varnames$C), varnames,
                          type = Yfamily, learners, "ID", bounded)
      }
    }
    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        varnames$Y,
                        c(varnames$M, varnames$Z,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:nrow(tr_vals)) {
      mu[folds[[v]]$validation_set, glue("mu(t{tr_vals$tt[jj]},r{tr_vals$r[jj]},m,z,c)")] <-  preds[, jj]
    }
  }
  mu_mzc <- mu
  mu_mzc
}


muM.zc_t0 <- function(r, mu_mzc, h_mstarstar, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  # h_mstarstar <- h.m.starstar(r, rstarstar, tt, h_zm, r_tc, r_mtc)
  # h_mstarstar <- h.zstar.mstarstar(r, r, rstarstar, tt, h_zm, r_tc, r_mtc, r_ztc)
  tt <- 0
  data_in[["mu*hm"]] <- mu_mzc[, glue("mu(t0,r{r},m,z,c)")] * h_mstarstar
  muM_zc <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muM_zc) <- glue("muM(z,c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- r
      valid_v[, c(varnames$tt)] <- 0
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu*hm",
                        c(varnames$R, varnames$Z,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu*hm",
                        c(varnames$Z, varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:1) {
      muM_zc[folds[[v]]$validation_set, glue("muM(z,c)")] <- 1 - preds[, jj]
      muM_zc[folds[[v]]$validation_set, glue("muM(z,c)")] <- preds[, jj]
    }
  }

  muM_zc
}

muZ.mc_t0 <- function(r, mu_mzc, h_zstar, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  # h_zstar <- h.z.star(r, rstar, tt, h_zm, r_tc, r_ztc)
  # h_zstar <- h.zstar.mstarstar(r, rstar, r, tt, h_zm, r_tc, r_mtc, r_ztc)
  tt <- 0
  data_in[["mu*hz"]] <- mu_mzc[, glue("mu(t0,r{r},m,z,c)")] * h_zstar

  # tvals <- expand.grid(tt=c(0,1))
  # muZ_mc <- matrix(nrow = nrow(data_in), ncol = nrow(tvals))
  # colnames(muZ_mc) <- glue("muM(t{tvals$tt},r{r},z,c)")
  muZ_mc <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muZ_mc) <- glue("muZ(m,c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- r
      valid_v[, c(varnames$tt)] <- 0
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })
    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu*hz",
                        c(varnames$R, varnames$M,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu*hz",
                        c(varnames$M, varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (jj in 1:1) {
      muZ_mc[folds[[v]]$validation_set, glue("muZ(m,c)")] <- 1 - preds[, jj]
      muZ_mc[folds[[v]]$validation_set, glue("muZ(m,c)")] <- preds[, jj]
    }
  }

  muZ_mc
}


muZM.c_t0 <- function(rstarstar, muZ_mc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  tt <- 0
  muZM_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muZM_c) <- glue("muZM(c)")

  data_in[["muZ(m,c)"]] <- muZ_mc[, "muZ(m,c)"]
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    valid_list <- lapply(1:1, function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstarstar
      valid_v[, c(varnames$tt)] <- 0
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "muZ(m,c)",
                        c(varnames$R,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }

    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        "muZ(m,c)",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    muZM_c[folds[[v]]$validation_set, glue("muZM(c)")] <-  preds[, 1]
  }
  muZM_c
}


muMZ.c_t0 <- function(rstar, muM_zc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  tt <- 0
  muMZ_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(muMZ_c) <- glue("muMZ(c)")

  data_in[["muM(z,c)"]] <- muM_zc[, "muM(z,c)"]
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    valid_list <- lapply(1:1, function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstar
      valid_v[, c(varnames$tt)] <- 0
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "muM(z,c)",
                        c(varnames$R,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    if (full.sample == FALSE) {
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))
      # predict_subset <- ( (train[[varnames$tt]] == tt) & (train[[varnames$R]] == rstar) )
      alist <- crossfit(train[predict_subset, ], valid_list,
                        "muM(z,c)",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    muMZ_c[folds[[v]]$validation_set, glue("muMZ(c)")] <-  preds[, 1]
  }
  muMZ_c
}


mu.joint.c_t0 <- function(r, rstar, mu_mzc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = FALSE) {
  tt <- 0
  data_in[["mu"]] <- mu_mzc[, glue("mu(t0,r{r},m,z,c)")]

  # tvals <- expand.grid(tt=c(0,1))
  mu_joint_c <- matrix(nrow = nrow(data_in), ncol = 1)
  colnames(mu_joint_c) <- glue("mu_joint(c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:1, function(jj=1) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- rstar
      valid_v[, c(varnames$tt)] <- 0
      # update the interactions
      # two-way
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      # three-way and four-way
      valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
      valid_v
    })
    if (full.sample == TRUE) {
      alist <- crossfit(train, valid_list,
                        "mu",
                        c(varnames$R,
                          varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }


    if (full.sample == FALSE) {
      # predict_subset <- ( (train[[varnames$tt]] == tt) & (train[[varnames$R]] == rstar) )
      predict_subset <- (train[[varnames$R]] == unique(valid_list[[1]][[varnames$R]]))

      alist <- crossfit(train[predict_subset, ], valid_list,
                        "mu",
                        c(varnames$C),
                        varnames,
                        type = "gaussian", learners, "ID", bounded = FALSE)
    }


    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    mu_joint_c[folds[[v]]$validation_set, glue("mu_joint(c)")] <- preds[, 1]

  }

  mu_joint_c
}
