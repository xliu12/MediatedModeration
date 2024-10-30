
# subgroup
r.c <- function(data_in, varnames, folds, learners, bounded = TRUE) {
  r_c <- matrix(nrow = nrow(data_in), ncol = 2)
  colnames(r_c) <- c("r(0|c)", "r(1|c)")
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    alist <- crossfit(train, list(valid),
                      varnames$R,
                      c(varnames$C), varnames,
                      type = c("binomial"), learners, "ID", bounded)
    preds <- alist$preds
    r_c[folds[[v]]$validation_set, "r(0|c)"] <- 1 - preds[, 1]
    r_c[folds[[v]]$validation_set, "r(1|c)"] <- preds[, 1]
  }
  r_c
}


r.tc <- function(data_in, varnames, folds, learners, bounded = TRUE, no_tt = FALSE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  if (no_tt==TRUE) {
    tr_vals <- expand.grid(tt=c(0), r=c(0,1))
  }
  r_tc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(r_tc) <- (glue("r({tr_vals$r}|t{tr_vals$tt},c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    tval <- c(0,1)
    valid_list <- lapply(tval, function(tt=0) {
      valid_v <- valid
      valid_v[, c(varnames$tt)] <- tt
      valid_v
    })

    alist <- crossfit(train, valid_list,
                      varnames$R,
                      c(varnames$tt, varnames$C), varnames,
                      type = c("binomial"), learners, "ID", bounded)
    preds <- alist$preds
    for (jj in 1:length(val)) {
      r_tc[folds[[v]]$validation_set, glue("r(0|t{val[jj]},c)")] <- 1 - preds[, jj]
      r_tc[folds[[v]]$validation_set, glue("r(1|t{val[jj]},c)")] <- preds[, jj]
    }

  }
  r_tc
}



r.ztc <- function(data_in, varnames, interact, folds, learners, bounded = TRUE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  r_ztc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(r_ztc) <- (glue("r({tr_vals$r}|z,t{tr_vals$tt},c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    val <- c(0, 1)
    valid_list <- lapply(val, function(tt=0) {
      valid_v <- valid
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
    if(interact == TRUE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$Z, varnames$tt,
                          varnames$ttZ,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    if(interact == FALSE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$Z, varnames$tt,
                          # varnames$ttZ,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    preds <- alist$preds
    for (jj in 1:length(val)) {
      r_ztc[folds[[v]]$validation_set, glue("r(0|z,t{val[jj]},c)")] <- 1 - preds[, jj]
      r_ztc[folds[[v]]$validation_set, glue("r(1|z,t{val[jj]},c)")] <- preds[, jj]
    }

  }
  r_ztc
}


r.mtc <- function(data_in, varnames, interact, folds, learners, bounded = TRUE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  r_mtc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(r_mtc) <- (glue("r({tr_vals$r}|m,t{tr_vals$tt},c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    val <- c(0, 1)
    valid_list <- lapply(val, function(tt=0) {
      valid_v <- valid
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
    if(interact == TRUE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$M, varnames$tt,
                          varnames$ttM,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    if(interact == FALSE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$M, varnames$tt,
                          # varnames$ttM,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    preds <- alist$preds
    for (jj in 1:length(val)) {
      r_mtc[folds[[v]]$validation_set, glue("r(0|m,t{val[jj]},c)")] <- 1 - preds[, jj]
      r_mtc[folds[[v]]$validation_set, glue("r(1|m,t{val[jj]},c)")] <- preds[, jj]
    }

  }
  r_mtc
}



r.zmtc <- function(data_in, varnames, interact, folds, learners, bounded = TRUE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  r_zmtc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(r_zmtc) <- (glue("r({tr_vals$r}|z,m,t{tr_vals$tt},c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])
    val <- c(0, 1)
    valid_list <- lapply(val, function(tt=0) {
      valid_v <- valid
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
    if (interact == TRUE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$Z, varnames$M, varnames$tt,
                          varnames$ttM, varnames$ttZ, varnames$MZ,
                          varnames$ttMZ,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    if (interact == FALSE) {
      alist <- crossfit(train, valid_list,
                        varnames$R,
                        c(varnames$Z, varnames$M, varnames$tt,
                          # varnames$ttM, varnames$ttZ, varnames$MZ,
                          # varnames$ttMZ,
                          varnames$C), varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    preds <- alist$preds
    for (jj in 1:length(val)) {
      r_zmtc[folds[[v]]$validation_set, glue("r(0|z,m,t{val[jj]},c)")] <- 1 - preds[, jj]
      r_zmtc[folds[[v]]$validation_set, glue("r(1|z,m,t{val[jj]},c)")] <- preds[, jj]
    }

  }
  r_zmtc
}


# treatment
t.rc <- function(data_in, varnames, folds, learners, bounded = TRUE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  t_rc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(t_rc) <- (glue("t({tr_vals$tt}|r{tr_vals$r},c)"))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    val <- c(0, 1)
    valid_list <- lapply(val, function(r=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- r
      valid_v
    })
    alist <- crossfit(train, valid_list,
                      varnames$tt,
                      c(varnames$R, varnames$C), varnames,
                      type = c("binomial"), learners, "ID", bounded)
    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    for (rr in 1:length(val)) {
      t_rc[folds[[v]]$validation_set, glue("t(0|r{val[rr]},c)")] <- 1 - preds[, rr]
      t_rc[folds[[v]]$validation_set, glue("t(1|r{val[rr]},c)")] <- preds[, rr]
    }
  }
  t_rc
}


# t.rzmc <- function(data_in, varnames, folds, learners, bounded = TRUE) {
#   tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
#   t_rzmc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
#   colnames(t_rzmc) <- glue("t({tr_vals$tt}|r{tr_vals$r},z,m,c)")
#   v <- 1
#   for (v in seq_along(folds)) {
#     train <- origami::training(data_in, folds[[v]])
#     valid <- origami::validation(data_in, folds[[v]])
#     val <- c(0, 1)
#     valid_list <- lapply(val, function(r=0) {
#       valid_v <- valid
#       valid_v[, c(varnames$R)] <- r
#       # update the interactions
#       # two-way
#       valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
#       valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
#       valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
#       valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
#       valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
#       valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
#       # three-way and four-way
#       valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
#       valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
#       valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
#       valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
#       valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
#       valid_v
#     })
#     if (interact == TRUE) {
#       alist <- crossfit(train, valid_list,
#                         varnames$tt,
#                         c(varnames$Z, varnames$M, varnames$R,
#                           varnames$RM, varnames$RZ, varnames$MZ,
#                           varnames$RMZ,
#                           varnames$C), varnames,
#                         type = c("binomial"), learners, "ID", bounded)
#     }
#     if (interact == FALSE) {
#       alist <- crossfit(train, valid_list,
#                         varnames$tt,
#                         c(varnames$Z, varnames$M, varnames$R,
#                           # varnames$ttM, varnames$ttZ, varnames$MZ,
#                           # varnames$ttMZ,
#                           varnames$C), varnames,
#                         type = c("binomial"), learners, "ID", bounded)
#     }
#     preds <- alist$preds
#     for (rr in 1:length(val)) {
#       t_rzmc[folds[[v]]$validation_set, glue("t(0|r{val[rr]},z,m,c)")] <- 1 - preds[, rr]
#       t_rzmc[folds[[v]]$validation_set, glue("t(1|r{val[rr]},z,m,c)")] <- preds[, rr]
#     }
#   }
#   t_rzmc
# }



# p(T, R| C)
tr.c <- function(r_c, t_rc) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  tr_c <- matrix(nrow = nrow(r_c), ncol = nrow(tr_vals))
  colnames(tr_c) <- glue("p(t{tr_vals$tt},r{tr_vals$r}|c)")
  for (jj in 1:nrow(tr_vals)) {
    tt <- tr_vals$tt[jj]
    r <- tr_vals$r[jj]
    tr_c[, glue("p(t{tt},r{r}|c)")] <-  r_c[, glue("r({r}|c)")] * t_rc[, glue("t({tt}|r{r},c)")]
  }
  tr_c
}

