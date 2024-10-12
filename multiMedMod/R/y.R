
# conditional outcome -------------------
# E(y|m1,m2,t,r,c) for intervened
y.mzc <- function(data_in, varnames, Yfamily = "gaussian", interact = TRUE, folds, learners, bounded = FALSE) {
  mztr_vals <- expand.grid(m = c(0,1), z = c(0,1), tt=c(0,1), r=c(0,1))
  mtr_vals <- expand.grid(m = c(0,1), z = NA, tt=c(0,1), r=c(0,1))
  ztr_vals <- expand.grid(m=NA, z = c(0,1), tt=c(0,1), r=c(0,1))
  tr_vals <- expand.grid(m=NA, z=NA, tt=c(0,1), r=c(0,1))
  NA_vals <- expand.grid(m=NA, z=NA, tt=NA, r=NA)

  namevec <- c(glue("mu(z{mztr_vals$z},m{mztr_vals$m},t{mztr_vals$tt},r{mztr_vals$r},c)"),
               glue("mu(Z,m{mtr_vals$m},t{mtr_vals$tt},r{mtr_vals$r},c)"),
               glue("mu(z{ztr_vals$z},M,t{ztr_vals$tt},r{ztr_vals$r},c)"),
               glue("mu(Z,M,t{tr_vals$tt},r{tr_vals$r},c)"),
               glue("mu(Z,M,T,R,c)")
  )

  mu <- matrix(nrow = nrow(data_in), ncol = length(namevec))
  colnames(mu) <- namevec
  vals <- do.call(rbind, list(mztr_vals, mtr_vals, ztr_vals, tr_vals, NA_vals))

  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(vals), function(jj=1) {
      valid_v <- valid
      if ( !is.na(vals$r[jj]) ) { valid_v[, c(varnames$R)] <- vals$r[jj] }
      if ( !is.na(vals$tt[jj]) ) { valid_v[, c(varnames$tt)] <- vals$tt[jj] }
      if ( !is.na(vals$m[jj]) ) { valid_v[, c(varnames$M)] <- vals$m[jj] }
      if ( !is.na(vals$z[jj]) ) { valid_v[, c(varnames$Z)] <- vals$z[jj] }

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

    if (interact==TRUE) {
      alist <- crossfit(
        train, valid_list,
        varnames$Y,
        c(varnames$tt, varnames$R, varnames$M, varnames$Z,
          # two-way
          varnames$Rtt, varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
          # three-way and four-way
          varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ
          ,varnames$C), varnames,
        type = Yfamily, learners, "ID", bounded)
    }
    if (interact==FALSE) {
      alist <- crossfit(
        train, valid_list,
        varnames$Y,
        c(varnames$tt, varnames$R, varnames$M, varnames$Z,
          # only Rtt
          # varnames$Rtt,
          varnames$C), varnames,
        type = Yfamily, learners, "ID", bounded)
    }

    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    mu[folds[[v]]$validation_set, ]  <- preds

  }
  y_mzc <- mu
  y_mzc
}


# E(y|t,r,c) for natural (i.e., no mediator intervention)
y.c <- function(data_in, varnames, Yfamily = "gaussian", folds, learners, bounded = FALSE) {

  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  NA_vals <- expand.grid(tt=NA, r=NA)

  namevec <- c(glue("mu(t{tr_vals$tt},r{tr_vals$r},c)"),
               glue("mu(T,R,c)")
  )

  mu <- matrix(nrow = nrow(data_in), ncol = length(namevec))
  colnames(mu) <- namevec
  vals <- do.call(rbind, list(tr_vals, NA_vals))

  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(vals), function(jj=1) {
      valid_v <- valid
      if ( !is.na(vals$r[jj]) ) { valid_v[, c(varnames$R)] <- vals$r[jj] }
      if ( !is.na(vals$tt[jj]) ) { valid_v[, c(varnames$tt)] <- vals$tt[jj] }

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


    alist <- crossfit(
      train, valid_list,
      varnames$Y,
      c(varnames$tt, varnames$R, varnames$Rtt
        ,varnames$C), varnames,
      type = Yfamily, learners, "ID", bounded)

    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    mu[folds[[v]]$validation_set, ]  <- preds

  }
  y_c <- mu
  y_c
}


# partially marginalized outcome ----

# integrate over p(m2 | t,r**,c)

y.zc <- function(y_mzc, tt, r0, m_c, r2) {
  mvals <- c(0,1)
  y_zc <- rowSums( y_mzc[, glue("mu(Z,m{mvals},t{tt},r{r0},c)")] *
    m_c[, glue("m({mvals}|t{tt},r{r2},c)")] )

  y_zc
}


# integrate over p(m1 | t,r*,c)

y.mc <- function(y_mzc, tt, r0, z_c, r1) {
  zvals <- c(0,1)
  y_mc <- rowSums( y_mzc[, glue("mu(z{zvals},M,t{tt},r{r0},c)")] *
                     z_c[, glue("z({zvals}|t{tt},r{r1},c)")] )

  y_mc
}


# fully marginalized outcome ----

# integrate over p(m1 | t,r*,c)p(m2 | t,r**,c)

y.c_12 <- function(y_mzc, tt, r0, z_c, r1, m_c, r2) {
  zmvals <- expand.grid(z = c(0, 1), m = c(0, 1))
  y_c_12 <- rowSums( y_mzc[, glue("mu(z{zmvals$z},m{zmvals$m},t{tt},r{r0},c)")] *
                     z_c[, glue("z({zmvals$z}|t{tt},r{r1},c)")] *
                     m_c[, glue("m({zmvals$m}|t{tt},r{r2},c)")]  )

  y_c_12
}

# integrate over p(m1, m2 | t,rjo,c)

y.c_jo <- function(y_mzc, tt, r0, mz_c, rjo) {
  zmvals <- expand.grid(z = c(0, 1), m = c(0, 1))
  y_c_jo <- rowSums( y_mzc[, glue("mu(z{zmvals$z},m{zmvals$m},t{tt},r{r0},c)")] *
                    mz_c[, glue("mz(m{zmvals$m},z{zmvals$z}|t{tt},r{rjo},c)")]   )

  y_c_jo
}


# ...........--------------------------------------------
# no_tt = TRUE ----------------------------------
y.c_t0 <- function(data_in, varnames, Yfamily = "gaussian", folds, learners, bounded = FALSE) {

  tr_vals <- expand.grid(tt=c(0), r=c(0,1))
  NA_vals <- expand.grid(tt=NA, r=NA)

  namevec <- c(glue("mu(t{tr_vals$tt},r{tr_vals$r},c)"),
               glue("mu(T,R,c)")
  )

  mu <- matrix(nrow = nrow(data_in), ncol = length(namevec))
  colnames(mu) <- namevec
  vals <- do.call(rbind, list(tr_vals, NA_vals))

  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(vals), function(jj=1) {
      valid_v <- valid
      if ( !is.na(vals$r[jj]) ) { valid_v[, c(varnames$R)] <- vals$r[jj] }
      if ( !is.na(vals$tt[jj]) ) { valid_v[, c(varnames$tt)] <- vals$tt[jj] }

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


    alist <- crossfit(
      train, valid_list,
      varnames$Y,
      c(varnames$R
        ,varnames$C), varnames,
      type = Yfamily, learners, "ID", bounded)

    # alist$fit$fitLibrary$SL.glm_All$object$coefficients
    preds <- alist$preds
    mu[folds[[v]]$validation_set, ]  <- preds

  }
  y_c <- mu
  y_c
}
