# Mediators (Z, M)
z.c <- function(data_in, varnames, interact = TRUE, folds, learners, bounded = FALSE) {
  ztr_vals <- expand.grid(z = c(0,1), tt=c(0,1), r=c(0,1))
  z_c <- matrix(nrow = nrow(data_in), ncol = nrow(ztr_vals))
  colnames(z_c) <- glue("z({ztr_vals$z}|t{ztr_vals$tt},r{ztr_vals$r},c)")

  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(tr_vals), function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- tr_vals$r[jj]
      valid_v[, c(varnames$tt)] <- tr_vals$tt[jj]
      # update the interactions
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v
    })

    if (interact == TRUE) {
      alist <- crossfit(train, valid_list,
                        varnames$Z,
                        c(varnames$R, varnames$tt,
                          varnames$Rtt,
                          varnames$C),
                        varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    if (interact == FALSE) {
      alist <- crossfit(train, valid_list,
                        varnames$Z,
                        c(varnames$R, varnames$tt,
                          # varnames$Rtt,
                          varnames$C),
                        varnames,
                        type = c("binomial"), learners, "ID", bounded)
    }
    preds <- alist$preds

    for (jj in 1:nrow(tr_vals)) {
      tt <- tr_vals$tt[jj]
      r <- tr_vals$r[jj]
      z_c[folds[[v]]$validation_set, glue("z(0|t{tt},r{r},c)")] <- 1 - preds[, jj]
      z_c[folds[[v]]$validation_set, glue("z(1|t{tt},r{r},c)")]  <- preds[, jj]
    }
  }
  z_c
}


m.zc <- function(data_in, varnames, interact = TRUE, folds, learners, bounded = FALSE) {
  mztr_vals <- expand.grid(m = c(0,1), z = c(0,1), tt=c(0,1), r=c(0,1))
  m_zc <- matrix(nrow = nrow(data_in), ncol = nrow(mztr_vals))
  colnames(m_zc) <- glue("m({mztr_vals$m}|z{mztr_vals$z},t{mztr_vals$tt},r{mztr_vals$r},c)")

  ztr_vals <- expand.grid(z=c(0,1), tt=c(0,1), r=c(0,1))
  v <- 1
  for (v in seq_along(folds)) {
    train <- origami::training(data_in, folds[[v]])
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(ztr_vals), function(jj=0) {
      valid_v <- valid
      valid_v[, c(varnames$R)] <- ztr_vals$r[jj]
      valid_v[, c(varnames$tt)] <- ztr_vals$tt[jj]
      valid_v[, c(varnames$Z)] <- ztr_vals$z[jj]
      # update the interactions involved
      grep("R", varnames, value = T) ; grep("tt", varnames, value = T)
      valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
      valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
      valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
      valid_v
    })

    if (interact==TRUE) {
      alist <- crossfit(
        train, valid_list,
        varnames$M,
        c(varnames$R, varnames$tt, varnames$Z,
          # two-way and three-way
          varnames$Rtt,
          varnames$RZ, varnames$ttZ, varnames$ttRZ,
          varnames$C),
        varnames,
        type = c("binomial"), learners, "ID", bounded)
    }
    if (interact==FALSE) {
      alist <- crossfit(
        train, valid_list,
        varnames$M,
        c(varnames$R, varnames$tt, varnames$Z,
          # only Rtt
          # varnames$Rtt,
          varnames$C),
        varnames,
        type = c("binomial"), learners, "ID", bounded)
    }

    preds <- alist$preds

    for (jj in 1:nrow(ztr_vals)) {
      tt <- ztr_vals$tt[jj]
      r <- ztr_vals$r[jj]
      z <- ztr_vals$z[jj]
      m_zc[folds[[v]]$validation_set, glue("m(0|z{z},t{tt},r{r},c)")] <- 1 - preds[, jj]
      m_zc[folds[[v]]$validation_set, glue("m(1|z{z},t{tt},r{r},c)")]  <- preds[, jj]
    }
  }
  m_zc
}

# p(m,z|t,r,c) joint mediator distribution
mz.c <- function(z_c, m_zc) {
  mztr_vals <- expand.grid(m = c(0,1), z = c(0,1), tt=c(0,1), r=c(0,1))
  mz_c <- matrix(nrow = nrow(z_c), ncol = nrow(mztr_vals))
  colnames(mz_c) <- glue("mz(m{mztr_vals$m},z{mztr_vals$z}|t{mztr_vals$tt},r{mztr_vals$r},c)")

  for(i in 1:nrow(mztr_vals)) {
    m <- mztr_vals$m[i]
    z <- mztr_vals$z[i]
    tt <- mztr_vals$tt[i]
    r <- mztr_vals$r[i]
    mz_c[, glue("mz(m{m},z{z}|t{tt},r{r},c)")] <-
      m_zc[, glue("m({m}|z{z},t{tt},r{r},c)")] * z_c[, glue("z({z}|t{tt},r{r},c)")]
  }

  mz_c
}
# p(m|t,r,c) marginal mediator distribution
m.c <- function(mz_c) {
  mtr_vals <- expand.grid(m=c(0,1),tt=c(0,1), r=c(0,1))
  m_c <- matrix(nrow = nrow(mz_c), ncol = nrow(mtr_vals))

  colnames(m_c) <- glue("m({mtr_vals$m}|t{mtr_vals$tt},r{mtr_vals$r},c)")

  z_vec <- c(0, 1)
  for (i in 1:nrow(mtr_vals)) {
    m <- mtr_vals$m[i]
    tt <- mtr_vals$tt[i]
    r <- mtr_vals$r[i]
    m_c[, glue("m({m}|t{tt},r{r},c)")] <-
      rowSums( mz_c[, glue("mz(m{m},z{z_vec}|t{tt},r{r},c)")] )
  }
  m_c
}

