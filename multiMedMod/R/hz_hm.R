
h.zm <- function(data_in, varnames, Mfamily, interact, folds, learners, bounded = FALSE) {
  tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
  h_zm <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
  colnames(h_zm) <- (glue("h_zm(t{tr_vals$tt},r{tr_vals$r})"))
  v <- 1
  for (v in seq_along(folds)) {
    # train <- origami::training(data_in, folds[[v]])
    # valid <- origami::validation(data_in, folds[[v]])
    # train <- stack_data(origami::training(data_in, folds[[v]]), varnames$Mediators)
    train <- stack_data(origami::training(data_in, folds[[v]]), varnames$M, Mfamily, varnames, learners, bounded, no_tt = FALSE)
    valid <- origami::validation(data_in, folds[[v]])

    valid_list <- lapply(1:nrow(tr_vals), function(jj=0) {
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
    if (interact == TRUE) {
      p_zmc <- crossfit(train, valid_list,
                        "Lambda",
                        c(varnames$tt, varnames$R,
                          varnames$M,
                          varnames$Z,
                          # two-way
                          varnames$Rtt,
                          varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
                          # three-way and four-way
                          varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ,
                          varnames$C), varnames,
                        type = "binomial", learners,
                        IDname = "tmp_id",
                        bounded = TRUE)
      p_mc <- crossfit(train, valid_list,
                       "Lambda",
                       c(varnames$tt, varnames$R,
                         varnames$M,
                         # two-way
                         varnames$Rtt,
                         varnames$RM, varnames$ttM,
                         # three-way and four-way
                         varnames$ttRM,
                         varnames$C), varnames,
                       type = "binomial", learners,
                       IDname = "tmp_id",
                       bounded = TRUE)
    }

    if (interact == FALSE) {
      p_zmc <- crossfit(train, valid_list,
                        "Lambda",
                        c(varnames$tt, varnames$R,
                          varnames$M,
                          varnames$Z,
                          # # two-way
                          # varnames$Rtt,
                          # varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
                          # # three-way and four-way
                          # varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ,
                          varnames$C), varnames,
                        type = "binomial", learners,
                        IDname = "tmp_id",
                        bounded = TRUE)
      p_mc <- crossfit(train, valid_list,
                       "Lambda",
                       c(varnames$tt, varnames$R,
                         varnames$M,
                         # # two-way
                         # varnames$Rtt,
                         # varnames$RM, varnames$ttM,
                         # # three-way and four-way
                         # varnames$ttRM,
                         varnames$C), varnames,
                       type = "binomial", learners,
                       IDname = "tmp_id",
                       bounded = TRUE)
    }
    for (jj in 1:nrow(tr_vals)) {
      h_zm[folds[[v]]$validation_set, glue("h_zm(t{tr_vals$tt[jj]},r{tr_vals$r[jj]})")] <- p_zmc$preds[, jj] * (1 - p_mc$preds[, jj]) / ( (1 - p_zmc$preds[, jj])*p_mc$preds[, jj] )
    }
  }
  h_zm
}

h.zstar.mstarstar <- function(r, rstar, rstarstar, tt, h_zm, tr_c, r_mtc, r_ztc) {

  iorw_m <- tr_c[, glue("p(t{tt},r{r}|c)")]*r_mtc[, glue("r({rstarstar}|m,t{tt},c)")] /
    bound(tr_c[, glue("p(t{tt},r{rstarstar}|c)")]*r_mtc[, glue("r({r}|m,t{tt},c)")])

  iorw_z <- tr_c[, glue("p(t{tt},r{r}|c)")]*r_ztc[, glue("r({rstar}|z,t{tt},c)")] /
    bound(tr_c[, glue("p(t{tt},r{rstar}|c)")]*r_ztc[, glue("r({r}|z,t{tt},c)")])

  iorw_mz <- (tr_c[, glue("p(t{tt},r{r}|c)")]*r_mtc[, glue("r({rstarstar}|m,t{tt},c)")] * tr_c[, glue("p(t{tt},r{r}|c)")]*r_ztc[, glue("r({rstar}|z,t{tt},c)")]) /
    bound( (tr_c[, glue("p(t{tt},r{rstarstar}|c)")]*r_mtc[, glue("r({r}|m,t{tt},c)")])* (tr_c[, glue("p(t{tt},r{rstar}|c)")]*r_ztc[, glue("r({r}|z,t{tt},c)")]) )

  if ((rstar==r) & (rstarstar!=r)) {
    iorw_mz <- iorw_m
  }
  if ((rstar!=r) & (rstarstar==r)) {
    iorw_mz <- iorw_z
  }
  if ((rstar!=r) & (rstarstar!=r)) {
    iorw_mz <- iorw_mz
  }
  if ((rstar==r) & (rstarstar==r)) {
    iorw_mz <- 1
  }
  h_zstar_mstarstar <- h_zm[, glue("h_zm(t{tt},r{r})")] * iorw_mz
  h_zstar_mstarstar
}

h.m.starstar <- function(r, rstarstar, tt, h_zm, r_tc, r_mtc) {

  iorw_m <- r_tc[, glue("r({r}|t{tt},c)")]*r_mtc[, glue("r({rstarstar}|m,t{tt},c)")] /
    (r_tc[, glue("r({rstarstar}|t{tt},c)")]*r_mtc[, glue("r({r}|m,t{tt},c)")])

  h_mstarstar <- h_zm[, glue("h_zm(t{tt},r{r})")] * iorw_m
  h_mstarstar
}

h.z.star <- function(r, rstar, tt, h_zm, r_tc, r_ztc) {

  iorw_z <- r_tc[, glue("r({r}|t{tt},c)")]*r_ztc[, glue("r({rstar}|z,t{tt},c)")] /
    (r_tc[, glue("r({rstar}|t{tt},c)")]*r_ztc[, glue("r({r}|z,t{tt},c)")])

  h_zstar <- h_zm[, glue("h_zm(t{tt},r{r})")] * iorw_z
  h_zstar
}


h.zm.joint <- function(r, rstar, tt, tr_c, r_zmtc) {
  iorw_joint <- bound( tr_c[, glue("p(t{tt},r{r}|c)")]*r_zmtc[, glue("r({rstar}|z,m,t{tt},c)")] ) /
    bound(tr_c[, glue("p(t{tt},r{rstar}|c)")]*r_zmtc[, glue("r({r}|z,m,t{tt},c)")])

  h_zm_joint <- iorw_joint
  h_zm_joint
}

sample_M <- function(M) { # univariate
  set.seed(12345)
  vals_M <- unique(M)
  vals_M[sample.int(length(vals_M), size = length(M), replace = TRUE)]
}



stack_data <- function(data_in, Mediators, Mfamily, varnames, learners = "SL.glm", bounded=FALSE, no_tt=FALSE) {
  delta <- data_in
  for(jj in 1:length(Mediators)) {
    delta[, Mediators[jj]] <- sample_M( data_in[[Mediators[jj]]] )
  }

  out <- rbind(data_in, delta)
  out[["Lambda"]] <- rep(c(0, 1), each = nrow(data_in))
  out[["tmp_id"]] <- rep(1:nrow(data_in), times = 2)
  out
}

