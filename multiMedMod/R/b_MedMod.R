


# Estimation with binary mediators --------------------------------------
b.MedMod <- function(
  data_in,
  Yname, M1names, M2names, ttname, Rname,
  Cnames,
  learners,
  num_folds,
  Yfamily,
  fity.interact = TRUE,
  fitm.interact = TRUE,
  TotMod_dr = FALSE,
  out_eif = FALSE
) {

  Znames <- M1names #Mediators[[1]]
  Mnames <- M2names #Mediators[[2]]
  data_in <- data_in %>% mutate(
    Rtt = data_in[[Rname]] * data_in[[ttname]],
    ttM = data_in[[ttname]] * data_in[, Mnames],
    RM = data_in[[Rname]] * data_in[, Mnames],
    ttRM = data_in[[Rname]] * data_in[[ttname]] * data_in[, Mnames]
  )
  data_in <- data_in %>% mutate(
    RZ = data_in[[Rname]] * data_in[, Znames],
    ttZ = data_in[[ttname]] * data_in[, Znames],
    MZ = data_in[, Znames] * data_in[, Mnames]
  )
  data_in <- data_in %>% mutate(
    ttRZ = data_in[[Rname]] * data_in[[ttname]] * data_in[, Znames],
    ttMZ = data_in[[ttname]] * data_in[, Znames] * data_in[, Mnames],
    RMZ = data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames],
    ttRMZ = data_in[[ttname]] * data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames]
  )
  data_in <- do.call(data.frame, data_in)

  varnames <- list("R" = Rname, "tt" = ttname, "M" = Mnames, "Y" = Yname,
                   "Z" = Znames, "C" = Cnames,
                   "Rtt" = paste0("Rtt"), "RM" = paste0("RM"), "ttM" = paste0("ttM"), "ttRM" = paste0("ttRM"),
                   "RZ" = paste0("RZ"), "ttZ" = paste0("ttZ"), "MZ" = paste0("MZ"),
                   "ttRZ" = "ttRZ", "ttMZ" = "ttMZ", "RMZ" = "RMZ", "ttRMZ" = "ttRMZ")

  folds <- make_folds(n = nrow(data_in), fold_fun = folds_vfold, V = num_folds)
  if (num_folds==1) { folds[[1]]$training_set <- folds[[1]]$validation_set }


  r_c <- r.c(data_in, varnames, folds, learners, bounded = TRUE)
  t_rc <- t.rc(data_in, varnames, folds, learners, bounded = TRUE)
  tr_c <- tr.c(r_c, t_rc)

  z_c <- z.c(data_in, varnames, interact = fitm.interact, folds, learners, bounded = TRUE)
  m_zc <- m.zc(data_in, varnames, interact = fitm.interact, folds, learners, bounded = TRUE)
  mz_c <- mz.c(z_c, m_zc)
  m_c <- m.c(mz_c)

  y_mzc <- y.mzc(data_in, varnames, Yfamily = Yfamily, interact = fity.interact, folds, learners, bounded = FALSE)

  TT <- data_in[[varnames$tt]]
  R <- data_in[[varnames$R]]
  Z <- data_in[[varnames$Z]]
  M <- data_in[[varnames$M]]
  Y <- data_in[[varnames$Y]]

  vals <- rbind(
    expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(0,1), r2 = c(0), rjo = NA),
    expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(1), r2 = c(0,1), rjo = NA),
    expand.grid(tt = c(0, 1), r0 = c(1), r1 = NA, r2 = NA, rjo = c(0,1)),
    expand.grid(tt = c(0, 1), r0 = c(0), r1 = NA, r2 = NA, rjo = c(0))
  )

  theta <- eif <- reg <- rmpw <- list()

  for (j in 1:nrow(vals)) {
    tt <- vals$tt[j]
    r0 <- vals$r0[j]
    r1 <- vals$r1[j]
    r2 <- vals$r2[j]
    rjo <- vals$rjo[j]

    mu <- y_mzc[, glue("mu(Z,M,T,R,c)")]
    # observed (joint) mediator distribution
    p_Mjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{r0},c)")] +
      (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{r0},c)")] +
      (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{r0},c)")] +
      (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{r0},c)")]

    if (!is.na(r1)) {
      # intervened mediator distribution
      p_M1 <- (Z==1)*z_c[, glue("z(1|t{tt},r{r1},c)")] +
        (Z==0)*z_c[, glue("z(0|t{tt},r{r1},c)")]

      p_M2 <- (M==1)*m_c[, glue("m(1|t{tt},r{r2},c)")] +
        (M==0)*m_c[, glue("m(0|t{tt},r{r2},c)")]

      # partially marginalized outcome over p(M2|t,r2,c) (so dependent on M1)
      y_m1c <- y.zc(y_mzc, tt, r0, m_c, r2)
      # partially marginalized outcome over p(M1|t,r1,c) (so dependent on M2)
      y_m2c <- y.mc(y_mzc, tt, r0, z_c, r1)
      # fully marginalized outcome over p(M1|t,r1,c)p(M2|t,r2,c)
      y_c_12 <- y.c_12(y_mzc, tt, r0, z_c, r1, m_c, r2)


      rm1m2ipw_y <- 1*(TT==tt)*(R==r0)*(p_M1 * p_M2) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")]*p_Mjo)
      eify <- rm1m2ipw_y / mean(rm1m2ipw_y) * (Y - mu)

      ipw_m1 <- 1*(TT==tt)*(R==r1) / bound(tr_c[, glue("p(t{tt},r{r1}|c)")])
      eifm1 <- ipw_m1 / mean(ipw_m1) * (y_m1c - y_c_12)

      ipw_m2 <- 1*(TT==tt)*(R==r2) / bound(tr_c[, glue("p(t{tt},r{r2}|c)")])
      eifm2 <- ipw_m2 / mean(ipw_m2) * (y_m2c - y_c_12)

      eif_mar <- eify + eifm1 + eifm2 + y_c_12

      eif[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- eif_mar
      theta[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean(eif_mar)

      rmpw[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( rm1m2ipw_y / mean(rm1m2ipw_y) * Y )
      reg[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( y_c_12 )

    }

    if (!is.na(rjo)) {
      p_Mjo_rjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{rjo},c)")] +
        (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{rjo},c)")] +
        (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{rjo},c)")] +
        (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{rjo},c)")]

      y_c_jo <- y.c_jo(y_mzc, tt, r0, mz_c, rjo)

      rmjo_ipw_y <- 1*(TT==tt)*(R==r0)*(p_Mjo_rjo) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")]*p_Mjo)
      eify_ajo <- rmjo_ipw_y / mean(rmjo_ipw_y) * (Y - mu)

      y_m1m2c_r0 <- y_mzc[, glue("mu(Z,M,t{tt},r{r0},c)")]

      ipw_m1m2 <- 1*(TT==tt)*(R==rjo) / bound(tr_c[, glue("p(t{tt},r{rjo}|c)")])

      eifm1m2 <- ipw_m1m2 / mean(ipw_m1m2) * (y_m1m2c_r0 - y_c_jo)

      eif_jo <- eify_ajo + eifm1m2 + y_c_jo

      eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]] <- eif_jo
      theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean(eif_jo)

      rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( rmjo_ipw_y / mean(rmjo_ipw_y) * (Y) )
      reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( y_c_jo )

      ## when r0=rjo,
      if(r0==rjo) {
        if(TotMod_dr==TRUE) {
          y_c <- y.c(data_in, varnames, Yfamily = Yfamily, folds, learners, bounded = FALSE)
          mu_c <- y_c[, glue("mu(t{tt},r{r0},c)")]
          ipw_y <- 1*(TT==tt)*(R==r0) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")])

          eify_c <- ipw_y / mean(ipw_y) * (Y - mu_c) + mu_c

          eif[[ glue("delta(t{tt},r{r0})") ]] <- eify_c
          theta[[ glue("delta(t{tt},r{r0})")  ]] <- mean(eify_c)

          rmpw[[ glue("delta(t{tt},r{r0})")  ]] <- mean( ipw_y / mean(ipw_y) * (Y) )
          reg[[ glue("delta(t{tt},r{r0})")  ]] <- mean( mu_c )
        }
        if (TotMod_dr==FALSE) {
          eif[[ glue("delta(t{tt},r{r0})") ]] <- eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]]
          theta[[ glue("delta(t{tt},r{r0})")  ]] <- theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
          rmpw[[ glue("delta(t{tt},r{r0})")  ]] <- rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
          reg[[ glue("delta(t{tt},r{r0})")  ]] <- reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
        }
      }

    }
  }

  # multiply-robust
  eifs_effs <- effect(eif)
  thetas_effs <- sapply(eifs_effs, mean)
  se_effs <- sapply(eifs_effs, function(s) {
    sqrt(var(s) / nrow(data_in))
  })
  ci_effs <- sapply(eifs_effs, function(s) {
    mean(s) + c(-1, 1) * qnorm(0.975) * sqrt(var(s) / nrow(data_in))
  })

  # regression
  regs_effs <- unlist(effect(reg))
  # weighting
  rmpw_effs <- unlist(effect(rmpw))

  estimates <- cbind(
    est_multi = thetas_effs,
    se_multi = se_effs,
    ci1_multi = ci_effs[1, ],
    ci2_multi = ci_effs[2, ]
    ,est_reg = regs_effs,
    est_rmpw = rmpw_effs
    )

  if (out_eif == TRUE) {
    out <- list(eifs_effs = eifs_effs, estimates = estimates)
  }
  if (out_eif == FALSE) {
    out <- estimates
  }
  return(out)
}







effect <- function(theta) {

  theta$`TotMod`  <- with(theta, {
    `delta(t1,r1)`-`delta(t1,r0)`-`delta(t0,r1)`+`delta(t0,r0)`
  })
  # the (remaining) moderated treatment effect while fixing the joint mediator distribution at the reference subgroup; or, the treatment effect on the IDE of the subgroup
  theta$`RemainMod` <- with(theta, {
    `theta(t1,r1,rjo0)` - `delta(t1,r0)` -
      `theta(t0,r1,rjo0)` + `delta(t0,r0)`
  })
  # the (joint mediated) moderated treatment effect mediated by the difference in  the (joint distribution of) potential mediators between subgroups; or, the treatment effect on the IIE_jo of the subgroup
  theta$`MedMod_jo`  <- with(theta, {
    `delta(t1,r1)` - `theta(t1,r1,rjo0)` -
      `delta(t0,r1)` + `theta(t0,r1,rjo0)`
  })
  # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M1 between subgroups; or, the treatment effect on the IIE_M1 of the subgroup
  theta$`MedMod_M1`  <- with(theta, {
    `theta(t1,r1,r1,r0)` - `theta(t1,r1,r0,r0)` -
      `theta(t0,r1,r1,r0)` + `theta(t0,r1,r0,r0)`
  })
  # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M2 between subgroups; or, the treatment effect on the IIE_M2 of the subgroup
  theta$`MedMod_M2`  <- with(theta, {
    `theta(t1,r1,r1,r1)` - `theta(t1,r1,r1,r0)` -
      `theta(t0,r1,r1,r1)` + `theta(t0,r1,r1,r0)`
  })
  # the (mediated) moderated treatment effect mediated by the difference in the mutual dependence of potential mediators M1 and M2 between subgroups; or, the treatment effect on the IIE_mu of the subgroup
  theta$`MedMod_mu`  <- with(theta, {
    `delta(t1,r1)`-`theta(t1,r1,r1,r1)`-`theta(t1,r1,rjo0)`+`theta(t1,r1,r0,r0)` -
      (`delta(t0,r1)`-`theta(t0,r1,r1,r1)`-`theta(t0,r1,rjo0)`+`theta(t0,r1,r0,r0)`)
  })



  theta
}

