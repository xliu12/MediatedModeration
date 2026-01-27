
# Estimation with non-binary mediators --------------------------------------
h.MedMod <- function(
  medmod_data,
  fity.interact = FALSE,
  fitm.interact = FALSE,
  TotMod_dr = FALSE,
  full.sample = TRUE
) {
  list2env(medmod_data, envir = environment())

  r_c <- r.c(data_in, varnames, folds, learners, bounded = TRUE)
  t_rc <- t.rc(data_in, varnames, folds, learners, bounded = TRUE)
  tr_c <- tr.c(r_c, t_rc)

  r_ztc <- r.ztc(data_in, varnames, fitm.interact, folds, learners, bounded = TRUE)
  if (!is.null(varnames$M)) {
    r_mtc <- r.mtc(data_in, varnames, fitm.interact, folds, learners, bounded = TRUE)
    r_zmtc <- r.zmtc(data_in, varnames, fitm.interact, folds, learners, bounded = TRUE)
    h_zm <- h.zm(data_in, varnames, Mfamily, fitm.interact, folds, learners_h, bounded = TRUE, full.sample = FALSE) # stratify by tt, R
  }
  if (is.null(varnames$M)) {
    r_zmtc <- r_ztc
    # colnames(r_zmtc)
  }

  mu_mzc <- mu.mzc(r=NA, tt=NA, data_in, varnames, Yfamily = Yfamily, fity.interact, folds, learners, bounded = FALSE, full.sample = TRUE)

  TT <- data_in[[varnames$tt]]
  R <- data_in[[varnames$R]]
  Z <- data_in[, varnames$Z]
  if (!is.null(varnames$M)) {
    M <- data_in[, varnames$M]
  }
  Y <- data_in[[varnames$Y]]

  vals <- rbind(
    expand.grid(tt = c(0, 1), r0 = c(1), r1 = NA, r2 = NA, rjo = c(0,1)),
    expand.grid(tt = c(0, 1), r0 = c(0), r1 = NA, r2 = NA, rjo = c(0))
  )
  if (!is.null(varnames$M)) {
    vals <- rbind(vals,
                  expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(0,1), r2 = c(0), rjo = NA),
                  expand.grid(tt = c(0, 1), r0 = c(1), r1 = c(1), r2 = c(0,1), rjo = NA)
    )
  }

  theta <- eif <- reg <- rmpw <- list()

  j <- 1
  for (j in 1:nrow(vals)) {
    tt <- vals$tt[j]
    r0 <- vals$r0[j]
    r1 <- vals$r1[j]
    r2 <- vals$r2[j]
    rjo <- vals$rjo[j]

    mu <- mu_mzc[, glue("mu(t{tt},r{r0},m,z,c)")]
    ipw_y <- 1*(TT==tt)*(R==r0) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")])

    if (!is.na(r1)) {

      h_zstar_mstarstar <- h.zstar.mstarstar(r0, r1, r2, tt, h_zm, tr_c, r_mtc, r_ztc)
      h_zstar <- h.zstar.mstarstar(r0, rstar=r1, r0, tt, h_zm, tr_c, r_mtc, r_ztc)
      h_mstarstar <- h.zstar.mstarstar(r0, r0, rstarstar=r2, tt, h_zm, tr_c, r_mtc, r_ztc)

      muM_zc <- muM.zc(r0, tt, mu_mzc, h_mstarstar, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
      muZ_mc <- muZ.mc(r0, tt, mu_mzc, h_zstar, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)

      muMZ_c <- muMZ.c(rstar=r1, tt, muM_zc, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
      muZM_c <- muZM.c(rstarstar=r2, tt, muZ_mc, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)

      rm1m2ipw_y <- 1*(TT==tt)*(R==r0)* ipw_y * h_zstar_mstarstar
      eify <- rm1m2ipw_y * (Y - mu)  / mean(rm1m2ipw_y)

      ipw_m1 <- 1*(TT==tt)*(R==r1) / bound(tr_c[, glue("p(t{tt},r{r1}|c)")])
      # eifm1 <- ipw_m1 * (muM_zc - muMZ_c) / mean(ipw_m1)
      eifm1 <- ipw_m1 * (muM_zc - muZM_c) / mean(ipw_m1)

      ipw_m2 <- 1*(TT==tt)*(R==r2) / bound(tr_c[, glue("p(t{tt},r{r2}|c)")])
      eifm2 <- ipw_m2 * (muZ_mc - muZM_c) / mean(ipw_m2)

      eif_mar <- eify + eifm1 + eifm2 + muZM_c

      eif[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- eif_mar
      theta[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean(eif_mar)

      rmpw[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( rm1m2ipw_y / mean(rm1m2ipw_y) * Y )
      reg[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( muZM_c )

    }

    if (!is.na(rjo)) {

      h_zm_joint <- h.zm.joint(r0, rstar=rjo, tt, tr_c, r_zmtc)
      mu_joint_c <- mu.joint.c(r0, rstar=rjo, tt, mu_mzc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = full.sample)

      rmjo_ipw_y <- 1*(TT==tt)*(R==r0) * ipw_y * h_zm_joint
      eify_ajo <- rmjo_ipw_y * (Y - mu) / mean(rmjo_ipw_y)

      ipw_m1m2 <- 1*(TT==tt)*(R==rjo) / bound(tr_c[, glue("p(t{tt},r{rjo}|c)")])

      eifm1m2 <- ipw_m1m2 * (mu - mu_joint_c) / mean(ipw_m1m2)

      eif_jo <- eify_ajo + eifm1m2 + mu_joint_c

      eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]] <- eif_jo
      theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean(eif_jo)

      rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( rmjo_ipw_y / mean(rmjo_ipw_y) * (Y) )
      reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( mu_joint_c )

      ## when r0=rjo,
      if(r0==rjo) {
        if(TotMod_dr==TRUE) {
          y_c <- y.c(data_in, varnames, Yfamily = Yfamily, folds, learners, bounded = FALSE)
          mu_c <- y_c[, glue("mu(t{tt},r{r0},c)")]

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

  out <- mget(ls(), envir = environment())

  return(out)
}







# effect <- function(theta) {
#
#   theta$`TotMod`  <- with(theta, {
#     `delta(t1,r1)`-`delta(t1,r0)`-`delta(t0,r1)`+`delta(t0,r0)`
#   })
#   # the (remaining) moderated treatment effect while fixing the joint mediator distribution at the reference subgroup; or, the treatment effect on the IDE of the subgroup
#   theta$`RemainMod` <- with(theta, {
#     `theta(t1,r1,rjo0)` - `delta(t1,r0)` -
#       `theta(t0,r1,rjo0)` + `delta(t0,r0)`
#   })
#   # the (joint mediated) moderated treatment effect mediated by the difference in  the (joint distribution of) potential mediators between subgroups; or, the treatment effect on the IIE_jo of the subgroup
#   theta$`MedMod_jo`  <- with(theta, {
#     `delta(t1,r1)` - `theta(t1,r1,rjo0)` -
#       `delta(t0,r1)` + `theta(t0,r1,rjo0)`
#   })
#   # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M1 between subgroups; or, the treatment effect on the IIE_M1 of the subgroup
#   theta$`MedMod_M1`  <- with(theta, {
#     `theta(t1,r1,r1,r0)` - `theta(t1,r1,r0,r0)` -
#       `theta(t0,r1,r1,r0)` + `theta(t0,r1,r0,r0)`
#   })
#   # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M2 between subgroups; or, the treatment effect on the IIE_M2 of the subgroup
#   theta$`MedMod_M2`  <- with(theta, {
#     `theta(t1,r1,r1,r1)` - `theta(t1,r1,r1,r0)` -
#       `theta(t0,r1,r1,r1)` + `theta(t0,r1,r1,r0)`
#   })
#   # the (mediated) moderated treatment effect mediated by the difference in the mutual dependence of potential mediators M1 and M2 between subgroups; or, the treatment effect on the IIE_mu of the subgroup
#   theta$`MedMod_mu`  <- with(theta, {
#     `delta(t1,r1)`-`theta(t1,r1,rjo0)`-`theta(t1,r1,r1,r1)`+`theta(t1,r1,r0,r0)` -
#       (`delta(t0,r1)`-`theta(t0,r1,rjo0)`-`theta(t0,r1,r1,r1)`+`theta(t0,r1,r0,r0)`)
#   })
#
#
#   # theta$`MedMod_M1`+theta$`MedMod_M2`+theta$`MedMod_mu` - theta$`MedMod_jo`
#
#   theta
# }


