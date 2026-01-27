


# Estimation with binary mediators --------------------------------------
b.MedMod <- function(
  medmod_data,
  fity.interact = TRUE,
  fitm.interact = TRUE,
  TotMod_dr = FALSE
) {
  list2env(medmod_data, envir = environment())

  r_c <- r.c(data_in, varnames, folds, learners, bounded = TRUE)
  t_rc <- t.rc(data_in, varnames, folds, learners, bounded = TRUE)
  tr_c <- tr.c(r_c, t_rc)

  z_c <- z.c(data_in, varnames, interact = fitm.interact, folds, learners, bounded = TRUE) # Z is prior to M
  if (!is.null(varnames$M)) {
    m_zc <- m.zc(data_in, varnames, interact = fitm.interact, folds, learners, bounded = TRUE)
    mz_c <- mz.c(z_c, m_zc)
    m_c <- m.c(mz_c)
  }
  if (is.null(Mnames)) {
    mz_c <- z_c
  }

  y_mzc <- y.mzc(data_in, varnames, Yfamily = Yfamily, interact = fity.interact, folds, learners, bounded = FALSE)

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

    mu <- y_mzc[, glue("mu(Z,M,T,R,c)")]
    # observed (joint) mediator distribution
    if (!is.null(varnames$M)) {
      p_Mjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{r0},c)")] +
        (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{r0},c)")] +
        (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{r0},c)")] +
        (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{r0},c)")]
    }
    if (is.null(varnames$M)) {
      p_Mjo <- (Z==1)*z_c[, glue("z(1|t{tt},r{r0},c)")] +
        (Z==0)*z_c[, glue("z(0|t{tt},r{r0},c)")]
    }


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
      if (!is.null(varnames$M)) {
        p_Mjo_rjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{rjo},c)")] +
          (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{rjo},c)")] +
          (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{rjo},c)")] +
          (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{rjo},c)")]
      }
      if (is.null(varnames$M)) {
        p_Mjo_rjo <- (Z==1)*z_c[, glue("z(1|t{tt},r{rjo},c)")] +
          (Z==0)*z_c[, glue("z(0|t{tt},r{rjo},c)")]
      }

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

  out <- mget(ls(), envir = environment())


  return(out)
}






