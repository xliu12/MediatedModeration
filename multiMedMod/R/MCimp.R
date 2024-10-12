

par.MedMod <- function(data_in,
                       Yname,
                       Rname,
                       ttname,
                       Mediators,
                       Cnames,
                       Yfamily = "gaussian", linear.Y = TRUE,
                       TotMod_dr = FALSE) {

  Znames <- Mediators[[1]]
  Mnames <- Mediators[[2]]
  # add interactions among the key variables ----
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

  # ttMnames <- paste0("ttM.", Mnames)
  # RMnames <- paste0("RM.", Mnames)
  # ttRMnames <- paste0("ttRM.", Mnames)
  # colnames(data_in)[(ncol(data_in)+1-3*length(Mnames)): ncol(data_in)] <- c(ttMnames, RMnames, ttRMnames)

  varnames <- list("R" = Rname, "tt" = ttname, "M" = Mnames, "Y" = Yname,
                   "Z" = Znames, "C" = Cnames,
                   "Rtt" = paste0("Rtt"), "RM" = paste0("RM"), "ttM" = paste0("ttM"), "ttRM" = paste0("ttRM"),
                   "RZ" = paste0("RZ"), "ttZ" = paste0("ttZ"), "MZ" = paste0("MZ"),
                   "ttRZ" = "ttRZ", "ttMZ" = "ttMZ", "RMZ" = "RMZ", "ttRMZ" = "ttRMZ")

  # fit models ----
  train <- data_in
  valid.list <- list(data_in)

  z_mar <- glmfitting(train, valid.list,
                    varnames$Z,
                    c(varnames$R, varnames$tt,
                      varnames$Rtt,
                      varnames$C),
                    varnames,
                    type = "gaussian", bounded=FALSE)

  m_mar <- glmfitting(train, valid.list,
                      varnames$M,
                      c(varnames$R, varnames$tt,
                        varnames$Rtt,
                        varnames$C),
                      varnames,
                      type = "gaussian", bounded=FALSE)

  m_zc <- glmfitting(train, valid.list,
                      varnames$M,
                      c(varnames$R, varnames$tt, varnames$Z,
                        varnames$Rtt, varnames$RZ, varnames$ttZ, varnames$ttRZ,
                        varnames$C),
                      varnames,
                      type = "gaussian", bounded=FALSE)

  y_mzc <- glmfitting(train, valid.list,
                      varnames$Y,
                      c(varnames$tt, varnames$R, varnames$M, varnames$Z,
                        # two-way
                        varnames$Rtt, varnames$RM, varnames$ttM, varnames$RZ, varnames$ttZ, varnames$MZ,
                        # three-way and four-way
                        varnames$ttRM, varnames$ttRZ, varnames$ttMZ, varnames$RMZ, varnames$ttRMZ
                        ,varnames$C),
                      varnames,
                      type = Yfamily, bounded=FALSE)
  r_c <- glmfitting(train, valid.list,
                    varnames$R,
                    c(varnames$C),
                    varnames,
                    type = "binomial", bounded=TRUE)
  t_rc <- glmfitting(train, valid.list,
                     varnames$tt,
                     c(varnames$R, varnames$C),
                     varnames,
                     type = "binomial", bounded=TRUE)


  # compute -----
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
  # theta(tt,1,1) = delta(tt, 1)

  theta <- eif <- reg <- rmpw <- list()

  j <- 1
  for (j in 1:nrow(vals)) {
    tt <- vals$tt[j]
    r0 <- vals$r0[j]
    r1 <- vals$r1[j]
    r2 <- vals$r2[j]
    rjo <- vals$rjo[j]

    valid_obs <- data_in
    mu <- predict(y_mzc$fit, valid_obs)
    # observed (joint) mediator distribution ----
    pZM <- dnorm(M, m_zc$preds, sigma(m_zc$fit)) * dnorm(Z, z_mar$preds, sigma(z_mar$fit))

    # p_Mjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{r0},c)")] +
    #   (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{r0},c)")] +
    #   (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{r0},c)")] +
    #   (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{r0},c)")]
    pTR <- (TT==1)*(R==1) * t_rc$preds * r_c$preds + (TT==1)*(R==0) * t_rc$preds * (1 - r_c$preds) +
      (TT==0)*(R==1) * (1 - t_rc$preds) * r_c$preds + (TT==0)*(R==0) * (1 - t_rc$preds) * (1 - r_c$preds)
    ipw_y <- 1*(TT==tt)*(R==r0) / bound(pTR)

    if (!is.na(r1)) {
      # intervened mediator distribution
      valid_r1 <- valid_obs;
      valid_r1[, c(varnames$R)] <- r1
      valid_r1[, c(varnames$tt)] <- tt
      valid_r1 <- update.interaction(valid_r1, varnames)
      EZ_mar <- predict(z_mar$fit, valid_r1, type = "response")
      pZ_mar <- dnorm(Z, EZ_mar, sigma(z_mar$fit))

      valid_r2 <- valid_obs;
      valid_r2[, c(varnames$R)] <- r2
      valid_r2[, c(varnames$tt)] <- tt
      valid_r2 <- update.interaction(valid_r2, varnames)
      EM_mar <- predict(m_mar$fit, valid_r2, type = "response")
      pM_mar <- dnorm(M, EM_mar, sigma(m_mar$fit))

      # p_M1 <- (Z==1)*z_c[, glue("z(1|t{tt},r{r1},c)")] + (Z==0)*z_c[, glue("z(0|t{tt},r{r1},c)")]

      # p_M2 <- (M==1)*m_c[, glue("m(1|t{tt},r{r2},c)")] + (M==0)*m_c[, glue("m(0|t{tt},r{r2},c)")]

      # partially marginalized outcome -----
      if (linear.Y) {
        valid_r0EZ <- valid_obs;
        valid_r0EZ[, c(varnames$R)] <- r0
        valid_r0EZ[, c(varnames$tt)] <- tt
        valid_r0EZ[, c(varnames$Z)] <- EZ_mar
        valid_r0EZ <- update.interaction(valid_r0EZ, varnames)
        yZ_mc <- predict(y_mzc$fit, newdata = valid_r0EZ, type = "response")


        valid_r0EM <- valid_obs;
        valid_r0EM[, c(varnames$R)] <- r0
        valid_r0EM[, c(varnames$tt)] <- tt
        valid_r0EM[, c(varnames$M)] <- EM_mar
        valid_r0EM <- update.interaction(valid_r0EM, varnames)
        yM_zc <- predict(y_mzc$fit, newdata = valid_r0EM, type = "response")


        valid_r0EZEM <- valid_obs;
        valid_r0EZEM[, c(varnames$R)] <- r0
        valid_r0EZEM[, c(varnames$tt)] <- tt
        valid_r0EZEM[, c(varnames$Z)] <- EZ_mar
        valid_r0EZEM[, c(varnames$M)] <- EM_mar
        valid_r0EZEM <- update.interaction(valid_r0EZEM, varnames)
        yZM_c <- predict(y_mzc$fit, newdata = valid_r0EZEM, type = "response")
      }

      if (!linear.Y) {
        EZ_mar_jo <- predict(z_mar$fit, valid_jo, type = "response")

        one.k <- function(k=1) {
          set.seed(k)
          valid_r0Z <- valid_obs;
          valid_r0Z[, c(varnames$R)] <- r0
          valid_r0Z[, c(varnames$tt)] <- tt
          valid_r0Z[, c(varnames$Z)] <- rnorm(length(EZ_mar), EZ_mar, sigma(z_mar$fit))
          valid_r0Z <- update.interaction(valid_r0Z, varnames)

          yZ_mc.k <- predict(y_mzc$fit, newdata = valid_r0Z, type = "response")
          yZ_mc.k
        }
        yZ_mc <- rowMeans( do.call(cbind, lapply(1:500, one.k)) )

        one.k <- function(k=1) {
          set.seed(k)
          valid_r0M <- valid_obs;
          valid_r0M[, c(varnames$R)] <- r0
          valid_r0M[, c(varnames$tt)] <- tt
          valid_r0M[, c(varnames$M)] <- rnorm(length(EM_mar), EM_mar, sigma(m_mar$fit))
          valid_r0M <- update.interaction(valid_r0M, varnames)
          yM_zc.k <- predict(y_mzc$fit, newdata = valid_r0M, type = "response")
          yM_zc.k
        }
        yM_zc <- rowMeans( do.call(cbind, lapply(1:500, one.k)) )

      }
      # fot the effect via M1 or M2 alone

      rzmipw_y <- ipw_y*1*(TT==tt)*(R==r0)* (pZ_mar * pM_mar)/bound(pZM)
      eify <- rzmipw_y / mean(rzmipw_y) * (Y - mu)

      ipw_z <- 1*(TT==tt)*(R==r1) / bound(pTR)
      eifz <- ipw_z / mean(ipw_z) * (yM_zc - yZM_c)

      ipw_m <- 1*(TT==tt)*(R==r2) / bound(pTR)
      eifm <- ipw_m / mean(ipw_m) * (yZ_mc - yZM_c)

      # eif with marginal mediator distributions
      eif_mar <- eify + eifz + eifm + yZM_c

      eif[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- eif_mar
      theta[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean(eif_mar)

      rmpw[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( rzmipw_y / mean(rzmipw_y) * Y )
      reg[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( yZM_c )

    }

    if (!is.na(rjo)) {
      valid_jo <- valid_obs
      valid_jo[, c(varnames$R)] <- rjo
      valid_jo[, c(varnames$tt)] <- tt
      valid_jo <- update.interaction(valid_jo, varnames)

      # observed (joint) mediator distribution ----
      pZM_rjo <- dnorm(M, predict(m_zc$fit, valid_jo, type = "response"), sigma(m_zc$fit)) *
        dnorm(Z, predict(z_mar$fit, valid_jo, type = "response"), sigma(z_mar$fit))

      # p_Mjo_rjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{rjo},c)")] +
      #   (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{rjo},c)")] +
      #   (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{rjo},c)")] +
      #   (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{rjo},c)")]
      if (linear.Y) {
        EZ_mar_jo <- predict(z_mar$fit, valid_jo, type = "response")
        valid_joEZ <- valid_jo
        valid_joEZ[, varnames$Z] <- EZ_mar_jo
        valid_joEZ <- update.interaction(valid_joEZ, varnames)
        EM_zc_jo <- predict(m_zc$fit, valid_joEZ, type = "response")

        valid_r0EZM <- valid_obs;
        valid_r0EZM[, c(varnames$R)] <- r0
        valid_r0EZM[, c(varnames$tt)] <- tt
        valid_r0EZM[, c(varnames$Z)] <- EZ_mar_jo
        valid_r0EZM[, c(varnames$M)] <- EM_zc_jo
        valid_r0EZM <- update.interaction(valid_r0EZM, varnames)
        yZM_c_jo <- predict(y_mzc$fit, newdata = valid_r0EZM, type = "response")
      }


      # h_Mjo <- p_Mjo_rjo / p_Mjo
      # eify_ajo <- ipw_y*h_Mjo / mean(ipw_y*h_Mjo) * (Y - mu)
      rmjo_ipw_y <- ipw_y*1*(TT==tt)*(R==r0) * (pZM_rjo) / bound(pZM)
      eify_ajo <- rmjo_ipw_y / mean(rmjo_ipw_y) * (Y - mu)

      # fixing R = r0 for the outcome theta(tt,r0,rjo)
      valid_r0 <- valid_obs;
      valid_r0[, c(varnames$R)] <- r0
      valid_r0[, c(varnames$tt)] <- tt
      valid_r0 <- update.interaction(valid_r0, varnames)
      y_mzc_r0 <- predict(y_mzc$fit, newdata = valid_r0, type = "response")

      ipw_mz <- 1*(TT==tt)*(R==rjo) / bound(pZM)
      eifmz <- ipw_mz / mean(ipw_mz) * (y_mzc_r0 - yZM_c_jo)

      # eif with the joint mediator distribution
      eif_jo <- eify_ajo + eifmz + yZM_c_jo

      eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]] <- eif_jo
      theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean(eif_jo)

      rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( rmjo_ipw_y / mean(rmjo_ipw_y) * (Y) )
      reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( yZM_c_jo )

      ## when r0=rjo, ----
      # can also direct estimating E[E(Y|t,r,c)p(c)] without mediator intervention
      if(r0==rjo) {
        if(TotMod_dr==TRUE) {
          y_c <- glmfitting(train, valid.list,
                     varnames$Y,
                     c(varnames$tt, varnames$R,
                       # two-way
                       varnames$Rtt
                       ,varnames$C),
                     varnames,
                     type = Yfamily, bounded=FALSE)

          valid_r0 <- valid_obs;
          valid_r0[, c(varnames$R)] <- r0
          valid_r0[, c(varnames$tt)] <- tt
          valid_r0 <- update.interaction(valid_r0, varnames)
          mu_c <- predict(y_c$fit, newdata = valid_r0, type = "response")

          ipw_y <- 1*(TT==tt)*(R==r0) / bound(pTR)

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
    ci2_multi = ci_effs[2, ],
    est_reg = regs_effs,
    est_rmpw = rmpw_effs)


  return(estimates)
}





update.interaction <- function(valid_v, varnames) {
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
}
