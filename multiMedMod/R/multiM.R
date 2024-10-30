#
#
# # Estimating the effects
# multiM <- function(
#     data_in,
#     Yname, Rname, ttname,
#     Mediators,
#     Cnames,
#     learners = c("SL.glm"),
#     learners_h,
#     learners_mu,
#     Yfamily = "binomial",
#     Mfamily,
#     fity.interact = F,
#     fitm.interact = F,
#     num_folds = 1,
#     Tot_dr = FALSE,
#     full.sample = TRUE,
#     no_tt = TRUE
# ) {
#
#   Znames <- Mediators[[1]]
#   Mnames <- Mediators[[2]]
#   # add interactions among the key variables
#   data_in <- data_in %>% mutate(
#     Rtt = data_in[[Rname]] * data_in[[ttname]]
#     , ttM = data_in[[ttname]] * data_in[, Mnames],
#     RM = data_in[[Rname]] * data_in[, Mnames],
#     ttRM = data_in[[Rname]] * data_in[[ttname]] * data_in[, Mnames]
#   )
#   data_in <- data_in %>% mutate(
#     RZ = data_in[[Rname]] * data_in[, Znames],
#     ttZ = data_in[[ttname]] * data_in[, Znames],
#     MZ = data_in[, Znames] * data_in[, Mnames]
#   )
#   data_in <- data_in %>% mutate(
#     ttRZ = data_in[[Rname]] * data_in[[ttname]] * data_in[, Znames],
#     ttMZ = data_in[[ttname]] * data_in[, Znames] * data_in[, Mnames],
#     RMZ = data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames],
#     ttRMZ = data_in[[ttname]] * data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames]
#   )
#   data_in <- do.call(data.frame, data_in)
#
#   # ttMnames <- paste0("ttM.", Mnames)
#   # RMnames <- paste0("RM.", Mnames)
#   # ttRMnames <- paste0("ttRM.", Mnames)
#   # colnames(data_in)[(ncol(data_in)+1-3*length(Mnames)): ncol(data_in)] <- c(ttMnames, RMnames, ttRMnames)
#
#   varnames <- list("R" = Rname, "tt" = ttname, "M" = Mnames, "Y" = Yname,
#                    "Z" = Znames, "C" = Cnames,
#                    "Rtt" = paste0("Rtt"), "RM" = paste0("RM"), "ttM" = paste0("ttM"), "ttRM" = paste0("ttRM"),
#                    "RZ" = paste0("RZ"), "ttZ" = paste0("ttZ"), "MZ" = paste0("MZ"),
#                    "ttRZ" = "ttRZ", "ttMZ" = "ttMZ", "RMZ" = "RMZ", "ttRMZ" = "ttRMZ")
#
#   folds <- make_folds(n = nrow(data_in), fold_fun = folds_vfold, V = num_folds)
#   if (num_folds==1) { folds[[1]]$training_set <- folds[[1]]$validation_set }
#
#
#   r_c <- r.c(data_in, varnames, folds, learners, bounded = TRUE) # later in the denominator
#   # t_rc <- t.rc(data_in, varnames, folds, learners, bounded = TRUE)
#   # tr_c <- r_c
#
#   # r_tc <- r.tc(data_in, varnames, folds, learners, bounded = TRUE)
#   r_zc <- r.zc_t0(data_in, varnames, folds, learners, bounded = TRUE)
#   r_mc <- r.mc_t0(data_in, varnames, folds, learners, bounded = TRUE)
#   r_zmc <- r.zmc_t0(data_in, varnames, fitm.interact, folds, learners, bounded = TRUE)
#
#   h_zm <- h.zm_t0(data_in, varnames, Mfamily, fitm.interact, folds, learners_h, bounded = TRUE)
#
#   # y_mzc <- y.mzc(data_in, varnames, Yfamily = Yfamily, interact = fity.interact, folds, learners, bounded = FALSE) # the same as y.mzc() used with zm.R
#   mu_mzc <- mu.mzc_t0(r=NA, data_in, varnames, Yfamily = Yfamily, fity.interact, folds, learners, bounded = FALSE, full.sample = TRUE)
#
#   # tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
#   # mu_mzc <- matrix(nrow = nrow(data_in), ncol = nrow(tr_vals))
#   # colnames(mu_mzc) <- (glue("mu(t{tr_vals$tt},r{tr_vals$r},m,z,c)"))
#   # for (jj in 1:nrow(tr_vals)) {
#   #   tt <- tr_vals$tt[jj]
#   #   r <- tr_vals$r[jj]
#   #   mu_mzc[, glue("mu(t{tt},r{r},m,z,c)")] <- mu.mzc(r, tt, data_in, varnames, Yfamily = "gaussian", folds, learners, bounded = FALSE, full.sample = FALSE)
#   # }
#
#   # TT <- data_in[[varnames$tt]]
#   R <- data_in[[varnames$R]]
#   Z <- data_in[, varnames$Z]
#   M <- data_in[, varnames$M]
#   Y <- data_in[[varnames$Y]]
#
#   vals <- rbind(
#     expand.grid(tt = c(0), r0 = c(1), r1 = c(0,1), r2 = c(0), rjo = NA),
#     expand.grid(tt = c(0), r0 = c(1), r1 = c(1), r2 = c(0,1), rjo = NA),
#     expand.grid(tt = c(0), r0 = c(1), r1 = NA, r2 = NA, rjo = c(0,1)),
#     expand.grid(tt = c(0), r0 = c(0), r1 = NA, r2 = NA, rjo = c(0))
#   )
#   # theta(tt,1,1) = delta(tt, 1)
#
#   theta <- eif <- reg <- rmpw <- list()
#
#
#   j <- 1
#   for (j in 1:nrow(vals)) {
#     tt <- 0 # vals$tt[j]
#     r0 <- vals$r0[j]
#     r1 <- vals$r1[j]
#     r2 <- vals$r2[j]
#     rjo <- vals$rjo[j]
#
#     mu <- mu_mzc[, glue("mu(t{tt},r{r0},m,z,c)")]
#     # mu_mzc <- mu.mzc(r0, tt, data_in, varnames, Yfamily = Yfamily, folds, learners, bounded = FALSE, full.sample = full.sample)
#     # mu <- mu_mzc
#
#     ipw_y <- 1*(R==r0) / bound(r_c[, glue("r({r0}|c)")])
#
#     # mu <- y_mzc[, glue("mu(Z,M,T,R,c)")]
#     # observed (joint) mediator distribution
#     # p_Mjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{r0},c)")] +
#     #   (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{r0},c)")] +
#     #   (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{r0},c)")] +
#     #   (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{r0},c)")]
#
#     if (!is.na(r1)) {
#       # intervened mediator distribution
#
#       # h_zstar_mstarstar <- h.zstar.mstarstar(r0, r1, r2, tt, h_zm, r_tc, r_mtc, r_ztc)
#       # h_zstar <- h.zstar.mstarstar(r0, rstar=r1, r0, tt, h_zm, r_tc, r_mtc, r_ztc)
#       # h_mstarstar <- h.zstar.mstarstar(r0, r0, rstarstar=r2, tt, h_zm, r_tc, r_mtc, r_ztc)
#
#       h_zstar_mstarstar <- h.zstar.mstarstar_t0(r0, r1, r2, h_zm, r_c, r_mc, r_zc)
#       h_zstar <- h.zstar.mstarstar_t0(r0, rstar=r1, r0, h_zm, r_c, r_mc, r_zc)
#       h_mstarstar <- h.zstar.mstarstar_t0(r0, r0, rstarstar=r2, h_zm, r_c, r_mc, r_zc)
#
#       muM_zc <- muM.zc_t0(r0, mu_mzc, h_mstarstar, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
#       muZ_mc <- muZ.mc_t0(r0, mu_mzc, h_zstar, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
#
#       muMZ_c <- muMZ.c_t0(rstar=r1, muM_zc, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
#       muZM_c <- muZM.c_t0(rstarstar=r2, muZ_mc, data_in, varnames, folds, learners_mu, bounded = FALSE, full.sample = full.sample)
#       # muMZ_c <- muZM_c
#       # partially marginalized outcome over p(M2|t,r2,c) (so dependent on M1)
#       # y_m1c <- y.zc(y_mzc, tt, r0, m_c, r2)
#       # partially marginalized outcome over p(M1|t,r1,c) (so dependent on M2)
#       # y_m2c <- y.mc(y_mzc, tt, r0, z_c, r1)
#       # fully marginalized outcome over p(M1|t,r1,c)p(M2|t,r2,c)
#       # y_c_12 <- y.c_12(y_mzc, tt, r0, z_c, r1, m_c, r2)
#
#       # fot the effect via M1 or M2 alone
#       # h_M1M2 <- p_M1 * p_M2 / p_Mjo
#       # eify <- ipw_y*h_M1M2 / mean(ipw_y*h_M1M2) * (Y - mu)
#
#       # rm1m2ipw_y <- 1*(TT==tt)*(R==r0)*(p_M1 * p_M2) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")]*p_Mjo)
#       rm1m2ipw_y <- 1*(R==r0)* ipw_y * h_zstar_mstarstar
#       eify <- rm1m2ipw_y * (Y - mu)  / mean(rm1m2ipw_y)
#
#       ipw_m1 <- 1*(R==r1) / bound(r_c[, glue("r({r1}|c)")])
#       # eifm1 <- ipw_m1 * (y_m1c - y_c_12) / mean(ipw_m1)
#       eifm1 <- ipw_m1 * (muM_zc - muMZ_c) / mean(ipw_m1)
#
#       ipw_m2 <- 1*(R==r2) / bound(r_c[, glue("r({r2}|c)")])
#       # eifm2 <- ipw_m2 * (y_m2c - y_c_12) / mean(ipw_m2)
#       eifm2 <- ipw_m2 * (muZ_mc - muZM_c) / mean(ipw_m2)
#
#       # eif with marginal mediator distributions
#       eif_mar <- eify + eifm1 + eifm2 + muZM_c # (muZM_c+muMZ_c)/2
#
#       eif[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- eif_mar
#       theta[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean(eif_mar)
#
#       rmpw[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( rm1m2ipw_y / mean(rm1m2ipw_y) * Y )
#       reg[[ glue("theta(t{tt},r{r0},r{r1},r{r2})") ]] <- mean( muMZ_c )
#
#     }
#
#     if (!is.na(rjo)) {
#       # p_Mjo_rjo <- (M==1)*(Z==1)*mz_c[, glue("mz(m1,z1|t{tt},r{rjo},c)")] +
#       #   (M==0)*(Z==1)*mz_c[, glue("mz(m0,z1|t{tt},r{rjo},c)")] +
#       #   (M==1)*(Z==0)*mz_c[, glue("mz(m1,z0|t{tt},r{rjo},c)")] +
#       #   (M==0)*(Z==0)*mz_c[, glue("mz(m0,z0|t{tt},r{rjo},c)")]
#
#       # h_zm_joint <- h.zm.joint(r0, rstar=rjo, tt, r_tc, r_zmtc)
#       h_zm_joint <- h.zm.joint_t0(r0, rstar=rjo, r_c, r_zmc)
#       mu_joint_c <- mu.joint.c_t0(r0, rstar=rjo, mu_mzc, data_in, varnames, folds, learners, bounded = FALSE, full.sample = full.sample)
#       # y_c_jo <- y.c_jo(y_mzc, tt, r0, mz_c, rjo)
#
#       # h_Mjo <- p_Mjo_rjo / p_Mjo
#       # eify_ajo <- ipw_y*h_Mjo / mean(ipw_y*h_Mjo) * (Y - mu)
#       # rmjo_ipw_y <- 1*(TT==tt)*(R==r0)*(p_Mjo_rjo) / bound(tr_c[, glue("p(t{tt},r{r0}|c)")]*p_Mjo)
#       rmjo_ipw_y <- 1*1*(R==r0) * ipw_y * h_zm_joint
#       eify_ajo <- rmjo_ipw_y * (Y - mu) / mean(rmjo_ipw_y)
#
#       # fixing R = r0 for the outcome theta(tt,r0,rjo)
#       # y_m1m2c_r0 <- y_mzc[, glue("mu(Z,M,t{tt},r{r0},c)")]
#
#       # already run mu above
#       # mu <- mu_mzc[, glue("mu(t{tt},r{r0},m,z,c)")]
#
#       ipw_m1m2 <- 1*1*(R==rjo) / bound(r_c[, glue("r({rjo}|c)")])
#
#       # eifm1m2 <- ipw_m1m2 * (y_m1m2c_r0 - y_c_jo) / mean(ipw_m1m2)
#       eifm1m2 <- ipw_m1m2 * (mu - mu_joint_c) / mean(ipw_m1m2)
#
#       # eif with the joint mediator distribution
#       eif_jo <- eify_ajo + eifm1m2 + mu_joint_c
#
#       eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]] <- eif_jo
#       theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean(eif_jo)
#
#       rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( rmjo_ipw_y / mean(rmjo_ipw_y) * (Y) )
#       reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]] <- mean( mu_joint_c )
#
#       ## when r0=rjo, ----
#       # can also direct estimating E[E(Y|t,r,c)p(c)] without mediator intervention
#       if(r0==rjo) {
#         if(Tot_dr==TRUE) {
#           y_c <- y.c_t0(data_in, varnames, Yfamily = Yfamily, folds, learners, bounded = FALSE)
#           mu_c <- y_c[, glue("mu(t{tt},r{r0},c)")]
#           # ipw_y <- 1*(TT==tt)*(R==r0) / (tr_c[, glue("p(t{tt},r{r0}|c)")])
#
#           eify_c <- ipw_y / mean(ipw_y) * (Y - mu_c) + mu_c
#
#           eif[[ glue("delta(t{tt},r{r0})") ]] <- eify_c
#           theta[[ glue("delta(t{tt},r{r0})")  ]] <- mean(eify_c)
#
#           rmpw[[ glue("delta(t{tt},r{r0})")  ]] <- mean( ipw_y / mean(ipw_y) * (Y) )
#           reg[[ glue("delta(t{tt},r{r0})")  ]] <- mean( mu_c )
#         }
#         if (Tot_dr==FALSE) {
#           eif[[ glue("delta(t{tt},r{r0})") ]] <- eif[[ glue("theta(t{tt},r{r0},rjo{rjo})") ]]
#           theta[[ glue("delta(t{tt},r{r0})")  ]] <- theta[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
#           rmpw[[ glue("delta(t{tt},r{r0})")  ]] <- rmpw[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
#           reg[[ glue("delta(t{tt},r{r0})")  ]] <- reg[[ glue("theta(t{tt},r{r0},rjo{rjo})")  ]]
#         }
#       }
#
#     }
#   }
#
#   # multiply-robust
#   eifs_effs <- effect_nott(eif)
#   thetas_effs <- sapply(eifs_effs, mean)
#   se_effs <- sapply(eifs_effs, function(s) {
#     sqrt(var(s) / nrow(data_in))
#   })
#   ci_effs <- sapply(eifs_effs, function(s) {
#     mean(s) + c(-1, 1) * qnorm(0.975) * c(sqrt(var(s) / nrow(data_in)))
#   })
#
#   # regression
#   regs_effs <- unlist(effect_nott(reg))
#   # weighting
#   rmpw_effs <- unlist(effect_nott(rmpw))
#
#   estimates <- cbind(
#     est_multi = thetas_effs,
#     se_multi = se_effs,
#     ci1_multi = ci_effs[1, ],
#     ci2_multi = ci_effs[2, ],
#     est_reg = regs_effs,
#     est_rmpw = rmpw_effs)
#
#
#   return(estimates)
# }
#
#
#
#
#
#
#
# # effect <- function(theta) {
# #
# #   theta$`TotMod`  <- with(theta, {
# #     `delta(t1,r1)`-`delta(t1,r0)`-`delta(t0,r1)`+`delta(t0,r0)`
# #   })
# #   # the (remaining) moderated treatment effect while fixing the joint mediator distribution at the reference subgroup; or, the treatment effect on the IDE of the subgroup
# #   theta$`RemainMod` <- with(theta, {
# #     `theta(t1,r1,rjo0)` - `delta(t1,r0)` -
# #       `theta(t0,r1,rjo0)` + `delta(t0,r0)`
# #   })
# #   # the (joint mediated) moderated treatment effect mediated by the difference in  the (joint distribution of) potential mediators between subgroups; or, the treatment effect on the IIE_jo of the subgroup
# #   theta$`MedMod_jo`  <- with(theta, {
# #     `delta(t1,r1)` - `theta(t1,r1,rjo0)` -
# #       `delta(t0,r1)` + `theta(t0,r1,rjo0)`
# #   })
# #   # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M1 between subgroups; or, the treatment effect on the IIE_M1 of the subgroup
# #   theta$`MedMod_M1`  <- with(theta, {
# #     `theta(t1,r1,r1,r0)` - `theta(t1,r1,r0,r0)` -
# #       `theta(t0,r1,r1,r0)` + `theta(t0,r1,r0,r0)`
# #   })
# #   # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M2 between subgroups; or, the treatment effect on the IIE_M2 of the subgroup
# #   theta$`MedMod_M2`  <- with(theta, {
# #     `theta(t1,r1,r1,r1)` - `theta(t1,r1,r1,r0)` -
# #       `theta(t0,r1,r1,r1)` + `theta(t0,r1,r1,r0)`
# #   })
# #   # the (mediated) moderated treatment effect mediated by the difference in the mutual dependence of potential mediators M1 and M2 between subgroups; or, the treatment effect on the IIE_mu of the subgroup
# #   theta$`MedMod_mu`  <- with(theta, {
# #     `delta(t1,r1)`-`theta(t1,r1,rjo0)`-`theta(t1,r1,r1,r1)`+`theta(t1,r1,r0,r0)` -
# #       (`delta(t0,r1)`-`theta(t0,r1,rjo0)`-`theta(t0,r1,r1,r1)`+`theta(t0,r1,r0,r0)`)
# #   })
# #
# #
# #   # theta$`MedMod_M1`+theta$`MedMod_M2`+theta$`MedMod_mu` - theta$`MedMod_jo`
# #
# #   theta
# # }
#
#
#
# # update.interaction <- function(valid_v, varnames) {
# #   # update the interactions
# #   # two-way
# #   valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
# #   valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]
# #   valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
# #   valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
# #   valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
# #   valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
# #   # three-way and four-way
# #   valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
# #   valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
# #   valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
# #   valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
# #   valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
# #
# #   valid_v
# # }
