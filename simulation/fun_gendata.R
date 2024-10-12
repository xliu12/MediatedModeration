pr <- function(r, r_given, gen_r) {
  prob1 <- pnorm(r_given, mean = 0, sd = 1)
  r * prob1 + (1 - r) * (1 - prob1)
}
ptt <- function(tt, tt_given, gen_tt) {
  prob1 <- pnorm(tt_given, mean = 0, sd = 1)
  tt * prob1 + (1 - tt) * (1 - prob1)
}

pz <- function(tt, r, z_given, gen_z, if_binary) {
  latent <- gen_z[["z_r"]] * r + gen_z[["z_tt"]] * tt + gen_z[["z_rtt"]] * r * tt + z_given
  # prob1 <- pnorm(latent, mean = 0, sd = 1)
  # z * prob1 + (1 - z) * (1 - prob1)
  if (if_binary) {
    cond_mean <- pnorm(latent, mean = 0, sd = 1)
  }
  if (!if_binary) {
    cond_mean <- latent
  }

  cond_mean
}

pm <- function(z, tt, r, m_given, gen_m, if_binary) {
  latent <- m_given + gen_m[["m_z"]] * z + gen_m[["m_tt"]] * tt + gen_m[["m_r"]] * r +
    gen_m[["m_rtt"]] * r * tt + gen_m[["m_rz"]] * r * z + gen_m[["m_ttz"]] * tt * z +
    gen_m[["m_zrtt"]] * z * r * tt

  # prob1 <- pnorm(latent, mean = 0, sd = 1)
  # m * prob1 + (1 - m) * (1 - prob1)
  if (if_binary) {
    cond_mean <- pnorm(latent, mean = 0, sd = 1)
  }
  if (!if_binary) {
    cond_mean <- latent
  }

  cond_mean
}



pm_intz <- function(tt, r, m_given, z_given, gen_m, gen_z, if_binary) {
  if (if_binary) {
    pm_mar <- pm(z = 1, tt, r, m_given, gen_m, if_binary = if_binary) * pz(tt, r, z_given, gen_z, if_binary = if_binary) +
      pm(z = 0, tt, r, m_given, gen_m, if_binary = if_binary) * (1 - pz(tt, r, z_given, gen_z, if_binary = if_binary))
  }
  if (!if_binary) {
    condZ <- pz(tt, r, z_given, gen_z, if_binary)

    pm_mar <- pm(z = condZ, tt, r, m_given, gen_m, if_binary)
  }
  pm_mar
}


my <- function(m, z, tt, r, y_given, gen_y, if_binary) {
  latent <- y_given + gen_y[["y_m"]] * m + gen_y[["y_tt"]] * tt + gen_y[["y_r"]] * r +gen_y[["y_z"]] * z +
    gen_y[["y_rtt"]] * r * tt + gen_y[["y_rm"]]*r*m + gen_y[["y_ttm"]] * tt * m +
    gen_y[["y_rz"]] * r * z + gen_y[["y_ttz"]] *tt*z + gen_y[["y_zm"]] * z*m +
    gen_y[["y_rttm"]] * r * tt * m + gen_y[["y_rttz"]] * r * tt * z + gen_y[["y_rmz"]] * r * m * z + gen_y[["y_ttmz"]] * tt * m * z + gen_y[["y_ttrzm"]] * tt * r * m * z


  if (if_binary) {
    cond_mean <- pnorm(latent, mean = 0, sd = 1)
  }
  if (!if_binary) {
    cond_mean <- latent
  }

  cond_mean
}

# marginalized outcome ----

true.val <- function(tt, r0, r1, r2, rjo, y_given, gen_y, if_binaryZM, if_binaryY, m_given, z_given, gen_m, gen_z, linear=FALSE) {

  # binomial Mfamily ----
  if (if_binaryZM) {
    vals <- expand.grid(m1 = c(0, 1), m2 = c(0, 1))

    intm1m2 <- lapply(1:nrow(vals), function(i=1) {
      m1 <- vals$m1[i]
      m2 <- vals$m2[i]

      # glue("theta(t{tt},r{r0},r{r1},r{r2})")
      if (!is.na(r1)) {
        p_marginal <- ((m1==1) * pz(tt, r1, z_given, gen_z, if_binaryZM) + (m1==0) * (1-pz(tt, r1, z_given, gen_z, if_binaryZM))) *
          ((m2==1) * pm_intz(tt, r2, m_given, z_given, gen_m, gen_z, if_binaryZM) + (m2==0) * (1-pm_intz(tt, r2, m_given, z_given, gen_m, gen_z, if_binaryZM)))

        p_m1m2 <- p_marginal
      }
      # glue("theta(t{tt},r{r0},rjo{rjo})")
      if (!is.na(rjo)) {
        p_joint <- ((m1==1) * pz(tt, rjo, z_given, gen_z, if_binaryZM) + (m1==0) * (1-pz(tt, rjo, z_given, gen_z, if_binaryZM))) *
          ((m2==1) * pm(z=m1, tt, rjo, m_given, gen_m, if_binaryZM) + (m2==0) * (1-pm(z=m1, tt, rjo, m_given, gen_m, if_binaryZM)))

        p_m1m2 <- p_joint
      }

      p_m1m2 * my(m=m2, z=m1, tt, r0, y_given, gen_y, if_binaryY)

    } )

    int_m1m2 <- reduce(intm1m2, `+`) #rowSums(intm1m2)

    trueval <- mean(int_m1m2)
  }

  # gaussian Mfamily ----
  if (!if_binaryZM) {
    if (linear==TRUE) {
      if (!is.na(r1)) {
        Em_mar <- pm_intz(tt, r2, m_given, z_given, gen_m, gen_z, if_binaryZM)
        Ez_mar <- pz(tt, r1, z_given, gen_z, if_binaryZM)
        EYtr0r1r2 <- my(m=Em_mar, z=Ez_mar, tt, r0, y_given, gen_y, if_binaryY)
        trueval <- mean(EYtr0r1r2)

      }
      if (!is.na(rjo)) {
        Ez_mar <- pz(tt, rjo, z_given, gen_z, if_binaryZM)

        Em_z <- pm(Ez_mar, tt, rjo, m_given, gen_m, if_binaryZM)
        EYtr0rjo <- my(m=Em_z, z=Ez_mar, tt, r0, y_given, gen_y, if_binaryY)
        trueval <- mean(EYtr0rjo)

      }
    }
    if (linear==FALSE) {
      n <- length(m_given)

      if (!is.na(r1)) {
        one.draw <- function(i) {
          m_mar <- pm_intz(tt, r2, m_given, z_given, gen_m, gen_z, if_binaryZM) + rnorm(n)
          z_mar <- pz(tt, r1, z_given, gen_z, if_binaryZM) + rnorm(n)
          Ytr0r1r2 <- my(m=m_mar, z=z_mar, tt, r0, y_given, gen_y, if_binaryY)
          Ytr0r1r2
        }
        trueval <- rowMeans(do.call(cbind, lapply(1:3, one.draw))) %>% mean(.)

      }
      if (!is.na(rjo)) {
        one.draw <- function(i) {
          z_mar <- pz(tt, rjo, z_given, gen_z, if_binaryZM) + rnorm(n)
          m_z <- pm(z_mar, tt, rjo, m_given, gen_m, if_binaryZM) + rnorm(n)
          Ytr0rjo <- my(m=m_z, z=z_mar, tt, r0, y_given, gen_y, if_binaryY)
          Ytr0rjo
        }
        trueval <- rowMeans(do.call(cbind, lapply(1:3, one.draw))) %>% mean(.)
      }
    }
  }
  trueval
}

