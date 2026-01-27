
crossfit <- function(train, valid.list, yname, xnames, varnames,
                     type, learners, IDname = NULL, bounded = FALSE) {
  family <- ifelse(type == "binomial", binomial(), gaussian())
  df_lm <- data.frame(Y = train[[yname]],
                      train[, c(xnames), drop = FALSE] )

  ID <- train[[IDname]]
  set.seed(12345)
  fit <- SuperLearner::SuperLearner(
    df_lm$Y,
    df_lm[, c(xnames), drop = FALSE],
    family = family[[1]],
    id = ID,
    SL.library = learners
  )

  preds <- sapply(valid.list, function(validX) {
    newX <- data.frame(validX[, c(xnames), drop = FALSE])
    preds <- predict(fit, newX[, fit$varNames])$pred
    if (!bounded) {
      return(preds)
    }
    bound(preds)
  }, simplify = TRUE)

  out <- list(fit = fit, preds = preds)
  return(out)
}


update.interactions <- function(valid_v, varnames) {
  # update the interactions
  valid_v[, c(varnames$Rtt)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]
  valid_v[, c(varnames$RZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]
  valid_v[, c(varnames$ttZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$Z)]
  valid_v[, c(varnames$ttRZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$Z)]

  # if two mediators
  if (!is.null(varnames$M)) {
    valid_v[, c(varnames$RM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)] # if no M, then this is a data frame with 0 columns and nrow(valid_v) rows
    valid_v[, c(varnames$ttM)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
    valid_v[, c(varnames$ttRM)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]
    valid_v[, c(varnames$MZ)] <- valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
    valid_v[, c(varnames$ttMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
    valid_v[, c(varnames$RMZ)] <- valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
    valid_v[, c(varnames$ttRMZ)] <- valid_v[, c(varnames$tt)]*valid_v[, c(varnames$R)]*valid_v[, c(varnames$M)]*valid_v[, c(varnames$Z)]
  }

  valid_v
}


bound <- function(vals, tol = 0.01) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}
