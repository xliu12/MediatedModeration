
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



bound <- function(vals, tol = 0.01) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}
