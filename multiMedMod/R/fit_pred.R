
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


glmfitting <- function(train, valid.list, yname, xnames, varnames,
                     type, bounded = FALSE) {
  family <- ifelse(type == "binomial", binomial(), gaussian())
  df_lm <- data.frame(Y = train[[yname]],
                      train[, c(xnames), drop = FALSE] )

  fit <- glm(Y ~ ., data = df_lm,
    family = family[[1]]
  )

  preds <- sapply(valid.list, function(validX) {
    newX <- data.frame(validX[, c(xnames), drop = FALSE])
    preds <- predict(fit, newdata = newX, type = "response")
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

#' @export
SL.lightgbm <- function(Y, X, newX, family, obsWeights, id, nrounds = 1000, verbose = -1,
                        learning_rate = 0.1, min_data_in_leaf = 10, max_depth = -1, ...) {
  if (!requireNamespace("lightgbm", quietly = FALSE)) {
    stop("loading required package (lightgbm) failed", call. = FALSE)
  }

  if (family$family == "gaussian") {
    objective <- "regression"
    evalu <- ""
  }

  if (family$family == "binomial") {
    objective <- "binary"
    evalu <- "binary_logloss"
  }

  if (!is.matrix(X)) {
    X <- model.matrix(~. - 1, X)
  }

  lgb_data <- try(
    lightgbm::lgb.Dataset(
      data = X,
      label = as.numeric(Y)
    ), silent = TRUE
  )

  try(lightgbm::set_field(lgb_data, "weight", as.numeric(obsWeights)), silent = TRUE)

  params <- list(
    min_data_in_leaf = min_data_in_leaf,
    learning_rate = learning_rate,
    max_depth = max_depth
  )

  model <- lightgbm::lgb.train(params, data = lgb_data, obj = objective, eval = evalu,
                               nrounds = nrounds, verbose = verbose)

  if (!is.matrix(newX)) {
    newX <- model.matrix(~. - 1, newX)
  }

  pred <- predict(model, newX)
  fit <- list(object = model)
  class(fit) <- c("SL.lightgbm")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
predict.SL.lightgbm <- function(object, newdata, family, ...) {
  if (!requireNamespace("lightgbm", quietly = FALSE)) {
    stop("loading required package (lightgbm) failed", call. = FALSE)
  }

  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~. - 1, newdata)
  }
  pred <- predict(object$object, newdata)
  return(pred)
}


SL.hal.modified <- function(...) {
  SL.hal9001(...,
             # tuning
             max_degree = 1, smoothness_orders = 1,
             num_knots = c(5), #num_knots =c(25, 10, 5),
             # reduce_basis = NULL,
             X_unpenalized = NULL
  )
}
