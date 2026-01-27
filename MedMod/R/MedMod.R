
#' Mediated moderation analysis with the debiased machine learning method
#' @param data a \code{data.frame} containing variables involved in the analysis.
#' @param outcome a character string of the name of the column in "data" that correspond to the outcome variable (Y).
#' @param mediators a character vector of the names of the columns in "data" that correspond to the mediator M1. The definition for the mediated moderation via mediator M1 (i.e., MedMod_M1) considers shifting the levels of M1, while fixing the levels of mediator M2 at its levels among the reference subgroup.
#' @param treatment a character string of the name of the column in "data" that correspond to the treatment variable (T). The treatment variable should be dummy coded (e.g., 1 for the intervention condition, 0 for the control condition).
#' @param subgroup a character string of the name of the column in "data" that correspond to the moderator subgroup status (R). The subgroup status should be dummy coded, with 1 for the focal subgroup and 0 for the reference subgroup.
#' @param covariates a character vector of the names of the columns in "data" that correspond to pretreatment (i.e., baseline) covariates (C).
#' @param learners a character vector specifying the methods for estimating the nuisance models. Currently supported options include the list of methods included in the "SuperLearner" package (<https://cran.r-project.org/package=SuperLearner>), which can be listed via function listWrappers() of the "SuperLearner" package.
#' @param num_folds the number of folds used for the cross-fitting procedure.
#'
#'
#'
#' @return A data frame containing the results for the mediated moderation estimands and the expected potential outcomes in the estimands. The results include the estimates, standard error estimates, and 0.95 confidence intervals.
#'
#' @examples
#'

#' library(tidyverse)
#' library(glue)
#' library(origami)
#' library(mvtnorm)
#' library(SuperLearner)
#' library(ranger)
#' library(nnet)
#'
#' # Run ----------
#' data <- MedMod::data
#' head(data)
#'
#' set.seed(12345)
#' out <- MedMod::MedMod(
#'   data = data, # data containing all variables
#'   outcome = "Outcome", # name of outcome
#'   mediators = c("M1", "M2"), # name of mediators. When more than two mediators are included, the first mediator is designated as M1, and all remaining mediators are collectively treated as M2.
#'   treatment = "Intervention", # name of treatment variable
#'   subgroup = "Male", # name of moderator subgroup variable
#'   covariates = c("C.1", "C.2", "C.3"), # names of baseline covariates
#'   learners = c("SL.mean", "SL.glm", "SL.ranger", "SL.nnet"),
#'   # methods from the SuperLearner package to include in the super learner ensemble (intercept-only model, generalized linear model, random forest, neural network).
#'   # To see all available methods, run: SuperLearner::listWrappers()
#'   num_folds = 4, # number of folds for cross-fitting
#'   ci.level = 0.95 # default: 95% confidence interval
#' )
#' # Extract results for mediated moderation
#' out %>%
#'   filter(Estimand %in% c("TotMod", "MedMod", "RemainMod", # Total moderation, Mediated moderation, Remaining  moderation
#'                          "MedMod_M1", "MedMod_M2", "MedMod_mu")) # When more than one mediators, the output also include: Mediated moderation via M1, Mediated moderation via M2, Mediated moderation due to mediators' mutual dependence
#'
#'
#' @export
#'


#'
#'
MedMod <- function(
    data,
    outcome,
    mediators,
    treatment,
    subgroup,
    covariates,
    learners = c("SL.glm"),
    learners_h = NULL,
    learners_mu = NULL,
    num_folds = 5,
    ci.level = 0.95
) {

  # inputs -----------
  M1names <- mediators[1]

  if (length(mediators) > 1) {
    M2names <- mediators[-1]
  } else {
    M2names <- NULL
  }


  Yname <- outcome
  ttname <- treatment
  Rname <- subgroup
  Cnames <- covariates

  # Mediators <- list(M1names, M2names)
  Yfamily <- ifelse(length(unique(data[[Yname]])) > 2, "gaussian", "binomial")

  Mfamily <- "h" # one or more continuous mediators
  if((nrow(unique(data[, M1names, drop=FALSE])) == 2) & (nrow(unique(data[, M2names, drop=FALSE])) <= 2)) {
    Mfamily <- "b"
  }

  data_in <- data
  Znames <- M1names #Mediators[[1]]
  Mnames <- M2names #Mediators[[2]]
  data_in <- data_in %>% mutate(
    Rtt = data_in[[Rname]] * data_in[[ttname]]
  )
  if (!is.null(Znames)) {
    data_in <- data_in %>% mutate(
      RZ = data_in[[Rname]] * data_in[, Znames],
      ttZ = data_in[[ttname]] * data_in[, Znames],
      ttRZ = data_in[[ttname]] * data_in[[Rname]] * data_in[, Znames]
    )
  }
  if (!is.null(Mnames)) {
    data_in <- data_in %>% mutate(
      RM = data_in[[Rname]] * data_in[, Mnames],
      ttM = data_in[[ttname]] * data_in[, Mnames],
      ttRM = data_in[[Rname]] * data_in[[ttname]] * data_in[, Mnames]
    )
  }

  if (!is.null(Znames) & !is.null(Mnames)) {
    data_in <- data_in %>% mutate(
      MZ = data_in[, Znames] * data_in[, Mnames],
      ttMZ = data_in[[ttname]] * data_in[, Znames] * data_in[, Mnames],
      RMZ = data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames],
      ttRMZ = data_in[[ttname]] * data_in[[Rname]] * data_in[, Znames] * data_in[, Mnames]
    )
  }

  data_in <- do.call(data.frame, data_in)

  varnames <- list("R" = Rname, "tt" = ttname, "M" = Mnames, "Y" = Yname,
                   "Z" = Znames, "C" = Cnames,
                   "Rtt" = paste0("Rtt"), "RM" = paste0("RM"), "ttM" = paste0("ttM"), "ttRM" = paste0("ttRM"),
                   "RZ" = paste0("RZ"), "ttZ" = paste0("ttZ"), "MZ" = paste0("MZ"),
                   "ttRZ" = "ttRZ", "ttMZ" = "ttMZ", "RMZ" = "RMZ", "ttRMZ" = "ttRMZ")
  if (is.null(Mnames)) {
    varnames[grep("M", names(varnames), value = TRUE)] <- NULL
  }
  data_in$ID <- 1:nrow(data_in)

  folds <- origami::make_folds(n = nrow(data_in), V = num_folds)
  if (num_folds==1) {
    folds[[1]]$training_set <- folds[[1]]$validation_set
  }
  if (num_folds > 1) {
    folds <- vector("list", length = num_folds)

    tr_vals <- expand.grid(tt=c(0,1), r=c(0,1))
    j <- 1
    for(j in 1:nrow(tr_vals)) {
      id_tr <- (data_in[[varnames$tt]] == tr_vals$tt[j]) &
        (data_in[[varnames$R]] == tr_vals$r[j])

      z <- origami::make_folds(data_in[id_tr, ], V = num_folds)

      k <- 1
      for (k in 1:num_folds) {
        folds[[k]]$v <- k
        folds[[k]]$training_set <- c(folds[[k]]$training_set, data_in[id_tr, "ID"][z[[k]]$training_set])
        folds[[k]]$validation_set <- c(folds[[k]]$validation_set, data_in[id_tr, "ID"][z[[k]]$validation_set])
        attr(folds[[k]], "class")<-"fold"
      }
    }


  }

  if (is.null(learners_h)) {
    learners_h <- learners
  }
  if (is.null(learners_mu)) {
    learners_mu <- learners
  }
  medmod_data <- mget(ls(), envir = environment())

  if (Mfamily=="b") {
    out <- b.MedMod(
      medmod_data = medmod_data,
      fity.interact = TRUE,
      fitm.interact = TRUE,
      TotMod_dr = FALSE
    )
  }

  if (Mfamily=="h") {

    out <- h.MedMod(
      medmod_data = medmod_data,
      fity.interact = TRUE,
      fitm.interact = TRUE,
      TotMod_dr = FALSE,
      full.sample = TRUE
    )
  }

  list2env(out, envir = environment())

  # estimates and inferences --------------------
  eifs_effs <- effect(eif)
  thetas_effs <- sapply(eifs_effs, mean)
  se_effs <- sapply(eifs_effs, function(s) {
    sqrt(var(s) / length(s))
  })
  ci_effs <- sapply(eifs_effs, function(s) {
    mean(s) + c(-1, 1) * qnorm(1-(1-ci.level)/2) * sqrt(var(s) / length(s))
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
  out <- list(estimates = estimates, eifs = eifs_effs)

  out$estimates <- data.frame(out$estimates) %>%
    rownames_to_column(var = "Estimand") %>%
    mutate(Estimand = ifelse(Estimand == "MedMod_jo", "MedMod", Estimand)) %>%
    rename(Estimate=est_multi, std.error=se_multi,
           `CI.lower`=ci1_multi, `CI.upper`=ci2_multi) %>%
    select(!c(est_rmpw,est_reg))

  # out
  out$estimates
}
