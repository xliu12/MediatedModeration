
#' Mediated moderation for multiple mediators
#' @param data_in a \code{data.frame} containing variables involved in the analysis.
#' @param Yname a character string of the name of the column in "data_in" that correspond to the outcome variable (Y).
#' @param M1names a character vector of the names of the columns in "data_in" that correspond to the mediator M1. The definition for the mediated moderation via mediator M1 (i.e., MedMod_M1) considers shifting the levels of M1, while fixing the levels of mediator M2 at its levels among the reference subgroup.
#' @param M2names a character vector of the names of the columns in "data_in" that correspond to the mediator M2. The definition for the mediated moderation via mediator M2 (i.e., MedMod_M2) considers shifting the levels of M2, while fixing the levels of mediator M1 at its levels among the focal subgroup.
#' @param ttname a character string of the name of the column in "data_in" that correspond to the treatment variable (T). The treatment variable should be dummy coded (e.g., 1 for the intervention condition, 0 for the control condition).
#' @param Rname a character string of the name of the column in "data_in" that correspond to the moderator subgroup status (R). The subgroup status should be dummy coded, with 1 for the focal subgroup and 0 for the reference subgroup.
#' @param Cnames a character vector of the names of the columns in "data_in" that correspond to pretreatment (i.e., baseline) covariates (C).
#' @param learners a character vector specifying the methods for estimating the nuisance models. Currently supported options include the list of methods included in the "SuperLearner" package (<https://cran.r-project.org/package=SuperLearner>), which can be listed via function listWrappers() of the "SuperLearner" package.
#' @param num_folds the number of folds used for the cross-fitting procedure.
#' @param Yfamily a character specifying the link function for a generalized linear model for the outcome. Currently supported options include "gaussian" and "binomial".
#'
#' @return A data frame containing the results for the mediated moderation estimands and the expected potential outcomes in the estimands. The results include the estimates, standard error estimates, and 0.95 confidence intervals.
#'
#'  The column "Estimand" in the returned data frame contains the names of the estimands. Specifically:
#'  "MedMod_M1": The mediated moderation via mediator M1.
#'  "MedMod_M2": The mediated moderation via mediator M2.
#'  "MedMod_mu": The mediated moderation due to mutual dependence of the mediators M1 and M2.
#'  "MedMod_jo": The mediated moderation via both mediators M1 and M2 jointly. "MedMod_jo" is the sum of "MedMod_M1", "MedMod_M2", and "MedMod_mu".
#'  "RemainMod": The remaining moderation (i.e., not via (M1, M2)).
#'  "TotMod": The total moderation in the treatment effect on outcome between subgroups, controlling for the pretreatment covariates. "TotMod" is the sum of "RemainMod" and "MedMod_jo".


#'
#'
#' @export
#'


#'
#'
MedMod <- function(
    data_in,
    Yname,
    M1names,
    M2names,
    ttname,
    Rname,
    Cnames,
    learners = c("SL.glm"),
    num_folds = 1,
    Yfamily = NULL
) {
  Mediators <- list(M1names, M2names)
  if (is.null(Yfamily)) {
    Yfamily <- ifelse(length(unique(data_in[[Yname]])) > 2, "gaussian", "binomial")
  }

  Mfamily <- "h"
  if ( length(M1names)+length(M2names) == 2 )  {
    if (max( map_dbl(Mediators, ~length(unique(data_in[[.]]))) ) == 2) {
      Mfamily <- "b"
    }
  }

  out_eif <- FALSE
  if (Mfamily=="b") {
    out <- b.MedMod(
      data_in = data_in,
      Yname = Yname,
      M1names = M1names,
      M2names = M2names,
      ttname = ttname,
      Rname = Rname,
      Cnames = Cnames,
      learners = learners,
      num_folds = num_folds,
      Yfamily = Yfamily,
      fity.interact = TRUE,
      fitm.interact = TRUE,
      TotMod_dr = FALSE,
      out_eif = out_eif
    )
  }

  if (Mfamily=="h") {
    out <- h.MedMod(
      data_in = data_in,
      Yname = Yname,
      M1names = M1names,
      M2names = M2names,
      ttname = ttname,
      Rname = Rname,
      Cnames = Cnames,
      learners = learners,
      learners_h = learners,
      learners_mu = learners,
      num_folds = num_folds,
      Yfamily = Yfamily,
      fity.interact = FALSE,
      fitm.interact = FALSE,
      TotMod_dr = FALSE,
      full.sample = TRUE,
      out_eif = out_eif
    )
  }

  if (out_eif == TRUE) {
    out$estimates <- data.frame(out$estimates) %>%
      rownames_to_column(var = "Estimand") %>%
      rename(Estimate=est_multi, SE=se_multi,
             `95% CI.lower`=ci1_multi, `95% CI.upper`=ci2_multi) %>%
      select(!c(est_rmpw,est_reg))
  }
  if (out_eif == FALSE) {
    out <- data.frame(out) %>%
      rownames_to_column(var = "Estimand") %>%
      rename(Estimate=est_multi, SE=se_multi,
             `95% CI.lower`=ci1_multi, `95% CI.upper`=ci2_multi) %>%
      select(!c(est_rmpw,est_reg))
  }
  out
}
