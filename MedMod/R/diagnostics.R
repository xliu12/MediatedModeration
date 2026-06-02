
# Resolve the fit-inputs bundle for a positivity diagnostic, accepting either the object
# returned by MedMod() (which carries the "fit_inputs" attribute) or that fit-inputs list
# passed directly. Shared by MedMod_overlap() and MedMod_mediator_positivity() so the
# dual-input handling stays identical between them. Internal; not exported.
.medmod_resolve_fit_inputs <- function(object) {
  fi <- MedMod_fit_inputs(object)
  if (is.null(fi)) {
    if (is.list(object) && !is.null(object$data_in) && !is.null(object$Mfamily)) {
      fi <- object
    } else {
      stop("`object` must be the result of MedMod() (carrying a \"fit_inputs\" attribute) ",
           "or the fit-inputs list returned by MedMod_fit_inputs().")
    }
  }
  fi
}


#' Positivity / common-support overlap diagnostic for a MedMod fit
#'
#' Diagnoses the \emph{first} component of the positivity assumption underlying
#' \code{\link{MedMod}}: \eqn{p(T = t, R = r \mid C = c) > 0} for all values. The joint
#' probability \eqn{p(T, R \mid C)} is the "propensity score" whose inverse forms the
#' inverse-probability weights. Practically, the subgroups being compared must have sufficient
#' overlap in the covariates \eqn{C}: if an estimated probability concentrates near 0, that
#' (treatment, subgroup) cell has a practical positivity violation and the corresponding
#' estimates can be unstable. This is a controlled descriptive comparison across subgroups, so
#' the diagnostic contrasts the subgroups within each treatment condition rather than a single
#' treated-vs-control propensity. (The \emph{second} positivity component, the mediator
#' condition \eqn{p(M = m \mid T, R, C) > 0}, is diagnosed by
#' \code{\link{MedMod_mediator_positivity}}.)
#'
#' The joint propensity \eqn{p(T, R \mid C)} and subgroup probability \eqn{p(R \mid C)} are
#' \strong{re-fit with} \code{bounded = FALSE} (on the original cross-fitting folds carried by
#' the fit) so that the probability truncation used during estimation -- \code{bound()} with
#' tolerance \code{bound_tol} -- cannot hide violations behind the truncation. \code{bound_tol}
#' is then used purely as the flag threshold and as the dashed reference lines, not as a clamp.
#'
#' The figure has three panels sharing a common subgroup legend: one panel per treatment
#' condition showing the density of the joint propensity \eqn{p(T = t, R = r \mid C = c)} for
#' each subgroup, plus a panel showing the subgroup-membership probability \eqn{p(R \mid C)}.
#' Dashed vertical lines mark \code{bound_tol} (and \code{1 - bound_tol}); density mass beyond
#' these lines flags practical positivity concerns.
#'
#' @param object either the object returned by \code{\link{MedMod}} (a data frame carrying the
#'   \code{"fit_inputs"} attribute) or the fit-inputs list itself (as returned by
#'   \code{\link{MedMod_fit_inputs}}).
#' @param subgroup_labels length-2 character vector labelling the reference (R = 0) and focal
#'   (R = 1) subgroups, in that order.
#' @param treatment_labels length-2 character vector labelling the control (T = 0) and
#'   treatment (T = 1) conditions, in that order.
#' @param bound_tol the probability bounding tolerance used by the estimator (default 0.01);
#'   drawn as dashed reference lines and used to compute the near-0 / near-1 shares.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{plot}}{a \code{ggplot} object: overlap densities of the estimated
#'       probabilities, faceted by treatment condition (plus the subgroup-propensity panel) and
#'       coloured by subgroup.}
#'     \item{\code{summary}}{a data frame with one row per estimated-probability column
#'       (\code{quantity}, \code{treatment}, \code{subgroup}, \code{n}, \code{min},
#'       \code{median}, \code{max}, \code{pct_near0}, \code{pct_near1}). The common-support
#'       range of the joint propensity across the four (T, R) cells is attached as
#'       \code{attr(summary, "common_support")}.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(12345)
#' out <- MedMod::MedMod(
#'   data = MedMod::data, outcome = "Outcome",
#'   mediators = c("M1", "M2"), treatment = "Intervention",
#'   subgroup = "Male", covariates = c("C.1", "C.2", "C.3"),
#'   learners = c("SL.mean", "SL.glm"), num_folds = 4
#' )
#' res <- MedMod_overlap(out)
#' res$plot
#' res$summary
#' attr(res$summary, "common_support")
#' }
#'
#' @export
MedMod_overlap <- function(object,
                           subgroup_labels  = c("Reference (R=0)", "Focal (R=1)"),
                           treatment_labels = c("Control (T=0)", "Intervention (T=1)"),
                           bound_tol = 0.01) {

  # Resolve the fit-inputs bundle, whether given the MedMod() data frame or the list directly.
  fi <- .medmod_resolve_fit_inputs(object)

  # Re-fit the joint propensity p(T, R | C) and subgroup probability p(R | C) with
  # bounded = FALSE on the original cross-fitting folds. The estimator runs these with
  # bounded = TRUE (truncating to [bound_tol, 1 - bound_tol]); using those truncated values
  # here would clamp positivity violations out of view, which is exactly what this diagnostic
  # is meant to surface. bound_tol is therefore only a flag threshold / reference line below.
  r_c          <- r.c (fi$data_in, fi$varnames, fi$folds, fi$learners, bounded = FALSE)
  t_rc         <- t.rc(fi$data_in, fi$varnames, fi$folds, fi$learners, bounded = FALSE)
  joint_tr_c   <- tr.c(r_c, t_rc)      # p(T, R | C): cols p(t{tt},r{r}|c)
  subgroup_r_c <- r_c                  # p(R | C):    cols r({r}|c)

  sub_lab <- function(r) factor(subgroup_labels[r + 1], levels = subgroup_labels)

  # --- Long data frame for the joint-propensity panels (one per treatment condition) ---
  jt <- colnames(joint_tr_c)
  jt_tt <- as.integer(sub("p\\(t([01]),r[01]\\|c\\)", "\\1", jt))
  jt_r  <- as.integer(sub("p\\(t[01],r([01])\\|c\\)", "\\1", jt))
  joint_long <- do.call(rbind, lapply(seq_along(jt), function(j) {
    data.frame(
      panel    = treatment_labels[jt_tt[j] + 1],
      subgroup = sub_lab(jt_r[j]),
      prob     = joint_tr_c[, j],
      stringsAsFactors = FALSE
    )
  }))

  # --- Long data frame for the subgroup-propensity panel p(R | C) ---
  sc <- colnames(subgroup_r_c)
  sc_r <- as.integer(sub("r\\(([01])\\|c\\)", "\\1", sc))
  sub_long <- do.call(rbind, lapply(seq_along(sc), function(j) {
    data.frame(
      panel    = "Subgroup p(R|C)",
      subgroup = sub_lab(sc_r[j]),
      prob     = subgroup_r_c[, j],
      stringsAsFactors = FALSE
    )
  }))

  plot_df <- rbind(joint_long, sub_long)
  plot_df$panel <- factor(plot_df$panel,
                          levels = c(treatment_labels, "Subgroup p(R|C)"))

  cols <- scales::hue_pal()(2)
  names(cols) <- subgroup_labels

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = prob,
                                    fill = subgroup,
                                    colour = subgroup)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = c(bound_tol, 1 - bound_tol),
                        linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual("Subgroup", values = cols) +
    ggplot2::scale_colour_manual("Subgroup", values = cols) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(
      x = "Estimated probability",
      y = "Density",
      title = "Positivity / overlap diagnostic",
      subtitle = paste0("Joint propensity p(T,R|C) by treatment condition, and subgroup ",
                        "propensity p(R|C)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12),
      legend.position = "bottom"
    )

  # --- Numeric common-support summary ---
  summ_one <- function(quantity, treatment, subgroup, x) {
    data.frame(
      quantity  = quantity,
      treatment = treatment,
      subgroup  = subgroup,
      n         = length(x),
      min       = min(x),
      median    = stats::median(x),
      max       = max(x),
      pct_near0 = mean(x < bound_tol),
      pct_near1 = mean(x > 1 - bound_tol),
      stringsAsFactors = FALSE
    )
  }
  summ <- rbind(
    do.call(rbind, lapply(seq_along(jt), function(j) {
      summ_one("joint p(T,R|C)", treatment_labels[jt_tt[j] + 1],
               subgroup_labels[jt_r[j] + 1], joint_tr_c[, j])
    })),
    do.call(rbind, lapply(seq_along(sc), function(j) {
      summ_one("subgroup p(R|C)", NA_character_,
               subgroup_labels[sc_r[j] + 1], subgroup_r_c[, j])
    }))
  )
  rownames(summ) <- NULL

  # Common support across the four joint (T,R) cells: [max of mins, min of maxes].
  joint_mins <- apply(joint_tr_c, 2, min)
  joint_maxs <- apply(joint_tr_c, 2, max)
  attr(summ, "common_support") <- c(lower = max(joint_mins), upper = min(joint_maxs))

  list(plot = p, summary = summ)
}


#' Mediator positivity diagnostic for a MedMod fit
#'
#' Diagnoses the second component of the positivity assumption underlying
#' \code{\link{MedMod}}. The assumption has two parts -- \eqn{p(T = t, R = r \mid C = c) > 0}
#' (the joint propensity, checked by \code{\link{MedMod_overlap}}) and
#' \eqn{p(M = m \mid T = t, R = r, C = c) > 0} (the mediator condition, checked here). The
#' quantity examined depends on the mediator family of the fit, because the estimator handles
#' the two cases differently:
#'
#' \describe{
#'   \item{Binary mediators (\code{Mfamily = "b"})}{The mediator condition involves a genuine
#'     joint probability \eqn{p(M \mid T, R, C)} (built by \code{mz.c()} and consumed
#'     directly in \code{b.MedMod()}). The bound is literally \eqn{[0, 1]}; the substantive
#'     positivity violation is a mediator-cell probability \emph{near 0}.}
#'   \item{Continuous mediators (\code{Mfamily = "h"})}{Here \eqn{p(M = m \mid \cdot)} is a
#'     density, which the estimator never touches: the mediator density ratio is reparameterized
#'     by Bayes (see \code{h.zm.joint()}) so the mediator-side
#'   positivity depends on \eqn{p(R \mid M, T, C)} (built by \code{r.zmtc()}). Because
#'     this term enters the inverse-odds-ratio weight as a denominator under both subgroup
#'     orderings, values \emph{near 0 or near 1} are where the weights destabilize.}
#' }
#'
#' Both quantities are probabilities in \eqn{[0, 1]} that can be plotted directly; for the
#' continuous case the plotted probability is the reparameterized proxy the estimator actually
#' uses in place of the uncheckable density. The relevant nuisances are \strong{re-fit with}
#' \code{bounded = FALSE} (on the original cross-fitting folds carried by the fit) so that the
#' probability truncation used during estimation -- \code{bound()} with tolerance
#' \code{bound_tol} -- cannot hide violations behind the truncation. \code{bound_tol} is then
#' used purely as the flag threshold and as the dashed reference lines, not as a clamp.
#'
#' @param object either the object returned by \code{\link{MedMod}} (a data frame carrying the
#'   \code{"fit_inputs"} attribute) or the fit-inputs list itself (as returned by
#'   \code{\link{MedMod_fit_inputs}}).
#' @param subgroup_labels length-2 character vector labelling the reference (R = 0) and focal
#'   (R = 1) subgroups, in that order.
#' @param treatment_labels length-2 character vector labelling the control (T = 0) and
#'   treatment (T = 1) conditions, in that order.
#' @param bound_tol the probability bounding tolerance used by the estimator (default 0.01);
#'   drawn as dashed reference lines and used to compute the near-0 / near-1 shares.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{plot}}{a \code{ggplot} object: densities of the 
#'       mediator-side probabilities. For continuous mediators, faceted by treatment condition
#'       and coloured by the subgroup \eqn{p(R \mid M, T, C)}; for binary mediators,
#'       faceted by treatment-by-subgroup and coloured by mediator cell \eqn{p(M \mid T, R, C)}.}
#'     \item{\code{summary}}{a data frame with one row per estimated-probability column
#'       (\code{quantity}, \code{treatment}, \code{subgroup}, \code{mediator_cell}, \code{n},
#'       \code{min}, \code{median}, \code{max}, \code{pct_near0}, \code{pct_near1}). The
#'       common-support range across the plotted columns is attached as
#'       \code{attr(summary, "common_support")}.}
#'   }
#'
#' @seealso \code{\link{MedMod_overlap}} for the joint-propensity \eqn{p(T, R \mid C)} component.
#'
#' @examples
#' \dontrun{
#' set.seed(12345)
#' out <- MedMod::MedMod(
#'   data = MedMod::data, outcome = "Outcome",
#'   mediators = c("M1", "M2"), treatment = "Intervention",
#'   subgroup = "Male", covariates = c("C.1", "C.2", "C.3"),
#'   learners = c("SL.mean", "SL.glm"), num_folds = 4
#' )
#' res <- MedMod_mediator_positivity(out)
#' res$plot
#' res$summary
#' attr(res$summary, "common_support")
#' }
#'
#' @export
MedMod_mediator_positivity <- function(object,
                                       subgroup_labels  = c("Reference (R=0)", "Focal (R=1)"),
                                       treatment_labels = c("Control (T=0)", "Intervention (T=1)"),
                                       bound_tol = 0.01) {

  # Resolve the fit-inputs bundle, whether given the MedMod() data frame or the list directly.
  fi <- .medmod_resolve_fit_inputs(object)
  data_in       <- fi$data_in
  varnames      <- fi$varnames
  folds         <- fi$folds
  learners      <- fi$learners
  fitm.interact <- fi$fitm.interact
  Mfamily       <- fi$Mfamily

  has_M <- !is.null(varnames$M)

  # --- Re-fit the relevant mediator-side probabilities with bounded = FALSE -------------------
  # Truncation (bounded = TRUE) during estimation would clamp violations into [bound_tol,
  # 1 - bound_tol] and hide them; we deliberately want the un-truncated estimates here.
  if (Mfamily == "b") {
    # Binary mediators: joint mediator probability p(M, Z | T, R, C).
    z_c <- z.c(data_in, varnames, interact = fitm.interact,
               folds = folds, learners = learners, bounded = FALSE)
    if (has_M) {
      m_zc <- m.zc(data_in, varnames, interact = fitm.interact,
                   folds = folds, learners = learners, bounded = FALSE)
      mat <- mz.c(z_c, m_zc)
    } else {
      mat <- z_c  # only a single binary mediator (Z); p(Z | T, R, C)
    }
  } else {
    # Continuous mediators: reparameterized density ratios in terms of p(R | M, T, C).
    if (has_M) {
      mat <- r.zmtc(data_in, varnames, interact = fitm.interact,
                    folds = folds, learners = learners, bounded = FALSE)
    } else {
      mat <- r.ztc(data_in, varnames, interact = fitm.interact,
                   folds = folds, learners = learners, bounded = FALSE)
    }
  }

  sub_lab <- function(r) factor(subgroup_labels[r + 1], levels = subgroup_labels)

  # --- One-row-per-column summary (shared across families) ------------------------------
  summ_one <- function(quantity, treatment, subgroup, mediator_cell, x) {
    data.frame(
      quantity      = quantity,
      treatment     = treatment,
      subgroup      = subgroup,
      mediator_cell = mediator_cell,
      n             = length(x),
      min           = min(x),
      median        = stats::median(x),
      max           = max(x),
      pct_near0     = mean(x < bound_tol),
      pct_near1     = mean(x > 1 - bound_tol),
      stringsAsFactors = FALSE
    )
  }

  cols <- colnames(mat)

  if (Mfamily == "h") {
    # Columns: r({r}|z,m,t{tt},c) or r({r}|z,t{tt},c). Parse the predicted-R digit and T digit.
    col_r  <- as.integer(sub("^r\\(([01])\\|.*$", "\\1", cols))
    col_tt <- as.integer(sub("^.*t([01]),c\\)$",  "\\1", cols))
    # quantity_lab <- if (has_M) "P(R|Z,M,T,C)" else "P(R|Z,T,C)"
    # in the manuscript, M = (M1, M2) if two mediators, and there is no Z
    quantity_lab <- "P(R|M,T,C)"

    plot_df <- do.call(rbind, lapply(seq_along(cols), function(j) {
      data.frame(
        panel    = treatment_labels[col_tt[j] + 1],
        subgroup = sub_lab(col_r[j]),
        prob     = mat[, j],
        stringsAsFactors = FALSE
      )
    }))
    plot_df$panel <- factor(plot_df$panel, levels = treatment_labels)

    pal <- scales::hue_pal()(2)
    names(pal) <- subgroup_labels

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = prob, fill = subgroup, colour = subgroup)) +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::geom_vline(xintercept = c(bound_tol, 1 - bound_tol),
                          linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free_y") +
      ggplot2::scale_fill_manual("Subgroup", values = pal) +
      ggplot2::scale_colour_manual("Subgroup", values = pal) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(
        x = "Estimated probability",
        y = "Density",
        title = "Mediator positivity diagnostic",
        subtitle = paste0(quantity_lab, " in the reparameterized mediator density ratio by treatment condition")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 12),
        legend.position = "bottom"
      )

    summ <- do.call(rbind, lapply(seq_along(cols), function(j) {
      summ_one(quantity_lab, treatment_labels[col_tt[j] + 1],
               subgroup_labels[col_r[j] + 1], NA_character_, mat[, j])
    }))

  } else {
    # Binary. Columns: mz(m{m},z{z}|t{tt},r{r},c) (with M) or z({z}|t{tt},r{r},c) (no M).
    col_tt <- as.integer(sub("^.*\\|t([01]),r[01],c\\)$", "\\1", cols))
    col_r  <- as.integer(sub("^.*,r([01]),c\\)$",         "\\1", cols))
    if (has_M) {
      col_m <- as.integer(sub("^mz\\(m([01]),.*$",       "\\1", cols))
      col_z <- as.integer(sub("^mz\\(m[01],z([01])\\|.*$", "\\1", cols))
      cell_lab     <- paste0("M2=", col_m, ", M1=", col_z)
      quantity_lab <- "P(M2,M1|T,R,C)"
    } else {
      col_m <- rep(NA_integer_, length(cols))
      col_z <- as.integer(sub("^z\\(([01])\\|.*$", "\\1", cols))
      cell_lab     <- paste0("M=", col_z)
      quantity_lab <- "P(M|T,R,C)"
    }
    cell_levels <- unique(cell_lab)

    plot_df <- do.call(rbind, lapply(seq_along(cols), function(j) {
      data.frame(
        treatment = factor(treatment_labels[col_tt[j] + 1], levels = treatment_labels),
        subgroup  = sub_lab(col_r[j]),
        cell      = factor(cell_lab[j], levels = cell_levels),
        prob      = mat[, j],
        stringsAsFactors = FALSE
      )
    }))

    pal <- scales::hue_pal()(length(cell_levels))
    names(pal) <- cell_levels

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = prob, fill = cell, colour = cell)) +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::geom_vline(xintercept = c(bound_tol, 1 - bound_tol),
                          linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      ggplot2::facet_grid(subgroup ~ treatment, scales = "free_y") +
      ggplot2::scale_fill_manual("Mediator cell", values = pal) +
      ggplot2::scale_colour_manual("Mediator cell", values = pal) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(
        x = "Estimated probability",
        y = "Density",
        title = "Mediator positivity diagnostic",
        subtitle = paste0("Joint mediator probability ", quantity_lab,
                          " by treatment x subgroup")
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 12),
        legend.position = "bottom"
      )

    summ <- do.call(rbind, lapply(seq_along(cols), function(j) {
      summ_one(quantity_lab, treatment_labels[col_tt[j] + 1],
               subgroup_labels[col_r[j] + 1], cell_lab[j], mat[, j])
    }))
  }

  rownames(summ) <- NULL

  # Common support across the plotted columns: [max of mins, min of maxes].
  col_mins <- apply(mat, 2, min)
  col_maxs <- apply(mat, 2, max)
  attr(summ, "common_support") <- c(lower = max(col_mins), upper = min(col_maxs))

  list(plot = p, summary = summ)
}
