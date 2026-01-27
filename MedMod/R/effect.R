

effect <- function(theta) {

  theta$`TotMod`  <- with(theta, {
    `delta(t1,r1)`-`delta(t1,r0)`-`delta(t0,r1)`+`delta(t0,r0)`
  })
  # the (remaining) moderated treatment effect while fixing the joint mediator distribution at the reference subgroup; or, the treatment effect on the IDE of the subgroup
  theta$`RemainMod` <- with(theta, {
    `theta(t1,r1,rjo0)` - `delta(t1,r0)` -
      `theta(t0,r1,rjo0)` + `delta(t0,r0)`
  })
  # the (joint mediated) moderated treatment effect mediated by the difference in  the (joint distribution of) potential mediators between subgroups; or, the treatment effect on the IIE_jo of the subgroup
  theta$`MedMod_jo`  <- with(theta, {
    `delta(t1,r1)` - `theta(t1,r1,rjo0)` -
      `delta(t0,r1)` + `theta(t0,r1,rjo0)`
  })
  if (any(grepl("t1,r1,r1,r0", names(theta)))) {
    # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M1 between subgroups; or, the treatment effect on the IIE_M1 of the subgroup
    theta$`MedMod_M1`  <- with(theta, {
      `theta(t1,r1,r1,r0)` - `theta(t1,r1,r0,r0)` -
        `theta(t0,r1,r1,r0)` + `theta(t0,r1,r0,r0)`
    })
    # the (mediated) moderated treatment effect mediated by the difference in the (marginal distribution of) potential mediator M2 between subgroups; or, the treatment effect on the IIE_M2 of the subgroup
    theta$`MedMod_M2`  <- with(theta, {
      `theta(t1,r1,r1,r1)` - `theta(t1,r1,r1,r0)` -
        `theta(t0,r1,r1,r1)` + `theta(t0,r1,r1,r0)`
    })
    # the (mediated) moderated treatment effect mediated by the difference in the mutual dependence of potential mediators M1 and M2 between subgroups; or, the treatment effect on the IIE_mu of the subgroup
    theta$`MedMod_mu`  <- with(theta, {
      `delta(t1,r1)`-`theta(t1,r1,r1,r1)`-`theta(t1,r1,rjo0)`+`theta(t1,r1,r0,r0)` -
        (`delta(t0,r1)`-`theta(t0,r1,r1,r1)`-`theta(t0,r1,rjo0)`+`theta(t0,r1,r0,r0)`)
    })
  }

  theta
}
