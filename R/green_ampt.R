# pit_recharge.R


#' Green-Ampt flow time
#'
#' @param VWC_0 Soil volumetric water content prior to event
#' @param n Soil porosity
#' @param Fcum cumulative infiltration in mm water equivalent
#' @param Ksat Saturated hydraulic conductivity
#' @param h_b Hydraulic head at soil surface boundary
#' @param h_0 Hydraulic head in soil prior to event
#' @export
#' @description
#' This function models the Green-Ampt equation, with a saturated
#' surface flowing vertically downward into a soil profile. The
#' assumptions are:
#'
#' 1. Homogeneous soil with initial water content (VWC_0)
#' 2. Constant pressure head (h_0) at the wetting front
#' 3. Saturated soil above above wetting front
#' 4. Continuous supply of water with constant head (h_b) at the soil surface boundary
#'
#' With these assumptions, the amount of time for a particular infiltration depth (Fcum)
#' can be calculated using Darcy's law as:
#'
#' $$t = 1/Ksat * (Fcum - (n - VWC_0)(h_b - h_0) * ln(1 + Fcum / ((n - VWC_0)(h_b - h_))))$
#'
#' @returns Returns the time at which a cumulative amount of
#' infiltration occurs.
#' @examples
#'
#' library(units)
#' VWC_0 <- 0.2 # unitless
#' n <- 0.35 # unitless
#' Fcum <- set_units(1:20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' get_greenampt_time(VWC_0, n, Fcum, Ksat, h_b, h_0)
get_greenampt_time <- function(VWC_0, n, Fcum, Ksat, h_b, h_0) {
  dVWC <- n - VWC_0
  d_h <- h_b - h_0

  depth_over_dVWC <- Fcum / dVWC / d_h
  if (length(units(depth_over_dVWC)$numerator) != 0) {
    stop(" Fcum / dVWC / d_h not yielding unitless result")
  }

  units(depth_over_dVWC) <- NULL

  t <- 1/Ksat * (Fcum - dVWC * d_h * log(1 + depth_over_dVWC))

  return(t)
}



#' Green-Ampt flow
#'
#' @inheritParams get_greenampt_time
#' @param times Times to calculate total infiltration
#' @export
#' @description
#' This function models the Green-Ampt equation, with a saturated
#' surface flowing vertically downward into a soil profile. The
#' assumptions are:
#'
#' 1. Homogeneous soil with initial water content (VWC_0)
#' 2. Constant pressure head (h_0) at the wetting front
#' 3. Saturated soil above above wetting front
#' 4. Continuous supply of water with constant head (h_b) at the soil surface boundary
#'
#' With these assumptions, the amount of time for a particular infiltration depth (Fcum)
#' can be calculated using Darcy's law as:
#'
#' $$t = 1/Ksat * (Fcum - (n - VWC_0)(h_b - h_0) * ln(1 + Fcum / ((n - VWC_0)(h_b - h_))))$$
#'
#' The cumulative infiltration is then found using a numerical root solver,
#' `uniroot`.
#'
#' @returns Returns the cumulative infiltration at times given by `times`.
#' @export
#' @examples
#'
#' library(units)
#' VWC_0 <- 0.2 # unitless
#' n <- 0.35 # unitless
#' Fcum <- set_units(1:20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- set_units(seq(5, 60, by= 5), "min")
#' fcum <- get_greenampt_flow_numerical(VWC_0, n, Ksat, h_b, h_0, times)
#' # Double check that we get the original times back with the Green-Ampt equation
#' set_units(get_greenampt_time(VWC_0, n, fcum, Ksat, h_b, h_0),"min")
get_greenampt_flow_numerical <- function(VWC_0, n, Ksat, h_b, h_0, times) {
  Fcum_vert <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_time",
                                     VWC_0 = VWC_0, n = n, Ksat = Ksat, h_b = h_b, h_0 = h_0)

  return(Fcum_vert)
}



#' Green-Ampt horizontal flow
#'
#'
#' @inheritParams get_greenampt_flow_numerical
#' @export
#' @description
#' This function models saturated horizontal flow from a saturated
#' surface into a soil profile. The assumptions are the same as
#' in the Green-Ampt equation, except that there is no effect of
#' gravity because all flow is assumed to be horizontal. The equation is:
#'
#' $$F^2$$
#' $$F^2$$
#' $$F^2$$
#' @returns Returns the cumulative infiltration at the specified times
#' @examples
#'
#' library(units)
#' VWC_0 <- 0.2 # unitless
#' n <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- set_units(seq(5, 60, by= 5), "min")
#' Fcum <- get_greenampt_horiz_flow(VWC_0, n, Ksat, h_b, h_0, times)
get_greenampt_horiz_flow <- function(VWC_0, n, Ksat, h_b, h_0, times) {

  dVWC <- n - VWC_0

  units::units_options(set_units_mode = "standard")

  h_b <- set_units(h_b, units(Ksat)$numerator)
  h_0 <- set_units(h_0, units(Ksat)$numerator)

  Fcum <- sqrt(2 * dVWC * Ksat * (h_b - h_0) * times)

  return(Fcum)
}


#' Green-Ampt horizontal flow time
#'
#'
#' @inheritParams get_greenampt_time
#' @export
#' @description
#' This function models saturated horizontal flow from a saturated
#' surface into a soil profile. The assumptions are the same as
#' in the Green-Ampt equation, except that there is no effect of
#' gravity because all flow is assumed to be horizontal. The equation is:
#'
#' $$F^2$$
#' $$F^2$$
#' $$F^2$$
#' @returns Returns the time at which a cumulative amount of
#' infiltration occurs.
#' @examples
#'
#' library(units)
#' VWC_0 <- 0.2 # unitless
#' n <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_horiz_time(VWC_0, n, Fcum, Ksat, h_b, h_0)
get_greenampt_horiz_time <- function(VWC_0, n, Fcum, Ksat, h_b, h_0) {
  dVWC = n - VWC_0

  t = 1/2 * Fcum^2 / (dVWC * Ksat * (h_b - h_0))

  return(t)
}



#' Horizontal Green-Ampt flow integrated
#'
#'
#' @inheritParams get_greenampt_time
#' @param thickness Depth over which to be integrated. If NULL, set to h_b
#' @export
#' @description
#' This function models saturated horizontal flow from a saturated
#' surface into a soil profile. The assumptions are the same as
#' in the Green-Ampt equation, except that there is no effect of
#' gravity because all flow is assumed to be horizontal. The function
#' further integrates over the vertical depth of a pit such that
#' the hydraulic head at the boundary goes from (h_b - thickness) to h_b.
#' The equation is:
#'
#' $$F^2$$
#' $$F^2$$
#' $$F^2$$
#' @returns Returns the time at which a cumulative amount of
#' infiltration occurs.
#' @examples
#'
#' library(units)
#' VWC_0 <- 0.2 # unitless
#' n <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' Fv <- get_greenampt_horiz_flow_integrated(VWC_0, n, Ksat, h_b, h_0, t = set_units(1,"hr"))
#' Fv <- get_greenampt_horiz_flow_integrated(VWC_0, n, Ksat, h_b, h_0, t = set_units(1,"hr"), thickness = set_units(1, "ft"))
get_greenampt_horiz_flow_integrated <- function(VWC_0, n, Ksat, h_b, h_0, t, thickness = NULL) {

  if (is.null(thickness)) {
    thickness <- h_b
  }

  dVWC = n - VWC_0
  # h_diff <- h_b - h_0

  Ksat2 <- Ksat * set_units(1, units(Ksat)$numerator)
  h_b_bottom2 <- h_b * set_units(1, units(h_b)$numerator)
  h_b_top2 <- h_b * set_units(1, units(h_b)$numerator) - thickness * set_units(1, units(thickness)$numerator)
  h_02 <- h_0 * set_units(1, units(h_0)$numerator)

  Fv2 <- sqrt(8/9 * dVWC * Ksat2 * t) * ((h_b_bottom2 - h_02)^(3/2) - (h_b_top2-h_02)^(3/2))
  Fv <- Fv2 / set_units(1, units(h_b)$numerator) / set_units(1, units(Ksat)$numerator)

  return(Fv)
}



