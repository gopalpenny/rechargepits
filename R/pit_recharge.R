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


#' Horizontal Green-Ampt flow time
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



