# pit_recharge.R

#' Horizontal Green-Ampt flow time
#'
#'
#' @param VWC_0 Soil volumetric water content prior to event
#' @param n Soil porosity
#' @param Fcum cumulative infiltration in mm water equivalent
#' @param Ksat Saturated hydraulic conductivity
#' @param h_b Hydraulic head at soil surface boundary
#' @param h_0 Hydraulic head in soil prior to event
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
#' VWC_0 <- 0.2
#' n <- 0.35
#' Fcum <- 20 # mm
#' Ksat <- 5 # UNITS??
#' h_b <- 6
#' h_0 <- -20
#' get_greenampt_time(VWC_0, n, Fcum, Ksat, h_b, h_0)
get_greenampt_time <- function(VWC_0, n, Fcum, Ksat, h_b, h_0) {
  dVWC = n - VWC_0

  t = 1/2 * Fcum^2 / (dVWC * Ksat * (h_b - h_0))

  return(t)
}

