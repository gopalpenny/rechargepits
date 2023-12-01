# pit_recharge.R


#' Green-Ampt flow time
#'
#' @param theta_0 Soil volumetric water content prior to event
#' @param theta_s Soil porosity
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
#' 1. Homogeneous soil with initial water content (theta_0)
#' 2. Constant pressure head (h_0) at the wetting front
#' 3. Saturated soil above above wetting front
#' 4. Continuous supply of water with constant head (h_b) at the soil surface boundary
#'
#' With these assumptions, the amount of time for a particular infiltration depth (Fcum)
#' can be calculated using Darcy's law as:
#'
#' $$t = 1/Ksat * (Fcum - (theta_s - theta_0)(h_b - h_0) * ln(1 + Fcum / ((theta_s - theta_0)(h_b - h_))))$
#'
#' @returns Returns the time at which a cumulative amount of
#' infiltration occurs.
#' @examples
#'
#' library(units)
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Fcum <- set_units(1:20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' get_greenampt_time(theta_0, theta_s, Fcum, Ksat, h_b, h_0)
get_greenampt_time <- function(theta_0, theta_s, Fcum, Ksat, h_b, h_0) {
  dVWC <- theta_s - theta_0
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
#' 1. Homogeneous soil with initial water content (theta_0)
#' 2. Constant pressure head (h_0) at the wetting front
#' 3. Saturated soil above above wetting front
#' 4. Continuous supply of water with constant head (h_b) at the soil surface boundary
#'
#' With these assumptions, the amount of time for a particular infiltration depth (Fcum)
#' can be calculated using Darcy's law as:
#'
#' $$t = 1/Ksat * (Fcum - (theta_s - theta_0)(h_b - h_0) * ln(1 + Fcum / ((theta_s - theta_0)(h_b - h_))))$$
#'
#' The cumulative infiltration is then found using a numerical root solver,
#' `uniroot`.
#'
#' @returns Returns the cumulative infiltration at times given by `times`.
#' @export
#' @examples
#'
#' library(units)
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Fcum <- set_units(1:20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- set_units(seq(5, 60, by= 5), "min")
#' fcum <- get_greenampt_flow_numerical(theta_0, theta_s, Ksat, h_b, h_0, times)
#' # Double check that we get the original times back with the Green-Ampt equation
#' set_units(get_greenampt_time(theta_0, theta_s, fcum, Ksat, h_b, h_0),"min")
get_greenampt_flow_numerical <- function(theta_0, theta_s, Ksat, h_b, h_0, times) {
  Fcum_vert <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_time",
                                     theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0)

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
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- set_units(seq(5, 60, by= 5), "min")
#' Fcum <- get_greenampt_horiz_flow(theta_0, theta_s, Ksat, h_b, h_0, times)
get_greenampt_horiz_flow <- function(theta_0, theta_s, Ksat, h_b, h_0, times) {

  dVWC <- theta_s - theta_0

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
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_horiz_time(theta_0, theta_s, Fcum, Ksat, h_b, h_0)
get_greenampt_horiz_time <- function(theta_0, theta_s, Fcum, Ksat, h_b, h_0) {
  dVWC = theta_s - theta_0

  t = 1/2 * Fcum^2 / (dVWC * Ksat * (h_b - h_0))

  return(t)
}



#' Horizontal Green-Ampt flow integrated
#'
#' @inheritParams get_greenampt_flow_numerical
#' @param d Depth over which to be integrated. If NULL, set to h_b
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
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#'
#' Fv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr"), d = NULL)
#'
#' # Get the infiltration over a 1 mm differential depth to compare with the point infiltration
#' Fv_dv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr"), d = set_units(1, "mm"))
#' Fv_point <- Fv_dv/ set_units(1, "mm")
#' Fv_point
#' # Get the point infiltration
#' F_point <- get_greenampt_horiz_flow(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr")) %>% set_units("mm")
#' F_point
#' # percent error:
#' (Fv_point - F_point) / F_point * 100
get_greenampt_horiz_flow_integrated <- function(theta_0, theta_s, Ksat, h_b, h_0, times, d = NULL) {

  units::units_options(set_units_mode = "standard")

  if (is.null(d)) {
    d <- h_b
  }

  units_obj <- d^2

  dVWC <- theta_s - theta_0
  # h_diff <- h_b - h_0

  t_num <- as.numeric(set_units(times, "h"))

  Ksat_num <- as.numeric(set_units(Ksat, "cm/h"))
  h_b_bottom_num <- as.numeric(set_units(h_b, "cm"))
  h_b_top_num <- as.numeric(set_units(h_b, "cm") - set_units(d, "cm"))
  h_0_num <- as.numeric(set_units(h_0, "cm"))

  Fv_num <- sqrt(8/9 * dVWC * Ksat_num * t_num) * ((h_b_bottom_num - h_0_num)^(3/2) - (h_b_top_num-h_0_num)^(3/2))
  Fv_cm2 <- set_units(Fv_num, "cm^2")
  Fv <- set_units(Fv_cm2, units(units_obj))

  return(Fv)
}



#' Green-Ampt cylindrical horizontal flow time
#'
#'
#' @inheritParams get_greenampt_time
#' @param r_b radius from the centroid to the free water--soil boundary
#' @param F_r cumulative radial infiltration through the cylinder (units of L^2)
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
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' F_r <- set_units(10, "ft^2") # units of length^2
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_cyl_horiz_time(theta_0, theta_s, F_r, Ksat, h_b, h_0, r_b)
get_greenampt_cyl_horiz_time <- function(theta_0, theta_s, F_r, Ksat, h_b, h_0, r_b) {
  dVWC <- theta_s - theta_0

  r_f <- sqrt(F_r * dVWC / pi + r_b^2)

  t <- dVWC / (4 * Ksat * (h_b - h_0)) * (r_b^2 + 2 * r_f^2 * as.numeric(log(r_f/r_b)) - r_f^2)
  # t <- dVWC / (4 * Ksat * (h_b - h_0)) * (r_b^2 + r_f^2 * (2 * as.numeric(log(r_f/r_b)) - 1))

  return(t)
}



