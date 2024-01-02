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
#' \deqn{t = \frac{1}{K_{sat}}\Big[F - \Delta \theta (h_b - h_0) \ln (1 + \frac{F}{\Delta \theta (h_b - h_0)}) \Big]}
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
#' \deqn{t = \frac{1}{K_{sat}}\Big[F - \Delta \theta (h_b - h_0) \ln (1 + \frac{F}{\Delta \theta (h_b - h_0)}) \Big]}
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
  Fcum_vert <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function_name = "get_greenampt_time",
                                     theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0)

  return(Fcum_vert)
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
#' \deqn{t = \frac{F^2}{2 \Delta \theta K_{sat} (h_b - h_0)}}
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
#' \deqn{F = \sqrt{2 \Delta \theta K_{sat} (h_b - h_0) t}}
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

  h_b <- units::set_units(h_b, units(Ksat)$numerator)
  h_0 <- units::set_units(h_0, units(Ksat)$numerator)

  Fcum <- sqrt(2 * dVWC * Ksat * (h_b - h_0) * times)

  return(Fcum)
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
#' \deqn{F_v = \sqrt{\frac{8}{9} \Delta \theta K_{sat} t} \Big[(h_b - h_0)^{3/2} - (h_b - d - h_0)^{3/2} \Big]}
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
#' Fv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s,
#'         Ksat, h_b, h_0, times = set_units(1,"hr"), d = NULL)
#'
#' # Get infiltration over a 1 mm differential depth to compare with the point infiltration
#' Fv_dv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s,
#'         Ksat, h_b, h_0, times = set_units(1,"hr"), d = set_units(1, "mm"))
#' Fv_point <- Fv_dv/ set_units(1, "mm")
#' Fv_point
#' # Get the point infiltration
#' F_point <- get_greenampt_horiz_flow(theta_0, theta_s, Ksat, h_b, h_0,
#'         times = set_units(1,"hr")) %>% set_units("mm")
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

  t_num <- as.numeric(units::set_units(times, "h"))

  Ksat_num <- as.numeric(units::set_units(Ksat, "cm/h"))
  h_b_bottom_num <- as.numeric(units::set_units(h_b, "cm"))
  h_b_top_num <- as.numeric(units::set_units(h_b, "cm") - units::set_units(d, "cm"))
  h_0_num <- as.numeric(units::set_units(h_0, "cm"))

  Fv_num <- sqrt(8/9 * dVWC * Ksat_num * t_num) * ((h_b_bottom_num - h_0_num)^(3/2) - (h_b_top_num-h_0_num)^(3/2))
  Fv_cm2 <- units::set_units(Fv_num, "cm^2")
  Fv <- units::set_units(Fv_cm2, units(units_obj))

  return(Fv)
}



#' Green-Ampt cylindrical horizontal flow time
#'
#'
#' @inheritParams get_greenampt_time
#' @param r_b radius from the centroid to the free water--soil boundary
#' @param F_c cumulative radial infiltration through the cylinder (units of area or L^2)
#' @export
#' @description
#' This function models saturated outward radial flow from a circular saturated
#' surface horizontally into a soil profile. The assumptions are the same as
#' in the Green-Ampt equation, except that there is no effect of
#' gravity because all flow is assumed to be horizontal. The equation includes
#' \eqn{r_f}, which is the radius to the edge of the front from the center of the
#' cylinder:
#'
#' \deqn{t = \frac{\Delta \theta}{4 (h_b - h_0)} \Big[r_b^2 + r_f^2 (2  \ln \frac{r_f}{r_b} - 1) \Big]}
#'
#' It is straightforward to convert between cumulative infiltration \eqn{F_c}
#' and \eqn{r_f}:
#'
#' \deqn{F_c = \pi (r_f^2 - r_b^2) \Delta \theta}
#'
#' and the reverse:
#'
#' \deqn{r_f = \sqrt{ \frac{F_c}{\pi \Delta \theta} + r_b^2}}
#'
#' @returns Returns the time at which a cumulative amount of
#' infiltration occurs.
#' @examples
#'
#' library(units)
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' F_c <- set_units(10, "ft^2") # units of length^2
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_cyl_horiz_time(theta_0, theta_s, F_c, Ksat, h_b, h_0, r_b)
get_greenampt_cyl_horiz_time <- function(theta_0, theta_s, F_c, Ksat, h_b, h_0, r_b) {
  dVWC <- theta_s - theta_0

  r_f <- sqrt(F_c * dVWC / pi + r_b^2)

  t <- dVWC / (4 * Ksat * (h_b - h_0)) * (r_b^2 + 2 * r_f^2 * as.numeric(log(r_f/r_b)) - r_f^2)
  # t <- dVWC / (4 * Ksat * (h_b - h_0)) * (r_b^2 + r_f^2 * (2 * as.numeric(log(r_f/r_b)) - 1))

  return(t)
}

#' Green-Ampt cylindrical horizontal flow
#'
#'
#' @inheritParams get_greenampt_cyl_horiz_time
#' @param times times at which to calculate cumulative infiltration
#' @param F_units character indicating the areal (L^2) `units` for the output
#' @export
#' @description
#' This function models saturated outward radial flow from a circular saturated
#' surface horizontally into a soil profile. The assumptions are the same as
#' in the Green-Ampt equation, except that there is no effect of
#' gravity because all flow is assumed to be horizontal. The equation includes
#' \eqn{r_f}, which is the radius to the edge of the front from the center of the
#' cylinder:
#'
#' \deqn{t = \frac{\Delta \theta}{4 (h_b - h_0)} \Big[r_b^2 + r_f^2 (2  \ln \frac{r_f}{r_b} - 1) \Big]}
#'
#' It is straightforward to convert between cumulative infiltratioj \eqn{F_c}
#' and \eqn{r_f}:
#'
#' \deqn{F_c = \pi (r_f^2 - r_b^2) \Delta \theta}
#'
#' and the reverse:
#'
#' \deqn{r_f = \sqrt{ \frac{F_c}{\pi \Delta \theta} + r_b^2}}
#'
#' @returns Returns the cumulative amount of
#' infiltration after time `times`.
#' @examples
#'
#' library(units)
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' F_c <- set_units(c(1, 5, 10, 20), "ft^2") # units of length^2
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_cyl_horiz_time(theta_0, theta_s, F_c, Ksat, h_b, h_0, r_b)
#' F_c_calculated <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b, h_0,
#'                                                     r_b, times, F_units = "ft^2")
get_greenampt_cyl_horiz_numerical <- function(theta_0, theta_s, Ksat, h_b, h_0, r_b, times, F_units = "ft^2") {
  units::units_options(set_units_mode = "standard")

  units_check <- units::set_units(1, F_units)
  units_check <- units::set_units(units_check, "m^2") # check to ensure units are L^2

  F_c_cum <- get_greenampt_x_roots(times = times, x_units = F_units, green_ampt_function_name = "get_greenampt_cyl_horiz_time",
                                   theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0, r_b = r_b)

  return(F_c_cum)
}

#' Get greenampt cylindrical flow integrated
#'
#' @inheritParams get_greenampt_cyl_horiz_numerical
#' @param d Depth over which to be integrated. If NULL, set to h_b
#' @param num_sections Number of discrete bins over which to calculate recharge
#' @export
#' @description
#' This function calculates green-ampt cylindrical flow (using
#' `get_greenampt_cyl_horiz_numerical`) over a number of discrete vertical sections. This
#' allows calculation of volumetric flow moving horizontally outward from a vertical,
#' cylindrical pit filled with water. For instance, consider an idealized example of
#' a cylindrical recharge well with a radius of 2 ft (`r_b`), a depth of 4 ft (`h_b`), screened over
#' the bottom 3 ft (`d`). In this case, we might specify `num_sections = 3`, in which case
#' this function would estimate planar radial flow (\eqn{F_r}) using `get_greenampt_cyl_horiz_numerical` with
#' \eqn{h_b \in \{1.5, 2.5, 3.5 ft\}}. Then the volumetric flow would be estimated as:
#'
#' \deqn{F_{c,total} = 1 ft \ times (F_r |_{h_b = 1.5} + F_r |_{h_b = 2.5} + F_r |_{h_b = 3.5})}
#'
#' The \eqn{1 ft} at the front represents the width of each section.
#'
#' @return
#' This function returns cumulative volumetric infiltration.
#' @examples
#' library(units)
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' times <- set_units(c(0, 0.5, 1, 5), "hr")
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(4, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' d <- set_units(3, "ft")
#' num_sections <- 3
#' F_v_cum <- get_greenampt_cyl_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, r_b, times,
#'                                              F_units = "ft^3", num_sections = 3, d = d)
#'
#' # This is equivalent to a discretized version of `get_greenampt_cyl_horiz_numerical`:
#' Fc_1_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat,
#'                 h_b = set_units(1.5,"ft"), h_0, r_b, times, F_units = "ft^2")
#' Fc_2_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat,
#'                 h_b = set_units(2.5,"ft"), h_0, r_b, times, F_units = "ft^2")
#' Fc_3_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat,
#'                 h_b = set_units(3.5,"ft"), h_0, r_b, times, F_units = "ft^2")
#' Fc_sum <- d / num_sections * (Fc_1_5 + Fc_2_5 + Fc_3_5)
#' F_v_cum
#' Fc_sum
get_greenampt_cyl_flow_integrated <- function(theta_0, theta_s, Ksat, h_b, h_0, r_b, times, F_units = "ft^3", num_sections = 5, d = NULL) {
  if (is.null(d)) {
    d <- h_b
  }
  num_iter <- length(times) * num_sections
  if (num_iter > 20) {
    pbtrue <- TRUE
    pb <- utils::txtProgressBar(max = num_sections)
  } else{
    pbtrue <- FALSE
  }
  section_width <- d / num_sections
  section_centers <- (h_b - d - section_width/2) + (section_width) * 1:num_sections
  for (i in 1:length(section_centers)) {
    # i <- 1
    h_b_i <- section_centers[i]
    F_v_i <- section_width * get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b_i, h_0, r_b, times, F_units = "ft^2")

    if (i == 1) {
      F_v_cum <- F_v_i
    } else {
      F_v_cum <- F_v_cum + F_v_i
    }
    if (pbtrue) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  F_v_cum <- units::set_units(F_v_cum, F_units)
  return(F_v_cum)
}


#' Green-Ampt half-spherical pressure flow time
#'
#'
#' @inheritParams get_greenampt_cyl_horiz_time
#' @param F_s cumulative radial infiltration through the cylinder (units of volume or L^3)
#' @export
#' @description
#'
#' This function models saturated flow from a wetted half-spherical
#' surface into a soil profile. Flow is driven by pressure only and gravity
#' is ignored. This might mimic a situation of subsurface recharge where the
#' pressure head gradient is significantly larger than the elevation head
#' gradient. Otherwise, the assumptions are the same as in the Green-Ampt
#' equation. The equation includes
#' \eqn{r_f}, which is the radius to the edge of the front from the center of the
#' sphere. The equation to calculate the time required to achieve
#' some quantum of recharge is:
#'
#' \deqn{t = \frac{\Delta \theta}{K_{sat} (h_b - h_0)}  \Bigg[ \frac{r_f^3}{3 r_b}  - \frac{r_f^2}{2} + \frac{r_b^2}{6} \Bigg]}
#'
#' It is straightforward to convert between cumulative infiltration \eqn{F_c}
#' and \eqn{r_f}:
#'
#' \deqn{F_s = \Delta \theta \pi \frac{2}{3}\Big( r_f^3 - r_b^3 \Big)}
#'
#' and
#'
#' \deqn{r_f = \Big( \frac{3 F_s}{2 \Delta \theta \pi} + r_b^3 \Big)^{1/3}}
#'
#' @returns Returns the time at which a cumulative amount of volumetric
#' infiltration occurs through the half sphere.
#' @examples
#'
#' library(units)
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' F_s <- set_units(10, "ft^3") # units of length^2
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_hsphere_time(theta_0, theta_s, F_s, Ksat, h_b, h_0, r_b)
get_greenampt_hsphere_time <- function(theta_0, theta_s, F_s, Ksat, h_b, h_0, r_b) {
  dVWC <- theta_s - theta_0

  r_f <- (3 * F_s / (2 * dVWC * pi) + r_b^3)^(1/3)

  t <- dVWC / (Ksat * (h_b - h_0)) * (r_f^3/(3*r_b) - r_f^2/2 + r_b^2/6)
  # t <- dVWC / (4 * Ksat * (h_b - h_0)) * (r_b^2 + r_f^2 * (2 * as.numeric(log(r_f/r_b)) - 1))

  return(t)
}



#' Green-Ampt half-spherical pressure flow
#'
#'
#' @inheritParams get_greenampt_hsphere_time
#' @param times times at which to calculate cumulative infiltration
#' @param F_units character indicating the volumetric (L^3) `units` for the output
#' @export
#' @description
#' This function models saturated flow from a wetted half-spherical
#' surface into a soil profile. Flow is driven by pressure only and gravity
#' is ignored. This might mimic a situation of subsurface recharge where the
#' pressure head gradient is significantly larger than the elevation head
#' gradient. Otherwise, the assumptions are the same as in the Green-Ampt
#' equation. The equation includes
#' \eqn{r_f}, which is the radius to the edge of the front from the center of the
#' sphere. The equation to calculate the time required to achieve
#' some quantum of recharge is:
#'
#' \deqn{t = \frac{\Delta \theta}{K_{sat} (h_b - h_0)}  \Bigg[ \frac{r_f^3}{3 r_b}  - \frac{r_f^2}{2} + \frac{r_b^2}{6} \Bigg]}
#'
#' It is straightforward to convert between cumulative infiltration \eqn{F_c}
#' and \eqn{r_f}:
#'
#' \deqn{F_s = \Delta \theta \pi \frac{2}{3}\Big( r_f^3 - r_b^3 \Big)}
#'
#' and
#'
#' \deqn{r_f = \Big( \frac{3 F_s}{2 \Delta \theta \pi} + r_b^3 \Big)^{1/3}}
#'
#' @returns Returns the cumulative amount of volumetric
#' infiltration occurs through the half sphere.
#' @examples
#'
#' library(units)
#' r_b <- set_units(2, "ft") # length
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' F_s <- set_units(c(1, 5, 10, 20), "ft^3") # units of length^2
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_hsphere_time(theta_0, theta_s, F_s, Ksat, h_b, h_0, r_b)
#' F_s_calc <- get_greenampt_hsphere_numerical(theta_0, theta_s, Ksat, h_b, h_0, r_b, times)
#' F_s_calc
get_greenampt_hsphere_numerical <- function(theta_0, theta_s, Ksat, h_b, h_0, r_b, times, F_units = "ft^3") {
  units::units_options(set_units_mode = "standard")
  units_check <- units::set_units(1, F_units)
  units_check <- units::set_units(units_check, "m^3") # check to ensure units are L^3
  F_s_cum <- get_greenampt_x_roots(times = times, x_units = F_units, green_ampt_function_name = "get_greenampt_hsphere_time",
                                     theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0, r_b = r_b)

  return(F_s_cum)
}

