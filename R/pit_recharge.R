# pit_recharge.R


#' Rectangular pit flow
#'
#' Horizontal and vertical parallel flow from rectangular pit
#' @inheritParams get_greenampt_flow_numerical
#' @param height vertical height (depth) of the pit
#' @param length horizontal length of the pit
#' @param width horizontal width of the pit
#' @export
#' @details
#' This function calculates infiltration from a rectangular pit,
#' using the assumptions in the Green Ampt equation. The infiltration
#' is the sum of the following
#'
#' horizontal infiltration * (width + length) * 2
#' vertical infiltration * width * length
#'
#' The function assumed the pit is full for the duration of
#' the infiltration experiment, meaning the surface has a
#' pressure head of zero and the bottom has a pressure head
#' of length.
#' @examples
#'
#' library(units)
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- set_units(1,"hr")
#' length <- set_units(4, "ft")
#' width <- set_units(4, "ft")
#' height <- set_units(6, "ft")
get_pit_recharge_rectangular <- function(height, length, width, theta_s, Ksat, theta_0, h_0, times) {
  h_b <- height

  # get horizontal flow:
  Fv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, times = times)
  flow_horizontal <- Fv * (width * 2 + length * 2)

  fcum_vert <- get_greenampt_flow_numerical(theta_0, theta_s, Ksat, h_b, h_0, times)
  flow_vertical <- fcum_vert * width * length
  flow_vertical <- units::set_units(flow_vertical, units(flow_horizontal))

  flow <- flow_horizontal + flow_vertical


  return(flow)
}
