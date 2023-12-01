#' Get green ampt roots
#'
#' @inheritParams get_greenampt_flow_numerical
#' @param x_units Units for output variable
#' @export
#' @examples
#' library(units)
#' theta_0 <- 0.2 # unitless
#' theta_s <- 0.35 # unitless
#' Fcum <- set_units(20, "mm") # depth
#' Ksat <- set_units(0.2, "cm/h") # length / time
#' h_b <- set_units(6, "ft") # hydraulic head (length)
#' h_0 <- set_units(-10, "cm") # hydraulic head (length)
#' times <- get_greenampt_horiz_time(theta_0, theta_s, Fcum, Ksat, h_b, h_0)
#' times <- set_units(1:60, "min")
#' Fcum_horiz <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_horiz_time",
#'    theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0)
#' Fcum_vert <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_time",
#'    theta_0 = theta_0, theta_s = theta_s, Ksat = Ksat, h_b = h_b, h_0 = h_0)
get_greenampt_x_roots <- function(times, x_units, green_ampt_function_name, ...) {

  units::units_options(set_units_mode = "standard")

  for (i in 1:length(times)) {

    time <- times[i]

    green_ampt_function_i <- function(x) {
      x <- set_units(x, x_units)
      if (green_ampt_function_name=="get_greenampt_time") {
        time_root <- get_greenampt_time(Fcum = x, ...) - time
      } else if (green_ampt_function_name=="get_greenampt_horiz_time") {
        time_root <- get_greenampt_horiz_time(Fcum = x, ...) - time
      } else {
        stop("green_ampt_function_name does not match acceptable names")
      }
      return(as.numeric(time_root))
    }

    max_x <- ifelse(i == 1, 1e100, max(x * 2))
    x_i <- set_units(uniroot(green_ampt_function_i, c(0, 1e10))$root, x_units)

    if (i == 1) {
      x <- x_i
    } else {
      x <- c(x, x_i)
    }
  }
  return(x)
}
