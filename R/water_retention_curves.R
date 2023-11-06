# water_retention_curves


#' van Genuchten Water Retention Curve
#'
#' @param theta volumetric water content
#' @param theta_r residual water content
#' @param theta_s saturated water content
#' @param alpha scaling parameter (reciprocal of air entry pressure, units of 1/length)
#' @param n curve shape, dimensionless
#' @export
#' @details
#' See:
#' Pan, T., Hou, S., Liu, Y., & Tan, Q. (2019).
#' Comparison of three models fitting the soil water retention
#' curves in a degraded alpine meadow region. Scientific Reports, 9(1), 18407.
#'
#' @examples
#' # example code
#' library(units)
#' theta <- 0.24
#' theta_r <- 0.15
#' theta_s <- 0.25
#' alpha <- set_units(0.008, "1/cm") # should have units of 1/length
#' n <- 10.4
#' get_vg_hydraulic_head(theta, theta_r, theta_s, alpha, n)
get_vg_hydraulic_head <- function(theta, theta_r, theta_s, alpha, n) {
  S <- (theta - theta_r) / (theta_s - theta_r)

  m <- 1 - 1/n

  h <- -1/alpha * (S^(-1/m) - 1)^(1/n)

  return(h)
}


#' Brooks-Corey Water Retention Curve
#'
#' @inheritParams get_vg_hydraulic_head
#' @param h_d air entry pressure (units of length)
#' @param lambda curve shape, dimensionless
#' @export
#' @details
#' See:
#' Pan, T., Hou, S., Liu, Y., & Tan, Q. (2019).
#' Comparison of three models fitting the soil water retention
#' curves in a degraded alpine meadow region. Scientific Reports, 9(1), 18407.
#'
#' @examples
#' # example code
#' library(units)
#'
#' sand_rawls <- get_rawls_soil(soil_name = "Sandy loam")
#' theta <- seq(0.15, sand_rawls$porosity, length.out = 100)
#' theta_r <- sand_rawls$theta_r
#' theta_s <- sand_rawls$porosity
#' h_d <- set_units(sand_rawls$bubbling_pressure_arithmetic_cm, "cm") # should have units of length
#' lambda <- sand_rawls$pore_size_distrib_arithmetic
#' hydraulic_head <- get_bc_hydraulic_head(theta, theta_r, theta_s, h_d, lambda)
#'
#' library(tidyr)
#' rel_theta <- data.frame(rel_theta = seq(0.01, 1, length.out = 100))
#' df <- crossing(rawls_soils, rel_theta)
#' df$theta <- df$rel_theta * (df$porosity - df$theta_r) + df$theta_r
#' df$hydraulic_head <- get_bc_hydraulic_head(
#'   theta = theta,
#'   theta_r = df$theta_r,
#'   theta_s = df$porosity,
#'   h_d = set_units(df$bubbling_pressure_arithmetic_cm, "cm"),
#'   lambda = df$pore_size_distrib_arithmetic)
#' df$hydraulic_head_cm <- as.numeric(df$hydraulic_head)
#'
#' library(ggplot2)
#' ggplot(df) +
#'   geom_line(aes(theta, hydraulic_head_cm, color = texture_class, linetype = texture_class)) +
#'   ylim(c(-1000, 0)) +
#'   scale_linetype_manual(values = rep(c("solid","dashed","dotted"),4))
get_bc_hydraulic_head <- function(theta, theta_r, theta_s, h_d, lambda) {
  S <- (theta - theta_r) / (theta_s - theta_r)

  units::units_options(set_units_mode = "standard")

  h <- -h_d * S^(-1/lambda)
  h_units <- units(h_d)$numerator
  h <- replace(h, theta >= theta_s, units::set_units(0, units(h_d))) # adjust for saturated soils

  return(h)
}


