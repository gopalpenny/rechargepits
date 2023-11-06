# data

#' Rawls soil data
#'
#' @details
#' Soil properties from Rawls et al (1982)
#'
#' Rawls, W. J., Brakensiek, D. L., & Saxtonn, K. E. (1982).
#' Estimation of Soil Water Properties. Transactions of the
#' ASAE. https://doi.org/10.13031/2013.33720
"rawls_soils"


#' Get rawls soil params
#'
#' @param soil_name A soil name matching rawls_soils$texture_class
#' @export
#' @details
#' Returns soil properties from Rawls et al (1982)
#'
#' @examples
#'
#' print(rawls_soils$texture_class)
#' sand_rawls <- get_rawls_soil(soil_name = "Sandy loam")
#' get_rawls_soil(soil_name = c("Sand","Sandy loam"))
get_rawls_soil <- function(soil_name) {
  all_soils <- rawls_soils
  soil_rawls <- all_soils[all_soils$texture_class %in% soil_name,]
  if (nrow(soil_rawls) <= 1) {
    stop("Invalid soil name: ", soil_name)
  }
  return(soil_rawls)
}
