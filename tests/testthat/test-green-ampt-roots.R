
library(units)
VWC_0 <- 0.2 # unitless
n <- 0.35 # unitless
# Fcum <- set_units(20, "mm") # depth
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
# times <- get_greenampt_horiz_time(VWC_0, n, Fcum, Ksat, h_b, h_0)
times <- set_units(c(0, 15, 30, 60), "min")

## get_greenampt_horiz_time
Fcum_horiz <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_horiz_time",
   VWC_0 = VWC_0, n = n, Ksat = Ksat, h_b = h_b, h_0 = h_0)
# paste(round(Fcum_horiz,4), collapse = ",")
test_that("get_greenampt_horiz_time works as roots", {
  expect_identical(round(Fcum_horiz,4), set_units(c(0,17.0094,24.0549,34.0188), "mm"))
})


## get_greenampt_time
Fcum_vert <- get_greenampt_x_roots(times = times, x_units = "mm", green_ampt_function = "get_greenampt_time",
                                   VWC_0 = VWC_0, n = n, Ksat = Ksat, h_b = h_b, h_0 = h_0)
# paste(round(Fcum_vert,4), collapse = ",")

test_that("get_greenampt_horiz_time works as roots", {
  expect_identical(round(Fcum_vert,4), set_units(c(0,17.3444,24.7262,35.365), "mm"))
})
