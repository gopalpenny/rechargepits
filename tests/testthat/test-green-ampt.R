

library(units)
VWC_0 <- 0.2 # unitless
n <- 0.35 # unitless
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- set_units(c(0, 15, 30, 60), "min")

#### Traditional Green Ampt

Fcum <- get_greenampt_flow_numerical(VWC_0, n, Ksat, h_b, h_0, times)
# paste(round(Fcum,4), collapse = ",")

test_that("get_greenampt_flow_numerical works", {
  expect_identical(round(Fcum,4), set_units(c(0,17.3444,24.7262,35.365), "mm"))
})

# Double check that we get the original times back with the Green-Ampt equation
times_calc <- set_units(get_greenampt_time(VWC_0, n, Fcum, Ksat, h_b, h_0),"min")
# paste(round(times_calc,4), collapse = ",")

test_that("get_greenampt_time works", {
  expect_identical(round(times_calc,4), set_units(c(0,15,30,60), "min"))
})

# manually written equation
# 1/set_units(Ksat, "cm/min") * (Fcum - (n - VWC_0) * (h_b - h_0) * log(1 + as.numeric(Fcum / (n - VWC_0)/(h_b-h_0))))


