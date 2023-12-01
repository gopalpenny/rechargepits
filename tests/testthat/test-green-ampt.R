

library(units)
theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- set_units(c(0, 15, 30, 60), "min")

#### Traditional Green Ampt

Fcum <- get_greenampt_flow_numerical(theta_0, theta_s, Ksat, h_b, h_0, times)
# paste(round(Fcum,4), collapse = ",")

test_that("get_greenampt_flow_numerical works", {
  expect_identical(round(Fcum,4), set_units(c(0,17.3444,24.7262,35.365), "mm"))
})

# Double check that we get the original times back with the Green-Ampt equation
times_calc <- set_units(get_greenampt_time(theta_0, theta_s, Fcum, Ksat, h_b, h_0),"min")
# paste(round(times_calc,4), collapse = ",")

test_that("get_greenampt_time works", {
  expect_identical(round(times_calc,4), set_units(c(0,15,30,60), "min"))
})

# manually written equation
# 1/set_units(Ksat, "cm/min") * (Fcum - (theta_s - theta_0) * (h_b - h_0) * log(1 + as.numeric(Fcum / (theta_s - theta_0)/(h_b-h_0))))


