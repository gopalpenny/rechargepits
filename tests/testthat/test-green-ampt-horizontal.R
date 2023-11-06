
#### Horizontal green ampt

library(units)
VWC_0 <- 0.2 # unitless
n <- 0.35 # unitless
Fcum <- set_units(20, "mm") # depth
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- set_units(c(0, 15, 30, 60), "min")
Fcum <- get_greenampt_horiz_flow(VWC_0, n, Ksat, h_b, h_0, times)
paste(round(Fcum,4), collapse = ",")

test_that("get_greenampt_horiz_flow works", {
  expect_equal(round(Fcum, 4), set_units(c(0,1.7009,2.4055,3.4019), "cm"))
})


# Double check that we get the original times back with the Green-Ampt equation
times_calc <- set_units(get_greenampt_horiz_time(VWC_0, n, Fcum, Ksat, h_b, h_0),"min")
# paste(round(times_calc,4), collapse = ",")

test_that("get_greenampt_horiz_time works", {
  expect_identical(round(times_calc,4), set_units(c(0,15,30,60), "min"))
})

# manually written equation
# set_units(Fcum^2 / 2 / (n - VWC_0) / Ksat / (h_b - h_0), "min")
