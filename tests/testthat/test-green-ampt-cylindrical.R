

library(units)
r_b <- set_units(2, "ft") # length
theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
F_c <- set_units(c(0, 1, 5, 10, 20), "ft^2") # units of length^2
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- get_greenampt_cyl_horiz_time(theta_0, theta_s, F_c, Ksat, h_b, h_0, r_b)
F_c_calculated <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b, h_0, r_b, times)


paste(round(times, 4), collapse = ", ")
paste(round(F_c_calculated, 5), collapse = ", ")


test_that("get_greenampt_cyl_horiz_time works", {
  expect_equal(round(times, 4), set_units(c(0, 3e-04, 0.0063, 0.0248, 0.0956), "h"))
})



test_that("get_greenampt_cyl_horiz_numerical works", {
  expect_equal(round(F_c_calculated, 5), set_units(c(0, 1, 5, 10, 20), "ft^2"))
})
