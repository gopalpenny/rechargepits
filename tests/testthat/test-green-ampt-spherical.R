

library(units)
r_b <- set_units(2, "ft") # length
theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
F_s <- set_units(c(0, 1, 5, 10, 20), "ft^3") # units of length^2
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- get_greenampt_hsphere_time(theta_0, theta_s, F_s, Ksat, h_b, h_0, r_b)
F_s_calculated <- get_greenampt_hsphere_numerical(theta_0, theta_s, Ksat, h_b, h_0, r_b, times)


paste(round(F_s_calculated, 5), collapse = ", ")


test_that("get_greenampt_hsphere_time works", {
  expect_equal(round(times, 4), set_units(c(0, 0.1088, 1.8141, 5.3234, 14.3918), "h"))
})



test_that("get_greenampt_hsphere_numerical works", {
  expect_equal(round(F_s_calculated, 5), set_units(c(0, 1, 4.99999, 9.99998, 20.00001), "ft^3"))
})
