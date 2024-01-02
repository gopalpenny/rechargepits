

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


# paste(round(times, 4), collapse = ", ")
# paste(round(F_c_calculated, 5), collapse = ", ")


test_that("get_greenampt_cyl_horiz_time works", {
  expect_equal(round(times, 4), set_units(c(0, 0.4367, 7.5107, 22.7765, 64.5078), "h"))
})



test_that("get_greenampt_cyl_horiz_numerical works", {
  expect_equal(round(F_c_calculated, 5), set_units(c(0, 1, 5, 10, 20), "ft^2"))
})

r_b <- set_units(2, "ft") # length
theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
times <- set_units(c(0, 0.5, 1, 5), "hr")
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(4, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
d <- set_units(3, "ft")
num_sections <- 3
F_v_cum <- get_greenampt_cyl_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, r_b, times, F_units = "ft^3", num_sections = 3, d = d)

# paste(round(F_v_cum, 5), collapse = ", ")

test_that("get_greenampt_cyl_flow_integrated works", {
  expect_equal(round(F_v_cum, 5), set_units(c(0, 2.08077, 3.0067, 7.30069), "ft^3"))
})

# # This is equivalent to a discretized version of `get_greenampt_cyl_horiz_numerical`:
# Fc_1_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b = set_units(1.5,"ft"), h_0, r_b, times, F_units = "ft^2")
# Fc_2_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b = set_units(2.5,"ft"), h_0, r_b, times, F_units = "ft^2")
# Fc_3_5 <- get_greenampt_cyl_horiz_numerical(theta_0, theta_s, Ksat, h_b = set_units(3.5,"ft"), h_0, r_b, times, F_units = "ft^2")
# Fc_sum <- d / num_sections * (Fc_1_5 + Fc_2_5 + Fc_3_5)
# F_v_cum
# Fc_sum
