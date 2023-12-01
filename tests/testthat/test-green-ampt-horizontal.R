
#### Horizontal green ampt

library(units)
theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
Fcum <- set_units(20, "mm") # depth
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)
times <- set_units(c(0, 15, 30, 60), "min")
Fcum <- get_greenampt_horiz_flow(theta_0, theta_s, Ksat, h_b, h_0, times)
paste(round(Fcum,4), collapse = ",")

test_that("get_greenampt_horiz_flow works", {
  expect_equal(round(Fcum, 4), set_units(c(0,1.7009,2.4055,3.4019), "cm"))
})


# Double check that we get the original times back with the Green-Ampt equation
times_calc <- set_units(get_greenampt_horiz_time(theta_0, theta_s, Fcum, Ksat, h_b, h_0),"min")
# paste(round(times_calc,4), collapse = ",")

test_that("get_greenampt_horiz_time works", {
  expect_identical(round(times_calc,4), set_units(c(0,15,30,60), "min"))
})

# manually written equation
# set_units(Fcum^2 / 2 / (theta_s - theta_0) / Ksat / (h_b - h_0), "min")

### Integrated version

theta_0 <- 0.2 # unitless
theta_s <- 0.35 # unitless
Ksat <- set_units(0.2, "cm/h") # length / time
h_b <- set_units(6, "ft") # hydraulic head (length)
h_0 <- set_units(-10, "cm") # hydraulic head (length)

Fv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr"), d = NULL)

test_that("get_greenampt_horiz_flow_integrated works", {
  expect_identical(round(Fv,4), set_units(0.4653, "ft^2"))
})

# Get the infiltration over a 1 mm differential depth to compare with the point infiltration
Fv_dv <- get_greenampt_horiz_flow_integrated(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr"), d = set_units(1, "mm"))
Fv_point <- Fv_dv/ set_units(1, "mm")
# Fv_point


# Get the point infiltration
F_point <- get_greenampt_horiz_flow(theta_0, theta_s, Ksat, h_b, h_0, times = set_units(1,"hr")) %>% set_units("mm")
# F_point

test_that("get_greenampt_horiz_flow_integrated differential matches get_greenampt_horiz_flow", {
  expect_identical(as.numeric(round((Fv_point - F_point) / F_point * 100, 6)), -0.012963)
})

