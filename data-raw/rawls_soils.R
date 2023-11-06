## code to prepare `rawls_soils` dataset goes here

conv_numeric <- function(x) {
  tryCatch(
    {
      x <- as.numeric(x)
      return(x)
    },
    error = function(e) e <- e,
    warning = function(w) w <- w,
    finally={
      return(x)
    }
  )
}

rawls_raw <- readr::read_csv("data-raw/rawls_et_al_table.csv")

# x <- "0.23"
# conv_numeric(x)

rawls_soils <- rawls_raw %>%
  mutate(across(everything(), function(x) gsub("(\\(.*\\))|(\\\n)","",x)),
         across(everything(), conv_numeric)) %>%
  dplyr::select(texture_class:pore_size_distrib_geometric)

usethis::use_data(rawls_soils, overwrite = TRUE)
