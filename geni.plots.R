
library(dplyr)
library(geni.plots)

geni.plots::fig_region_stack(
  data = geni.plots::geni_test_stack_region$assoc,
  traits = c("Interleukin-6 levels", "Interleukin-6 receptor levels"),
  corr = geni.plots::geni_test_stack_region$corr,
  build = 37,
  highlights = "rs11265611",
  title_center = TRUE
)
