library(tidyverse)
library(here)

out_pth <- here("data/exp2/detection/range_estim/")
if(!dir.exists(out_pth)) {
  dir.create(out_pth,recursive = T)
}
p <- create_protocol_exp2_detection_range_estimation()
write_csv2(p, path = file.path(out_pth,"design_range_estim.csv"))

psychopy_pth <- here("psychopy/exp2/detection/range_estim")
file.copy( file.path(out_pth,"design_range_estim.csv"), 
           file.path(psychopy_pth,"design_range_estim.csv"),overwrite = T)
