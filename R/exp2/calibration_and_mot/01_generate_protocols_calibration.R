library(tidyverse)
library(here)

source(here("R/utils.R"))

out_pth <- here("data/exp2/detection/calibration/")
if(!dir.exists(out_pth)) {
  dir.create(out_pth,recursive = T)
}
p <- create_protocol_exp2_detection_calibration()
write_csv2(p, path = file.path(out_pth,"design_calibration.csv"))

psychopy_pth <- here("psychopy/exp2/detection/calibration")
file.copy( file.path(out_pth,"design_calibration.csv"), 
           file.path(psychopy_pth,"design_calibration.csv"),overwrite = T)
