# create gabors for MOT pilot
library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R", "utils.R"))

gabor_MOT_dir   <- here("data", "gabor_MOT_pilot")
if(!dir.exists(gabor_MOT_dir)) { dir.create(gabor_MOT_dir, recursive = T) }

# values of FD
beta <- 1.51
cT <- 0.04

dprimes <- c(4,8,12,16)
contrasts_for_mot <- dprime_fun(dprimes, cT, beta)

# generate gabors
set.seed(1902246)

gabor <- gabor_patch()
for (i in 1:length(dprimes)) {
  gaborx <- 128 + round(gabor * contrasts_for_mot[i] * 128)  
  png::writePNG(gaborx, file.path(gabor_MOT_dir, sprintf("gabor_contr_id_%d.png", i)))
}
