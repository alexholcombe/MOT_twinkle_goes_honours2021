rm(list = ls())
set.seed(190224)

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R","utils.R"))

stim_dir <- here("stimuli")
noise_dir   <- file.path(stim_dir, "noise")
gabor_dir   <- file.path(stim_dir, "gabor")

if(!dir.exists(noise_dir)) {
  dir.create(noise_dir, recursive = T)
}

if(!dir.exists(gabor_dir)) {
  dir.create(gabor_dir, recursive = T)
}

noise <- pink_noise_2d()
noise <- rescale_noise(noise, 0.1)
noise[1,1] <- 255
noise[1,2] <- 0
noise2 <- round(noise)

imager::save.image(imager::as.cimg(noise), file = file.path(noise_dir, "noise01.png"))
png::writePNG(noise2, file.path(noise_dir, "noise02.png"))

t_contr <- seq(0.1,1,by = 0.1)
for (i in 1:length(t_contr)) {
  gabor <- 128 + round(gabor_patch() * t_contr[i] * 128)  
  png::writePNG(gabor, file.path(gabor_dir, sprintf("gabor_contr%.1f.png", t_contr[i])))
}

