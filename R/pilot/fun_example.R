# this scripts generate trajectories for the experiment
rm(list = ls())
set.seed(190224)

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R","utils.R"))

noise <- pink_noise_2d()
max(noise)
min(noise)
noise <- rescale_noise(noise, 0.1)
noise[1,1] <- 255
noise[1,2] <- 0

fanda <- imager::load.image("stimuli_faces/fanda.png")
fanda <- imager::rm.alpha(fanda)

alpha_fanda <-attr(fanda,"alpha")%>% imager::resize_halfXY() %>% imager::imresize(1/10)

fanda <- fanda %>% imager::grayscale() %>% imager::resize_halfXY() %>% imager::imresize(1/10)

fanda <- 127 * (fanda - mean(fanda))
fanda <- fanda*alpha_fanda

sett_generate <-
  new_settings(xlim = c(-7, 7), ylim = c(-7, 7), min_distance = 2,
               arena_shape = "circle")
# they can move in a circular arena, bounce off borders and other objects
sett_move <-
  new_settings(speed = 5, xlim = c(-9, 9), ylim = c(-9, 9),
               bounce_off_square = F,
               bounce_off_circle = T, circle_bounce_jitter = pi / 6)
# here you can separately adjust settings for visualization

sett_show2 <-
  new_settings(show_labels = F, fill_target = NA, fill_object = NA, border_object = NA, border_target = "green", background = noise, target_img = fanda)

# set seed for replicability
set.seed(1001)

position <- generate_positions_random(8, sett_generate)
timescale <- seq(0, 8, by = 0.1)
trajectory_v <- 
  make_random_trajectory(position, timescale, sett_move, 
                         step_vonmises, kappa = 16)
animation::ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.8-Q16/ffmpeg")

render_trajectory_video_additive("trajectory_MOF.mp4", trajectory_v, 
                        sett_show2, targets = 1:4, targets_cue_only = T
)