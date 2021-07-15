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

gabor <- imager::load.image(system.file("img", "Rlogo.png", package="png")) %>% imager::rm.alpha() 
alpha_gabor <-attr (gabor,"alpha") %>% imager::resize_halfXY()
gabor <- gabor %>% imager::grayscale() %>% imager::resize_halfXY()
gabor <- 127 * (gabor - mean(gabor))
gabor <- gabor*alpha_gabor

sett_generate <-
  new_settings(xlim = c(-7, 7), ylim = c(-7, 7), min_distance = 2,
               arena_shape = "circle")
# they can move in a circular arena, bounce off borders and other objects
sett_move <-
  new_settings(speed = 5, xlim = c(-9, 9), ylim = c(-9, 9),
               bounce_off_square = F,
               bounce_off_circle = T, circle_bounce_jitter = pi / 6)
# here you can separately adjust settings for visualization
sett_show <-
  new_settings(show_labels = F, background = noise, target_img = gabor)

sett_show2 <-
  new_settings(show_labels = F, fill_target = NA, fill_object = NA, border_object = NA, border_target = "green", background = noise, target_img = gabor)

# set seed for replicability
set.seed(1001)

position <- generate_positions_random(8, sett_generate)

plot_position_image(position,sett_show2,targets = 1:4)

# our starting positions
plot_position(position)
plot_position(position, sett_show)

# Make objects move for 8 seconds (0 to 8), adjust direction every 100 ms
timescale <- seq(0, 8, by = 0.1)

# Objects only bounce
#trajectory_d <- 
#  make_random_trajectory(position, timescale, sett_move, 
#                         step_direct)

# Object change direction randomly every 0.5-1.5 seconds
#trajectory_z <- 
#  make_random_trajectory(position, timescale, sett_move, 
#                         step_zigzag, ttt = c(.5, 1.5), syncstart = F)

# Object change direction smoothly
trajectory_v <- 
  make_random_trajectory(position, timescale, sett_move, 
                         step_vonmises, kappa = 16)


plot_trajectory(trajectory_v, sett_show)

animation::ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.8-Q16/ffmpeg")

#render_trajectory_video("trajectory_v.mp4", trajectory_v, 
#                        sett_show, targets = 1:4
#)

#render_trajectory_video_additive("trajectory_v.mp4", trajectory_v, 
#                        sett_show2, targets = 1:4
#)

gabor <- 127 * gabor_patch()
#gabor <- 127 * (gabor - mean(gabor))
sett_show2 <-
  new_settings(show_labels = F, fill_target = NA, fill_object = NA, border_object = NA, border_target = "green", background = noise, target_img = gabor)
render_trajectory_video_additive("trajectory_gabors.mp4", trajectory_v, 
                                 sett_show2, targets = 1:4, targets_cue_only = T
)

gabor <- 127 * gabor_patch(bandwidth = 0.8) * 0.2

sett_show2 <-
  new_settings(show_labels = F, fill_target = NA, fill_object = NA, border_object = NA, border_target = "green", background = noise, target_img = gabor)
render_trajectory_video_additive("trajectory_gabors_band0_8_contr0_2.mp4", trajectory_v, 
                                 sett_show2, targets = 1:4, targets_cue_only = T
)
