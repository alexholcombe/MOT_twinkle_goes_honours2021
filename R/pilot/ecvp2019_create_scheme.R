# create stimuli for poster


rm(list = ls())
set.seed(190224)

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)



source(here("R","utils.R"))

noise <- pink_noise_2d()
noise <- rescale_noise(noise, 0.1)
noise[1,1] <- 255
noise[1,2] <- 0

gabor <- 127 * gabor_patch(contr = 0.4,ppd = 76)

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
  new_settings(show_labels = F, target_img = gabor)

sett_show2 <-
  new_settings(show_labels = F, fill_target = NA, fill_object = NA, border_object = NA, border_target = "green", background = noise, target_img = gabor)

# set seed for replicability
set.seed(1001)

position <- generate_positions_random(8, sett_generate)

p_back    <- plot_position_image(position,sett_show2)
p_no_back <- plot_position(position,sett_show)

ggsave(here("plots", "scheme_noise.svg"), p_back)
ggsave(here("plots", "scheme_unif.svg"), p_no_back)
