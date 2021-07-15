# this scripts generate trajectories for the experiment
rm(list = ls())
set.seed(190327)

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R", "utils.R"))

traj_MOT_dir <- here("data","trajectories")
traj_plot_dir <- here("plots","trajectories")
if(!dir.exists(traj_MOT_dir)) { dir.create(traj_MOT_dir, recursive = T) }
if(!dir.exists(traj_plot_dir)) { dir.create(traj_plot_dir, recursive = T) }

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
  new_settings(show_labels = T)

# Make objects move for 8 seconds (0 to 8), adjust direction every 100 ms
timescale <- seq(0, 8, by = 0.1)

for (i in 1:40) {
  
  position <- generate_positions_random(8, sett_generate)
  
  # Object change direction smoothly
  trajectory_v <- 
    make_random_trajectory(position, timescale, sett_move, 
                           step_vonmises, kappa = 10)
  save_trajectory(trajectory_v,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  p <- plot_trajectory(trajectory_v, targets = 1:4)
  ggsave(file.path(traj_plot_dir,sprintf("T%03d.png",i)),p)
                  
}
