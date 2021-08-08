#Trajectories for 4 objects in 2 cages.
#Monitor size (1440, 900) i.e., x: (-720, 720) y: (-450, 450)
#Object size: Radius: 28 pix

setwd("~/Documents/OneDrive - The University of Sydney (Students)/Empirical Thesis/Code/MOT_twinkle_goes_honours2021-main")

rm(list = ls())
set.seed(4)

library(tidyverse)

library(dplyr)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R", "utils.R"))

traj_MOT_dir <- here("exp2_JJC","mac13_trajectories_2cages_4objects")
if(!dir.exists(traj_MOT_dir)) { dir.create(traj_MOT_dir, recursive = T) }

no_objects = 4

no_practice_trials = 8

no_real_trials = 1000

no_trajectories = no_practice_trials + no_real_trials

timeGrain <- .01666666666666666666667 #Speed, direction are specified by R every timeGrain seconds, which gets
#turned into a trajectory file that specifies position every timeGrain seconds

sett_generate <-
  new_settings(xlim = c(-620, 620), ylim = c(-400, 400), min_distance = 150,
               arena_shape = "square")

# they can move in a square arena, bounce off borders and other objects
sett_move <-
  new_settings(speed = 400, xlim = c(-620, 620), ylim = c(-400, 400),
               bounce_off_square = T,
               bounce_off_circle = F, circle_bounce_jitter = pi / 6, bounce_off_others = TRUE, min_distance = 64)
# here you can separately adjust settings for visualization
sett_show = sett_move

sett_show$show_labels = T

##  VERTICAL MIDLINE TRAJECTORIES ##

#custom cage starting coordinates
random_coords_in_square <- function(n, xmin, xmax, ymin, ymax) {
  hemifield_order <- seq(1,2)
  hemifield_order = sample(hemifield_order)
  t <- tibble()
  #targets represent first four items, need to ensure one in each cage
  #first run through loop for targets, second for distractors
  for (target_vs_distractors in 1:2) {
    for (i in 1:2) {
      h = hemifield_order[i]
      if (h == 1) {
        x = stats::runif(1, min = -620, max = -50)
        y = stats::runif(1, min = -400, max = 400)
      }
      else if (h == 2){
        x = stats::runif(1, min = 50, max = 620)
        y = stats::runif(1, min = -400, max = 400)
      }
      this_object = tibble(x = x, y = y)
      t = rbind(t, this_object)
    }
  }
  t
}

generate_positions_cages <- function(n, settings, check_distance = T, border_distance = 0) {
  xlim <- settings$xlim
  ylim <- settings$ylim
  stopifnot(!is.null(xlim))
  stopifnot(!is.null(ylim))
  shapes <- c("square", "circle", "donut")
  shape <- pmatch(settings$arena_shape, shapes)
  if (is.na(shape)) {
    stop("Arena shape unknown")
  }
  while (T) {
    p <- switch(
      shapes[shape],
      square = random_coords_in_square( ),
      circle = random_coords_in_circle(
        n,
        mean(xlim), mean(ylim),
        sum(c(diff(xlim), diff(ylim)) / 4) - border_distance
      ),
      donut = random_coords_in_donut(
        n,
        mean(xlim), mean(ylim),
        sum(c(diff(xlim), diff(ylim)) / 4) - border_distance,
        settings$arena_inside_radius + border_distance
      )
    )
    p <- p %>%
      dplyr::mutate(object = 1:n) %>%
      dplyr::select(.data$object, dplyr::everything())
    if (check_distance) {
      #if (is_distance_at_least(p, min_distance = settings$min_distance)) {
      if (is_distance_at_least(p, min_distance = 150)) {
        return(p)
      }
    } else {
      return(p)
    }
  }
}

position <- generate_positions_cages(no_objects, sett_generate)

# our starting positions
plot_position(position, sett_move)

extrapolate_moment <- function(moment, timestep, cur_time, new_time = NULL) {
  s <- cur_speed(moment, cur_time)
  if (is.null(new_time)) new_time <- cur_time + timestep
  moment_next <- moment %>%
    dplyr::mutate(
      x = .data$x + cos(.data$direction) * s * timestep,
      y = .data$y + sin(.data$direction) * s * timestep,
      time = new_time
    )
  moment_next
}

bounce_off_square <- function(moment, time_next, settings) {
  time_now <- moment$time[1]
  timestep <- time_next - time_now
  # extrapolate future
  moment_next <-
    extrapolate_moment(moment, timestep, time_now, NA)
  
  for (obj in seq(1,length(moment_next$x))) {
    if (moment$x[obj] < 0) {
      cage_xlim = c(-620, -50)
      cage_ylim = c(-400, 400)
    }
    else if (moment$x[obj] > 0) {
      cage_xlim = c(50, 620)
      cage_ylim = c(-400, 400)
    }
    #print(cat(obj,cage_xlim))
    
    # check sides
    #Check whether each of the objects have exceeded their boundaries
    beyond_left <- moment_next$x[obj] < min(cage_xlim)
    beyond_right <- moment_next$x[obj] > max(cage_xlim)
    beyond_top <- moment_next$y[obj] > max(cage_ylim)
    beyond_bottom <- moment_next$y[obj] < min(cage_ylim)
    # if more than one, then just reverse
    
    # Do the bounce for those who need to be bounced
    beyond_more <- (beyond_left + beyond_right + beyond_top + beyond_bottom) > 1
    if (beyond_more) {
      moment$direction[obj] <-
        (moment$direction[obj] + pi) %% (2 * pi)
    } else if  (beyond_left | beyond_right) {#bounce from left/right
      moment$direction[obj] <-  (pi - moment$direction[obj]) %% (2 * pi)
    } else if (beyond_top | beyond_bottom) { # bounce from top/bottom
      moment$direction[obj] <- (2 * pi - moment$direction[obj]) %% (2 * pi)
    }
  } #end iterating through all of the objects
  moment
}

#temporary to debug below function line by line
#start <- position
#settings <- sett_move
#step_function <- step_zigzag
#    my_trajectory <- make_random_trajectoryJJC(position, times, sett_move, step_zigzag, ttt = c(.5, 1.5), syncstart = F)

make_random_trajectoryJJC <- function(start, times, settings, step_function, ...) {
  # take start position
  # are direction/speed present? - add if not
  # make step
  had_a_bounce_in_prohibited_interval <- FALSE
  changed_direction_in_final_interval <- FALSE
  
  trajectory_duration <- tail(times,1)
  
  
  #Determining length of trajectory (multipling by number of objects)
  trajectoryLength <- ((length(times) - 1) * no_objects) + 1
  
  duration_of_final_interval_when_bounces_are_prohibited <- 0.4 #Don't allow any bounces off walls, corners, or other balls in the last interval of the trial
  stopifnot(length(times) > 1)
  if ("direction" %in% names(start)) {
    moment <- start
  } else {
    moment <- start %>% add_random_direction()
  }
  if (!"speed" %in% names(moment)) {
    moment <- moment %>% dplyr::mutate(speed = settings$speed)
  }
  moment <- moment %>% dplyr::mutate(time = times[1])
  
  # moment contains position, direction, and speed
  # here we have tibble of moments with two columns:
  #   - time (from times)
  #   - position (embedded moment tibble for each time)
  # Later we expand tibbles into trajectory tibble
  moment_tbl <- tibble::tibble(time = times, position = list(tibble::tibble))
  moment_tbl$position[[1]] <- moment
  for (i in 2:length(times)) {
    
    if (settings$bounce_off_square) {
      bounced_off_wall <- FALSE #Bug where if bouncing of others, will not obey beyond more, this work around will lead to some eventual overlap between objects in corners
      old_moment = moment
      moment <- bounce_off_square(moment, times[i], settings)
      if (identical(old_moment$direction, moment$direction) == FALSE) {
        bounced_off_wall = TRUE
      }
    }
    
    if (settings$bounce_off_circle) {
      moment <- bounce_off_circle(moment, times[i], settings)
    }
    if (settings$bounce_off_inside) {
      moment <- bounce_off_inside(moment, times[i], settings)
    }
    if (!bounced_off_wall && settings$bounce_off_others) { #only accounts for 
      bounced_off_others <- FALSE
      old_moment = moment
      moment <- bounce_off_others(moment, times[i], settings)
      if (identical(old_moment$direction, moment$direction) == FALSE) {
        bounced_off_others <- TRUE
      }
    }
    
    #Only continue if any bounces that have occurred were not in the last no_bounces_in_last_interval. Otherwise, exit with error.
    # trajectory_duration <- tail(times,1)
    # if (times[i] >= trajectory_duration - duration_of_final_interval_when_bounces_are_prohibited) { #We are now in the interval during which bounces are prohibited
    #   if (bounced_off_wall | bounced_off_others) {
    #     had_a_bounce_in_prohibited_interval <- TRUE
    #     #print("It's a real bounce!")
    #   }
    # }
    
    moment_next <- step_function(moment, times[i]) #AH temporarily truncate off the following so can debug line-by-line:                settings, ...)
    
    moment_tbl$position[[i]] <- moment_next #moment_tbl is the concatenation of all individual moments (individual time steps)
    moment <- moment_next
  }
  my_moment_tbl <- moment_tbl %>%
    dplyr::select(-.data$time) %>% #This gets rid of the time variable - don't know why we want to do that, though
    tidyr::unnest(cols = c("position"))
  #return (list("moment_tbl" = my_moment_tbl, "had_a_bounce_in_prohibited_interval" = had_a_bounce_in_prohibited_interval))
  return (list("moment_tbl" = my_moment_tbl))
}


plot_trajectory <- function(trajectory,
                            settings = default_settings(),
                            targets = NULL,
                            ...)
{
  start_time <- min(trajectory$time)
  fig <- plot_position(trajectory %>% dplyr::filter(.data$time == start_time),
                       settings = settings,
                       targets = targets,
                       ...
  )
  fig <- fig +
    geom_path(
      aes_string(
        x = "x", y = "y", group = "object",
        fill = NULL, colour = "object"
      ),
      data = trajectory %>% dplyr::arrange(.data$time),
      #colour = "blue"
    ) +
    ggforce::geom_circle(aes_string(r = "settings$r"))
  
  if (settings$show_labels) {
    fig <-
      fig +
      ggplot2::geom_text(
        ggplot2::aes_string(x = "x", y = "y", label = "object"),
        colour = I("red")
      )
  }
  fig
}

step_zigzag <- function(moment, time_next, settings,
                        ttt = c(1, 1.5), syncstart = F) {
  time_now <- moment$time[1]
  timestep <- time_next - time_now
  # ttt is time-to-turn - when object should turn randomly
  n <- nrow(moment)
  if (!"ttt" %in% names(moment)) {
    times <- stats::runif(n, min = ttt[1], max = ttt[2])
    if (!syncstart) {
      passed <- stats::runif(n, min = 0, max = ttt[1])
      times <- times - passed
    }
    moment <- moment %>% dplyr::mutate(ttt = times)
  }
  which_turn <- moment$ttt < time_next
  if (any(which_turn)) {
    moment$direction[which_turn] <-
      stats::runif(sum(which_turn), 0, 2 * pi)
    moment$ttt[which_turn] <-
      moment$ttt[which_turn] +
      stats::runif(sum(which_turn), min = ttt[1], max = ttt[2])
    #print("time-to-turn reached!")
  }
  moment_next <-
    extrapolate_moment(moment, timestep, time_now, time_next)
  moment_next
}

for (i in 1:(no_trajectories/2)) {
  
  #print(i)
  
  position <- generate_positions_cages(no_objects, sett_generate)
  
  # Make objects move for X secs, adjust direction every timeGrain ms
  times <- seq(0, runif(1, 1.2016, 3.016), by = timeGrain)
  
  # Object change direction randomly every 0.5-1.5 seconds
  #But each time we calculate a trajectory, we check whether it had a bounce in the last X seconds and reject it if so and calculate a different one
  had_a_bounce_in_prohibited_interval <- TRUE
  changed_direction_in_final_interval <- TRUE
  attempts = 0; max_attempts <- 50
  while (attempts<max_attempts && changed_direction_in_final_interval == TRUE) {
    my_trajectory <-
      make_random_trajectoryJJC(position, times, sett_move,
                                step_zigzag, ttt = c(1, 1.5), syncstart = F)
    
    total_time_steps_needed <- length(times)
    
    duration_of_final_interval_when_bounces_are_prohibited <- 0.4
    steps_in_prohibited_interval <- duration_of_final_interval_when_bounces_are_prohibited/timeGrain
    #print(steps_in_prohibited_interval)
    if (round(steps_in_prohibited_interval) - steps_in_prohibited_interval != 0){
      print("timeGrain does not go evenly into the interval you are trying to set, so you better change timeGrain or duration_of_each_segment_of_last_portion_of_trial")
    }
    # print("Steps in prohibited interval: ")
    # print(steps_in_prohibited_interval)
    
    #had_a_bounce_in_prohibited_interval <- my_trajectory$had_a_bounce_in_prohibited_interval
    directions_the_same_all_times <- TRUE
    
    trajectory_duration <- tail(times,1)
    
    #print("trajectory duration:")
    #print(trajectory_duration)
    
    for (b in 1:1) { #step through the how_many_objects moving objects, check whether the direction is always the same for each
      #check whether there is more than one direction for this object
      final_index_this_object <- b + no_objects*total_time_steps_needed - 2*no_objects
      #print(final_index_this_object)
      start_time_step_of_prohibited_interval <- b + length(times)*no_objects - steps_in_prohibited_interval*no_objects - 2*no_objects
      idxs_this_obj_all_times_in_prohibited_interval <- seq(start_time_step_of_prohibited_interval, final_index_this_object, by=no_objects)
      #print("b:")
      #print(b)
      
      z <- my_trajectory$moment_tbl
      
      #all_directions <- z$direction[idxs_this_obj_all_times_in_prohibited_interval]
      #print("All directions:")
      #print(all_directions)
      
      directions <- unique(z$direction[idxs_this_obj_all_times_in_prohibited_interval])
      
      #print("Unique directions:")
      #print(directions)
      
      if (length(directions) > 1) {
        directions_the_same_all_times <- FALSE
        #print("Bounce detected, directions =")
        #print(length(directions))
      }
    }
    changed_direction_in_final_interval <- !directions_the_same_all_times
    attempts <- attempts + 1
    #print("attempts")
    #print(attempts)
  } #end while loop running to get a trajectory that doesnt have a bounce in the final portions
  trajectory_z <- my_trajectory$moment_tbl
  trajectory_last <- trajectory_z[start_time_step_of_prohibited_interval : nrow(trajectory_z), ]
  
  if (attempts==max_attempts){
    print("Sorry, could not generate a trajectory without a bounce in the prohibited final interval on this trial. Trial number =")
    print(i) } else {
      print(i)
    }
  
  save_trajectory(trajectory_z,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  p <- plot_trajectory(trajectory_z, targets = 1:4, sett_show)
  #p_last <- plot_trajectory(trajectory_last, targets = 1:4, sett_show)
  
  # # Object change direction smoothly
  # trajectory_d <-
  #   make_random_trajectoryJJC(position, times, sett_move,
  #                             step_direct)
  # 
  # save_trajectory(trajectory_d,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  # p <- plot_trajectory(trajectory_d, targets = 1:4, sett_show)
  
  # # Object change direction smoothly
  # trajectory_v <-
  #   make_random_trajectoryJJC(position, times, sett_move,
  #                             step_vonmises, kappa = 10)
  
  # save_trajectory(trajectory_v,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  # p <- plot_trajectory(trajectory_v, targets = 1:4, sett_show)
  #ggsave(file.path(traj_plot_dir,sprintf("T%03d.png",i)),p)
  
}

p

##    HORIZONTAL MIDLINE TRAJECTORIES   ##

random_coords_in_square <- function(n, xmin, xmax, ymin, ymax) {
  hemifield_order <- seq(1,2)
  hemifield_order = sample(hemifield_order)
  t <- tibble()
  #targets represent first four items, need to ensure one in each cage
  #first run through loop for targets, second for distractors
  for (target_vs_distractors in 1:2) {
    for (i in 1:2) {
      h = hemifield_order[i]
      if (h == 1) {
        x = stats::runif(1, min = -620, max = 620)
        y = stats::runif(1, min = -400, max = -50)
      }
      else if (h == 2){
        x = stats::runif(1, min = -620, max = 620)
        y = stats::runif(1, min = 50, max = 400)
      }
      this_object = tibble(x = x, y = y)
      t = rbind(t, this_object)
    }
  }
  t
}

generate_positions_cages <- function(n, settings, check_distance = T, border_distance = 0) {
  xlim <- settings$xlim
  ylim <- settings$ylim
  stopifnot(!is.null(xlim))
  stopifnot(!is.null(ylim))
  shapes <- c("square", "circle", "donut")
  shape <- pmatch(settings$arena_shape, shapes)
  if (is.na(shape)) {
    stop("Arena shape unknown")
  }
  while (T) {
    p <- switch(
      shapes[shape],
      square = random_coords_in_square( ),
      circle = random_coords_in_circle(
        n,
        mean(xlim), mean(ylim),
        sum(c(diff(xlim), diff(ylim)) / 4) - border_distance
      ),
      donut = random_coords_in_donut(
        n,
        mean(xlim), mean(ylim),
        sum(c(diff(xlim), diff(ylim)) / 4) - border_distance,
        settings$arena_inside_radius + border_distance
      )
    )
    p <- p %>%
      dplyr::mutate(object = 1:n) %>%
      dplyr::select(.data$object, dplyr::everything())
    if (check_distance) {
      #if (is_distance_at_least(p, min_distance = settings$min_distance)) {
      if (is_distance_at_least(p, min_distance = 150)) {
        return(p)
      }
    } else {
      return(p)
    }
  }
}

position <- generate_positions_cages(no_objects, sett_generate)

# our starting positions
plot_position(position, sett_move)

extrapolate_moment <- function(moment, timestep, cur_time, new_time = NULL) {
  s <- cur_speed(moment, cur_time)
  if (is.null(new_time)) new_time <- cur_time + timestep
  moment_next <- moment %>%
    dplyr::mutate(
      x = .data$x + cos(.data$direction) * s * timestep,
      y = .data$y + sin(.data$direction) * s * timestep,
      time = new_time
    )
  moment_next
}

bounce_off_square <- function(moment, time_next, settings) {
  time_now <- moment$time[1]
  timestep <- time_next - time_now
  # extrapolate future
  moment_next <-
    extrapolate_moment(moment, timestep, time_now, NA)
  
  for (obj in seq(1,length(moment_next$x))) {
    if (moment$y[obj] < 0) {
      cage_xlim = c(-620, 620)
      cage_ylim = c(-400, -50)
    }
    else if (moment$y[obj] > 0) {
      cage_xlim = c(-620, 620)
      cage_ylim = c(50, 400)
    }
    #print(cat(obj,cage_xlim))
    
    # check sides
    #Check whether each of the objects have exceeded their boundaries
    beyond_left <- moment_next$x[obj] < min(cage_xlim)
    beyond_right <- moment_next$x[obj] > max(cage_xlim)
    beyond_top <- moment_next$y[obj] > max(cage_ylim)
    beyond_bottom <- moment_next$y[obj] < min(cage_ylim)
    # if more than one, then just reverse
    
    # Do the bounce for those who need to be bounced
    beyond_more <- (beyond_left + beyond_right + beyond_top + beyond_bottom) > 1
    if (beyond_more) {
      moment$direction[obj] <-
        (moment$direction[obj] + pi) %% (2 * pi)
    } else if  (beyond_left | beyond_right) {#bounce from left/right
      moment$direction[obj] <-  (pi - moment$direction[obj]) %% (2 * pi)
    } else if (beyond_top | beyond_bottom) { # bounce from top/bottom
      moment$direction[obj] <- (2 * pi - moment$direction[obj]) %% (2 * pi)
    }
  } #end iterating through all of the objects
  moment
}

#temporary to debug below function line by line
#start <- position
#settings <- sett_move
#step_function <- step_zigzag
#    my_trajectory <- make_random_trajectoryJJC(position, times, sett_move, step_zigzag, ttt = c(.5, 1.5), syncstart = F)

make_random_trajectoryJJC <- function(start, times, settings, step_function, ...) {
  # take start position
  # are direction/speed present? - add if not
  # make step
  had_a_bounce_in_prohibited_interval <- FALSE
  changed_direction_in_final_interval <- FALSE
  
  trajectory_duration <- tail(times,1)
  
  
  #Determining length of trajectory (multipling by number of objects)
  trajectoryLength <- ((length(times) - 1) * no_objects) + 1
  
  duration_of_final_interval_when_bounces_are_prohibited <- 0.4 #Don't allow any bounces off walls, corners, or other balls in the last interval of the trial
  stopifnot(length(times) > 1)
  if ("direction" %in% names(start)) {
    moment <- start
  } else {
    moment <- start %>% add_random_direction()
  }
  if (!"speed" %in% names(moment)) {
    moment <- moment %>% dplyr::mutate(speed = settings$speed)
  }
  moment <- moment %>% dplyr::mutate(time = times[1])
  
  # moment contains position, direction, and speed
  # here we have tibble of moments with two columns:
  #   - time (from times)
  #   - position (embedded moment tibble for each time)
  # Later we expand tibbles into trajectory tibble
  moment_tbl <- tibble::tibble(time = times, position = list(tibble::tibble))
  moment_tbl$position[[1]] <- moment
  for (i in 2:length(times)) {
    
    if (settings$bounce_off_square) {
      bounced_off_wall <- FALSE #Bug where if bouncing of others, will not obey beyond more, this work around will lead to some eventual overlap between objects in corners
      old_moment = moment
      moment <- bounce_off_square(moment, times[i], settings)
      if (identical(old_moment$direction, moment$direction) == FALSE) {
        bounced_off_wall = TRUE
      }
    }
    
    if (settings$bounce_off_circle) {
      moment <- bounce_off_circle(moment, times[i], settings)
    }
    if (settings$bounce_off_inside) {
      moment <- bounce_off_inside(moment, times[i], settings)
    }
    if (!bounced_off_wall && settings$bounce_off_others) { #only accounts for 
      bounced_off_others <- FALSE
      old_moment = moment
      moment <- bounce_off_others(moment, times[i], settings)
      if (identical(old_moment$direction, moment$direction) == FALSE) {
        bounced_off_others <- TRUE
      }
    }
    
    #Only continue if any bounces that have occurred were not in the last no_bounces_in_last_interval. Otherwise, exit with error.
    # trajectory_duration <- tail(times,1)
    # if (times[i] >= trajectory_duration - duration_of_final_interval_when_bounces_are_prohibited) { #We are now in the interval during which bounces are prohibited
    #   if (bounced_off_wall | bounced_off_others) {
    #     had_a_bounce_in_prohibited_interval <- TRUE
    #     #print("It's a real bounce!")
    #   }
    # }
    
    moment_next <- step_function(moment, times[i]) #AH temporarily truncate off the following so can debug line-by-line:                settings, ...)
    
    moment_tbl$position[[i]] <- moment_next #moment_tbl is the concatenation of all individual moments (individual time steps)
    moment <- moment_next
  }
  my_moment_tbl <- moment_tbl %>%
    dplyr::select(-.data$time) %>% #This gets rid of the time variable - don't know why we want to do that, though
    tidyr::unnest(cols = c("position"))
  #return (list("moment_tbl" = my_moment_tbl, "had_a_bounce_in_prohibited_interval" = had_a_bounce_in_prohibited_interval))
  return (list("moment_tbl" = my_moment_tbl))
}


plot_trajectory <- function(trajectory,
                            settings = default_settings(),
                            targets = NULL,
                            ...)
{
  start_time <- min(trajectory$time)
  fig <- plot_position(trajectory %>% dplyr::filter(.data$time == start_time),
                       settings = settings,
                       targets = targets,
                       ...
  )
  fig <- fig +
    geom_path(
      aes_string(
        x = "x", y = "y", group = "object",
        fill = NULL, colour = "object"
      ),
      data = trajectory %>% dplyr::arrange(.data$time),
      #colour = "blue"
    ) +
    ggforce::geom_circle(aes_string(r = "settings$r"))
  
  if (settings$show_labels) {
    fig <-
      fig +
      ggplot2::geom_text(
        ggplot2::aes_string(x = "x", y = "y", label = "object"),
        colour = I("red")
      )
  }
  fig
}

step_zigzag <- function(moment, time_next, settings,
                        ttt = c(1, 1.5), syncstart = F) {
  time_now <- moment$time[1]
  timestep <- time_next - time_now
  # ttt is time-to-turn - when object should turn randomly
  n <- nrow(moment)
  if (!"ttt" %in% names(moment)) {
    times <- stats::runif(n, min = ttt[1], max = ttt[2])
    if (!syncstart) {
      passed <- stats::runif(n, min = 0, max = ttt[1])
      times <- times - passed
    }
    moment <- moment %>% dplyr::mutate(ttt = times)
  }
  which_turn <- moment$ttt < time_next
  if (any(which_turn)) {
    moment$direction[which_turn] <-
      stats::runif(sum(which_turn), 0, 2 * pi)
    moment$ttt[which_turn] <-
      moment$ttt[which_turn] +
      stats::runif(sum(which_turn), min = ttt[1], max = ttt[2])
    #print("time-to-turn reached!")
  }
  moment_next <-
    extrapolate_moment(moment, timestep, time_now, time_next)
  moment_next
}

failed_trajectories <- c()

for (i in (no_trajectories/2 + 1):no_trajectories) {
  
  #print(i)
  
  position <- generate_positions_cages(no_objects, sett_generate)
  
  # Make objects move for X secs, adjust direction every timeGrain ms
  times <- seq(0, runif(1, 1.2016, 3.016), by = timeGrain)
  
  # Object change direction randomly every 0.5-1.5 seconds
  #But each time we calculate a trajectory, we check whether it had a bounce in the last X seconds and reject it if so and calculate a different one
  had_a_bounce_in_prohibited_interval <- TRUE
  changed_direction_in_final_interval <- TRUE
  attempts = 0; max_attempts <- 50
  while (attempts<max_attempts && changed_direction_in_final_interval == TRUE) {
    my_trajectory <-
      make_random_trajectoryJJC(position, times, sett_move,
                                step_zigzag, ttt = c(1, 1.5), syncstart = F)
    
    total_time_steps_needed <- length(times)
    
    duration_of_final_interval_when_bounces_are_prohibited <- 0.4
    steps_in_prohibited_interval <- duration_of_final_interval_when_bounces_are_prohibited/timeGrain
    #print(steps_in_prohibited_interval)
    if (round(steps_in_prohibited_interval) - steps_in_prohibited_interval != 0){
      print("timeGrain does not go evenly into the interval you are trying to set, so you better change timeGrain or duration_of_each_segment_of_last_portion_of_trial")
    }
    # print("Steps in prohibited interval: ")
    # print(steps_in_prohibited_interval)
    
    #had_a_bounce_in_prohibited_interval <- my_trajectory$had_a_bounce_in_prohibited_interval
    directions_the_same_all_times <- TRUE
    
    trajectory_duration <- tail(times,1)
    
    #print("trajectory duration:")
    #print(trajectory_duration)
    
    for (b in 1:1) { #step through the how_many_objects moving objects, check whether the direction is always the same for each
      #check whether there is more than one direction for this object
      final_index_this_object <- b + no_objects*total_time_steps_needed - 2*no_objects
      #print(final_index_this_object)
      start_time_step_of_prohibited_interval <- b + length(times)*no_objects - steps_in_prohibited_interval*no_objects - 2*no_objects
      idxs_this_obj_all_times_in_prohibited_interval <- seq(start_time_step_of_prohibited_interval, final_index_this_object, by=no_objects)
      #print("b:")
      #print(b)
      
      z <- my_trajectory$moment_tbl
      
      #all_directions <- z$direction[idxs_this_obj_all_times_in_prohibited_interval]
      #print("All directions:")
      #print(all_directions)
      
      directions <- unique(z$direction[idxs_this_obj_all_times_in_prohibited_interval])
      
      #print("Unique directions:")
      #print(directions)
      
      if (length(directions) > 1) {
        directions_the_same_all_times <- FALSE
        #print("Bounce detected, directions =")
        #print(length(directions))
      }
    }
    changed_direction_in_final_interval <- !directions_the_same_all_times
    attempts <- attempts + 1
    #print("attempts")
    #print(attempts)
  } #end while loop running to get a trajectory that doesnt have a bounce in the final portions
  trajectory_z <- my_trajectory$moment_tbl
  trajectory_last <- trajectory_z[start_time_step_of_prohibited_interval : nrow(trajectory_z), ]
  
  if (attempts==max_attempts){
    print("Sorry, could not generate a trajectory without a bounce in the prohibited final interval on this trial. Trial number =")
    print(i) } else {
      print(i)
    }
  
  save_trajectory(trajectory_z,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  p <- plot_trajectory(trajectory_z, targets = 1:4, sett_show)
  #p_last <- plot_trajectory(trajectory_last, targets = 1:4, sett_show)
  
  # # Object change direction smoothly
  # trajectory_d <-
  #   make_random_trajectoryJJC(position, times, sett_move,
  #                             step_direct)
  # 
  # save_trajectory(trajectory_d,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  # p <- plot_trajectory(trajectory_d, targets = 1:4, sett_show)
  
  # # Object change direction smoothly
  # trajectory_v <-
  #   make_random_trajectoryJJC(position, times, sett_move,
  #                             step_vonmises, kappa = 10)
  
  # save_trajectory(trajectory_v,filename = file.path(traj_MOT_dir,sprintf("T%03d.csv",i)))
  # p <- plot_trajectory(trajectory_v, targets = 1:4, sett_show)
  #ggsave(file.path(traj_plot_dir,sprintf("T%03d.png",i)),p)
  
}

p
