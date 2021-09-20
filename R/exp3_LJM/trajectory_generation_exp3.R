# Install devtools package if necessary
if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")

# Install the stable verion from GitHub
#devtools::install_github("jirilukavsky/motrack")
library(motrack)
library(here)
library(tidyverse)
library(gtools)

timeGrain<- .004166666666666666666 #Speed, direction are specified by R every timeGrain seconds, which gets
#turned into a trajectory file that specifies position every timeGrain seconds

#set number of objects in trajectory here
how_many_objects <- 1

#starting location settings
sett_generate <-
  new_settings(xlim = c(-600.0000000000001, -600), ylim = c(150, 150.0000000000001), min_distance = 120,
               arena_shape = "circle")

sett_move <-
  new_settings(
    xlim = c(-640, 640), ylim = c(-400, 400),
    bounce_off_square = TRUE,
    bounce_off_circle = FALSE, circle_bounce_jitter = pi / 6,
    bounce_off_others = FALSE
  )

step_direct <- function(moment, time_next, settings) {
  # moment is position + direction + speed
  # just going in the same direction, possibly bouncing
  time_now <- moment$time[1]
  timestep <- time_next - time_now

  moment_next <-
    extrapolate_moment(moment, timestep, time_now, time_next)
  moment_next
}

make_random_trajectoryLJM <- function(start, times, settings, step_function, ...) {
  # make step
  stopifnot(length(times) > 1)

  num_objects <- max(start$object)
  start_time_zero <- start[1:num_objects,]
  if ("direction" %in% names(start)) {
    moment <- start_time_zero
  } else {
    moment <- start_time_zero %>% add_random_direction()
  }

  if (!"speed" %in% names(moment)) {
    moment <- moment %>% dplyr::mutate(speed = settings$speed)
  }
  moment <- moment %>% dplyr::mutate(time = times[1])

  # moment is position + direction + speed
  # here we have tibble of moments with two columns:
  #   - time (from times)
  #   - position (embedded moment tibble for each time)
  # Later we expand tibbles into trajectory tibble
  moment_tbl <- tibble::tibble(time = times, position = list(tibble::tibble))
  moment_tbl$position[[1]] <- moment #Initial position is specified by start, but not the subsequent position
  for (i in 2:length(times)) {

    #Check whether start actually specified speed/direction not just for the first time step but for the current time as well
    if (nrow(start) > i*num_objects) {
      #Set the speeds and directons to what is specified by start for this time
      moment$speed <- start[(i*num_objects+1):(i*num_objects+num_objects),]$speed
      moment$direction <- start[(i*num_objects+1):(i*num_objects+num_objects),]$direction
    }

    if (settings$bounce_off_square) {
      #error when this function is called
      moment <- bounce_off_square(moment, times[i], settings)
    }
    if (settings$bounce_off_circle) {
      moment <- bounce_off_circle(moment, times[i], settings)
    }
    if (settings$bounce_off_inside) {
      moment <- bounce_off_inside(moment, times[i], settings)
    }

    if (settings$bounce_off_others) {
      moment <- bounce_off_others(moment, times[i], settings)
        }

    moment_next <- step_function(moment, times[i], settings, ...)

    moment_tbl$position[[i]] <- moment_next
    moment <- moment_next
  }

  moment_tbl %>%
    dplyr::select(-.data$time) %>%
    tidyr::unnest(cols = c("position"))
}

#setting up path for trajectory files #LJM
traj_MOT_dir <- here("trajectories")
if(!dir.exists(traj_MOT_dir)) { dir.create(traj_MOT_dir, recursive = T) }

#pilot with 2 speeds for first and final speeds
first_speeds <- c(500,500,1000,1000)
final_speeds <- c(500,1000,500,1000)
possible_speeds <- tibble(first_speed = first_speeds, final_speed = final_speeds)
#duplicate again to make noise and no_noise trials
duplicated_possible_speeds <- rbind(possible_speeds,possible_speeds)
noise_present <- c(rep("dynamic", nrow(possible_speeds)),
                   rep("static", nrow(possible_speeds)))
noise_column <- tibble(noise_present = noise_present)
possible_speeds_with_noise <- cbind(duplicated_possible_speeds,noise_column)

#step function with a factor of 2 (like Ryo speed prereg)
stepfunction <- c(-8,-4,-2,-1,0,1,2,4,8)

#generating trajectory files LJM
for (step in 1:length(stepfunction)) {
  thisStep <- stepfunction[step]
  for (i in 1:nrow(possible_speeds_with_noise)) {
  #times of first trajectory
  first_portion_speed <- possible_speeds_with_noise[i,1]
  final_portion_speed <- possible_speeds_with_noise[i,2]
  distance_to_center <- abs(sett_generate$xlim[2])
  duration_of_final_portion <- 0.1

  duration_of_first_portion <- ((distance_to_center - duration_of_final_portion*final_portion_speed) / first_portion_speed) + 
    thisStep*timeGrain*(1000/first_portion_speed)
  
  times <- seq(from = 0,to = duration_of_first_portion,by=timeGrain)
  
  first_traj_length <- tail(times, n=1)
  
  #Make the trajectory for time zero up to the last critical moments when the speed is set specially
  position <- generate_positions_random(how_many_objects, sett_generate)
  
  position$direction = 0
  position$timeStep <- 0
  position$speed <- first_portion_speed
  
  trajectory_firstPortionOfTheTrial <-
    make_random_trajectoryLJM(position, times, sett_move,
                              step_direct)
  
  #Determining length of first trajectory (multipling by number of objects)
  firstPortionLength <- length(times)
  
  #copying position for the first of the last 3 portions where the speed for each is set by the factorial design
  #setting position for final_trajectory to trajectory_firstPortionOfTheTrial last how_many_objects rows position values
  position <- trajectory_firstPortionOfTheTrial[(firstPortionLength):(firstPortionLength+(how_many_objects-1)),]
  
  #defining specific times at which next trajectory takes place
  finalTimeOfFirstTrajectory <- tail(times, n=1)
  times <- seq(finalTimeOfFirstTrajectory, finalTimeOfFirstTrajectory+duration_of_final_portion, by = timeGrain)
  
  beginning_of_final_trajectory <- position
  beginning_of_final_trajectory$timeStep <- 1
  beginning_of_final_trajectory$speed <- final_portion_speed
  final_trajectory <-
    make_random_trajectoryLJM(beginning_of_final_trajectory,times,sett_move,step_direct)
  finaltlength <- nrow(final_trajectory)*timeGrain
  
  #clipping copied rows of each trajectory to avoid duplicating rows at each stitch
  trajectory_firstPortionOfTheTrial <- trajectory_firstPortionOfTheTrial[-c((firstPortionLength):(firstPortionLength+(how_many_objects-1))),]
  trajectory_wholeTrial <- rbind(trajectory_firstPortionOfTheTrial, final_trajectory)
  
 trajectory_wholeTrialreversed <- trajectory_wholeTrial
 trajectory_wholeTrialreversed$x <- trajectory_wholeTrialreversed$x*-1
 trajectory_wholeTrialreversed$y <- trajectory_wholeTrialreversed$y*-1
 trajectory_wholeTrialreversed$object <- trajectory_wholeTrialreversed$object+1
 
 for (thisRow in 1:nrow(trajectory_wholeTrial)){
   if (thisRow == 1) {trajectory_wholeTrialcombined <- NA}
   trajectory_wholeTrialcombined <- rbind(trajectory_wholeTrialcombined,trajectory_wholeTrial[thisRow,],
                                          trajectory_wholeTrialreversed[thisRow,])
   if (thisRow == nrow(trajectory_wholeTrial)) {
     trajectory_wholeTrialcombined <- trajectory_wholeTrialcombined[-1,]
   }
 }
 p<- plot_trajectory(trajectory_wholeTrialcombined, sett_move)
 show(p)
 position_of_step <- match(stepfunction[step],stepfunction) - 5
 save_trajectory(trajectory_wholeTrialcombined,filename = file.path(traj_MOT_dir,
                                                                    sprintf("T%00d%00dx%00dx%s.csv",
                                                                            position_of_step,possible_speeds_with_noise[i,1],
                                                                            possible_speeds_with_noise[i,2],possible_speeds_with_noise[i,3])))
 trajectory_wholeTrialflipped <- trajectory_wholeTrialcombined
 trajectory_wholeTrialflipped$y <- trajectory_wholeTrialflipped$y*-1
 p<- plot_trajectory(trajectory_wholeTrialflipped, sett_move)
 show(p)
 save_trajectory(trajectory_wholeTrialflipped,filename = file.path(traj_MOT_dir,
                                                                   sprintf("T%00d%00dx%00dx%sf.csv",
                                                                           position_of_step,possible_speeds_with_noise[i,1],
                                                                           possible_speeds_with_noise[i,2],possible_speeds_with_noise[i,3])))
}
}
