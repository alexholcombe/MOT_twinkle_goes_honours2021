#create a tibble which loops through the 2(target no) x 2(noise vs. no noise) combinations as many times as necessary (e.g. 30 times for 180 trials)
#have trajectories 1:180 shuffled randomly
setwd("~/Documents/OneDrive - The University of Sydney (Students)/Empirical Thesis/Code/MOT_twinkle_goes_honours2021-main")

library(tidyverse)
library(here)

set.seed(123)

prot_dir <- here("exp2_JJC", "protocols_2cages_4objects")

targets = c(1, 2)

noise = c("noise", "no_noise")

line = c("vertical", "horizontal")

n_trials_per_block = 100

n_test_trials = 1000

n_blocks = n_test_trials/n_trials_per_block

n_practice_trials = 8

reps <- n_test_trials/(length(targets)*length(noise)*length(line))

reps_p <- n_practice_trials/(length(targets)*length(noise)*length(line))

#create practice protocols
combinations <- expand.grid(noise_present = noise, no_targets = targets, rep = seq(1, reps_p), line_orientation = line) #8 practice trials

trajectory_id <- c(1:(n_practice_trials/2), (n_practice_trials/2 + n_test_trials/2 + 1):(n_test_trials/2 + n_practice_trials))
combinations$trajectory_id <- trajectory_id

combinations <- slice(combinations, sample(1:n()))

combinations <- sample(combinations)

prot_id <- rep(0, n_practice_trials)
combinations$prot_id <- prot_id

trial_id <- seq(1, n_practice_trials)
combinations$trial_id <- trial_id

combinations$rep <- NULL

combinations <- combinations %>% relocate(trial_id)
combinations <- combinations %>% relocate(prot_id)

write_csv(combinations, file.path(prot_dir, sprintf("P_practice_noise.csv")))

#create test protocols
for (i in 1:9) {
  
  combinations <- expand.grid(noise_present = noise, no_targets = targets, rep = seq(1, reps), line_orientation = line)
  
  trajectory_id <- c((n_practice_trials/2 + 1):(n_practice_trials/2 + n_test_trials/2),(n_practice_trials + n_test_trials/2 + 1):(n_practice_trials + n_test_trials))
  
  combinations$trajectory_id <- trajectory_id
  
  combinations <- slice(combinations, sample(1:n()))
  
  combinations <- sample(combinations)
  
  prot_id <- rep(i, n_test_trials)
  combinations$prot_id <- prot_id
  
  #print("we got here!")
  
  trial_id <- seq(1, n_test_trials)
  combinations$trial_id <- trial_id
  
  #print("we got here!")

  combinations$rep <- NULL
  
  combinations <- combinations %>% relocate(trial_id)
  combinations <- combinations %>% relocate(prot_id)
  
  for (j in 1:n_blocks) {
    start_index = j*n_trials_per_block - (n_trials_per_block - 1)
    stop_index = j*n_trials_per_block
    combinations_j <- combinations[start_index: stop_index, ]
    prot_subset <- rep(j, n_trials_per_block)
    combinations_j$prot_subset <- prot_subset
    combinations_j <- combinations_j %>% relocate(prot_subset)
    combinations_j <- combinations_j %>% relocate(prot_id)
    write_csv(combinations_j, file.path(prot_dir, sprintf("P%d%02d_noise.csv", i, j)))
  }
}