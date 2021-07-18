#create a tibble which loops through the 3(target no) x 2(noise vs. no noise) combinations as many times as necessary (e.g. 30 times for 180 trials)
#have trajectories 1:180 shuffled randomly

library(tidyverse)
library(here)

set.seed(123)

prot_dir <- here("exp1", "protocols_n")

targets = c(1, 2, 4)

noise = c("noise", "no_noise")

n_test_trials = 180

n_practice_trials = 6

n_trials = n_test_trials + n_practice_trials

reps <- n_test_trials/(length(targets)*length(noise))

reps_p <- n_practice_trials/(length(targets)*length(noise))

#create practice protocols
combinations <- expand.grid(noise_present = noise, no_targets = targets, rep = seq(1, reps_p)) #6 practice trials

combinations <- slice(combinations, sample(1:n()))

combinations <- sample(combinations)

prot_id <- rep(0, n_practice_trials)
combinations$prot_id <- prot_id

trial_id <- seq(1, n_practice_trials)
combinations$trial_id <- trial_id

trajectory_id <- seq(1, n_practice_trials)
print(trajectory_id)
#trajectory_id <- sample(trajectory_id)
#trajectory_id <- trajectory_id[c(1:6)]
combinations$trajectory_id <- trajectory_id

# cages <- c(1, 2, 3, 4, 1, 2)
# combinations$cage_indicated <- cages

combinations$rep <- NULL

combinations <- combinations %>% relocate(trial_id)
combinations <- combinations %>% relocate(prot_id)

write_csv(combinations, file.path(prot_dir, sprintf("P_practice_noise.csv")))

#create test protocols
for (i in 1:36) {
  
  combinations <- expand.grid(noise_present = noise, no_targets = targets, rep = seq(1, reps))
  
  combinations <- slice(combinations, sample(1:n()))
  
  combinations <- sample(combinations)
  
  prot_id <- rep(i, n_test_trials)
  combinations$prot_id <- prot_id
  
  #print("we got here!")
  
  
  trial_id <- seq(1, n_test_trials)
  combinations$trial_id <- trial_id
  
  #print("we got here!")
  
  trajectory_id <- seq(n_practice_trials + 1, n_trials)
  trajectory_id <- sample(trajectory_id)
  combinations$trajectory_id <- trajectory_id
  
  combinations$rep <- NULL
  
  combinations <- combinations %>% relocate(trial_id)
  combinations <- combinations %>% relocate(prot_id)
  
  write_csv(combinations, file.path(prot_dir, sprintf("P%03d_noise.csv", i)))
  
}