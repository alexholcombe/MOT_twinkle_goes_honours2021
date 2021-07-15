# this scripts generate trajectories for the experiment
rm(list = ls())
set.seed(191104001)

# preparations ------------------------------------------------------------

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R", "utils.R"))

stim_detection_dir    <- here("stimuli","exp2","detection")
noise_detection_dir   <- file.path(stim_detection_dir, "noise")
gabor_detection_dir   <- file.path(stim_detection_dir, "gabor")

prot_dir              <- here("data","protocols","exp2")
prot_detection_dir    <- file.path(prot_dir, "detection")
prot_MOT_dir          <- file.path(prot_dir, "MOT")

if(!dir.exists(prot_detection_dir)) { dir.create(prot_detection_dir, recursive = T) }
if(!dir.exists(prot_MOT_dir)) { dir.create(prot_MOT_dir, recursive = T) } 
if(!dir.exists(noise_detection_dir)) { dir.create(noise_detection_dir, recursive = T) }
if(!dir.exists(gabor_detection_dir)) { dir.create(gabor_detection_dir, recursive = T) }

stim_MOT_dir <- here("stimuli","exp2","MOT")
noise_MOT_dir   <- file.path(stim_MOT_dir, "noise")
if(!dir.exists(noise_MOT_dir)) { dir.create(noise_MOT_dir, recursive = T) }


# detection ---------------------------------------------------------------

## generate noise for detection and MOT

generate_noise_exp2(seed_id = 2020204, n_noises = 500, out_dir = noise_detect_dir)

generate_noise_exp2(seed_id = 2020404, n_noises = 500, out_dir = noise_MOT_dir)

## generate gabors for detection

generate_gabor_exp2(out_dir = gabor_detection_dir, file_name = "gabor_detection_template.png")
gb <- gabor_patch()
write.table(gb,file.path(gabor_detection_dir,"gabor_template.csv"),row.names = F,col.names=F)

## create protocols for range estimation exp

p <- tibble(contr_level = rep(c(0.040,0.065,0.090,0.115,0.140,0.165,0.190,0.215),5))
write_csv(p,path = file.path(prot_detection_dir,"design_range_estim.csv"))
# MOT ---------------------------------------------------------------------

set.seed(191104005)

nProt <- 36
n_levels <- 3
t_contr_id_mot  <- 1:n_levels
nTrials_per_level_mot <- 30

t_contr_all <- perm(t_contr_id_mot)

mark_nomark_raw <- 
  expand.grid(c("mark","noMark"), c("mark","noMark"),c("mark","noMark")) %>% as.matrix()
mark_nomark_all <- matrix(NA, ncol = 6, nrow = 8)

for (i in 1:nrow(mark_nomark_raw)) {
  mnm <- c(mark_nomark_raw[i,])
  mark_nomark_raw_complement <- if_else(mnm == "mark", "noMark","mark")
  mark_nomark_all[i,c(1,3,5)] <- mnm
  mark_nomark_all[i,c(2,4,6)] <- mark_nomark_raw_complement
}
for (i in 1:nProt) {
  
  p <- create_protocol_mot_exp1(i, t_contr_all[((i-1) %% 6)+1, ], mark_nomark_all[((i-1) %% 8)+1, ],nTrials_per_level_mot)
  p_noiseMot <- p %>% filter(t_contr != 0)
  p_normalMot <- p %>% filter(t_contr == 0) %>% mutate(trial_id = 1:n()) %>% mutate(t_contr = "MOT") 
  write_csv(p_noiseMot, file.path(prot_MOT_dir, sprintf("P%03d_noise.csv", i)))
  write_csv(p_normalMot, file.path(prot_MOT_dir, sprintf("P%03d_normal.csv", i)))
  
}

## MOT part - practice
set.seed(191104006)

n_levels <- 3
t_contr_id_mot  <- 1:n_levels
nTrials_per_level_mot <- 30

mark_nomark_raw <- 
  expand.grid(c("mark","noMark"), c("mark","noMark"),c("mark","noMark")) %>% as.matrix()
mark_nomark_all <- matrix(NA, ncol = 6, nrow = 8)

for (i in 1:nrow(mark_nomark_raw)) {
  mnm <- c(mark_nomark_raw[i,])
  mark_nomark_raw_complement <- if_else(mnm == "mark", "noMark","mark")
  mark_nomark_all[i,c(1,3,5)] <- mnm
  mark_nomark_all[i,c(2,4,6)] <- mark_nomark_raw_complement
}

p_practice <- create_protocol_mot_exp1(0, t_contr_all[1,], mark_nomark_all[1,], nTrials_per_level_mot)
p_practice_noiseMOT <- p_practice %>% 
  filter(t_contr!=0) %>% 
  group_by(t_contr,mark_type) %>% 
  sample_n(2) %>% 
  ungroup() %>% 
  arrange(-t_contr) %>% 
  mutate(trial_id = 1:n()) 

p_practice_normalMOT <- p_practice %>% 
  filter(t_contr==0) %>% 
  group_by(t_contr,mark_type) %>% 
  sample_n(3) %>% 
  ungroup() %>% 
  arrange(-t_contr) %>% 
  mutate(t_contr = "MOT") %>% 
  mutate(trial_id = 1:n()) 

write_csv(p_practice_noiseMOT,  file.path(prot_MOT_dir, "P_practice_noise.csv"))
write_csv(p_practice_normalMOT, file.path(prot_MOT_dir, "P_practice_normal.csv"))

# generate gabors
set.seed(191104007)

# values found by script 02_fitdata_pilot
c1 <- 0.11679
c2 <- 0.14049
c3 <- 0.16675
t_contr_val <- c(c1,c2,c3)

gabor <- gabor_patch()
for (i in 1:n_levels) {
  gaborx <- 128 + round(gabor * t_contr_val[i] * 128)  
  png::writePNG(gaborx, file.path(gabors_MOT_dir, sprintf("gabor_contr_id_%d_contr_val%.3f.png", i, t_contr_val[i])))
}

# copy to psychopy folder

psychopy_dir    <- here("psychopy","exp2")
psychopy_gabors <- file.path(psychopy_dir, "gabors")
psychopy_noise  <- file.path(psychopy_dir, "noise")
psychopy_protocols  <- file.path(psychopy_dir, "protocols")

if(!dir.exists(psychopy_gabors)) { dir.create(psychopy_gabors, recursive = T) }
if(!dir.exists(psychopy_noise)) { dir.create(psychopy_noise, recursive = T) }
if(!dir.exists(psychopy_protocols)) { dir.create(psychopy_protocols, recursive = T) }

# copy gabors
gabors_files_from <- list.files(gabors_MOT_dir,full.names = T)
gabors_files_to   <- 
  gabors_files_from %>% 
  basename() %>% 
  file.path(psychopy_gabors,.) %>% 
  str_remove("_contr_val[-+]?[0-9]*\\.?[0-9]+")

file.copy(gabors_files_from, gabors_files_to,overwrite = T)


# copy noise
noise_files_from <- list.files(noise_MOT_dir,full.names = T)
noise_files_to   <- 
  noise_files_from %>% 
  basename() %>% 
  file.path(psychopy_noise,.)
file.copy(noise_files_from, noise_files_to,overwrite = T)

# copy protocols
protocols_files_from <- list.files(prot_MOT_dir,full.names = T)
protocols_files_to   <- 
  protocols_files_from %>% 
  basename() %>% 
  file.path(psychopy_protocols,.)
file.copy(protocols_files_from, protocols_files_to,overwrite = T)

