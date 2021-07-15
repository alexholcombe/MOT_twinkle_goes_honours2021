# this scripts generate trajectories for the experiment
rm(list = ls())
set.seed(1902241)

library(tidyverse)

library(here)
detach("package:motrack", unload=TRUE)
library(motrack)

source(here("R", "utils.R"))

stim_detection_dir    <- here("stimuli_detection")
noise_detection_dir   <- file.path(stim_detection_dir, "noise")
gabor_detection_dir   <- file.path(stim_detection_dir, "gabor")

prot_dir              <- here("data","protocols")
prot_detection_dir    <- file.path(prot_dir, "detection")
prot_MOT_dir          <- file.path(prot_dir, "MOT")

if(!dir.exists(prot_detection_dir)) { dir.create(prot_detection_dir, recursive = T) }
if(!dir.exists(prot_MOT_dir)) { dir.create(prot_MOT_dir, recursive = T) } 
if(!dir.exists(noise_detection_dir)) { dir.create(noise_detection_dir, recursive = T) }
if(!dir.exists(gabor_detection_dir)) { dir.create(gabor_detection_dir, recursive = T) }

stim_MOT_dir <- here("stimuli_MOT")
noise_MOT_dir   <- file.path(stim_MOT_dir, "noise")
if(!dir.exists(noise_MOT_dir)) { dir.create(noise_MOT_dir, recursive = T) }


# generate noise for detection
set.seed(1902242)

n_noises <- 500

tm <- FDhelpers::create.time.measure(n_noises)

for(i in 1:n_noises) {
  noise <- pink_noise_2d()
  noise <- rescale_noise(noise, 0.1)
  noise[1,1] <- 255
  noise[1,2] <- 0
  noise2 <- round(noise)
  png::writePNG(noise2, file.path(noise_detection_dir, sprintf("noise%03d.png",i)))
  tm <- FDhelpers::update.tm(tm)
  FDhelpers::print.tm(tm)
}

# generate noise for MOT
set.seed(1902243)

n_noises <- 500

tm <- FDhelpers::create.time.measure(n_noises)

for(i in 1:n_noises) {
  noise <- pink_noise_2d()
  noise <- rescale_noise(noise, 0.1)
  noise[1,1] <- 255
  noise[1,2] <- 0
  noise2 <- round(noise)
  png::writePNG(noise2, file.path(noise_MOT_dir, sprintf("noise%03d.png",i)))
  
  tm <- FDhelpers::update.tm(tm)
  FDhelpers::print.tm(tm)
}

# Generate protocols for detection
set.seed(1902244)

nProt <- 20
nTrials_per_level <- 50
# five contrast levels
t_contr_id  <- 1:5
t_contr_val <- c(0.06, 0.075, 0.1, 0.13, 0.16)
contr_values <- tibble(t_contr_id, t_contr_val)

p0_det <- create_protocol_detection(0, t_contr_id, nTrials_per_level)


for (i in 1:nProt) {
  p <- p0_det %>% sample_frac(1) %>% mutate(prot_id = i)
  p <- p %>% mutate(trial_id = 1:n())
  write_excel_csv(p, file.path(prot_detection_dir, sprintf("P%03d.csv", i)))
}

## MOT part
set.seed(1902245)

nProt <- 20
n_levels <- 4
t_contr_id_mot  <- 1:n_levels
nTrials_per_level_mot <- 40

for (i in 1:nProt) {
  p <- create_protocol_mot(i, t_contr_id_mot, nTrials_per_level_mot)
  write_csv(p, file.path(prot_MOT_dir, sprintf("P%03d.csv", i)))
  
}

## MOT part - practice
set.seed(19022451)

n_levels <- 4
t_contr_id_mot  <- 1:n_levels
nTrials_per_level_mot <- 40

p_practice <- create_protocol_mot(0, t_contr_id_mot, nTrials_per_level_mot)
p_practice <- p_practice %>% group_by(t_contr) %>% sample_n(4) %>% ungroup() %>% arrange(-t_contr)
p_practice <- p_practice %>% mutate(trial_id = 1:n()) 
write_csv(p_practice, file.path(prot_MOT_dir, "P_practice.csv"))

# generate gabors
set.seed(1902246)

gabor <- gabor_patch()
for (i in 1:n_levels) {
  gaborx <- 128 + round(gabor * t_contr_val[i] * 128)  
  png::writePNG(gaborx, file.path(gabor_detection_dir, sprintf("gabor_contr_id_%d_contr_val%.3f.png", i, t_contr_val[i])))
}

# generate images for instructions
set.seed(1902247)

gaborx <- round(gabor * t_contr_val[5] * 128)  
im2 <- blend_images(noise2, gaborx)
png::writePNG(noise2, file.path(gabor_detection_dir, "stim_detection_instruction_noise.png"))
png::writePNG(im2, file.path(gabor_detection_dir, "stim_detection_instruction_gabor_easiest.png"))
