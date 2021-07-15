# read data
library(tidyverse)
library(purrr)
library(here)

outpth <- here("data","cleandata")
if(!dir.exists(outpth)) {dir.create(outpth,recursive = T)}

data_path_with_mouse <- here::here("data","results","with_allmouse_coords")   # path to the data
files <- dir(data_path_with_mouse, pattern = "*.csv") # get file names

df_withmouse <- files %>%
  map(~ read_csv(file.path(data_path_with_mouse, .)) %>% 
        select(-(practice_trials.thisRepN:practice_trials.thisIndex)) %>% 
        select(-(trials.thisRepN:trials.thisIndex)) %>% 
        select(-(mouse.leftButton:mouse.rightButton)) %>% 
        select(-date,-frameRate,-X30) %>% 
        filter(prot_id != 0)) %>% 
  reduce(rbind)
df_withmouse


data_path_without_mouse <- here::here("data","results","without_allmouse_coords")   # path to the data
files <- dir(data_path_without_mouse, pattern = "*.csv") # get file names

df_withoutmouse <- files %>%
  map(~ read_csv(file.path(data_path_without_mouse, .))  %>% 
        select(-(practice_trials.thisRepN:practice_trials.thisIndex)) %>% 
        select(-(trials.thisRepN:trials.thisIndex)) %>% 
        select(-(mouse.leftButton:mouse.rightButton)) %>% 
        select(-date,-frameRate,-X27) %>% 
        filter(prot_id != 0)) %>% 
  reduce(rbind)
df_withoutmouse


df_all <- rbind(df_withmouse %>% select(-mouse_xall,-mouse_yall,-mouse_tall), df_withoutmouse)

saveRDS(df_all, file.path(outpth, "df_all_190328.rds"))

saveRDS(df_withmouse, file.path(outpth, "df_withmouse_190328.rds"))

data_trajectories <- here::here("data","trajectories")  
files <- dir(data_trajectories, pattern = "*.csv") # get file names
df_trajectories <-  files %>%
  map(~ read_csv(file.path(data_trajectories, .),col_names =F)) %>% 
  reduce(rbind)

df_trajectories$trajectory_id <- rep(1:40, each = 81)

colnames(df_trajectories)[1:17] <- c("t",
  expand.grid(c("X","y"),1:8) %>% tidyr::unite(cn, Var1, Var2,sep="") %>% pull(cn))

df_trajectories <- df_trajectories %>% select(trajectory_id,everything())

saveRDS(df_trajectories, file.path(outpth, "trajectories_190328.rds"))
