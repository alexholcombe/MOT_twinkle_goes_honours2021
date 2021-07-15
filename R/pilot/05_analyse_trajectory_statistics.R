# trajectory statistics on repeated presentations

library(tidyverse)
theme_set(theme_classic(16))
library(here)

source(here("R","utils.R"))

data_pth <- here("data","cleandata")
df_all <- readRDS(file.path(data_pth, "df_all_190328.rds"))

answ <- df_all$mouse.clicked_name %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_remove_all("o1_copy_") %>% 
  str_split(",", simplify = T)
corr_trials <- answ[,5]==""
df_all$nCorrect <- matrix(data = answ %in% c("0","1","2","3"), ncol = 5) %>% rowSums()
df_all$nCorrect[!corr_trials] <- NA
df_withmouse <- readRDS(file.path(data_pth, "df_withmouse_190328.rds"))

answ <- df_withmouse$mouse.clicked_name %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_remove_all("o1_copy_") %>% 
  str_split(",", simplify = T)
corr_trials <- answ[,5]==""

df_withmouse$nCorrect <- matrix(data = answ %in% c("0","1","2","3"), ncol = 5) %>% rowSums()
df_withmouse$nCorrect[!corr_trials] <- NA

df_trajectories <- readRDS(file.path(data_pth, "trajectories_190328.rds"))
#df_trajends <- df_withmouse %>% select(trajectory_id,start_time,) %>% mutate(end_time = start_time + 6)
#df_trajectories %>% left_join(df_trajends) %>% filter(t == end_time)

mouse_xall <- df_withmouse$mouse_xall %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_split(",", simplify = T)

mouse_xall <- matrix(as.numeric(mouse_xall), ncol = ncol(mouse_xall))

mouse_yall <- df_withmouse$mouse_yall %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_split(",", simplify = T)

mouse_yall <- matrix(as.numeric(mouse_yall), ncol = ncol(mouse_yall))

nclicks <- rowSums(!is.na(mouse_xall))
df_withmouse$nClicks <- nclicks

df_complete <- df_all %>% filter(participant > 15)


df_traj1 <- df_trajectories %>% reshape_trajectories()
df_traj2a <- df_traj1 %>% group_by(trajectory_id) %>% do(describe_trajectory(.))
df_traj2b <- df_traj1 %>% group_by(trajectory_id) %>% do(describe_trajectory_tgtdist(.))

df_traj_stata <- df_all %>% group_by(trajectory_id,t_contr) %>% summarize(nCorrect = mean(nCorrect, na.rm = T)) %>% left_join(df_traj2a, by = "trajectory_id") %>% 
  mutate(target = if_else(target, "target", "distractor"))

df_traj_statb <- df_all %>% group_by(trajectory_id,t_contr) %>% summarize(nCorrect = mean(nCorrect, na.rm = T)) %>% left_join(df_traj2b, by = "trajectory_id")

df_traj_stata %>% ggplot(aes(x = chull_sd, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(target~t_contr)+ ggtitle("SD of convex hull for targets and distractor",subtitle = "per each contrast level")
df_traj_stata %>% ggplot(aes(x = chull_mean, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(target~t_contr)+ ggtitle("Mean of convex hull for targets and distractor",subtitle = "per each contrast level")
df_traj_stata %>% ggplot(aes(x = dist_from_tgt_center_mean, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(target~t_contr)+ ggtitle("Minimum distance between target and distractor",subtitle = "per each contrast level")
df_traj_stata %>% ggplot(aes(x = dist_from_center, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(target~t_contr)+ ggtitle("Minimum distance between target and distractor",subtitle = "per each contrast level")

df_traj_statb %>% ggplot(aes(x = min_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(~t_contr) + ggtitle("Minimum distance between target and distractor",subtitle = "per each contrast level")
df_traj_statb %>% ggplot(aes(x = mean_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(~t_contr) + ggtitle("Average distance between target and distractor",subtitle = "per each contrast level")
df_traj_statb %>% ggplot(aes(x = lower_quar_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(~t_contr)+ ggtitle("Lower quartile of distances between target and distractor",subtitle = "per each contrast level")

df_traj_statb %>% ggplot(aes(x = mean_ch_inter, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ facet_grid(~t_contr)+ ggtitle("Area of intersection between convex hulls",subtitle = "per each contrast level")


df_traj_statb %>% group_by(trajectory_id) %>% summarise_all(mean) %>% ggplot(aes(x = min_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Minimum distance between target and distractor")
df_traj_statb %>% group_by(trajectory_id) %>% summarise_all(mean) %>% ggplot(aes(x = mean_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Average distance between target and distractor")
df_traj_statb %>% group_by(trajectory_id) %>% summarise_all(mean) %>% ggplot(aes(x = lower_quar_tgtdist, y = nCorrect)) + geom_point() + geom_smooth(method = "lm")+ ggtitle("Lower quartile of distances between target and distractor")

df_traj_stat <- df_traj_stata %>% left_join(df_traj_statb)

cor(df_traj_stat %>% filter(target == "target") %>% select(dist_from_center:med_ch_inter) %>% as.matrix()) %>% round(2)
# distance from target mean and area of convex hull are highly correlated

lm(nCorrect ~ t_contr+dist_from_tgt_center_mean+min_tgtdist + mean_ch_inter, df_traj_stat %>% filter(target == "target")) %>% summary()

# linear model is probably not correct as it does not have limited range, we need something like logistic regression
# problem with logistic is that is expects 0/1 values, but we have counts. One option is to normalize it (by dividing it by 4) and ue quasibinomial family

df_complete_withtraj <- df_complete %>% left_join(df_traj_stat %>% select(-nCorrect) %>% filter(target == "target"))

glm1 <- glm(nCorrect ~ t_contr+dist_from_tgt_center_mean+min_tgtdist+ mean_ch_inter, df_complete_withtraj %>% mutate(nCorrect = nCorrect / 4), family = "quasibinomial") 
exp(coef(glm1))

# values seems ok, minimum distance between target and distractors is the best predictor (in the terms of log odds)




