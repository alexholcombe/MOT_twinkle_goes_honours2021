# trajectory statistics on repeated presentations

library(tidyverse)
theme_set(theme_classic(16))
library(here)

source(here("R","utils.R"))

data_pth <- here("data","exp1")
df_exp1 <- readRDS(file.path(data_pth, "data_191129.rds"))

df_trajectories <- readRDS(file.path(data_pth, "trajectories_200406.rds")) %>% filter(trajectory_id %in% df_exp1$trajectory_id)
#df_trajends <- df_withmouse %>% select(trajectory_id,start_time,) %>% mutate(end_time = start_time + 6)
#df_trajectories %>% left_join(df_trajends) %>% filter(t == end_time)


df_traj1 <- df_trajectories %>% reshape_trajectories()
df_traj2a <- df_traj1 %>% group_by(trajectory_id) %>% do(describe_trajectory(.))
df_traj2b <- df_traj1 %>% group_by(trajectory_id) %>% do(describe_trajectory_tgtdist(.))

df_traj_stata <- df_exp1 %>% group_by(trajectory_id,t_contr) %>% summarize(nCorrect = mean(nCorrect, na.rm = T)) %>% left_join(df_traj2a, by = "trajectory_id") %>% 
  mutate(target = if_else(target, "target", "distractor"))

df_traj_statb <- df_exp1 %>% group_by(trajectory_id,t_contr) %>% summarize(nCorrect = mean(nCorrect, na.rm = T)) %>% left_join(df_traj2b, by = "trajectory_id")

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

df_complete_withtraj <- df_exp1 %>% left_join(df_traj_stat %>% select(-nCorrect) %>% filter(target == "target"))

glm1 <- glm(acc ~ t_contr+dist_from_tgt_center_mean+min_tgtdist+ mean_ch_inter, df_complete_withtraj, family = "quasibinomial") 
exp(coef(glm1))

glm1 %>% summary()
# values seems ok, minimum distance between target and distractors is the best predictor (in the terms of log odds)




