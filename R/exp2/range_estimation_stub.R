df <- read_csv("data/exp2/discrimination/range_estim/999_noiseMot_exp2_detection_range_estim_2020_Feb_11_1340.csv") 

df %>% 
  group_by(contr_level) %>% 
  mutate(correct = key_resp_2.corr) %>% 
  ggplot(aes(x = contr_level, y = correct)) + 
  stat_summary(fun.data = "mean_cl_boot") + 
  scale_x_log10()


fit_par <- df %>% 
  mutate(t_contr = contr_level,
         correct = key_resp_2.corr) %>% 
  mle_SDT()

plot_fit(fit_par,df %>% 
           mutate(t_contr = contr_level,
                  correct = key_resp_2.corr))
