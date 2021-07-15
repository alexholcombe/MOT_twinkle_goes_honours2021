
# analyze FD data on just staircase ---------------------------------------

df <- read_csv2("psychopy/exp2/detection/data/999_noiseMot_exp2_detection_2020_Feb_11_0933.csv") %>% 
  mutate(trials.intensity=as.numeric(trials.intensity))

df %>% 
  mutate(contr=round(trials.intensity,2)) %>% 
  ggplot(aes(x=trials.thisTrialN,y=contr,group=1)) + 
  geom_point() + 
  geom_path()

df %>% pull(trials.intensity)%>% tail(20) %>% mean()

