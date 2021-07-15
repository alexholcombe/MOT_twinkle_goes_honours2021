rm(list = ls())

df <- read_csv("data/results/99_demo_2019_Mar_27_1702_1.csv")

answ <- df$mouse.clicked_name %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_remove_all("o1_copy_") %>% 
  str_split(",", simplify = T)

df$nCorrect <- matrix(data = answ %in% c("0","1","2","3"), ncol = 4) %>% rowSums()

df %>% group_by(t_contr) %>% summarize(nCorrect = mean(nCorrect))

df %>% ggplot(aes(x = t_contr, y = nCorrect)) + stat_summary(fun.data = "mean_cl_boot") + ylim(0,4)


df2 <- read_csv("psychopy/data/999_demo_2019_Mar_27_1653.csv")

x <- df2$mouse_xall %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_split(",", simplify = T) %>% as.numeric()

y <- df2$mouse_yall %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_split(",", simplify = T) %>% as.numeric()

plot(x,y)
