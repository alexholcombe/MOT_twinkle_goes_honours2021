# analyze data and generate gabors
library(here)
library(tidyverse)

fpth <- ""

prot_id <- 1

subj_id <-1

diag_plot_dir <- here("data","plots","contrast_fits")
if(!dir.exists(diag_plot_dir)) {
  dir.create(diag_plot_dir, recursive = T)
}

gabor_dir <- here("psychopy","plots","contrast_fits")
if(!dir.exists(diag_plot_dir)) {
  dir.create(diag_plot_dir, recursive = T)
}

prot_dir <- here("data","protocols")
if(!dir.exists(prot_dir)) {
 stop("Protocols directory not found. Incorrect path?")
}

out_gabor_dir <- here("data","psychopy","gabors",sprintf("subj%02d",subj_id))
if(!dir.exists(out_gabor_dir)) {
  dir.create(out_gabor_dir, recursive = T)
}

px <- read_csv2(file.path(prot_dir, sprintf("P%03d.csv", prot_id)), 
                col_types = cols(
  prot_id = col_double(),
  trial_id = col_double(),
  trajectory_id = col_double(),
  t_contr = col_double(),
  start_time = col_double()
))

df <- read_csv(here("data","trial_example.csv")) %>% 
  rename(t_contr = Var1, response = Var2, tgt_present = Var3) %>% 
  mutate(correct = if_else((tgt_present == 1 & response == "r")|(tgt_present == 0 & response == "l"),1,0))


fit_par <- mle_SDT(df)

p <- plot_fit(fit_par,df)
ggsave(file.path(diag_plot_dir,sprintf("subj%02d_cT_%.2f_beta_%.2f.png", subj_id, fit_par[1], fit_par[2])), plot = p)


stopifnot(inv_p_corr_fun(p_corr_fun(0.05, fit_par[1], fit_par[2]), 
                         fit_par[1], fit_par[2]) - 0.05 < 0.0001)
#detection_accuracies <- c(0.7,0.8,0.9,0.999)
dprimes <- c(4,6,8,10)
contrasts_for_mot <- dprime_fun(dprimes, fit_par[1], fit_par[2])

px_calibrated <- change_contrasts_protocols(px, contrasts_for_mot)

write_csv2(px_calibrated, sprintf("P%03d_subj%02d.csv", prot_id, subj_id))

# generate gabor with default paramteres

gabor <- gabor_patch()

for (i in 1:length(contrasts_for_mot)) {
  gaborx <- 128 + round(gabor * contrasts_for_mot[i] * 128)
  png::writePNG(gaborx, file.path(gabor_dir, sprintf("gabor_subj%02d_contrid_%d.png", subj_id, i)))
   
}