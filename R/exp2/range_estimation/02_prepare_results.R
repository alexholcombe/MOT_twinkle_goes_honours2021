library(tidyverse)
theme_set(theme_classic(16))

library(osfr)
library(here)

source(here("R","utils.R"))

local_data_pth <- file.path(here("data"), "exp2", "detection","range_estim","results")
local_data_pth_participants <- file.path(here("data"), "exp2", "detection","range_estim")

should_download <- T

if(should_download) {
  osf_auth(token = read_lines(list.files("data/", pattern = "osf_token_write_*",full.names = T)))
  
  data_guid <- "9ht2f"
  
  noisemot_project <- osf_retrieve_node(data_guid)
  
  if(!dir.exists(local_data_pth)) {
    dir.create(local_data_pth,recursive = T) 
  }
  
  if(!dir.exists(local_data_pth_participants)) {
    dir.create(local_data_pth_participants,recursive = T) 
  }
  
  
  data_files <- 
    noisemot_project %>% 
    osf_ls_files() %>% 
    filter(name == "exp2") %>% 
    osf_ls_files() %>% 
    filter(name == "detection") %>% 
    osf_ls_files() %>%
    filter(name == "range_estim") %>%
    osf_ls_files() %>%
    filter(name == "results") %>% 
    osf_ls_files(n_max = Inf) 
    
    
  
  data_files %>% 
    do(download_files(.,local_data_pth))
  
  # download participant info
  
  data_files_participants <- 
    noisemot_project %>% 
    osf_ls_files() %>% 
    filter(name == "exp2") %>% 
    osf_ls_files() %>% 
    filter(name == "detection") %>% 
    osf_ls_files() %>%
    filter(name == "range_estim") %>%
    osf_ls_files() %>% 
    filter(name != "results")
  
  data_files_participants %>% 
    do(download_files(.,local_data_pth_participants))
  
}

files <- dir(local_data_pth, pattern = "*.csv") # get file names

df <- files %>%
  map(~ read_csv(file.path(local_data_pth, .), col_types = cols(
    .default = col_double(),
    where_gabor = col_character(),
    b_key_resp_2_.keys = col_character(),
    noise_pth = col_character(),
    expName = col_character(),
    date = col_character(),
    gabor_template_pth = col_character(),
    X22 = col_logical()
  ))) %>% 
  reduce(rbind) %>% as_tibble()

df <- df %>% 
  mutate(trial_id = 1:n()) %>% 
  select(subject_id = participant, 
        trial_id,
         noise1_id, 
         noise2_id,
         contr_level,
         where_gabor,
         response = b_key_resp_2_.keys,
         rt = b_key_resp_2_.rt,
         correct = b_key_resp_2_.corr) 


write_rds(df, here("data","exp2","detection","range_estim","data_200219.rds"))
