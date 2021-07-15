library(tidyverse)
theme_set(theme_classic(16))

library(here)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp1"), "data")
#local_data_pth <- file.path(here("..","exp1", "data"))
dir(local_data_pth)

outpth <- here("exp1")

files <- dir(local_data_pth, pattern = "*.csv") # get file names
files <- tibble(fname=files)

outputFname = "data.rds"

testWithTwoFiles<-T #hand-pick two files for testing purposes
if (testWithTwoFiles) {
  outputFname = "data_testAOH.rds"
}

if (testWithTwoFiles) {
  files <- files %>% filter(fname=="999_noiseMot_exp1_noise_2021_May_10_1229.csv" | 
                          fname=="999_noiseMot_exp1_noise_2021_May_10_1219.csv")
}
fnames <- as.array(files$fname)
msg<- paste0("Number of files = ",as.character(nrow(fnames)), ".")
print(msg)

# #Try reading an individual file
# testfile <- "999_noiseMot_exp1_noise_2021_May_10_1219.csv" # "999_noiseMot_exp1_noise_2021_Apr_26_1154.csv"
# df<- read_csv(testfile,
#          col_types = cols(.default = col_double(),
#                             t_contr = col_character(),
#                             mark_type = col_character(),
#                             mouse.x = col_character(),
#                             mouse.y = col_character(),
#                             mouse.leftButton = col_character(),
#                             mouse.midButton = col_character(),
#                             mouse.rightButton = col_character(),
#                             mouse.time = col_character(),
#                             mouse.clicked_name = col_character(),
#                             date = col_character(),
#                             motbox_path = col_character(),
#                             expName = col_character()
#                           )
#          )


# http://jenrichmond.rbind.io/post/use-map-to-read-many-csv-files/
df<- fnames %>%  #takes the filenames and reads each datafile into a tibble and then combines them all
  map_dfr(function(x) 
           read_csv(file.path(local_data_pth, x),
            col_types = cols(.default = col_double(),
                               t_contr = col_character(),
                               mark_type = col_character(),
                               mouse.x = col_character(),
                               mouse.y = col_character(),
                               mouse.leftButton = col_character(),
                               mouse.midButton = col_character(),
                               mouse.rightButton = col_character(),
                               mouse.time = col_character(),
                               mouse.clicked_name = col_character(),
                               date = col_character(),
                               motbox_path = col_character(),
                               expName = col_character()
                             )
             )
  )
#Practice trials have some different columns than real trials, which explains many of the "NA"s

#Be aware that if you quit the experiment early, the final trial will be almost entirely blank and filled in with "NA"
incompleteTrials<- df %>% filter(is.na(participant))
df<- df %>% filter(!is.na(participant))
msg=paste("Incomplete (because terminated early?) = ",as.character(nrow(incompleteTrials)), "trials, now deleted.")
print(msg)

print("Here is a table of participant number versus date/time of session:")
table(df$participant,df$date)
#Don't specify the column types for now and fix it later
# df<- fnames %>%  #takes the filenames and reads each datafile into a tibble and then combines them all
#   map_dfr(function(x) 
#             read_csv(file.path(local_data_pth, x),
#                 col_types = cols(.default= col_double()))
#           )

#Original Filip
# df <- fnames %>%
#   map(~ read_csv(file.path(local_data_pth, .), col_types = cols(.default = col_double(),
#                                                                 t_contr = col_character(),
#                                                                 mark_type = col_character(),
#                                                                 mouse.x = col_character(),
#                                                                 mouse.y = col_character(),
#                                                                 mouse.leftButton = col_character(),
#                                                                 mouse.midButton = col_character(),
#                                                                 mouse.rightButton = col_character(),
#                                                                 mouse.time = col_character(),
#                                                                 mouse.clicked_name = col_character(),
#                                                                 date = col_character(),
#                                                                 motbox_path = col_character(),
#                                                                 expName = col_character(),
#                                                                 X27 = col_logical()))) %>% 
#   reduce(rbind) %>% as_tibble()

#Eliminate all but the important columns and rename a few
dg <-df %>% 
  select(subject_id = participant, 
         protocol_id = prot_id,
         trial_id = trial_id,
         noise_id, trajectory_id, t_contr, mark_type,
         mouse_x = mouse.x,mouse_y = mouse.y, 
         mouse.clicked_name,
         rt = mouse.time) 
#mouse.clicked_name is a list of which objects were clicked with the mouse.
#Turn it into separate fields, with one row for each trial
answ <- dg$mouse.clicked_name %>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  str_remove_all("'") %>% 
  str_remove_all(" ") %>% 
  str_remove_all("o1_copy_") %>% 
  str_split(",", simplify = T)

#The targets are always object numbers 0, 1, 2, 3, so check how many of them the participant clicked on.
dg$nCorrect <- matrix(data = answ %in% c("0","1","2","3"), ncol = 4) %>% rowSums()
dg$acc <- dg$nCorrect/4
dg$mark_type[dg$t_contr == "MOT"] <- "noMark"

#delete the practice trials
dg <- dg %>% filter(protocol_id!=0) 

write_rds(dg, here("exp1","data_processed",outputFname))


# prepare trajectories ----------------------------------------------------

trajectories_pth <- file.path(here("exp1"),"trajectories")

files <- dir(trajectories_pth, pattern = "*.csv") # get file names

df_trajectories <- files %>% map_dfr(function(x) 
                     read_csv(file.path(trajectories_pth, x), col_names =F, col_types = cols(.default = "d")))

#Why are there 500 trajectory ids and 81 of each?
df_trajectories$trajectory_id <- rep(1:500, each = nrow(df_trajectories)/500)

colnames(df_trajectories)[1:17] <- c("t",
                                     expand.grid(c("X","y"),1:8) %>% tidyr::unite(cn, Var1, Var2,sep="") %>% pull(cn))

#This is just to reorder the columns
df_trajectories <- df_trajectories %>% select(trajectory_id,everything())

saveRDS(df_trajectories, file.path(outpth, "trajectories_200406.rds"))
