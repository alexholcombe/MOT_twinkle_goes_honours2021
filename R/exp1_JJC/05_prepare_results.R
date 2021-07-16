library(tidyverse)
theme_set(theme_classic(16))

library(here)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp1_JJC"), "data")
#dir(local_data_pth)

outpth <- here("exp1_JJC")

files <- dir(local_data_pth, pattern = "*.csv") # get file names
files <- tibble(fname=files)

outputFname = "data_exp1.rds"

testWithTwoFiles<-T #hand-pick two files for testing purposes
if (testWithTwoFiles) {
  outputFname = "data_testAOH.rds"
}

if (testWithTwoFiles) {
  files <- files %>% filter(fname=="999_noiseMot_exp1_noise_2021_Jul_15_1439.csv" | 
                          fname=="999_noiseMot_exp1_noise_2021_Jul_16_0935_1.csv")
}
fnames <- as.array(files$fname)
msg<- paste0("Number of files = ",as.character(nrow(fnames)), ".")
print(msg)

#Try reading an individual file
try_reading_one_file <- T
if (try_reading_one_file) {
  #testfile <- "999_noiseMot_exp1_noise_2021_May_10_1219.csv"
  testfile<- file.path(local_data_pth, fnames[1])
  df_try <- read_csv(testfile)
  # df_try <- read_csv(testfile,
  #                            col_types = cols(.default = col_double(),
  #                                mark_type = col_character(),
  #                                mouse.x = col_character(),
  #                                mouse.y = col_character(),
  #                                mouse.leftButton = col_character(),
  #                                mouse.midButton = col_character(),
  #                                mouse.rightButton = col_character(),
  #                                mouse.clicked_name = col_character(),
  #                                date = col_character(),
  #                                motbox_path = col_character(),
  #                                expName = col_character()
  #               )
  # )
  
  #If the last trial is not completed, then some columns such as prot_id will be "NA"
}  

# http://jenrichmond.rbind.io/post/use-map-to-read-many-csv-files/
df<- fnames %>%                  # read each file into a tibble and combine them all
  map_dfr(function(x) 
           read_csv(file.path(local_data_pth, x)
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
  select(participant, 
         protocol = prot_id,
         trial = trial_id,
         trajectory_id,
         practice_trials.thisN,
         trials.thisN,
         mouse.x.firstBadClick, mouse.y.firstBadClick,
         mouse.x, mouse.y,
         obj0finalX = P.objects[0].finalx0, 
         obj1finalY = P.objects[0].finaly0,
         timingHiccupsInLastFramesOfStimuli,
         win.nDroppedFrames,
         RT = mouse.time) 

#Investigate the breakdown of practice and real trials
print("Here is a table of participant number versus protocol (1 = real trials):")
table(dg$participant,dg$protocol) #Because prot_id = 0 means practice trials, 1=real trials

#delete the practice trials
dh <- dg %>% filter(protocol !=0) 

#Check number of timing hiccups
print("Here is a table of participant number versus timing hiccups in last frames per trial:")
table(dh$participant,dh$timingHiccupsInLastFramesOfStimuli) 


#Calculate distance between mouse.click and target
dh %>% mutate(xErr = obj0finalX - mouse.x,
              yErr = obj0finalY - mouse.y)

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
