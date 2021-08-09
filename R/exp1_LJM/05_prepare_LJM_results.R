library(tidyverse)
theme_set(theme_classic(16))

library(here)
library(Hmisc)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp1_LJM"), "data", "cjh_exp1")
#dir(local_data_pth)

outpth <- here("exp1_LJM")

#Get file list and sort by creation time, from most recent to oldest
#https://stackoverflow.com/a/13762544/302378
files_with_details = file.info(  list.files(local_data_pth,pattern="*.csv")        )
files_with_details = files_with_details[with(files_with_details,order(as.POSIXct(ctime))), ]
files = rownames(files_with_details) #sorted list of files
files <- tibble(files) %>% arrange(-row_number())

outputFname = "data_exp1.rds"

testWithTwoFiles<- FALSE #hand-pick two files for testing purposes
if (testWithTwoFiles) {
  outputFname = "data_test.rds"
}
if (testWithTwoFiles) {
  files <- files %>% filter(fname=="999_noiseMot_exp1_noise_2021_Jul_15_1439.csv" | 
                              fname=="999_noiseMot_exp1_noise_2021_Jul_16_0935_1.csv")
}
fnames <- as.array(files$files)
msg<- paste0("Number of files = ",as.character(nrow(fnames)), ".")
print(msg)

#Try reading an individual file
analyse_only_most_recent_file <- F
if (analyse_only_most_recent_file) {
  #testfile <- "999_noiseMot_exp1_noise_2021_May_10_1219.csv"
  testfile<- files$files[1]
  df <- read_csv(  file.path(local_data_pth, testfile) )
  
  #If the last trial is not completed, then some columns such as prot_id will be "NA"
}  else {
  # http://jenrichmond.rbind.io/post/use-map-to-read-many-csv-files/
  df<- fnames %>%                  # read each file into a tibble and combine them all
    map_dfr(function(x) 
      read_csv(file.path(local_data_pth, x)
      )
    )
}
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

#Eliminate all but the important columns and rename a few
dg <-df %>% 
  select(participant, 
         protocol = prot_id,
         trial = trial_id,
         trajectory_id,
         first_speed,
         final_speed, 
         noise_present,
         practice_trials.thisN,
         trials.thisN,
         mouse.x.firstBadClick, mouse.y.firstBadClick,
         mouse.x, mouse.y,
         obj0finalX = `P.objects[0].finalx0`, 
         obj0finalY = `P.objects[0].finaly0`,
         obj0penultimateX = `P.objects[0].penultimatex0`,
         obj0penultimateY = `P.objects[0].penultimatey0`,
         obj0preantepenultimateX = `P.objects[0].preantepenultimatex0`,
         obj0preantepenultimateY = `P.objects[0].preantepenultimatey0`,
         obj0antepenultimateX = `P.objects[0].antepenultimatex0`,
         obj0antepenultimateY = `P.objects[0].antepenultimatey0`,
         timingHiccupsInLastFramesOfStimuli,
         win.nDroppedFrames,
         finalStimFrameTimes,
         RT = mouse.time) 

#Investigate the breakdown of practice and real trials
print("Here is a table of participant number versus protocol (1 = real trials):")
table(dg$participant,dg$protocol) #Because prot_id = 0 means practice trials, 1=real trials

#delete the practice trials
dh <- dg %>% filter(protocol !=0) 

#Check number of timing hiccups
print("Here is a table of participant number versus num timing hiccups in last frames of trials:")
table(dh$participant,dh$timingHiccupsInLastFramesOfStimuli) 

dh <- dh %>% mutate(anyHiccupsLastFrames = timingHiccupsInLastFramesOfStimuli > 0)

print("For each participant num trials with hiccups in last frames, and avg missed frames in whole trial for those trials:")
dh %>% filter( anyHiccupsLastFrames > 0)  %>% 
       group_by(participant)  %>%
       summarise( trialsWithHiccupsLastFrames=  n(), avgMissedFramesInWholeTrial = mean(win.nDroppedFrames) )
#Finished timing checks

#Remove trials with timingHiccupsInLastFramesOfStimuli
dh <- dh %>% filter(timingHiccupsInLastFramesOfStimuli == 0)

#Create function that can parse array of finalStimFrameTimes
calcAvgOfListFromPsychopy <- function( listAsTextFromPsychopy ){
  #Psychopy lists, like the frametimes the program outputs, end up as character vectors (strings)
  #like this:
  # "[16.56737499 16.67500002 16.98479103 16.32583397\n
  # 16.76725002 16.52166597 16.80783404 16.22549997 ]"
  
  #Clean out brackets and whitespace including newlines
  cleaned <- listAsTextFromPsychopy %>% str_remove(fixed('[')) %>% str_remove(fixed(']')) %>% str_squish()
  #rely on the remaining whitespace to allow read.table to parse it
  parsed <- read.table(textConnection(cleaned))
  #print( round(t(parsed)-16) )
  avg = mean( t(parsed) ) #Transpose into a single column, then take mean
  return(avg)
}

#Try to apply that function to every row, calculating the mean refresh rate for the last trials
dhh<- dh %>% rowwise() %>% #If don't call rowwise(), have to vectorise
  mutate(avgDurLastFrames = calcAvgOfListFromPsychopy(finalStimFrameTimes) )

#Calculate final velocity from preantepenultimatex0
dhh <- dhh %>% mutate(dx = obj0finalX - obj0preantepenultimateX,
                      dy = obj0finalY - obj0preantepenultimateY)
dhh<- dhh %>% mutate(finalDirection = atan2(dy,dx)/pi*180)
dhh<- dhh %>% mutate(finalSpeed =    sqrt(dx*dx + dy*dy) / (3*avgDurLastFrames/1000))

#final direction
ggplot(dhh, aes(finalDirection)) + geom_histogram()
#final speed
ggplot(dhh, aes(finalSpeed)) + geom_histogram()

#Calculate distance between mouse.click and target
dhh <- dhh %>% mutate(xErr = mouse.x - obj0finalX,
                    yErr = mouse.y - obj0finalY)
ggplot(dhh, aes(xErr,yErr)) + geom_point()

#Calculate component of error relative to last direction
#Project error vector onto direction vector
# as.vector( (u %*% v) / (v %*% v) ) * v #https://stackoverflow.com/a/62495884/302378
# https://en.wikipedia.org/wiki/Vector_projection https://en.wikipedia.org/wiki/Dot_product
# np.dot(x, y) / np.linalg.norm(y)  https://stackoverflow.com/questions/55226011/using-python-to-calculate-vector-projection
lengthProjected <- function(u, v) {
  #Dot product of 
  dotproduct <- (u %*% v)
  lengthOfV <- sqrt( sum(v^2) )
  return( dotproduct / lengthOfV )
  #For whole vector can use the following, although I don't understand it: (as.vector( (u %*% v) / (v %*% v) ) * v)
}


#Calculate length projected, cos() and sin() are in radians
dhh <- dhh %>%
  mutate( amountExtrapolation =
                      lengthProjected( #We want to project the mouse response error vector onto the final direction vector
                      c(xErr,yErr), #mouse response error vector
                      c(cos(finalDirection/180*pi), sin(finalDirection/180*pi)))) #finalDirection vector

#EXCLUDE OUTLIERS
outlierCriterion <- 150
dhh <- dhh %>% mutate(isOutlier = sqrt(xErr^2+yErr^2) > outlierCriterion)

#plot amount of extrapolation
ggplot(dhh %>% filter(isOutlier==FALSE), 
       aes(amountExtrapolation)) + geom_histogram()

#plot amount of extrapolation for each possible speed for 3 final speeds
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(speed1~.)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(speed2~.)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(speed3~.)

#plot effect of noise on amount of extrapolation
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~.)

#plot interaction between noise and each possible speed for 3 final speeds
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~speed1)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~speed2)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~speed3)


gg<- ggplot(dhh %>% filter(isOutlier==FALSE),
      aes(x=noise_present,y=amountExtrapolation)) + geom_point() + 
      geom_hline(yintercept=0) + 
            stat_summary(fun.data = "mean_cl_boot", color="red", size=1) 
show(gg)

ggplot(dhh %>% filter(isOutlier==FALSE),
       aes(x=noise_present,y=amountExtrapolation)) + geom_point() + 
  geom_hline(yintercept=0) + 
  stat_summary(fun.data = "mean_cl_boot", color="red", size=1) +
  facet_grid(speed1~.)
               
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~no_targets)


#Plot error data

#Calculate various distance metrics

#Save data
write_rds(dg, here("exp1","data_processed",outputFname))


# prepare trajectories for analysis --------------------------------------

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
