if(!"quickpsy" %in% rownames(installed.packages())) install.packages("quickpsy")

library(tidyverse)
theme_set(theme_classic(16))
library(MPDiR)
library(dplyr)
library(quickpsy)
library(here)
library(Hmisc)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp3_LJM"), "data")
dir(local_data_pth)

outpth <- here("exp2_LJM")

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
analyse_only_most_recent_file <- T
if (analyse_only_most_recent_file) {
  #testfile <- "999_noiseMot_exp1_noise_2021_May_10_1219.csv"
  testfile<- files$files[1]
  df <- read_csv(  file.path(local_data_pth, testfile) )
  
  #If the last trial is not completed, then some columns such as prot_id will be "NA"
}  else {
  # http://jenrichmond.rbind.io/post/use-map-to-read-many-csv-files/
  df<- fnames %>%                  # read each file into a tibble and combine them all
    map_dfr(  function(x) {
                  df = read_csv(  file.path(local_data_pth, x) )
                  df
                 } 
            )
}

#The mysterious None case for direction
#dg<-read_csv(  file.path(local_data_pth, "cjh_noiseMot_exp1_noise_2021_Aug_06_0943.csv") ) #chr
#dh<- read_csv(  file.path(local_data_pth, "cjh_noiseMot_exp1_noise_2021_Aug_06_0959.csv") ) #dbl

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
         trial = trialNumber,
         condition,
         flipStatus,
         offset,
         #response,
         #wasUp,
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
         finalStimFrameTimes) 

#Investigate the breakdown of practice and real trials
print("Here is a table of participant number versus protocol (1 = real trials):")
table(dg$participant,dg$trial) #Should be 8 practice trials at the end of this clumsy table

#delete the practice trials
dh <- dg %>% filter(trial!='practice') 

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
#dh <- dh %>% filter(timingHiccupsInLastFramesOfStimuli == 0)

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

#Extract condition details from 'condition' column
dhh <- dhh %>% mutate(first_speed = if (grepl('500x1',condition) == TRUE) {500}
  else if (grepl('500x5',condition)==TRUE) {500} 
  else {1000})
dhh <- dhh %>% mutate(final_speed = if (grepl('0x500',condition) == TRUE) {500} 
                      else {1000})
dhh <- dhh %>% mutate(noise_status = if (grepl('dynamic',condition) == TRUE) {'dynamic'}
                      else {'static'})

#Calculate final velocity from preantepenultimatex0
dhh <- dhh %>% mutate(dx = obj0finalX - obj0preantepenultimateX,
                      dy = obj0finalY - obj0preantepenultimateY)
dhh<- dhh %>% mutate(finalSpeedCheck =    sqrt(dx*dx + dy*dy) / (3*avgDurLastFrames/1000))

#final speed
ggplot(dhh, aes(finalSpeedCheck)) + geom_histogram()

#create proportionUp, which is the mean of 'up's (1) and 'down's (0) for each offset in each condition
dhf <- dhh %>%
  group_by(first_speed,final_speed,noise_status, offset) %>%
  summarise(proportionUp = mean(wasUp))

#count how many times each offset is used in each condition
dhg <- dhh %>% group_by(first_speed,final_speed,noise_status) %>% count(offset, name = 'timesOffsetUsed')

#make a table with important data, including proportionUp and timesOffsetUsed
dhi <- cbind(dhf,dhg[,5])
dhi <- dhi %>% mutate(nUp = proportionUp*timesOffsetUsed)

#contrast code the interactions between the three levels (but not the three-way interaction)
dhi <- dhi %>% mutate(noise_firstspeed = if ((grepl('dynamic',noise_status) == TRUE) && (first_speed == 500)) {1}
                                    else if ((grepl('dynamic',noise_status) == TRUE) && (first_speed == 1000)) {0}
                                    else if ((grepl('static',noise_status) == TRUE) && (first_speed == 500)) {0}
                                    else if ((grepl('static',noise_status) == TRUE) && (first_speed == 1000)) {1})
dhi <- dhi %>% mutate(noise_finalspeed = if ((grepl('dynamic',noise_status) == TRUE) && (final_speed == 500)) {1}
                                    else if ((grepl('dynamic',noise_status) == TRUE) && (final_speed == 1000)) {0}
                                    else if ((grepl('static',noise_status) == TRUE) && (final_speed == 500)) {0}
                                    else if ((grepl('static',noise_status) == TRUE) && (final_speed == 1000)) {1})
dhi <- dhi %>% mutate(finalspeed_firstspeed = if ((first_speed == 500) && (final_speed == 500)) {1}
                                    else if ((first_speed == 500) && (final_speed == 1000)) {0}
                                    else if ((first_speed == 1000) && (final_speed == 500)) {0}
                                    else if ((first_speed == 1000) && (final_speed == 1000)) {1})

#put offset in terms of degrees
dhi <- dhi %>% mutate(offsetInDegrees = ifelse(offset == 0, 0, 0.1*(2^abs(offset))))
for (n in 1:nrow(dhi)) {
  if (dhi$offset[n] < 0)
    dhi$offsetInDegrees[n] <- -1*dhi$offsetInDegrees[n]
}

#run the quickpsy models
modelWithAllConditions <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status,first_speed,final_speed), B = 5)
modelWithJustNoise <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 5)
modelWithJustFirstSpeed <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 5)
modelWithJustFinalSpeed <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 5)
modelWithNoiseFirstSpeedInteraction <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 5)
modelWithNoiseFinalSpeedInteraction <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 5)
modelWithFirstSpeedFinalSpeedInteraction <- quickpsy(dhi,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 5)

modelWithAllConditions

#plot 'em
plotcurves(modelWithAllConditions)
plotcurves(modelWithJustNoise)
plotcurves(modelWithJustFirstSpeed)
plotcurves(modelWithJustFinalSpeed)
plotcurves(modelWithNoiseFirstSpeedInteraction)
plotcurves(modelWithNoiseFinalSpeedInteraction)
plotcurves(modelWithFirstSpeedFinalSpeedInteraction)

plotthresholds(modelWithAllConditions, x = final_speed, xpanel = first_speed, color = noise_status)
plotthresholds(modelWithJustNoise)
plotthresholds(modelWithJustFirstSpeed)
plotthresholds(modelWithJustFinalSpeed)
plotthresholds(modelWithNoiseFirstSpeedInteraction)
plotthresholds(modelWithNoiseFinalSpeedInteraction)
plotthresholds(modelWithFirstSpeedFinalSpeedInteraction)

#having a go at calculating what I think could be effect size for difference between noise conditions
threInfSupPooledStandardError <- function(model, N) {
  cumuCondSE <- 0
  numConds <- nrow(model$thresholds) 
  for (c in 1:numConds) {
    cumuCondSE <- cumuCondSE + ((model$thresholds$thresup[c] - model$thresholds$threinf[c])/3.92)
  }
  pooledSE <- cumuCondSE/numConds
  pooledSD <<- pooledSE * sqrt(N)
}
threInfSupPooledStandardError(modelWithJustNoise, nrow(dhh))
noiseThresholdMeanDiff <- modelWithJustNoise$thresholds$thre[1] - modelWithJustNoise$thresholds$thre[2]
noiseEffectSize <- noiseThresholdMeanDiff/pooledSD
noiseEffectSize

#effect size for difference between first_speed conditions
threInfSupPooledStandardError(modelWithJustFirstSpeed, nrow(dhh))
firstSpeedThresholdMeanDiff <- modelWithJustFirstSpeed$thresholds$thre[1] - modelWithJustFirstSpeed$thresholds$thre[2]
firstSpeedEffectSize <- firstSpeedThresholdMeanDiff/pooledSD
firstSpeedEffectSize

#effect size for difference between final_speed conditions
threInfSupPooledStandardError(modelWithJustFinalSpeed, nrow(dhh))
finalSpeedThresholdMeanDiff <- modelWithJustFinalSpeed$thresholds$thre[1] - modelWithJustFinalSpeed$thresholds$thre[2]
finalSpeedEffectSize <- finalSpeedThresholdMeanDiff/pooledSD
finalSpeedEffectSize

meanEffectSize <- (abs(noiseEffectSize) + abs(firstSpeedEffectSize) + abs(finalSpeedEffectSize))/3
meanEffectSize

#Save data
write_rds(dg, here("exp2","data_processed",outputFname))


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
