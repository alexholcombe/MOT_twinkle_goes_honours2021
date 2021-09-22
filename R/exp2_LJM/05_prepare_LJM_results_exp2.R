if(!"quickpsy" %in% rownames(installed.packages())) install.packages("quickpsy")
if(!"ggpubr" %in% rownames(installed.packages())) install.packages("ggpubr")

library(tidyverse)
theme_set(theme_classic(16))
library(MPDiR)
library(dplyr)
library(quickpsy)
library(here)
library(Hmisc)
library(ggpubr)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp2_LJM"), "data")
#dir(local_data_pth)

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
analyse_only_most_recent_file <- F
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
                  df$participant <- as.character(df$participant)
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
         response,
         wasUp,
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
dhh <- dhh %>% mutate(noise_status = if (grepl('dynamic',condition) == TRUE) {'Dynamic'}
                      else {'Static'})

#Calculate final velocity from preantepenultimatex0
dhh <- dhh %>% mutate(dx = obj0finalX - obj0preantepenultimateX,
                      dy = obj0finalY - obj0preantepenultimateY)
dhh<- dhh %>% mutate(finalSpeedCheck =    sqrt(dx*dx + dy*dy) / (3*avgDurLastFrames/1000))

#final speed
ggplot(dhh, aes(finalSpeedCheck)) + geom_histogram()

#anonymise participant results
dhh <- dhh %>% mutate(participant = if (grepl('999',participant)) {'A'}
                      else if (grepl('jc',participant)) {'B'}
                      else if (grepl('RN',participant)) {'C'}
                      else if (grepl('jvr',participant)) {'D'}
                      else if (grepl('OSAY',participant)) {'E'})

#Separate data by participant
datasetPA <- dhh %>% filter(participant=="A")
datasetPB <- dhh %>% filter(participant=="B")
datasetPC <- dhh %>% filter(participant=="C")
datasetPD <- dhh %>% filter(participant=="D")
datasetPE <- dhh %>% filter(participant=="E")

alldatasets <- list(datasetPA,datasetPB,datasetPC,datasetPD,datasetPE)

for (dataset in alldatasets) {
#create proportionUp, which is the mean of 'up's (1) and 'down's (0) for each offset in each condition
datasetUp <- dataset %>%
  group_by(participant, first_speed,final_speed,noise_status, offset) %>%
  summarise(proportionUp = mean(wasUp))

#contrast code the interactions between the three levels (but not the three-way interaction)
datasetUp <- datasetUp %>% mutate(noise_firstspeed = if ((grepl('Dynamic',noise_status) == TRUE) && (first_speed == 500)) {1}
                                  else if ((grepl('Dynamic',noise_status) == TRUE) && (first_speed == 1000)) {0}
                                  else if ((grepl('Static',noise_status) == TRUE) && (first_speed == 500)) {0}
                                  else if ((grepl('Static',noise_status) == TRUE) && (first_speed == 1000)) {1})
datasetUp <- datasetUp %>% mutate(noise_finalspeed = if ((grepl('Dynamic',noise_status) == TRUE) && (final_speed == 500)) {1}
           else if ((grepl('Dynamic',noise_status) == TRUE) && (final_speed == 1000)) {0}
           else if ((grepl('Static',noise_status) == TRUE) && (final_speed == 500)) {0}
           else if ((grepl('Static',noise_status) == TRUE) && (final_speed == 1000)) {1})
datasetUp <- datasetUp %>% mutate(finalspeed_firstspeed = if ((first_speed == 500) && (final_speed == 500)) {1}
           else if ((first_speed == 500) && (final_speed == 1000)) {0}
           else if ((first_speed == 1000) && (final_speed == 500)) {0}
           else if ((first_speed == 1000) && (final_speed == 1000)) {1})

datasetUp <- datasetUp %>% mutate(offsetInDegrees = ifelse(offset == 0, 0, 0.1*(2^abs(offset))))
for (n in 1:nrow(datasetUp)) {
  if (datasetUp$offset[n] < 0)
    datasetUp$offsetInDegrees[n] <- -1*datasetUp$offsetInDegrees[n]
}

#count how many times each offset is used in each condition
datasetOffset <- dataset %>% group_by(first_speed,final_speed,noise_status) %>% count(offset, name = 'timesOffsetUsed')

#make a table with important data, including proportionUp and timesOffsetUsed
#also run the quickpsy models
if (datasetUp$participant[1] == 'A') {
  datasetPA <- cbind(datasetUp,datasetOffset[,5])
  datasetPA <- datasetPA %>% mutate(nUp = proportionUp*timesOffsetUsed)
  PAmodelWithJustNoise <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 500)
  PAmodelWithJustFirstSpeed <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 500)
  PAmodelWithJustFinalSpeed <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 500)
  PAmodelWithNoiseFirstSpeedInteraction <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 500)
  PAmodelWithNoiseFinalSpeedInteraction <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 500)
  PAmodelWithFirstSpeedFinalSpeedInteraction <- quickpsy(datasetPA,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 500)
} else if (datasetUp$participant[1] == 'B') {
  datasetPB <- cbind(datasetUp,datasetOffset[,5])
  datasetPB <- datasetPB %>% mutate(nUp = proportionUp*timesOffsetUsed)
  PBmodelWithJustNoise <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 500)
  PBmodelWithJustFirstSpeed <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 500)
  PBmodelWithJustFinalSpeed <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 500)
  PBmodelWithNoiseFirstSpeedInteraction <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 500)
  PBmodelWithNoiseFinalSpeedInteraction <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 500)
  PBmodelWithFirstSpeedFinalSpeedInteraction <- quickpsy(datasetPB,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 500)
} else if (datasetUp$participant[1] == 'C') {
  datasetPC <- cbind(datasetUp,datasetOffset[,5])
  datasetPC <- datasetPC %>% mutate(nUp = proportionUp*timesOffsetUsed)
  PCmodelWithJustNoise <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 500)
  PCmodelWithJustFirstSpeed <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 500)
  PCmodelWithJustFinalSpeed <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 500)
  PCmodelWithNoiseFirstSpeedInteraction <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 500)
  PCmodelWithNoiseFinalSpeedInteraction <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 500)
  PCmodelWithFirstSpeedFinalSpeedInteraction <- quickpsy(datasetPC,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 500)
} else if (datasetUp$participant[1] == 'D') {
  datasetPD <- cbind(datasetUp,datasetOffset[,5])
  datasetPD <- datasetPD %>% mutate(nUp = proportionUp*timesOffsetUsed)
  PDmodelWithJustNoise <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 500)
  PDmodelWithJustFirstSpeed <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 500)
  PDmodelWithJustFinalSpeed <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 500)
  PDmodelWithNoiseFirstSpeedInteraction <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 500)
  PDmodelWithNoiseFinalSpeedInteraction <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 500)
  PDmodelWithFirstSpeedFinalSpeedInteraction <- quickpsy(datasetPD,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 500)
} else if (datasetUp$participant[1] == 'E') {
  datasetPE <- cbind(datasetUp,datasetOffset[,5])
  datasetPE <- datasetPE %>% mutate(nUp = proportionUp*timesOffsetUsed)
  PEmodelWithJustNoise <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_status), B = 500)
  PEmodelWithJustFirstSpeed <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(first_speed), B = 500)
  PEmodelWithJustFinalSpeed <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(final_speed), B = 500)
  PEmodelWithNoiseFirstSpeedInteraction <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_firstspeed), B = 500)
  PEmodelWithNoiseFinalSpeedInteraction <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(noise_finalspeed), B = 500)
  PEmodelWithFirstSpeedFinalSpeedInteraction <- quickpsy(datasetPE,offsetInDegrees,nUp,timesOffsetUsed,grouping = .(finalspeed_firstspeed), B = 500)
}
}

alldatasets <- list(datasetPA,datasetPB,datasetPC,datasetPD,datasetPE)

#bar graphs of noise threshold difference
PAnoiseplot <- plotthresholds(PAmodelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBnoiseplot <- plotthresholds(PBmodelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCnoiseplot <- plotthresholds(PCmodelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDnoiseplot <- plotthresholds(PDmodelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PEnoiseplot <- plotthresholds(PEmodelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PAnoiseplot,PBnoiseplot,PCnoiseplot,PDnoiseplot,PEnoiseplot, common.legend = TRUE)

#plot of noise curve differences
PAnoisecurves <- plotcurves(PAmodelWithJustNoise) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBnoisecurves <- plotcurves(PBmodelWithJustNoise) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCnoisecurves <- plotcurves(PCmodelWithJustNoise) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDnoisecurves <- plotcurves(PDmodelWithJustNoise) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PEnoisecurves <- plotcurves(PEmodelWithJustNoise) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("Dynamic" = "orange", "Static" = "gray")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PAnoisecurves,PBnoisecurves,PCnoisecurves,PDnoisecurves,PEnoisecurves, common.legend = TRUE)

#plot of initial speed curve differences
PAFirstSpeedcurves <- plotcurves(PAmodelWithJustFirstSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBFirstSpeedcurves <- plotcurves(PBmodelWithJustFirstSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCFirstSpeedcurves <- plotcurves(PCmodelWithJustFirstSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDFirstSpeedcurves <- plotcurves(PDmodelWithJustFirstSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PEFirstSpeedcurves <- plotcurves(PEmodelWithJustFirstSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PAFirstSpeedcurves,PBFirstSpeedcurves,PCFirstSpeedcurves,PDFirstSpeedcurves,PEFirstSpeedcurves, common.legend = TRUE)

#plot of final speed curve differences
PAFinalSpeedcurves <- plotcurves(PAmodelWithJustFinalSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBFinalSpeedcurves <- plotcurves(PBmodelWithJustFinalSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCFinalSpeedcurves <- plotcurves(PCmodelWithJustFinalSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDFinalSpeedcurves <- plotcurves(PDmodelWithJustFinalSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PEFinalSpeedcurves <- plotcurves(PEmodelWithJustFinalSpeed) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("500" = "blue", "1000" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PAFinalSpeedcurves,PBFinalSpeedcurves,PCFinalSpeedcurves,PDFinalSpeedcurves,PEFinalSpeedcurves, common.legend = TRUE)

#plot of noise first speed interaction differences
PANoiseFirstSpeedInteractioncurves <- plotcurves(PAmodelWithNoiseFirstSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBNoiseFirstSpeedInteractioncurves <- plotcurves(PBmodelWithNoiseFirstSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCNoiseFirstSpeedInteractioncurves <- plotcurves(PCmodelWithNoiseFirstSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDNoiseFirstSpeedInteractioncurves <- plotcurves(PDmodelWithNoiseFirstSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PENoiseFirstSpeedInteractioncurves <- plotcurves(PEmodelWithNoiseFirstSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PANoiseFirstSpeedInteractioncurves,PBNoiseFirstSpeedInteractioncurves,PCNoiseFirstSpeedInteractioncurves,
          PDNoiseFirstSpeedInteractioncurves,PENoiseFirstSpeedInteractioncurves, common.legend = TRUE)

#plot of noise final speed interaction differences
PANoiseFinalSpeedInteractioncurves <- plotcurves(PAmodelWithNoiseFinalSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("A")
PBNoiseFinalSpeedInteractioncurves <- plotcurves(PBmodelWithNoiseFinalSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("B")
PCNoiseFinalSpeedInteractioncurves <- plotcurves(PCmodelWithNoiseFinalSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) +
  theme(plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("C")
PDNoiseFinalSpeedInteractioncurves <- plotcurves(PDmodelWithNoiseFinalSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("D*")
PENoiseFinalSpeedInteractioncurves <- plotcurves(PEmodelWithNoiseFinalSpeedInteraction) + labs(x="Offset (deg)", y="Probability of an 'UP' response") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) + expand_limits(y = c(-2,2)) + 
  theme(legend.position = "none",plot.title = element_text(size = 32,hjust = .5, face = "bold")) + ggtitle("E*")
ggarrange(PANoiseFinalSpeedInteractioncurves,PBNoiseFinalSpeedInteractioncurves,PCNoiseFinalSpeedInteractioncurves,
          PDNoiseFinalSpeedInteractioncurves,PENoiseFinalSpeedInteractioncurves, common.legend = TRUE)

modelWithJustNoise
modelWithJustFirstSpeed
modelWithJustFinalSpeed
modelWithNoiseFirstSpeedInteraction
modelWithNoiseFinalSpeedInteraction
modelWithFirstSpeedFinalSpeedInteraction

#plot 'em
plotcurves(PAmodelWithJustNoise)
plotcurves(modelWithJustFirstSpeed)
plotcurves(modelWithJustFinalSpeed)
plotcurves(modelWithNoiseFirstSpeedInteraction)
plotcurves(modelWithNoiseFinalSpeedInteraction)
plotcurves(modelWithFirstSpeedFinalSpeedInteraction)

plotthresholds(modelWithJustNoise) + labs(x="Noise condition", y="Offset (deg)", fill ="Noise condition") +
  scale_fill_manual(values = c("Dynamic" = "orange", "Static" = "gray"))
plotthresholds(modelWithJustFirstSpeed)
plotthresholds(modelWithJustFinalSpeed)
plotthresholds(modelWithNoiseFirstSpeedInteraction)
plotthresholds(modelWithNoiseFinalSpeedInteraction)
plotthresholds(modelWithFirstSpeedFinalSpeedInteraction)

#t test of difference between noise conditions
turnGroupCIsIntoGroupDifferenceCIs <- function(model) {
  cumuCondSE <- 0
  numConds <- 2
  firstGroupDiffCI <- (model$thresholds$thresup[1] - model$thresholds$thre[1])/.95*.83
  secondGroupDiffCI <- ( model$thresholds$thresup[2] - model$thresholds$thre[2])/.95*.83
  groupDiffModel <- model
  groupDiffModel$thresholds$threinf[1] = (model$thresholds$thre[1] - firstGroupDiffCI)
  groupDiffModel$thresholds$thresup[1] = (model$thresholds$thre[1] + firstGroupDiffCI)
  groupDiffModel$thresholds$threinf[2] = (model$thresholds$thre[2] - firstGroupDiffCI)
  groupDiffModel$thresholds$thresup[2] = (model$thresholds$thre[2] + firstGroupDiffCI)
  groupDiffModel
}
turnGroupCIsIntoGroupDifferenceCIs(modelWithJustNoise)
turnGroupCIsIntoGroupDifferenceCIs(modelWithJustFirstSpeed)
turnGroupCIsIntoGroupDifferenceCIs(modelWithJustFinalSpeed)
turnGroupCIsIntoGroupDifferenceCIs(modelWithNoiseFirstSpeedInteraction)
turnGroupCIsIntoGroupDifferenceCIs(modelWithNoiseFinalSpeedInteraction)
turnGroupCIsIntoGroupDifferenceCIs(modelWithFirstSpeedFinalSpeedInteraction)

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

