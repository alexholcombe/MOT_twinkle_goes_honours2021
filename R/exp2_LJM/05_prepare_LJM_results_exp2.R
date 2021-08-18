library(tidyverse)
theme_set(theme_classic(16))

library(here)
library(Hmisc)

source(here("utils.R"))

local_data_pth <- file.path(here("..","exp2_LJM"), "data")
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
dhh<- dhh %>% mutate(finalSpeed =    sqrt(dx*dx + dy*dy) / (3*avgDurLastFrames/1000))

#final speed
ggplot(dhh, aes(finalSpeed)) + geom_histogram()

#log transformation of wasUp
dhh <- dhh %>% mutate(logwasUp = log(wasUp+1))



#Analyse response by offset and condition?
dii <- dhh %>%
  group_by(condition,offset) %>% 
  summarise(proportionUp = mean(wasUp), se = sd(wasUp)/sqrt(nrow(dhh)))

ggplot(dii, aes(x=offset, y=proportionUp)) + geom_point() + stat_smooth(method = "lm", formula = y ~ x, size = 1) + facet_grid(condition~.) 
                        
#EXCLUDE OUTLIERS
outlierCriterion <- 220 # 150
dhh <- dhh %>% mutate(isOutlier = sqrt(xErr^2+yErr^2) > outlierCriterion)

dhh %>% filter(isOutlier==FALSE) %>%
  group_by(noise_present, first_speed, final_speed) %>% 
  summarise(average = mean(amountExtrapolation), se = sd(amountExtrapolation)/sqrt(nrow(dhh)))

ggplot(dhh, aes(x=offset, y=logwasUp)) + geom_point() + facet_grid(condition~.) 

#get means and standard error

#plot amount of extrapolation
ggplot(dhh %>% filter(isOutlier==FALSE), 
       aes(amountExtrapolation)) + geom_histogram()

#https://datavizpyr.com/rain-cloud-plots-using-half-violin-plot-with-jittered-data-points-in-r/
#Load half violin plot: geom_flat_violin()
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

#plot amount of extrapolation for each possible speed for first and final speed
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(first_speed~.)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(final_speed~.)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(first_speed~final_speed)

#plot effect of noise on amount of extrapolation
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~.)

#plot interaction between noise and each possible speed for first and final speed
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~first_speed)
ggplot(dhh, aes(amountExtrapolation)) + geom_histogram() + facet_grid(noise_present~final_speed)
ggplot(dhh, aes(x=noise_present, y=amountExtrapolation)) + geom_point() + facet_grid(first_speed~final_speed) + stat_summary(fun.data = mean_cl_boot, fun.args=(conf.int=0.95), 
                                                                                                                             geom="errorbar", size=2, width=0.2, color='green4', alpha=0.84) 

gg<- ggplot(dhh %>% filter(isOutlier==FALSE),
       aes(x=noise_present,y=amountExtrapolation)) +
  geom_jitter(alpha=0.1, size=.5, width=0.15, height=0) + #geom_point() + 
  geom_hline(yintercept=0) + 
  #geom_flat_violin(fill="gray80",color="gray80") +
  stat_summary(fun.data = mean_cl_boot, fun.args=(conf.int=0.95), 
               geom="errorbar", size=2, width=0.2, color='green4', alpha=0.82) +
  stat_summary(fun.data = "mean_cl_boot", color="green", size=.5) +
  facet_grid(first_speed~.)
show(gg)
#ggsave( file.path("figures","MultipleStudiesPrevalencePerceptionDistributionsComplete.png"), width = 50, height = 30, units = "cm" )

#Twinkle goes t-test
nonoise = dhh %>% filter(isOutlier==FALSE , noise_present=="no_noise")
noise = dhh %>% filter(isOutlier==FALSE , noise_present=="noise")
t.test(nonoise$amountExtrapolation,noise$amountExtrapolation)

#Temporal integration t-test (although maybe a t-test is not appropriate, these conditions are meant to be similar rather than different)
slowinitialfastfinal = dhh %>% filter(isOutlier==FALSE , first_speed == 500, final_speed == 1000)
fastinitialslowfinal = dhh %>% filter(isOutlier==FALSE , first_speed == 1000, final_speed == 500)
t.test(slowinitialfastfinal$amountExtrapolation,fastinitialslowfinal$amountExtrapolation)

#Effect of speed t-test
slowinitialslowfinal = dhh %>% filter(isOutlier==FALSE , first_speed == 500, final_speed == 500)
fastinitialfastfinal = dhh %>% filter(isOutlier==FALSE , first_speed == 1000, final_speed == 1000)
t.test(slowinitialslowfinal$amountExtrapolation,fastinitialfastfinal$amountExtrapolation)

#Graph effect of eccentricity on extrapolation by noise
dhh <- dhh %>% mutate(eccentricity = sqrt(obj0finalX^2 + obj0finalY^2))
ggplot(dhh %>% filter(isOutlier==FALSE), aes(x=eccentricity, y=amountExtrapolation, color=noise_present)) +
  geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

#Calculate linear regression model 
# model <- lm(dhh, formula = amountExtrapolation ~ noise_present + first_speed + final_speed +
#               first_speed*final_speed + first_speed*noise_present + final_speed*noise_present +
#               noise_present*first_speed*final_speed)
# summary(model)
options(contrasts = c("contr.sum","contr.poly"))

finish_line <- lm(amountExtrapolation ~ first_speed + final_speed + noise_present, data = dhh)

drop1(finish_line, .~., test="F")

#Same model with eccentricity as a confound
options(contrasts = c("contr.sum","contr.poly"))

finish_line <- lm(amountExtrapolation ~ first_speed + final_speed + noise_present + eccentricity, data = dhh)

drop1(finish_line, .~., test="F")

#Linear regression model with only noise, eccentricity, and their interaction term
modeleccentricityinteraction <- lm(dhh, formula = amountExtrapolation ~ noise_present + 
                                     eccentricity + velocityToFovea + noise_present*velocityToFovea + 
                                     eccentricity*velocityToFovea)
summary(modeleccentricityinteraction)


modelveltofoveainteraction <- lm(dhh, formula = amountExtrapolation ~ noise_present + velocityToFovea + noise_present*velocityToFovea)
modelveltofoveainteraction <- lm(dhh, formula = amountExtrapolation ~ noise_present)
summary(modelveltofoveainteraction)

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
