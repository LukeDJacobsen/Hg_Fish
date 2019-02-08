################################################################################
#Dependencies
################################################################################
library(tidyverse)
library(lubridate)
library(tidyquant)
library(ranger)
library(readxl)
##################################################################

################################################################################
#Import data 
################################################################################
fish_clean <- read_excel("Documents/Hg_fish/Hg_Fish_Project_Git/raw_hg_fish.xlsx")
################################################################################


################################################################################
#Temporal investigation
################################################################################
fish_clean$DATEI <- as.integer(fish_clean$DATECOL)
plot(fish_clean$DATEI, fish_clean$HGPPM)
#A few improperly added to database
fish_clean %>% filter(DATEI < 10000000) %>% select(DATECOL)
#need to add 0 infront of 6 for the fish collected on 20110617
fish_clean$DATEI <- ifelse(fish_clean$DATEI < 10000000, 20110617, fish_clean$DATEI)

#Use lubridate to more easily manage data
fish_clean$DATE <- ymd(fish_clean$DATEI)
fish_hg_date <- fish_clean %>% select(DATE, HGPPM)

#some visualizations of temporal trends
ggplot(fish_hg_date, aes(DATE, HGPPM)) + geom_point() + geom_smooth()

ggplot(fish_hg_date, aes(DATE, HGPPM)) + geom_smooth(se = FALSE) + 
  coord_cartesian(ylim = c(.2,.5))

#is there within year seasonality? Doesn't really look like it. 
fish_hg <- fish_hg_date %>% mutate(Year = year(DATE))
fish_hg$DATE1 <- as.integer(fish_hg$DATE)
ggplot(fish_hg, aes(DATE1, HGPPM)) + geom_point()   + geom_smooth(se = F) +
  facet_wrap(~Year, scales = "free") +  theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
################################################################################
#Data cleaning: Fix Bad variable titles in fish_clean. Fix structure issues
################################################################################
names(fish_clean)
names(fish_clean)[8] <- "EcoRegionIII"
names(fish_clean)[9] <- "EcoRegionECS"
names(fish_clean)[13] <- "Fish_Group"
names(fish_clean)[19] <- "Length_Group"
names(fish_clean)[34] <- "all_sample_numbers"
names(fish_clean)[36] <- "TOTAL_ALK_ppm"
names(fish_clean)[37] <- "Lake_Size_A"
names(fish_clean)[39] <- "percent_littoral_area"
names(fish_clean)[40] <- "SD_m"
names(fish_clean)[49] <- "over_0.22"
names(fish_clean)[50] <- "Lake_Class"
names(fish_clean)[51] <- "LC_Group"

#CTY should be a factor
fish_clean$CTY <- as.factor(fish_clean$CTY)

fish_clean$Lake_Size_A <- as.numeric(fish_clean$Lake_Size_A)
fish_clean$MaxDepthft <- as.numeric(fish_clean$MaxDepthft)
fish_clean$SD_m <- as.numeric(fish_clean$SD_m)
fish_clean$TP <- as.numeric(fish_clean$TP)
fish_clean$WatershedA <- as.numeric(fish_clean$WatershedA)
fish_clean$CatchmentA <- as.numeric(fish_clean$CatchmentA)
fish_clean$WwetlandA <- as.numeric(fish_clean$WwetlandA)
fish_clean$CwetlandA <- as.numeric(fish_clean$CwetlandA)
fish_clean$WconiferA <- as.numeric(fish_clean$WconiferA)
fish_clean$CconiferA <- as.numeric(fish_clean$CconiferA)
fish_clean$Lake_Class <- as.factor(fish_clean$Lake_Class)
?save
save(fish_clean, file = "~/Documents/Hg_fish/fish_clean.csv")
