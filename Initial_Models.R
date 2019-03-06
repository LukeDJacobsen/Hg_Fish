#Importing fish datainstall.packages('forecast')
library(forecast)
library(readr)
library(nlme)
library(lme4)
library(tidyverse)

fish_clean <- read_csv("Documents/Hg_Fish/fish_clean.csv")
table(fish_clean$NOFISH)[1]/nrow(fish_clean) 
#80 % of rows are 1 fish

###########################################
#Models
###########################################
basic_mod <- lm(HGPPM ~ LGTHIN, data = fish_clean)
plot(basic_mod)
boxCox(basic_mod)
BoxCox.lambda(basic_mod, lower = 0, upper = 1)

basic_mod_ln <- lm(LN_HGPPM ~ LGTHIN, data = fish_clean)
plot(basic_mod_ln)

#not great fit, ln helps though, try by species
spec_mod <- lm(LN_HGPPM ~ LGTHIN +SPEC, data = fish_clean)
plot(spec_mod)

spec_mod_ln <- lm(LN_HGPPM ~ log(LGTHIN) +SPEC, data = fish_clean)
plot(spec_mod_ln)
#logging length of fish doesn't help. 

spec_reg_lm <- lm(LN_HGPPM ~ LGTHIN + SPEC + EcoRegionIII, data = fish_clean)
plot(spec_reg_lm)


#Some improvement when taking odd transformation
fish_clean$HGPPM25 <- fish_clean$HGPPM^(.25)
mod_25 <- lm(HGPPM25 ~ LGTHIN, data = fish_clean)
plot(mod_25)
boxCox(mod_25)



#looks good except stuff at lower tail.
#lets see if we can see what's going on
hist(fish_clean$LN_HGPPM, breaks = 100)
#does fitting gamma distribution make any sense?
hist(fish_clean$HGPPM, breaks = 1000)
ggplot(fish_clean) + geom_histogram(aes(HGPPM), bins = 4000) +
  coord_cartesian(x = c(0, .1))
#plot above indicates some, but not all labs rounded to 
#nearest hundreth place.

#####################################################
#Does anything change if we only look at the last 
#ten years
#####################################################

fish_clean_2008on <- fish_clean %>% filter(Year >= 2008)
spec_mod08 <- lm(LN_HGPPM ~ LGTHIN +SPEC, data = fish_clean_2008on)
plot(spec_mod08)

#not surprising, not better.

#####################################################
#Lets look at individual species
#####################################################
fish_clean_WE <- fish_clean %>% filter(SPEC == 'WE')
WE_mod <- lm(LN_HGPPM ~ LGTHIN, data = fish_clean_WE)
plot(WE_mod)

fish_clean_NP <- fish_clean %>% filter(SPEC == 'NP')
NP_mod <- lm(LN_HGPPM ~ LGTHIN, data = fish_clean_NP)
plot(NP_mod)

##not helpful

#some mixed effects models
mm1 <- lme(LN_HGPPM ~ LGTHIN, data = fish_clean,
           na.action = na.omit ,random = ~1|Year)
summary(mm1)
qqnorm(mm1)


mm2 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = fish_clean, na.action = na.omit ,random = ~1|Year)
summary(mm2)
qqnorm(mm2)
plot(mm2)
#this may follow assumptions the best, but still not great

mm3 <- lme(LN_HGPPM ~ LGTHIN, data = fish_clean, na.action = na.omit ,random = ~1|DOWID)
summary(mm3)
plot(mm3)
qqnorm(mm3)

mm4 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = fish_clean, na.action = na.omit, random = ~1|DOWID)
summary(mm4)
plot(mm4)
qqnorm(mm4)

###########################################################################
#individual lake models
###########################################################################
which(table(fish_clean$DOWID) > 100)

#elk lake
elk_lake <- fish_clean %>% filter(DOWID == 15001000)
elk_m <- lme(LN_HGPPM ~ LGTHIN + SPEC, data = elk_lake,  random = ~1|Year)
summary(elk_m)
plot(elk_m)
qqnorm(residuals(elk_m))

elk_mm6 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = elk_lake, na.action = na.omit, 
           random = ~1|Year)
plot(elk_mm6)
qqnorm(elk_mm6)

#
boulder <- fish_clean %>% filter(DOWID == 69037300)
boulder_m <- lme(LN_HGPPM ~ LGTHIN + SPEC, data = boulder,  random = ~1|Year)
summary(boulder_m)
plot(boulder_m)
qqnorm(boulder_m)

boulder_mm6 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = boulder, na.action = na.omit, 
               random = ~1|Year)
boulder_mm2 <- lm(LN_HGPPM ~ LGTHIN + Fish_Group, data = boulder)
plot(boulder_mm2)


# ################################################################################
# #looking at last ~10 years only
# ################################################################################
mm1_08 <- lme(LN_HGPPM ~ LGTHIN, data = fish_clean_2008on, na.action = na.omit, 
              random = ~1|Year)
qqnorm(mm1_08)

mm2_08 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = fish_clean_2008on, na.action = na.omit, 
           random = ~1|Year)
qqnorm(mm2_08)

mm3_08 <- lme(LN_HGPPM ~ LGTHIN , data = fish_clean_2008on, na.action = na.omit, 
              random = ~1|DOWID)
qqnorm(mm3_08)

mm4_08 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = fish_clean_2008on, na.action = na.omit, 
              random = ~1|DOWID)
qqnorm(mm4_08)
#Looking at last 10 years doesn't help

################################################################################
#Look by eco region
################################################################################
mm5 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group + EcoRegionIII, data = fish_clean, na.action = na.omit, 
              random = ~1|Year)
summary(mm5)
qqnorm(mm5)

mm6 <- lme(LN_HGPPM ~ LGTHIN + Fish_Group, data = fish_clean, na.action = na.omit, 
           random = ~1|EcoRegionIII)
plot(mm6)
qqnorm(mm6)

mm7 <- lme(LN_HGPPM ~ LGTHIN + EcoRegionIII, data = fish_clean, na.action = na.omit, 
           random = ~1|Year)
summary(mm7)
qqnorm(mm7)

#mm6 is a pretty good model. 