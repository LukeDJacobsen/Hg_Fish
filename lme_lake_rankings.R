################################################################################
#Data import and dependent packages
################################################################################

library(tidyverse)
library(readr)
library(lme4)
library(xtable)
library(lmerTest)
fish_clean <- read_csv("Documents/Hg_Fish/fish_clean.csv")

################################################################################
#Data formatting
################################################################################

#Getting rid of composite samples because they will bias
#variance estimates

no_comp <- fish_clean %>% filter(composite == "N")

#make a species variable that's fewer levels
library(forcats)
table(no_comp$SPEC)
no_comp$SPEC_clumped <- fct_lump(no_comp$SPEC, prop = .0014)

library(skimr)
skim(no_comp)

################################################################################
#Model stage one: stage one--computing lake random effects
################################################################################

no_comp$Year_F <- as.factor(no_comp$Year)
#just starting with straight forward model.
lme1 <- lmer(LN_HGPPM ~  SPEC_clumped*LGTHIN + 
       (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme1)
BIC(lme1)
lake_ranef <- ranef(lme1)

lme2 <- lmer(LN_HGPPM ~  SPEC_clumped + LGTHIN + 
               (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme2)
BIC(lme2)
summary(lme2)

lme3 <- lmer(LN_HGPPM ~  SPEC_clumped * LGTHIN + EcoRegionIII + 
               (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme3)
BIC(lme3)

lme4 <- lmer(LN_HGPPM ~  SPEC_clumped * LGTHIN * EcoRegionIII + 
               (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme4)
BIC(lme4)
#over parameterized

lme5 <- lmer(LN_HGPPM ~  SPEC_clumped * LGTHIN + EcoRegionIII*SPEC_clumped + 
               (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme5)
BIC(lme5)
#over parameterized

lme6 <- lmer(LN_HGPPM ~  SPEC_clumped * LGTHIN + EcoRegionIII*LGTHIN + 
               (1|Year_F) + (1|DOWID), data = no_comp)
AIC(lme6)
BIC(lme6)

AIC(lme1,lme2,lme3,lme4,lme5,lme6)
BIC(lme1,lme2,lme3,lme4,lme5,lme6)

#34821 obs
#having the interaction and main ecoRegion effect is helpful.


skim(no_comp %>% filter(complete.cases(SPEC_clumped) & complete.cases(LGTHIN) + 
                          complete.cases(Year_F) & complete.cases(DOWID)))
#adding CTY loses 1717 rows. Adding lake attributes loses about 10,000 cases
#for following models have to re-factor clump because when losing the ~10,000 cases
#we get many species with small samples
no_comp_reduced <- no_comp %>% filter(complete.cases(CatchmentA) & complete.cases(CwetlandA)
                                      & complete.cases(Lake_Size_A) & complete.cases(MaxDepthft) &
                                        complete.cases(percent_littoral_area) & complete.cases(TOTAL_ALK_ppm) &
                                        complete.cases(TP) & complete.cases(WatershedA) & 
                                        complete.cases(WconiferA) & complete.cases(WwetlandA) &
                                        complete.cases(LGTHIN) & complete.cases(EcoRegionIII) & 
                                        complete.cases(MaxDepthft) & complete.cases(Year_F) &
                                        complete.cases(DOWID))
dim(no_comp_reduced)
no_comp_reduced$SPEC_clumpedV2 <- fct_lump(no_comp_reduced$SPEC_clumped, prop = .001)
table(no_comp_reduced$SPEC_clumpedV2)

#note that lat and long still have 5000 missing
lme7 <- lmer(LN_HGPPM ~ SPEC_clumpedV2 * LGTHIN + Lake_Size_A + EcoRegionIII+
               MaxDepthft + percent_littoral_area + TOTAL_ALK_ppm + TP + WatershedA + 
               WconiferA + WwetlandA + (1|Year_F) + (1|DOWID), data = no_comp_reduced)

AIC(lme7)
BIC(lme7)
summary(lme7)
#watershed, wetland, and conifer have small t values so removed for lme4

lme8 <- lmer(LN_HGPPM ~ SPEC_clumpedV2 * LGTHIN + Lake_Size_A + EcoRegionIII+
               MaxDepthft + percent_littoral_area + TOTAL_ALK_ppm + TP + 
               (1|Year_F) + (1|DOWID), data = no_comp_reduced)
AIC(lme8)
BIC(lme8)
summary(lme8)
#AIC, BIC better

#now removing each variable with lowest t value one at a time. 

lme9 <- update(lme8, .~. - Lake_Size_A)
AIC(lme9)
BIC(lme9)
summary(lme9)

lme10 <- update(lme9, .~. - MaxDepthft)
AIC(lme10)
BIC(lme10)
summary(lme10)

lme11 <- update(lme10, .~. -percent_littoral_area)
AIC(lme11)
BIC(lme11)

lme12 <- update(lme11, .~. -TP)
summary(lme12)

lme13 <- update(lme12, .~. -TOTAL_ALK_ppm)

lme13

AIC(lme7, lme8,lme9,lme10,lme11, lme12, lme13)
BIC(lme7, lme8,lme9,lme10,lme11, lme12, lme13)


#24785 obs. 

#Conclusion: percent littoral area, total alk ppm, and total phosphorus improve model, however 
#using these variables results in a loss of a 11,039 cases (out of 36,310 total), so probably not 
#worth including in model (Is this correct Gary).


################################################################################
#Organize data to look at lake random effects. Using lme1 random effects
################################################################################
#goal: df with column for renef estimate, dowid, lat, and long

#making df with mean lat and long by lake.
fish_clean$lat <- ifelse(fish_clean$LAT == 0, NA, fish_clean$LAT)
fish_clean$long <- ifelse(fish_clean$LONG == 0, NA, fish_clean$LONG)


lake_extra <- fish_clean %>% filter(composite == "N") %>% select(DOWID, lat, long) %>%
  group_by(DOWID) %>% summarise(mean_LAT = mean(lat, na.rm =T), mean_LONG = mean(long, na.rm = T))

#get rid of duplicated rows
lake <- lake_extra[!duplicated(lake_extra),]

#awkward data munging stuff
a <- as.data.frame(lake_ranef[[1]])
names(a) <- c("Random_Lake_Effect")
DOWID <- rownames(a)
a$DOWID <- DOWID
a <- as_tibble(a)

lake_level_eff <- inner_join(a, lake, by = "DOWID")

#appears to be a data issue with lake 03003000. Changed to match what it should be
lake_level_eff$mean_LONG <- ifelse(lake_level_eff$DOWID == "03003000", 952649, lake_level_eff$mean_LONG)
lake_level_eff <- lake_level_eff %>% mutate(Latitude = mean_LAT/10000, Longitude = -mean_LONG/10000)

write.csv(lake_level_eff, file = 'documents/GitHub/Hg_Fish/lake_level_eff.csv')

################################################################################
#Spatial modeling with random lake effect as response.
################################################################################


#models
names(lake_level_eff)

mod_loc1 <- lm(Random_Lake_Effect ~ Latitude + Longitude, data = lake_level_eff)
AIC(mod_loc1)
BIC(mod_loc1)

mod_loc2 <- update(mod_loc1, .~. + I(Latitude^2) + I(Longitude^2))
AIC(mod_loc2)
BIC(mod_loc2)

mod_loc2.1 <- update(mod_loc2, .~. + Latitude:Longitude)
AIC(mod_loc2.1)
BIC(mod_loc2.1)


mod_loc3 <- update(mod_loc2, .~. + I(Latitude^3) + I(Longitude^3))
AIC(mod_loc3)
BIC(mod_loc3)
summary(mod_loc3)

mod_loc4 <- update(mod_loc3, .~. + I(Latitude^4) + I(Longitude^4))
AIC(mod_loc4)
BIC(mod_loc4)
summary(mod_loc4)

anova(mod_loc3, mod_loc4)

mod_loc5 <- update(mod_loc4, .~. + I(Latitude ^ 5) + I(Longitude ^ 5))
summary(mod_loc5)
mod_loc5$df.residual
#what is goin on? why doesn't this


##########################################################
#Plotting the distance by change in effect between two lakes
##########################################################
#distance between lake i and j matrix
distance <- matrix(nrow = nrow(lake_level_eff), ncol = nrow(lake_level_eff))
for (i in 1:nrow(lake_level_eff)){
  for (j in 1:nrow(lake_level_eff)){
    distance[i,j] <- sqrt((lake_level_eff$mean_LAT[i] - lake_level_eff$mean_LAT[j])^2 +
                            (lake_level_eff$mean_LONG[i] - lake_level_eff$mean_LONG[j])^2)
  }
}

#change in random effect estimate between lake i and j
change_eff <- matrix(nrow = nrow(lake_level_eff), ncol = nrow(lake_level_eff))
for (i in 1:nrow(lake_level_eff)){
  for (j in 1:nrow(lake_level_eff)){
    change_eff[i,j] <- abs(lake_level_eff$Random_Lake_Effect[i] - lake_level_eff$Random_Lake_Effect[j])
  }
}
plot(distance[1201,], change_eff[1201,])
distance <- ifelse(distance == 0, NA, distance)
dist_change <- tibble(as.vector(distance), as.vector(change_eff))
names(dist_change) <- c("distance", "change_effect")
dist_effect <- ggplot(dist_change, aes(distance, change_effect)) + geom_smooth()
dist_effect + coord_cartesian(xlim = c(0,50000))

################################################################################
#Extracting data on the lowest and highest random effect lakes
################################################################################

sort(lake_level_eff$Random_Lake_Effect)

good_lakes <- lake_level_eff %>% filter(Random_Lake_Effect < -1.956)
good_lakes_fish <- fish_clean %>% filter(DOWID %in% good_lakes$DOWID)
good_lakes_fish
for (i in 1:length(good_lakes$DOWID)){
  D <- good_lakes$DOWID[i]
  print(good_lakes_fish %>% filter(DOWID == D) %>% group_by(SPEC) %>% 
    summarise(mean_HG = mean(HGPPM), mean_length = mean(LGTHIN), count = n()))
}


bad_lakes <- lake_level_eff %>% filter(Random_Lake_Effect > 1.5)
bad_lakes_fish <- fish_clean %>% filter(DOWID %in% bad_lakes$DOWID)
bad_lakes_fish
tab_for_paper <- data.frame()
Lake <- c(bad_lakes$DOWID[1])
for (i in 1:length(bad_lakes$DOWID)){
  D <- bad_lakes$DOWID[i]
  a <- bad_lakes_fish %>% filter(DOWID == D) %>% group_by(SPEC) %>% 
          summarise("Mean HG (ppm)" = mean(HGPPM), 
                    "Mean Length (in)" = mean(LGTHIN), Count = n())
  tab_for_paper <- rbind(tab_for_paper, a)
  l_now <- nrow(tab_for_paper)
  Lake[l_now + 1] <- bad_lakes$DOWID[i+1] 
}
tab_for_paper
Lake <- Lake[1:22]
tab <- cbind(Lake, tab_for_paper)
names(tab)[1] <- "DOWID" 
names(tab)[2] <- "Species" 
names(tab)[3] <- "Mean Mercury (ppm)"
xtable(tab)


good_lakes <- lake_level_eff %>% filter(Random_Lake_Effect < -1.956)
good_lakes_fish <- fish_clean %>% filter(DOWID %in% good_lakes$DOWID)
tab_for_paper <- data.frame()
Lake <- c(good_lakes$DOWID[1])
for (i in 1:length(good_lakes$DOWID)){
  D <- good_lakes$DOWID[i]
  a <- good_lakes_fish %>% filter(DOWID == D) %>% group_by(SPEC) %>% 
    summarise("Mean HG (ppm)" = mean(HGPPM), 
              "Mean Length (in)" = mean(LGTHIN), Count = n())
  tab_for_paper <- rbind(tab_for_paper, a)
  l_now <- nrow(tab_for_paper)
  Lake[l_now + 1] <- good_lakes$DOWID[i+1] 
}
dim(tab_for_paper)
Lake
Lake <- Lake[1:46]
tab <- cbind(Lake, tab_for_paper)
names(tab)[1] <- "DOWID" 
names(tab)[2] <- "Species" 
names(tab)[3] <- "Mean Mercury (ppm)"
xtable(tab)


