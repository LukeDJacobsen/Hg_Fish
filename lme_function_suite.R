################################################################################
#Data import
################################################################################

library(tidyverse)
library(readr)
library('usmap')
fish_clean <- read_csv("Documents/Hg_Fish/fish_clean.csv")
devtools::install_github("dkahle/ggmap", ref = "tidyup")
library(ggmap)
################################################################################
#Data formatting
################################################################################

#Getting rid of composite samples because they will bias
#variance estimates

no_comp <- fish_clean %>% filter(composite == "N")

#make a species variable that's fewer levels
library(forcats)
table(no_comp$SPEC)
no_comp$SPEC_clumped <- fct_lump(no_comp$SPEC, n = 50)
no_comp$Year_F <- as.factor(no_comp$Year)

lake_suite <- function(model = LN_HGPPM ~ 1, data = no_comp){
  lme1 <- lmer(model, data = no_comp)
  lake_ranef <- ranef(lme1)
  #making df with mean lat and long by lake.
  fish_clean$lat <- ifelse(fish_clean$LAT == 0, NA, fish_clean$LAT)
  fish_clean$long <- ifelse(fish_clean$LONG == 0, NA, fish_clean$LONG)
  
  #appears to be a data issue with lake 03003000. Changed to match what it should be
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
  lake_level_eff$mean_LONG <- ifelse(lake_level_eff$DOWID == "03003000", 952649, lake_level_eff$mean_LONG)
  
  return(list(ggplot(lake_level_eff, aes(-mean_LONG, mean_LAT)) + geom_point(aes(color = Random_Lake_Effect)) +
    scale_colour_gradient2(low = "dark green", high = "dark red") + ggtitle(model), 
  summary(lme1)))
}

lake_suite(model = LN_HGPPM ~  SPEC_clumped + SPEC_clumped:LGTHIN + 
             (1|Year_F) + (1|DOWID), data = no_comp)
names(no_comp)
lake_suite(model = LN_HGPPM ~ 
             (1|Year_F) + (1|DOWID), data = no_comp)
lake_suite(model = LN_HGPPM ~ EcoRegionIII +SPEC_clumped + SPEC_clumped:LGTHIN +
             (1|Year_F) + (1|DOWID), data = no_comp)


########################################################
#Playing around with US map
########################################################

plot_usmap( include = c("MN"), lines = "navy"
) 
ggplot()+ geom_point(data = lake_level_eff,aes(mean_LONG, mean_LAT, color = Random_Lake_Effect))) +
  scale_colour_gradient2(low = "dark green", high = "dark red")
ggplot(lake_level_eff, aes(-mean_LONG, mean_LAT)) + geom_point(aes(color = Random_Lake_Effect)) +
  scale_colour_gradient2(low = "dark green", high = "dark red")
?get_googlemap
ggmap(get_googlemap(center = c(lon = 93, lat = 46)))
?ggmap
qmap("Minnesota")
?register_google
register_google(key = "AIzaSyACYMmUg-c4QQfeGHzmuKpYEBE3DyDkiYI")
