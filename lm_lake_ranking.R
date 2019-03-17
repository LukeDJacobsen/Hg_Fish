################################################################################
#Data import
################################################################################

library(tidyverse)
library(readr)
fish_clean <- read_csv("Documents/Hg_Fish/fish_clean.csv")

################################################################################
#Data formatting
################################################################################

#For consistency only use since 1990 and duplicate composite samples
#even though this adds variance bias
fish_clean_1990 <- fish_clean %>% filter(Year >= 1990)
dim(fish_clean_1990)
dim(fish_clean)
#lose ~6500 samples from '67-'89
rep_comp <- fish_clean_1990[rep(seq_len(nrow(fish_clean_1990)), 
                                      times=fish_clean_1990$NOFISH),]
dim(rep_comp)

#make a species variable that's fewer levels
library(forcats)
rep_comp$SPEC_clumped <- fct_lump(rep_comp$SPEC, prop = .0005)
table(rep_comp$SPEC_clumped)
library(skimr)
skim(rep_comp)

################################################################################
#model
################################################################################

#minimum reasonable model
fm1 <- lm(LN_HGPPM ~ DATEI + LGTHIN*SPEC_clumped + EcoRegionIII, data = rep_comp)
summary(fm1)

################################################################################
#Compute mean residual for a lake
################################################################################

lake_resid <- function(Lake_ID, model, dat){
  a <- dat %>% filter(DOWID == Lake_ID) %>% select(DATEI, SPEC, LGTHIN, EcoRegionIII, SPEC_clumped, HGPPM, LN_HGPPM)
  plus5_check <<- ifelse(sum(ifelse((a %>% group_by(SPEC) %>% summarise(count = length(SPEC)))$count >=5, 1, 0)) >0, 1, 0)
  resid <<-  a$HGPPM - exp(predict(model, a))
  resid_lns <<- a$LN_HGPPM - predict(model, a)
}

individual_lake_plus5_check <- list()
individual_lake_resid <- list()
individual_lake_resid_lns <- list()
for (i in 1:nrow(lake)){
  lake_resid(lake$DOWID[i], fm1, dat = rep_comp)
  individual_lake_resid[[i]] <- resid
  individual_lake_resid_lns[[i]] <- resid_lns
  individual_lake_plus5_check[[i]] <- plus5_check
}

#makes ordered dataframe with largest mean residuals at top
hg_compared_to_modeled <- data.frame(mean_resid = sapply(individual_lake_resid, mean, na.rm = T),
           mean_resid_ln = sapply(individual_lake_resid_lns, mean, na.rm = T),
           DOWID = lake$DOWID, sample_5plus = unlist(individual_lake_plus5_check))
ordered_df_ln <- hg_compared_to_modeled[order(hg_compared_to_modeled$mean_resid_ln, decreasing = T),]
ordered_df_ln[1:20,]
#checking out some of the high ranking lakes
rep_comp %>% filter(DOWID == '69087000') %>% select(HGPPM, LGTHIN, SPEC_clumped)

########################################################################
#other models to consider
########################################################################
#conversation to have with Gary
#example obviously not use only p-value
#how to know if over fitting, should I do some kind of CV to 
#choose a good model?  
#plots
fm9 <- lm(LN_HGPPM ~ poly(DATEI,9) + LGTHIN*SPEC_clumped +
            EcoRegionIII, data = rep_comp)




