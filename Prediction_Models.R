################################################################################
#Using, fish_clean, created in Hg_First_Upload, 
#now creating a walleye and a pike dataset
################################################################################

WE_only <- fish_clean %>% filter(SPEC == "WE")
NP_only <- fish_clean %>% filter(SPEC == "NP")
WE_NP_only <- rbind(WE_only, NP_only)

################################################################################
#Straightforward predictions Walleye Only first
################################################################################
#For starts making assumptions that HGPPM is stationary &
#all predictors have linear relationship with HGPPM (this will be easy to relax)

#First regression with all potentially helpful variables
full_WE_df <- WE_only %>% select(CTY, EcoRegionIII, LGTHIN,
                                 HGPPM, TOTAL_ALK_ppm, ANAT, 
                                 Lake_Size_A, MaxDepthft, percent_littoral_area, 
                                 SD_m, TP, WatershedA, CatchmentA, WwetlandA, 
                                 CwetlandA, WconiferA, CconiferA,
                                 NOFISH)

WE_no_NA <- full_WE_df[complete.cases(full_WE_df), ]

#I found there is only one walley sample from DA with 
#complete data. I choose to delete this sample. Same for CTY 55. 
#Without deleting will cause error because cannot predict for that case
WE_no_NA <- WE_no_NA %>% filter(EcoRegionIII != 'DA' | CTY != '55')

#Below performs LOOCV with large linear regression 
err <- c()
Sys.time()
for (i in 1:(nrow(WE_no_NA) - 1)){
  train_ind <-  intersect(c(0:(i-1), (i+1):nrow(WE_no_NA)), 1:nrow(WE_no_NA)) 
  
  train <- WE_no_NA[train_ind, ]
  test <- WE_no_NA[-train_ind, ]
  
  reg <- lm(HGPPM~.,data = train)
  err[i] <- predict(reg, test) - test$HGPPM
}

lm_Err <- err

#LOOCV Rmse is .2238 and mae is .152557 


#next step, ranger rf

for (i in 1:(nrow(WE_no_NA) - 1)){
  train_ind <-  intersect(c(0:(i-1), (i+1):nrow(WE_no_NA)), 1:nrow(WE_no_NA)) 
  
  train <- WE_no_NA[train_ind, ]
  test <- WE_no_NA[-train_ind, ]
  
  reg <- ranger(HGPPM~.,data = train, importance = 'permutation')
  err[i] <- predict(reg, test)$predictions - test$HGPPM
  var_importance[[i]] <- importance(reg)
}
rf_err <- err

plot(var_importance[[1]])

write_csv(as.data.frame(lm_err), path = "WE_lm_err")
write_csv(as.data.frame(rf_err), path = "WE_rf_err")
save(var_importance, file = "var_imp_rf_WE.rdata")

################################################################################
#investigate errors
################################################################################


ggplot(WE_rf_err, aes(x = rf_err)) + geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(-1, 5))
ggplot(WE_lm_err, aes(x = lm_err)) + geom_histogram(bins = 100) 

