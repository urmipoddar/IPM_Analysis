## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for cone production and serotiny
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")

#Data formatting---------------------------------------------------

#Converting area and site area to factors
pines_long$AREA <- factor(pines_long$AREA, 
                          levels = c("DW", "SCC", "RP"))
pines_long$SITEAREA <- factor(pines_long$SITEAREA, 
                         levels = c("DW1", "DW2A", "DW2B",
                                    "DW2C", "DW3", "SCC5", "RP1", 
                                    "RP2"))

#Data frame of cone counts
cones <- pines_long|>
  filter(STAT==1)|>
  filter(!is.na(CONE_Mature) & !is.na(HT))

#keeping data from census 15 only 
cone15 <- cones|>filter(CENSUS_NUM == 15)|>
  mutate(SEROTINY = ifelse(SEROTINY=="", NA, SEROTINY))|>
  #converting serotiny into a binary column
  mutate(SEROTINY_binary= case_when(SEROTINY=="O"~0,
                  SEROTINY == "S" ~1, 
                .default = NA))

#Modeling cone counts-----------------------------------------
#Only data from census 15 used because
# mature cones accumulate on trees, 
# so repeated censuses recount the same cones

#Fitting model:
# Log Height & AREA as predictors, 
# Negative Binomial Hurdle model
mCone_lht_nb <- glmmTMB(CONE_Mature ~
  log(HT)*AREA + (1|PLOTCODE),  
  ziformula = ~ log(HT)*AREA,
    family = truncated_nbinom2(), data = cone15)

#checking for convergence issues etc.
mCone_lht_nb$fit$convergence #no issues
summary(mCone_lht_nb)
VarCorr(mCone_lht_nb)

#checking residuals
hist(residuals(mCone_lht_nb), breaks = 50)
mf  <- model.frame(mCone_lht_nb) 
res <- simulateResiduals(mCone_lht_nb, plot = T)
plotResiduals(res, form = mf[["log(HT)"]])   # residuals vs height
plotResiduals(res, form = mf$AREA)           # residuals vs AREA
testDispersion(res)
testOutliers(res, type = "bootstrap")
hist(res) 
#no major issues

#Plotting shape of predicted curves
cone15$pred_cone <- predict(mCone_lht_nb, 
  newdata = cone15, type = "response", re.form = NA)

ggplot(cone15, aes(HT, CONE_Mature))+
  geom_point(aes(fill = AREA), pch = 21, alpha = 0.6)+
  theme_bw()+facet_wrap(vars(AREA))+
  geom_point( aes(HT, pred_cone), size =1)

#Modeling serotiny-----------------------------------------
mSer <- glm(SEROTINY_binary ~ AREA , 
  family = binomial, data = cone15)
mSer
summary(mSer)

#checking residuals
simulateResiduals(mSer, plot=T) #residuals look fine

#comparing predicted and observed proportions
cone15$pred_serotiny <- predict(mSer, 
    type = "response", newdata = cone15)
cone15|>
  filter(!is.na(SEROTINY_binary))|>
  group_by(AREA)|>
  summarise(NumSer = sum(SEROTINY_binary),
    NumObs=n(), PredictedProp = mean(pred_serotiny))|>
  mutate(ObservedProp = NumSer/NumObs) 
#predicted and observed proportions match
#Saving models---------------------------------
models <- list(mCone = mCone_lht_nb, mSer = mSer)
saveRDS(models,"Output/ConeProductionModels.rds")
