## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for growth
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(splines)
library(DHARMa)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")
plot_info <- read.csv("Data/plot_information.csv")


#Data formatting---------------------------------------------------
#adding plot information
plot_info <- plot_info|>
  select(!c(SUBPLOT, AREA, PLTAREA, SDLAREA, TMT))

#Converting census start date column to Date format
pines_long <- pines_long|>
  mutate(START_DATE = as.Date(START_DATE))

#Converting area and site area to factors
pines_long$AREA <- factor(pines_long$AREA, 
                          levels = c("DW", "SCC", "RP"))
pines_long$SITEAREA <- factor(pines_long$SITEAREA, 
                         levels = c("DW1", "DW2A", "DW2B",
                                    "DW2C", "DW3", "SCC5", "RP1", 
                                    "RP2"))

#Adding column for time since fires
fire_date <- as.Date("06/01/95",format= "%m/%d/%y")
pines_long <- pines_long|>
  mutate(TimeSinceFire = 
           time_length(difftime(START_DATE, fire_date),"years"))

#Creating dataframe for growth analysis
growth_dat <- pines_long|>
  group_by(IND_ID)|>
  #keeping only those censuses where height was measured
  filter(CENSUS_NUM %in% c(3, 5, 7, 9, 11:15))|>
  arrange(IND_ID, CENSUS_NUM)|>
  #pairing each census with values in the next census
  mutate(CENSUS_next = lead(CENSUS_NUM),
        STAT_next = lead(STAT),
         HT_next = lead(HT),
         Date_next = lead(START_DATE))|>
  #keeping only individuals that were alive in both pairs of censuses
  filter(STAT ==1 & STAT_next==1)|>
  filter(!is.na(HT) & !is.na(HT_next))|>
  #calculating growth increment
  mutate(deltaHT = HT_next - HT,   
    CensusInterval =year(Date_next)-year(START_DATE))|>
  #calculating annual growth increment &
  #and linearly interpolating height to annual timescales
  mutate(deltaHT_annual = deltaHT/CensusInterval)|>
  mutate(HT_next_annual = HT+deltaHT_annual)|>
  #converting plot code and individual id to factor
  mutate(IND_ID = factor(IND_ID),
      PLOTCODE = factor(PLOTCODE))

#checking whether consecutive censuses got correctly paired
#and whether census interval was calculated correctly
unique(growth_dat[,c("CENSUS_NUM", "CENSUS_next",
       "CensusInterval", "START_DATE",
        "Date_next")]) #looks correct

#Fitting models-----------------------------------
#Baseline model - random effects only
m_baseline <- lmer(deltaHT_annual ~ 
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height only model - linear effect of height
m_ht <- lmer(deltaHT_annual ~ HT+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height & AREA
m_ht_area <- lmer(deltaHT_annual ~ HT+AREA+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height & AREA - ht x area interaction
m_ht_area_htarea <- lmer(deltaHT_annual ~ HT*AREA+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire
m_ht_area_time <- lmer(deltaHT_annual ~ 
  HT+AREA+TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire, height x area interaction
m_ht_area_htarea_time <- lmer(deltaHT_annual ~ HT*AREA+
  TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire, linear + log height
m_htLL_area_time <- lmer(deltaHT_annual ~ 
  HT+log(HT)+AREA+TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire, height x area interaction, linear + log HT
m_htLL_area_htarea_time <- lmer(deltaHT_annual ~ 
  HT+log(HT)+AREA+HT:AREA + log(HT)*AREA+
  TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA, height x area interaction, linear + log HT
m_htLL_area_htarea <- lmer(deltaHT_annual ~ 
  HT+log(HT)+AREA+HT:AREA + log(HT)*AREA+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire, quadratic height
m_ht2_area_time <- lmer(deltaHT_annual ~ 
  poly(HT, 2, raw=T)+AREA+TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA & time since fire, height x area interaction, quadratic ht
m_ht2_area_htarea_time <- lmer(deltaHT_annual ~ 
  poly(HT, 2, raw=T)*AREA+TimeSinceFire+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#Height, AREA,  height x area interaction, quadratic HT
m_ht2_area_htarea <- lmer(deltaHT_annual ~ 
  poly(HT, 2, raw=T)*AREA+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat)

#model diagnostics-------------------------------------------------
#list of models
models <- list(
  m_baseline = m_baseline,
  m_ht = m_ht,
  m_ht_area = m_ht_area,
  m_ht_area_htarea = m_ht_area_htarea,
  m_ht_area_time = m_ht_area_time,
  m_ht_area_htarea_time = m_ht_area_htarea_time,
  m_ht2_area_time = m_ht2_area_time,
  m_ht2_area_htarea = m_ht2_area_htarea, 
  m_ht2_area_htarea_time = m_ht2_area_htarea_time,
  m_htLL_area_time = m_htLL_area_time,
  m_htLL_area_htarea =m_htLL_area_htarea, 
  m_htLL_area_htarea_time = m_htLL_area_htarea_time)

#checking for singular fit issues, convergence issues,
#  overdispersion and across-plots variation
model_diagnostics <- tibble::tibble(
  Name = names(models),
  Singular = sapply(models, isSingular),
  ConvergenceWarnings = sapply(models,
     \(x) x@optinfo$conv$lme4$messages),
  DispersionRatio = sapply(models, \(x) 
      check_overdispersion(x)$dispersion_ratio),
  AcrossPlotVar = sapply(models, \(x) 
          VarCorr(x)$PLOTCODE|>as.numeric()),
AcrossIndVar = sapply(models, \(x) 
          VarCorr(x)$IND_ID|>as.numeric()))
model_diagnostics

#QQ-plots
for(m in names(models)){
  model = models[[m]]
  qqnorm(residuals(model), main = m)
  qqline(residuals(model))
  dev.new()
}

#plotting DHARMa residuals
res <- lapply(models, simulateResiduals)
for (m in names(res)) {
  plot(res[[m]], title = m, quantreg = T)
  dev.new()
}
for (m in names(res)) {
  hist(res[[m]], main = m)
  dev.new()
}

#checking the shape of height vs height increment curves
newdat <- expand.grid( #generating new data for plotting
  HT = seq(
    min(growth_dat$HT, na.rm = TRUE),
    max(growth_dat$HT, na.rm = TRUE),
    length.out = 1000),
  AREA = levels(growth_dat$AREA))

#Plotting shape of curves
newdat$PredLinear <- predict(m_ht_area_htarea,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$PredLL <- predict(m_htLL_area_htarea,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$PredQuad <- predict(m_ht2_area_htarea,
          newdata = newdat, type = "response", 
        re.form = NA)

ggplot(newdat, aes(HT, PredLinear, col=AREA))+
  geom_line()+theme_bw()+xlab("Height")+
  ylab("Predicted Growth Increment")+
  ggtitle("Linear effect of height")
ggplot(newdat, aes(HT, PredLL, col=AREA))+
  geom_line()+theme_bw()+xlab("Height")+
  ylab("Predicted Growth Increment")+
  ggtitle("Linear + log effect of height")
ggplot(newdat, aes(HT, PredQuad, col=AREA))+
  geom_line()+theme_bw()+xlab("Height")+
  ylab("Predicted Growth Increment")+
  ggtitle("quadratic effect of height")

ggplot(growth_dat, aes(HT, deltaHT_annual, fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
   geom_line(data = newdat, aes(HT, PredQuad))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Quadratic")

ggplot(growth_dat, aes(HT, deltaHT_annual, fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
   geom_line(data = newdat, aes(HT, PredLinear))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Linear")

ggplot(growth_dat, aes(HT, deltaHT_annual, fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
   geom_line(data = newdat, aes(HT, PredLL))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Linear+Log")
#Model comparisons on full dataset-----------------------------------------
# Comparing models based on their fit to training data        
#comparing AICs
model_perform <- tibble::tibble(
  Name = names(models),
  AIC = sapply(models, AIC))
model_perform|>arrange(AIC)

#comparing mean squared error (MSE)
mse <- function(pred, obs){
  mean((obs-pred)^2)
}

model_perform$MSE <- sapply(models,
   function(x){mse(predict(x, re.form = NA),growth_dat$deltaHT_annual )})
model_perform|>
  mutate(MSE_by_var = MSE/var(growth_dat$deltaHT_annual))|>
  arrange(MSE)
