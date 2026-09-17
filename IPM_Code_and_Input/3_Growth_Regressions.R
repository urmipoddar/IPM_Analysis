## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for growth
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(glmmTMB)
library(gamlss)
library(performance)
library(splines)
library(DHARMa)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")
plot_info <- read.csv("Data/plot_information.csv")

#Functions--------------------------------------------------
#Mean Squared Error
mse <- function(pred, obs){
  mean((obs-pred)^2)
}

#Mean absolute error
mae <- function(pred, obs){
  mean(abs(obs-pred))
}
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
  #removing unnessary columns
  select(AREA, SITEAREA, PLOTCODE, 
    CENSUS_NUM ,IND_ID, START_DATE, 
    TimeSinceFire, STAT, HT)|>
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
    CensusInterval = 
      time_length(difftime(Date_next, START_DATE), "years"))|>
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

#Fitting models 1: Linear models-----------------------------------
# Based on prelim analysis, all models are fitted on log(HT_next)
# With log height as the predictor

#Baseline model - random effects & census interval only
mGrwt_baseline <- glmmTMB(log(HT_next) ~ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID),
  family = t_family(), 
  data = growth_dat)

#Height only model - log linear effect of height
mGrwt_ht <- glmmTMB(log(HT_next) ~ log(HT)+ CensusInterval+
   (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
  family = t_family(),
  dispformula = ~ log(HT))

#Height & AREA
mGrwt_ht_area <- glmmTMB(log(HT_next) ~ 
  log(HT)+AREA+ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height + AREA with quadratic height
mGrwt_ht2_area <- glmmTMB(log(HT_next) ~ 
  poly(log(HT),2)+AREA+ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height & AREA - ht x area interaction
mGrwt_ht_area_htarea <- glmmTMB(log(HT_next) ~ 
  log(HT)*AREA+ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height & AREA - ht x area interaction, quadratic height
mGrwt_ht2_area_htarea <- glmmTMB(log(HT_next) ~ 
  poly(log(HT),2)*AREA+ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), 
  data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height, AREA & time since fire
mGrwt_ht_area_time <- glmmTMB(log(HT_next) ~ 
  log(HT)+AREA + TimeSinceFire + CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height, AREA & time since fire, quadratic height
# mGrwt_ht2_area_time <- glmmTMB(log(HT_next) ~ 
#   poly(log(HT),2)+AREA + TimeSinceFire + CensusInterval+
#   (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
#     family = t_family(),
#   dispformula = ~ log(HT))

#Height, AREA & time since fire, height x area interaction
mGrwt_ht_area_htarea_time <- glmmTMB(log(HT_next) ~ 
  log(HT)*AREA + TimeSinceFire + CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

#Height, AREA & time since fire, height x area interaction, quadratic ht
mGrwt_ht2_area_htarea_time <- glmmTMB(log(HT_next) ~ 
  poly(log(HT),2)*AREA + TimeSinceFire + CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

  mGrwt_ht2_area_htarea_time_normalres <- glmmTMB(log(HT_next) ~ 
  poly(log(HT),2)*AREA + TimeSinceFire + CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
  dispformula = ~ log(HT))

m_skewt <- gamlss(
  log(HT_next) ~ poly(log(HT), 2) * AREA + 
    TimeSinceFire + CensusInterval  + 
    random(PLOTCODE) + random(IND_ID),
  sigma.formula = ~ log(HT),
  nu.formula = ~ 1,     
  family = ST3,
  data = growth_dat)

m_skewt3 <- gamlss(
  log(HT_next) ~ poly(log(HT), 2) * AREA + TimeSinceFire + CensusInterval +
    random(SITEAREA) + random(PLOTCODE) + random(IND_ID),
  sigma.formula = ~ log(HT) + I(log(HT)^2),
  nu.formula = ~ log(HT),
  tau.formula = ~ log(HT),
  family = ST3,
  data = growth_dat)

#model diagnostics-------------------------------------------------
#list of models
models <- list(
  mGrwt_baseline = mGrwt_baseline,
  mGrwt_ht = mGrwt_ht,
  mGrwt_ht_area = mGrwt_ht_area,
  Grwt_ht2_area = mGrwt_ht2_area,
  mGrwt_ht_area_htarea = mGrwt_ht_area_htarea,
  Grwt_ht2_area_htarea = mGrwt_ht2_area_htarea,
  mGrwt_ht_area_time = mGrwt_ht_area_time,
  #mGrwt_ht2_area_time = mGrwt_ht2_area_time,
  mGrwt_ht_area_htarea_time = mGrwt_ht_area_htarea_time,
  mGrwt_ht2_area_htarea_time = mGrwt_ht2_area_htarea_time,
  #m_skewt = m_skewt,
  mGrwt_ht2_area_htarea_time_normalres = 
  mGrwt_ht2_area_htarea_time_normalres)

#checking for singular fit issues, convergence issues,
#  overdispersion and across-plots variation
model_diagnostics <- tibble::tibble(
  Name = names(models),
  Singular = check_singularity(models),
  ConvergenceWarnings = sapply(models,
     \(x) x$fit$convergence),
  DispersionRatio = sapply(models, \(x) 
      check_overdispersion(x)$dispersion_ratio),
  AcrossPlotVar = sapply(models, \(x) 
          VarCorr(x)$cond$PLOTCODE|>as.numeric()|>sqrt()),
AcrossIndVar = sapply(models, \(x) 
          VarCorr(x)$cond$IND_ID|>as.numeric()|>sqrt()))
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

#checking shapes of predicted curves----------------------------
#Building newdat grid with ALL predictors
newdat <- expand.grid(
  HT = seq(
    min(growth_dat$HT, na.rm = TRUE),
    max(growth_dat$HT, na.rm = TRUE),
    length.out = 1000),
  AREA = levels(growth_dat$AREA),
  CensusInterval = 1,
  TimeSinceFire = mean(growth_dat$TimeSinceFire, na.rm = TRUE))

# Predictions on the log scale
newdat$PredLinear <- predict(
  mGrwt_ht_area_htarea,
  newdata = newdat,
  type = "response",
  re.form = NA)
newdat$PredQuad <- predict(
  mGrwt_ht2_area_htarea,
  newdata = newdat,
  type = "response",
  re.form = NA)

#Predictions on the raw scale - linear model
mu_marginal_linear <- predict(
  mGrwt_ht_area_htarea, 
  newdata = growth_dat,
  type = "response", 
  re.form = NA)
total_resid_linear <- log(growth_dat$HT_next) - mu_marginal_linear   # captures RE variance + residual + any skew
smear_factor <- mean(exp(total_resid_linear))
newdat$PredLinear_corrected <- 
      exp(newdat$PredLinear)*smear_factor

# Compute Duan's smearing factor from the QUADRATIC model's
# own marginal residuals (not the linear model's)
mu_marginal_quad <- predict(
  mGrwt_ht2_area_htarea,
  newdata = growth_dat,
  type = "response",
  re.form = NA)
total_resid_quad <- log(growth_dat$HT_next) - mu_marginal_quad
smear_factor_quad <- mean(exp(total_resid_quad))
newdat$PredQuad_corrected <- exp(newdat$PredQuad) * smear_factor_quad


#Log-scale plots
ggplot(growth_dat, aes(log(HT), log(HT_next), fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
   geom_line(data = newdat, aes(log(HT), (PredLinear)))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Linear")
ggplot(growth_dat, aes(log(HT), log(HT_next), fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
   geom_line(data = newdat, aes(log(HT), (PredQuad)))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Quadratic")

#Raw scale plots
ggplot(growth_dat, aes((HT), (HT_next), fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
  geom_line(data = newdat, aes((HT), 
     PredLinear_corrected))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Linear")

ggplot(growth_dat, aes((HT), (HT_next), fill = AREA))+
  geom_point(pch = 21, col = "white", alpha = 0.6)+
  geom_line(data = newdat, aes((HT), 
     PredQuad_corrected))+
  theme_bw()+facet_wrap(vars(AREA))+ggtitle("Quadratic")
#Model comparisons on full dataset-----------------------------------------
# Comparing models based on their fit to training data        
#comparing AICs
model_perform <- tibble::tibble(
  Name = names(models),
  AIC = sapply(models, AIC))
model_perform|>arrange(AIC)

model_perform$MSE <- sapply(models,
   function(x){mse(predict(x, re.form = NA),
     log(growth_dat$HT_next))})
model_perform|>
  mutate(MSE_by_var = 
    1-MSE/var(log(growth_dat$HT_next)))|>
  arrange(MSE)
