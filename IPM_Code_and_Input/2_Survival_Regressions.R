## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for survival
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(lme4)
library(performance)
library(splines)
library(DHARMa)
library(broom.mixed)

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

#creating dataframe for survival vs height analysis
survival_dat <- pines_long|>
  arrange(IND_ID, CENSUS_NUM)|>
  group_by(IND_ID)|>
  #creating column for death (1 if dead, 0 if alive)
  mutate(Dead = case_when(STAT == 0 ~ 1L, 
                    STAT==1 ~ 0L, 
                  .default = STAT))|>
  #pairing height & survival in census t with survival in census t+1
  mutate(STAT_next = lead(STAT),
        Dead_next = lead(Dead),
         CENSUS_next = lead(CENSUS_NUM),
         Date_next = lead(START_DATE))|>
  #only keeping trees that were alive in census t
  #and had height measurements in census t
  #and had their survival recorded in census t+1
  filter(STAT==1 & !is.na(STAT_next) & !is.na(HT))|>
  #calculating time interval b/w censuses, in years
  mutate(CensusInterval = 
           time_length(difftime(Date_next, START_DATE),"years"))|>
  mutate(CensusInterval = round(CensusInterval, digits = 1))

#ensuring that only consecutive censuses are paired
table(survival_dat$CENSUS_next - survival_dat$CENSUS_NUM) #everything is correct

#checking values of census interval
table(survival_dat$CensusInterval) #looks correct

#Exploratory model---------------------------------------
# fitting a model with binned height 
# to explore shape of height vs survival relationship

#creating height bins
survival_dat$HT_bin <- cut(survival_dat$HT,
   breaks = quantile(survival_dat$HT, 
    probs = seq(0,1,0.05), na.rm=TRUE), include.lowest = TRUE)

#checking sample size in each bin
survival_dat|>
  ungroup()|>
  count(HT_bin) #>500 points in each bin, i.e. sufficient sample size

#fitting model on binned height
m_ht_bins <- glmer(Dead_next ~ HT_bin + AREA + TimeSinceFire +
  offset(log(CensusInterval)) +
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat, family = binomial(link="cloglog"))

#extracting coefficients for plotting
coefs <- tidy(m_ht_bins, 
  effects = "fixed", conf.int = TRUE)
coefs <- coefs |> filter(str_detect(term, "^HT_bin"))
coefs <- coefs |>
  rename(cloglog_est = estimate, 
    cloglog_lo = conf.low, 
    cloglog_hi = conf.high)

#extracting intercept
intercept <- tidy(m_ht_bins, effects = "fixed") |>
  filter(term == "(Intercept)") |> pull(estimate)

#setting the first bin as reference
ref_bin <- levels(survival_dat$HT_bin)[1]
coefs <- bind_rows(
  tibble(term = paste0("HT_bin", ref_bin), 
  cloglog_est = 0, cloglog_lo = NA, cloglog_hi = NA),
  coefs)

#extracting mid point of each bin
bin_levels <- levels(survival_dat$HT_bin)
bin_bounds <- bin_levels |>
  # parsing the (a,b] interval labels into numeric midpoints
  str_remove_all("\\[|\\]|\\(|\\)") |>
  str_split(",", simplify = TRUE) |>
  apply(2, as.numeric)

bin_midpoints <- tibble(
  term = paste0("HT_bin", bin_levels),
  midpoint = rowMeans(bin_bounds)
)

coefs <- coefs |> left_join(bin_midpoints, by = "term")

#Plotting
ggplot(coefs, aes(x = midpoint, y = cloglog_est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = cloglog_lo,
     ymax = cloglog_hi), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Height (cm), bin midpoint",
       y = "CLog-log hazard relative to reference bin",
       title = "Nonparametric mortality-hazard shape by height bin") +
  theme_bw()

#Fitting models----------------------------------------------------
#Baseline model - random effects and time offset only
m_baseline <- glmer(Dead_next ~ 
  offset(log(CensusInterval))+
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height only model - linear effect of height
m_ht <- glmer(Dead_next ~ HT+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height & AREA
m_ht_area <- glmer(Dead_next ~ HT+AREA+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height & AREA - ht x area interaction
m_ht_area_htarea <- glmer(Dead_next ~ HT*AREA+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire
m_ht_area_time <- glmer(Dead_next ~ 
  HT+AREA+TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction
m_ht_area_htarea_time <- glmer(Dead_next ~ HT*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, quadratic height
m_ht2_area_time <- glmer(Dead_next ~ poly(HT, 2, raw=T)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, quadratic ht
m_ht2_area_htarea_time <- glmer(Dead_next ~ poly(HT, 2, raw=T)*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, linear plus log height
m_htLL_area_time <- glmer(Dead_next ~ HT+log(HT)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, linear + log ht
m_htLL_area_htarea_time <- glmer(Dead_next ~ HT+log(HT)+AREA+
  HT:AREA+log(HT):AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, log height
m_htLog_area_time <- glmer(Dead_next ~ log(HT)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, log ht
m_htLog_area_htarea_time <- glmer(Dead_next ~ log(HT)*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, spline for height
m_htS_area_time <- glmer(Dead_next ~ ns(HT, 3)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, spline for height
m_htS_area_htarea_time <- glmer(Dead_next ~ ns(HT, 3)*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|SITEAREA)+(1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

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
  m_ht2_area_htarea_time = m_ht2_area_htarea_time,
  m_htLL_area_time = m_htLL_area_time,
  m_htLL_area_htarea_time = m_htLL_area_htarea_time,
  m_htLog_area_time = m_htLog_area_time,
  m_htLog_area_htarea_time = m_htLog_area_htarea_time,
  m_htS_area_time = m_htS_area_time,
  m_htS_area_htarea_time = m_htS_area_htarea_time)

#checking for singular fit issues, convergence issues,
#  overdispersion and across-plots variation
model_diagnostics <- tibble::tibble(
  Name = names(models),
  Singular = sapply(models, isSingular),
  ConvergenceWarnings = sapply(models,
     \(x) x@optinfo$conv$lme4$messages),
  DispersionRatio = sapply(models, \(x) 
      check_overdispersion(x)$dispersion_ratio),
  AcrossPlotVar = sapply(models, \(x) VarCorr(x)$PLOTCODE|>as.numeric()),
  AcrossIndVar = sapply(models, \(x) VarCorr(x)$IND_ID|>as.numeric()))
model_diagnostics

#plotting DHARMa residuals
residuals <- lapply(models, simulateResiduals)
for (m in names(residuals)) {
  plot(residuals[[m]], title = m, quantreg = T)
  dev.new()
} #some of the models have significant overdispersion, but effect size is not too high

#checking for significant outliers
model_diagnostics$OutlierSig = NA
for(m in names(residuals)){
  print(m)
  outlier_test <- 
      testOutliers(residuals[[m]], type = "bootstrap")
  print(outlier_test) 

  model_diagnostics$OutlierSig[model_diagnostics$Name==m] <- 
    outlier_test$p.value < 0.05
}
View(model_diagnostics) #some models have significant outliers, but their number is small

#checking the shape of height vs survival curves
newdat <- expand.grid( #generating new data for plotting
  HT = seq(
    min(survival_dat$HT, na.rm = TRUE),
    max(survival_dat$HT, na.rm = TRUE),
    length.out = 1000),
  TimeSinceFire = mean(survival_dat$TimeSinceFire),
  AREA = levels(survival_dat$AREA),
  CensusInterval = 1)

#Plotting shape of curves
newdat$MortProb_linear <- predict(m_ht_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_spline <- predict(m_htS_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_LL <- predict(m_htLL_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_Log <- predict(m_htLog_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_quad <- predict(m_ht2_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)

ggplot(newdat, aes(HT, MortProb_linear, col=AREA))+
  geom_point()+theme_bw()+xlab("Height")+
  ylab("Predicted Mortality Probability")+
  ggtitle("Linear effect of height")
ggplot(newdat, aes(HT, MortProb_LL, col=AREA))+
  geom_point()+theme_bw()+xlab("Height")+
  ylab("Predicted Mortality Probability")+
  ggtitle("Linear + log effect of height")
ggplot(newdat, aes(HT, MortProb_Log, col=AREA))+
  geom_point()+theme_bw()+xlab("Height")+
  ylab("Predicted Mortality Probability")+
  ggtitle("Log effect of height")
ggplot(newdat, aes(HT, MortProb_quad, col=AREA))+
  geom_point()+theme_bw()+xlab("Height")+
  ylab("Predicted Mortality Probability")+
  ggtitle("quadratic effect of height")
ggplot(newdat, aes(HT, MortProb_spline, col=AREA))+
  geom_point()+theme_bw()+xlab("Height")+
  ylab("Predicted Mortality Probability")+
  ggtitle("Spline for height")


survival_dat|>filter(CENSUS_NUM==11)|>
  ggplot(aes(HT, Dead_next,fill = AREA))+
  geom_point(pch=21, size = 2, alpha =0.5)+theme_bw()+
  facet_wrap(~AREA)

#Model comparisons on full dataset-----------------------------------------
# Comparing models based on their fit to training data        
#comparing AICs
model_perform <- tibble::tibble(
  Name = names(models),
  AIC = sapply(models, AIC))
model_perform|>arrange(AIC)

#comparing brier scores
brier_score_train <- function(model, obs = survival_dat$Dead_next){
  pred <- predict(model, type="response",re.form = NA)
  return(mean((obs - pred)^2))
}
model_perform$BrierScore <- sapply(models, brier_score_train)
model_perform

#comparing log loss
logloss_train <- function(model, obs = survival_dat$Dead_next){
  eps <-   1e-15
  pred <- predict(model, type="response",re.form = NA)
  pred <- pmin(pmax(pred, eps), 1 - eps)
  return(-mean(obs * log(pred) +
      (1 - obs) * log(1 - pred)))
}

model_perform$LogLoss <- sapply(models, logloss_train)
model_perform


#Model comparisons 2: time-series test-train split-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data

model_perform$TT_Brier <- NA
model_perform$TT_LogLoss <- NA

# Test-train split - 
# data upto census 9 used for training, rest for testing
TrainDat <- survival_dat|>
  filter(CENSUS_NUM<=9)
TestDat <- survival_dat|>
  filter(CENSUS_NUM>9)

#Redefining functions for brier score and logloss
brier_score <- function(pred, obs){
  return(mean((obs - pred)^2))}

logloss <- function(pred, obs){
  eps <-   1e-15
  pred <- pmin(pmax(pred, eps), 1 - eps)
  return(-mean(obs * log(pred) +
      (1 - obs) * log(1 - pred)))}


#Re-fitting models and testing accuracy
for(m in names(models)){
  model <- models[[m]]
  new_model <- glmer(formula = formula(model),
    data = TrainDat,
    family = binomial(link="cloglog"))
  preds <- predict(new_model, newdata= TestDat, 
            type = "response", re.form = NA)
  model_perform$TT_Brier[model_perform$Name ==m] <-
    brier_score(preds, TestDat$Dead_next)
  model_perform$TT_LogLoss[model_perform$Name ==m] <-
    logloss(preds, TestDat$Dead_next)

  print(c(m, isSingular(new_model)))
  print(new_model@optinfo$conv$lme4$messages)
}
model_perform|>
  select(Name, AIC,TT_Brier,TT_LogLoss)|>
  arrange(TT_Brier)

#Model comparisons 3: time-series cross-validation-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data
# with 3 fold expanding window time-series cross validation (TSV)

#Adding columns for storing CV ccuracy scores
model_perform <- model_perform|>
  mutate(CV1_Brier = NA, CV2_Brier = NA, CV3_Brier =NA,
  CV1_LogLoss = NA, CV2_LogLoss = NA, CV3_LogLoss =NA)

folds <- 3 #number of cross-validation folds
min_training_size <- 8 #minimum number of censuses in training data
test_size <- 4 #number of census in test data

for( i in 1:folds){
  #censuses in training data
  train_censuses <- 1:(min_training_size+i-1)

  #censuses in test data
  test_censuses <- (min_training_size+i):(min_training_size+i+test_size-1)

  #test-train split
  TrainDat <- survival_dat|>
    filter(CENSUS_NUM %in% train_censuses)
  TestDat <- survival_dat|>
    filter(CENSUS_NUM %in% test_censuses)

  #Re-fitting models and testing accuracy
  for(m in names(models)){
    model <- models[[m]]
    
    #fitting model
    new_model <- glmer(formula = formula(model),
      data = TrainDat,
      family = binomial(link="cloglog"))
    
    #predicting heldout data
    preds <- predict(new_model, newdata= TestDat, 
            type = "response", re.form = NA)
    
    #calculating accuracy
    brier <-
      brier_score(preds, TestDat$Dead_next)
    LL <-
      logloss(preds, TestDat$Dead_next)
    
    #recording results
    col_brier <- paste0("CV", i, "_Brier")
    model_perform[which(model_perform$Name==m), col_brier] <- brier
    col_LL <- paste0("CV", i, "_LogLoss")
    model_perform[which(model_perform$Name==m), col_LL] <- LL

    #tracking progress
    print(c(i, m, isSingular(new_model)))
    print(new_model@optinfo$conv$lme4$messages)
  }
}

model_perform|>
  select(Name, AIC, starts_with("CV") & ends_with("Brier"))|>
  rowwise()|>
  mutate(mean_Brier = mean(c_across(starts_with("CV"))))|>
  arrange(mean_Brier)

model_perform|>
  select(Name, AIC, starts_with("CV") & ends_with("LogLoss"))|>
  rowwise()|>
  mutate(mean_LogLoss = mean(c_across(starts_with("CV"))))|>
  arrange(mean_LogLoss)

