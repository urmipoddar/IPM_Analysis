## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for survival
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(glmmTMB)
library(performance)
library(splines)
library(DHARMa)
library(broom.mixed)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")
plot_info <- read.csv("Data/plot_information.csv")
#Functions---------------------------------------------------
#Redefining functions for brier score and logloss
brier_score <- function(pred, obs){
  return(mean((obs - pred)^2))}

logloss <- function(pred, obs){
  eps <-   1e-15
  pred <- pmin(pmax(pred, eps), 1 - eps)
  return(-mean(obs * log(pred) +
      (1 - obs) * log(1 - pred)))}
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

#Fitting models----------------------------------------------------
#Baseline model - random effects and time offset only
mSurv_baseline <- glmmTMB(Dead_next ~ 
  offset(log(CensusInterval))+
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height only model - linear effect of height
mSurv_ht <- glmmTMB(Dead_next ~ HT+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height & AREA
mSurv_ht_area <- glmmTMB(Dead_next ~ HT+AREA+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height & AREA - ht x area interaction
mSurv_ht_area_htarea <- glmmTMB(Dead_next ~ HT*AREA+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire
mSurv_ht_area_time <- glmmTMB(Dead_next ~ 
  HT+AREA+TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction
mSurv_ht_area_htarea_time <- glmmTMB(Dead_next ~ HT*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, quadratic height
mSurv_ht2_area_time <- glmmTMB(Dead_next ~ 
  HT+ I(HT^2)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, quadratic ht
mSurv_ht2_area_htarea_time <- glmmTMB(Dead_next ~ 
  (HT+ I(HT^2))*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, linear plus log height
mSurv_htLL_area_time <- glmmTMB(Dead_next ~ HT+log(HT)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, linear + log ht
mSurv_htLL_area_htarea_time <- glmmTMB(Dead_next ~ HT+log(HT)+AREA+
  HT:AREA+log(HT):AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, log height
mSurv_htLog_area_time <- glmmTMB(Dead_next ~ log(HT)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
  family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, log ht
mSurv_htLog_area_htarea_time <- glmmTMB(Dead_next ~ log(HT)*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, spline for height
mSurv_htS_area_time <- glmmTMB(Dead_next ~ ns(HT, 3)+AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#Height, AREA & time since fire, height x area interaction, spline for height
mSurv_htS_area_htarea_time <- glmmTMB(Dead_next ~ ns(HT, 3)*AREA+
  TimeSinceFire+
  offset(log(CensusInterval))+ 
  (1|PLOTCODE)+(1|IND_ID),
  data = survival_dat,
 family = binomial(link="cloglog"))

#model diagnostics-------------------------------------------------
#list of models
models <- list(
  mSurv_baseline = mSurv_baseline,
  mSurv_ht = mSurv_ht,
  mSurv_ht_area = mSurv_ht_area,
  mSurv_ht_area_htarea = mSurv_ht_area_htarea,
  mSurv_ht_area_time = mSurv_ht_area_time,
  mSurv_ht_area_htarea_time = mSurv_ht_area_htarea_time,
  mSurv_ht2_area_time = mSurv_ht2_area_time,
  mSurv_ht2_area_htarea_time = mSurv_ht2_area_htarea_time,
  mSurv_htLL_area_time = mSurv_htLL_area_time,
  mSurv_htLL_area_htarea_time = mSurv_htLL_area_htarea_time,
  mSurv_htLog_area_time = mSurv_htLog_area_time,
  mSurv_htLog_area_htarea_time = mSurv_htLog_area_htarea_time,
  mSurv_htS_area_time = mSurv_htS_area_time,
  mSurv_htS_area_htarea_time = mSurv_htS_area_htarea_time)

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
View(model_diagnostics) #some models have significant outliers, but their num of outliers is small

#Shape of height vs survival curves-----------------------------
#Generating new data for plotting
newdat <- expand.grid( 
  HT = seq(
    min(survival_dat$HT, na.rm = TRUE),
    max(survival_dat$HT, na.rm = TRUE),
    length.out = 1000),
  TimeSinceFire = mean(survival_dat$TimeSinceFire),
  AREA = levels(survival_dat$AREA),
  CensusInterval = 1)

#Generating predictions from different models
newdat$MortProb_linear <- predict(mSurv_ht_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_spline <- predict(mSurv_htS_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_LL <- predict(mSurv_htLL_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_Log <- predict(mSurv_htLog_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)
newdat$MortProb_quad <- predict(mSurv_ht2_area_htarea_time,
          newdata = newdat, type = "response", 
        re.form = NA)

#Plotting
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

#Plotting observed data 
survival_dat|>filter(CENSUS_NUM==11)|>
  ggplot(aes(HT, Dead_next,fill = AREA))+
  geom_point(pch=21, size = 2, alpha =0.5)+theme_bw()+
  facet_wrap(~AREA)

#Model fit comparisons-----------------------------------------
# Comparing models based on their fit to training data        
#comparing AICs
model_perform <- tibble::tibble(
  Name = names(models),
  AIC = sapply(models, AIC))
model_perform|>arrange(AIC)

#comparing brier scores
model_perform$BrierScore <- sapply(models, 
          function(x){brier_score(
            predict(x, type="response",re.form = NA),
            survival_dat$Dead_next)})
model_perform

#comparing log loss
model_perform$LogLoss <- sapply(models, function(x){logloss(
            predict(x, type="response",re.form = NA),
            survival_dat$Dead_next)})
model_perform


#Model accuracy comparisons 1: time-series test-train split-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data

model_perform$TT_Brier <- NA
model_perform$TT_LogLoss <- NA

# Test-train split - 
# data upto census 9 used for training, rest for testing
TrainDat <- survival_dat|>
  filter(CENSUS_NUM<=9)
TestDat <- survival_dat|>
  filter(CENSUS_NUM>9)

#Re-fitting models and testing accuracy
for(m in names(models)){
  model <- models[[m]]
  new_model <- glmmTMB(formula = formula(model),
    data = TrainDat,
    family = binomial(link="cloglog"))
  preds <- predict(new_model, newdata= TestDat, 
            type = "response", re.form = NA)
  model_perform$TT_Brier[model_perform$Name ==m] <-
    brier_score(preds, TestDat$Dead_next)
  model_perform$TT_LogLoss[model_perform$Name ==m] <-
    logloss(preds, TestDat$Dead_next)

  print(c(m, check_singularity(new_model)))
  print(new_model$fit$convergence)
}
model_perform|>
  select(Name, AIC,TT_Brier,TT_LogLoss)|>
  arrange(TT_Brier)

#Model accuracy comparisons 2: time-series cross-validation-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data
# with 3 fold expanding window time-series cross validation (TSV)

#Adding columns for storing TSV ccuracy scores
model_perform <- model_perform|>
  mutate(TSV1_Brier = NA, TSV2_Brier = NA, TSV3_Brier =NA,
  TSV1_LogLoss = NA, TSV2_LogLoss = NA, TSV3_LogLoss =NA)

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
    new_model <- glmmTMB(formula = formula(model),
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
    col_brier <- paste0("TSV", i, "_Brier")
    model_perform[which(model_perform$Name==m), col_brier] <- brier
    col_LL <- paste0("TSV", i, "_LogLoss")
    model_perform[which(model_perform$Name==m), col_LL] <- LL
    
    # tracking progress
    print(c(i,  m, check_singularity(new_model)))
    print(model$fit$convergence)
  }
}

model_perform|>
  select(Name, AIC, starts_with("TSV") & ends_with("Brier"))|>
  rowwise()|>
  mutate(mean_Brier = mean(c_across(starts_with("TSV"))))|>
  arrange(mean_Brier)

model_perform|>
  select(Name, AIC, starts_with("TSV") & ends_with("LogLoss"))|>
  rowwise()|>
  mutate(mean_LogLoss = mean(c_across(starts_with("TSV"))))|>
  arrange(mean_LogLoss)

#Model accuracy comparison 3: Spatial cross validation------------------------------------------------
# Comparing models based on predictive accuracy on unseen SITEASREAs
# with leave-one-site-out (LOSO) cross validation
#
# SCC5 is excluded as a holdout target 
# because it is the only SITEAREA in SCC
# so holding it out would entirely drop that AREA level from training data
sites <- setdiff(levels(survival_dat$SITEAREA), "SCC5")
n_sites <- length(sites)   # 7

# Adding columns for storing spatial CV accuracy scores
model_perform <- model_perform |>
  mutate(!!!setNames(rep(list(NA_real_), n_sites*2),
                      c(paste0("SP", 1:n_sites, "_Brier"),
                        paste0("SP", 1:n_sites, "_LogLoss"))))

for (i in 1:n_sites) {

  held_out_site <- sites[i]

  # spatial train-test split: hold out one site 
  TrainDat <- survival_dat |> filter(SITEAREA != held_out_site)
  TestDat  <- survival_dat |> filter(SITEAREA == held_out_site)

  for (m in names(models)) {
    model <- models[[m]]

    # refitting model on all sites except the held-out one
    new_model <- glmmTMB(formula = formula(model),
                        data = TrainDat,
                        family = binomial(link = "cloglog"))

    # predicting the held-out data
    preds <- predict(new_model, newdata = TestDat, type = "response", re.form = NA)

    brier   <- brier_score(preds, TestDat$Dead_next)
    LL      <- logloss(preds, TestDat$Dead_next)

    # recording results
    col_brier <- paste0("SP", i, "_Brier")
    model_perform[which(model_perform$Name == m), col_brier] <- brier
    col_LL <- paste0("SP", i, "_LogLoss")
    model_perform[which(model_perform$Name == m), col_LL] <- LL

    # tracking progress
    print(c(i, held_out_site, m, check_singularity(new_model)))
    print(model$fit$convergence)
  }
}

model_perform |>
  select(Name, AIC, starts_with("SP") & ends_with("Brier")) |>
  rowwise() |>
  mutate(mean_Brier = mean(c_across(starts_with("SP")), na.rm = TRUE)) |>
  arrange(mean_Brier)

model_perform |>
  select(Name, AIC, starts_with("SP") & ends_with("LogLoss")) |>
  rowwise() |>
  mutate(mean_LogLoss = mean(c_across(starts_with("SP")), na.rm = TRUE)) |>
  arrange(mean_LogLoss)

#Choosing best models------------------------------------------------
# Ranking models based on AIC and CV accuracy
# and choosing top 4
 
n_top <- 4 #number of models to choose
 
model_perform_summary <- model_perform |>
  rowwise() |>
  mutate(mean_TSV_Brier = mean(c_across(starts_with("TSV") &
    ends_with("Brier"))),
         mean_SP_Brier = mean(c_across(starts_with("SP") &
                   ends_with("Brier")), na.rm = TRUE)) |>
  ungroup() |>
  mutate(rank_AIC = rank(AIC),
         rank_TSV = rank(mean_TSV_Brier),
         rank_SP  = rank(mean_SP_Brier),
         rank_sum = rank_AIC + rank_TSV + rank_SP) |>
  select(Name, AIC, mean_TSV_Brier, 
    mean_SP_Brier, rank_AIC, rank_TSV, rank_SP, rank_sum) |>
  arrange(rank_sum)
 
model_perform_summary|>View()
 
chosen_models <- model_perform_summary$Name[1:n_top]
chosen_models

#savingt top 4 models 
if (!dir.exists("Output")) {
  dir.create("Output", recursive = TRUE)
}
saveRDS(models[chosen_models], "Output/SurvivalModels.rds")
