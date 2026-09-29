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
library(qgam)
library(performance)
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

# function for translating glmmTMB-style random effects formula
#  into gamlss syntax
translate_to_gamlss_formula <- function(glmmTMB_formula){
  f_chr <- paste(deparse(glmmTMB_formula), collapse = " ") #extracting model formula
  f_chr <- gsub("\\s+", " ", f_chr) #collapsing all whitespaces to a single space
  f_chr <- gsub("\\(1 \\| SITEAREA\\)", "random(SITEAREA)", f_chr)
  f_chr <- gsub("\\(1 \\| PLOTCODE\\)", "random(PLOTCODE)", f_chr)
  f_chr <- gsub("\\(1 \\| IND_ID\\)", "random(IND_ID)", f_chr)
  as.formula(f_chr)
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

#Fitting models-------------------------------------
# Based on prelim analysis:
# All models are fitted on log(HT_next)
# with log(HT) as predictor
# Residuals are t-distributed rather than normal, to account for heavy tails
# and variance is modelled as a function of log(HT) to account for heteroskedascity

#Baseline model - random effects & census interval only
mGrwt_baseline <- glmmTMB(log(HT_next) ~ CensusInterval+
  (1|PLOTCODE)+(1|IND_ID),
  family = t_family(), 
  data = growth_dat)

#Height & AREA -log linear effect of height
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

#Height & AREA - ht x area interaction, linear height
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
mGrwt_ht2_area_time <- glmmTMB(log(HT_next) ~ 
  poly(log(HT), 2)+AREA + TimeSinceFire + CensusInterval+
  (1|PLOTCODE)+(1|IND_ID), data = growth_dat,
    family = t_family(),
  dispformula = ~ log(HT))

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

#model diagnostics----------------------------------------------
#list of models
models <- list(
  mGrwt_baseline = mGrwt_baseline,
  mGrwt_ht_area = mGrwt_ht_area,
  Grwt_ht2_area = mGrwt_ht2_area,
  mGrwt_ht_area_htarea = mGrwt_ht_area_htarea,
  Grwt_ht2_area_htarea = mGrwt_ht2_area_htarea,
  mGrwt_ht_area_time = mGrwt_ht_area_time,
  mGrwt_ht2_area_time = mGrwt_ht2_area_time,
  mGrwt_ht_area_htarea_time = mGrwt_ht_area_htarea_time,
  mGrwt_ht2_area_htarea_time = mGrwt_ht2_area_htarea_time)

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
model_diagnostics|>View() #no major issues

#histograms of residuals
for(m in names(models)){
  model = models[[m]]
  hist(residuals(model), breaks = 50, main = m)
  dev.new()
} #all residual distributions are heavy tailed and skewed

#QQ-plots
for(m in names(models)){
  model = models[[m]]
  qqnorm(residuals(model), main = m)
  qqline(residuals(model))
  dev.new()
} #this again shows the heavy tails

#DHARMa residual plots
res <- lapply(models, simulateResiduals)
for (m in names(res)) {
  plot(res[[m]], title = m, quantreg = T)
  dev.new()
} #KS test is siginificant for all models, some models also show heteroskedacity  
for (m in names(res)) {
  hist(res[[m]], main = m)
  dev.new()
} #skewed residual distributions

#checking shapes of predicted curves----------------------------
#Building new data grid with all predictors
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
total_resid_linear <- log(growth_dat$HT_next) - 
          mu_marginal_linear   # captures RE variance + residual + any skew
smear_factor <- mean(exp(total_resid_linear))#Daun's smearing factor
newdat$PredLinear_corrected <- 
      exp(newdat$PredLinear)*smear_factor

#Predictions on the raw scale - quadratic model
mu_marginal_quad <- predict(
  mGrwt_ht2_area_htarea,
  newdata = growth_dat,
  type = "response",
  re.form = NA)
total_resid_quad <- log(growth_dat$HT_next) - 
          mu_marginal_quad
smear_factor_quad <- mean(exp(total_resid_quad)) #Daun's smearing factor
newdat$PredQuad_corrected <- 
                exp(newdat$PredQuad) * smear_factor_quad


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
  mutate(Pseudo_R_squared = 
    1-MSE/var(log(growth_dat$HT_next)))|>
  arrange(MSE, AIC)


#Model accuracy comparisons 1: time-series test-train split-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data
 
model_perform$TT_MSE <- NA
model_perform$TT_MAE <- NA
 
# Test-train split:
# data upto census 9 used for training, rest for testing
TrainDat <- growth_dat|>
  filter(CENSUS_NUM<=9)
TestDat <- growth_dat|>
  filter(CENSUS_NUM>9)
 
#Re-fitting models and testing accuracy
for(m in names(models)){
  model <- models[[m]]
  new_model <- glmmTMB(formula = formula(model),
    data = TrainDat,
    family = t_family(),
    dispformula = model$modelInfo$allForm$dispformula)
  preds <- predict(new_model, newdata= TestDat, 
            type = "response", re.form = NA)
  model_perform$TT_MSE[model_perform$Name ==m] <-
    mse(preds, log(TestDat$HT_next))
  model_perform$TT_MAE[model_perform$Name ==m] <-
    mae(preds, log(TestDat$HT_next))
 
  print(c(m, check_singularity(new_model)))
  print(new_model$fit$convergence)
}

model_perform|>
  select(Name, AIC,TT_MSE,TT_MAE)|>
  mutate(TT_Psuedo_R2 =
     1- TT_MSE/var(log(TestDat$HT_next)))|>
  arrange(TT_MSE)

#Model accuracy comparisons 2: time-series cross-validation-----------------------------------------        
# Comparing models based on predictive accuracy on unseen data
# with 3 fold expanding window time-series cross validation (TSV)
#
# Size measurements only available for
#  censuses 3, 5, 7, 9, 11-15
# i.e. 8 distinct pairs of censuses

census_levels <- sort(unique(growth_dat$CENSUS_NUM))
 
#Adding columns for storing TSV accuracy scores
model_perform <- model_perform|>
  mutate(TSV1_MSE = NA, TSV2_MSE = NA, TSV3_MSE =NA,
  TSV1_MAE = NA, TSV2_MAE = NA, TSV3_MAE =NA)
 
folds <- 3 #number of cross-validation folds
min_training_size <- 5 #minimum number of census-points in training data
test_size <- 1 #number of census-points in test data -- 
              #kept at 1 so all 3 expanding-window folds fit within 
              #the 8 available #census points 
 
for( i in 1:folds){
  #census-points in training data
  train_censuses <- census_levels[1:(min_training_size+i-1)]
 
  #census-points in test data
  test_censuses <- 
    census_levels[(min_training_size+i):
      (min_training_size+i+test_size-1)]
 
  #test-train split
  TrainDat <- growth_dat|>
    filter(CENSUS_NUM %in% train_censuses)
  TestDat <- growth_dat|>
    filter(CENSUS_NUM %in% test_censuses)
 
  #Re-fitting models and testing accuracy
  for(m in names(models)){
    model <- models[[m]]
 
    #fitting model
    new_model <- glmmTMB(formula = formula(model),
      data = TrainDat,
      family = t_family(),
      dispformula = model$modelInfo$allForm$dispformula)
 
    #predicting heldout data
    preds <- predict(new_model, newdata= TestDat, 
            type = "response", re.form = NA)
 
    #calculating accuracy
    MSE_val <- mse(preds, log(TestDat$HT_next))
    MAE_val <- mae(preds, log(TestDat$HT_next))
 
    #recording results
    col_mse <- paste0("TSV", i, "_MSE")
    model_perform[which(model_perform$Name==m), col_mse] <- MSE_val
    col_mae <- paste0("TSV", i, "_MAE")
    model_perform[which(model_perform$Name==m), col_mae] <- MAE_val
 
    # tracking progress
    print(c(i,  m, check_singularity(new_model)))
    print(new_model$fit$convergence)
  }
}
 
model_perform|>
  select(Name, AIC, starts_with("TSV") & ends_with("MSE"))|>
  rowwise()|>
  mutate(mean_TSV_MSE = mean(c_across(starts_with("TSV"))),
          sd_TSV_MSE = sd(c_across(starts_with("TSV"))))|>
  arrange(mean_TSV_MSE)
 
model_perform|>
  select(Name, AIC, starts_with("TSV") & ends_with("MAE"))|>
  rowwise()|>
  mutate(mean_TSV_MAE = mean(c_across(starts_with("TSV"))),
          sd_TSV_MAE = sd(c_across(starts_with("TSV"))))|>
  arrange(mean_TSV_MAE)

#Model accuracy comparison 3: Spatial cross validation------------------------------------------------
# Comparing models based on predictive accuracy on unseen SITEAREAs
# with leave-one-site-out (LOSO) cross validation
#
# SCC5 is excluded as a holdout target 
# because it is the only SITEAREA in SCC
# so holding it out would entirely drop that AREA level from training data

sites <- setdiff(levels(growth_dat$SITEAREA), "SCC5")
n_sites <- length(sites)
 
# Adding columns for storing spatial CV accuracy scores
model_perform <- model_perform |>
  mutate(!!!setNames(rep(list(NA_real_), n_sites*2),
                      c(paste0("SP", 1:n_sites, "_MSE"),
                        paste0("SP", 1:n_sites, "_MAE"))))
 
for (i in 1:n_sites) {
 
  held_out_site <- sites[i]
 
  # spatial train-test split: hold out one site
  TrainDat <- growth_dat |> filter(SITEAREA != held_out_site)
  TestDat  <- growth_dat |> filter(SITEAREA == held_out_site)
 
  for (m in names(models)) {
    model <- models[[m]]
 
    # refitting model on all sites except the held-out one
    new_model <- glmmTMB(formula = formula(model),
                        data = TrainDat,
                        family = t_family(),
                        dispformula = model$modelInfo$allForm$dispformula)
 
    # predicting the held-out data
    preds <- predict(new_model, newdata = TestDat, type = "response", re.form = NA)
 
    MSE_val <- mse(preds, log(TestDat$HT_next))
    MAE_val <- mae(preds, log(TestDat$HT_next))
 
    # recording results
    col_mse <- paste0("SP", i, "_MSE")
    model_perform[which(model_perform$Name == m), col_mse] <- MSE_val
    col_mae <- paste0("SP", i, "_MAE")
    model_perform[which(model_perform$Name == m), col_mae] <- MAE_val
 
    # tracking progress
    print(c(i, held_out_site, m, check_singularity(new_model)))
    print(new_model$fit$convergence)
  }
}
 
model_perform |>
  select(Name, AIC, starts_with("SP") & ends_with("MSE")) |>
  rowwise() |>
  mutate(mean_SP_MSE = 
    mean(c_across(starts_with("SP")), na.rm = TRUE),
     sd_SP_MSE = sd(c_across(starts_with("SP")), na.rm = TRUE)) |>
  arrange(mean_SP_MSE)|>View()
 
model_perform |>
  select(Name, AIC, starts_with("SP") & ends_with("MAE")) |>
  rowwise() |>
  mutate(mean_SP_MAE = 
    mean(c_across(starts_with("SP")), na.rm = TRUE)) |>
  arrange(mean_SP_MAE)|>View()

#Choosing the best mean-structured models------------------------------------------------
# Ranking models based on AIC and CV accuracy
# and choosing top 4
 
n_top <- 4 #number of models to choose
 
model_perform_summary <- model_perform |>
  rowwise() |>
  mutate(mean_TSV_MSE = mean(c_across(starts_with("TSV") & ends_with("MSE"))),
         mean_SP_MSE  = mean(c_across(starts_with("SP") & ends_with("MSE")), na.rm = TRUE)) |>
  ungroup() |>
  mutate(rank_AIC = rank(AIC),
         rank_TSV = rank(mean_TSV_MSE),
         rank_SP  = rank(mean_SP_MSE),
         rank_sum = rank_AIC + rank_TSV + rank_SP) |>
  select(Name, AIC, mean_TSV_MSE, mean_SP_MSE, rank_AIC, rank_TSV, rank_SP, rank_sum) |>
  arrange(rank_sum)
 
model_perform_summary|>View()
 
chosen_models <- model_perform_summary$Name[1:n_top]
chosen_models


#Examining skewness/kurtosis relationship with predicted size------------------------------
# Following Miller & Ellner (2025, Ecology, "My, how you've grown"): 
# this section fits spline quantile regressions(qgam)  
# to scaled residuals of the best pilot models, 
# then computes their nonparametric (quantile-based) skewness and excess kurtosis
# as a function of predicted size.

# Selecting best pilot model
diag_model <- models[[chosen_models[1]]]
 
#  fitted mean (log scale) and dispersion (scale parameter) for every observation
growth_dat$fitted_mean <- predict(diag_model, newdata = growth_dat,
                                   type = "response", re.form = NA)
growth_dat$fitted_disp <- predict(diag_model, newdata = growth_dat,
                                   type = "disp", re.form = NA)
 
# extracting scaled residuals
growth_dat$scaled_resid <- (log(growth_dat$HT_next) - growth_dat$fitted_mean) /
                                                  growth_dat$fitted_disp
 
# fittingspline quantile regression at 5/10/25/50/75/90/95% quantiles
taus <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
 
qgam_fits <- lapply(taus, function(tau) {
  qgam(scaled_resid ~ s(fitted_mean, k = 4), data = growth_dat, qu = tau)
})
names(qgam_fits) <- paste0("q", taus)
 
# predicting each quantile across the observed range of fitted values
pred_grid <- data.frame(fitted_mean = seq(min(growth_dat$fitted_mean),
                                           max(growth_dat$fitted_mean),
                                           length.out = 200))
for (tau in taus) {
  pred_grid[[paste0("q", tau)]] <- predict(qgam_fits[[paste0("q", tau)]],
                                            newdata = pred_grid)}
 
# calculating NP Skewness (Eq. 1-2 in the paper, alpha = 0.1, "Bowley's skewness") and
# NP Excess Kurtosis (Eq. 3, alpha = 0.05, scaled relative to a Gaussian so
# that 0 = Gaussian-like tails and +-1 indicates an extreme departure)
gaussian_np_kurtosis <- (qnorm(0.95) - qnorm(0.05)) / (qnorm(0.75) - qnorm(0.25))
 
pred_grid$NP_skewness <- (pred_grid$q0.1 + pred_grid$q0.9 - 2*pred_grid$q0.5) /
  (pred_grid$q0.9 - pred_grid$q0.1)
pred_grid$NP_kurtosis <- (pred_grid$q0.95 - pred_grid$q0.05) /
  (pred_grid$q0.75 - pred_grid$q0.25)
pred_grid$NP_excess_kurtosis <- pred_grid$NP_kurtosis / gaussian_np_kurtosis - 1
 
# Plotting scaled residuals + quantile regression lines (background),
quantile_long <- pred_grid |>
  select(fitted_mean, starts_with("q0")) |>
  pivot_longer(-fitted_mean, names_to = "quantile", values_to = "value")


# Plotting NP skewness (blue) and NP excess kurtosis (red) 
p_resid <- ggplot() +
  geom_point(data = growth_dat, aes(fitted_mean, scaled_resid), alpha = 0.15) +
  geom_line(data = quantile_long, aes(fitted_mean, value, group = quantile)) +
  theme_bw() +
  labs(x = "Fitted log(HT_next)", y = "Scaled residual",
       title = "Scaled residuals and quantile trends")
p_resid
 
p_skew_kurt <- ggplot(pred_grid) +
  geom_hline(yintercept = 0, lty = "dashed", col = "grey50") +
  geom_line(aes(fitted_mean, NP_skewness), col = "blue", linewidth = 1) +
  geom_line(aes(fitted_mean, NP_excess_kurtosis), col = "red", linewidth = 1) +
  theme_bw() +
  labs(x = "Fitted log(HT_next)",
       y = "NP skewness (blue) / NP excess kurtosis (red)",
       title = "NP skewness and excess kurtosis vs. predicted size")
p_skew_kurt #both skewness and kurtosis vary with predicted size

#Fitting skewed-t (gamlss) models for the chosen mean structures------------------------------------------------
# Both skewness and kurtosis were just shown to vary with predicted size,
# so sigma (variance parameter), nu (skewness parameter), and 
# tau (kurtosis parameter) are all specified as functions of log(HT) directly

#specifying formula for sigma, nu and tau parameters
skewt_formula <- list(
  sigma.formula = ~ log(HT),
  nu.formula    = ~ log(HT),
  tau.formula   = ~ log(HT))

#list for storing gamlss models
gamlss_models <- list()

#fitting skew-t models
for (m in chosen_models) {
  print(m)
  mean_formula_gamlss <- 
    translate_to_gamlss_formula(formula(models[[m]]))

  fit <- tryCatch(
    gamlss(
      formula = mean_formula_gamlss,
      sigma.formula = skewt_formula$sigma.formula,
      nu.formula = skewt_formula$nu.formula,
      tau.formula = skewt_formula$tau.formula,
      family = ST3,
      data = growth_dat),
    error = function(e) { message(m, " failed to fit: ", e$message); NULL }
  )

  gamlss_models[[m]] <- fit
}

# dropping any models  that failed to fit
gamlss_models <- gamlss_models[!sapply(gamlss_models, is.null)]
names(gamlss_models)

# AIC comparison
 sapply(gamlss_models, GAIC)#AIC ranks are similar to those for original models
model_perform_summary|>
  filter(Name %in% chosen_models)|>
  select(Name, AIC) 
#Comparing different distributionw------------------------------------------------
# In te previous section, models were fit using ST3 distribution
# Here I comparing the best fit ST3 model
# with other types of skewed t-distributions

#selecting the best ST3 model
best_gamlss_name <- names(gamlss_models)[
  which.min(sapply(gamlss_models, GAIC))]
best_gamlss_model <- gamlss_models[[best_gamlss_name]]

best_mean_formula <- formula(best_gamlss_model)  
st_families <- c("ST1", "ST2", "ST3", "ST4", "ST5")

#Parallelizng
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores, type = "PSOCK")
clusterEvalQ(cl, library(gamlss))
clusterExport(cl, varlist = c("best_mean_formula",
 "skewt_formula", "growth_dat"))
 
#Fitting models
st_models_list <- parLapply(cl, st_families, function(fam) {
  tryCatch(
    gamlss(
      formula = best_mean_formula,
      sigma.formula = skewt_formula$sigma.formula,
      nu.formula = skewt_formula$nu.formula,
      tau.formula = skewt_formula$tau.formula,
      family = fam,
      data = growth_dat),
    error = function(e) NULL
  )
})
 
stopCluster(cl)

names(st_models_list) <- paste0(best_gamlss_name,
   "__", st_families)

#checking if any models failed to fit
failed_st <- st_families[sapply(st_models_list, is.null)]
if (length(failed_st) > 0) message("Failed to fit: ", paste(failed_st, collapse = ", "))
st_models <- st_models_list[!sapply(st_models_list, is.null)]
names(st_models)
 
# AIC comparison 
st_gaic <- tibble::tibble(
  Family = names(st_models),
  AIC = sapply(st_models, GAIC))
st_gaic |> arrange(AIC) #ST3 model is the best
 
#Model diagnostics for skew-t models------------------------------------------------
# DHARMa does not support gamlss objects directly, so diagnostics use
# gamlss's own quantile-residual tools instead (4-panel plot + worm plot)

for (m in names(gamlss_models)) {
  model <- gamlss_models[[m]]
  plot(model, main = m)          # prints the quantile-residual summary as a side effect
  dev.new()
  wp(model, ylim.all = 2)        # worm plot; ylim.all set generously to avoid clipped points
  title(main = m)
  dev.new()
}
