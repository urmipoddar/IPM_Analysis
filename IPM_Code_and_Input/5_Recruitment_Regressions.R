## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Fitting regressions for recruitment
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(MASS)
library(glmmTMB)
library(gamlss)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")
plot_info <- read.csv("Data/plot_information.csv")

#Data formatting---------------------------------------------------
#Converting area and site area to factors
pines_long$AREA <- factor(pines_long$AREA, 
                          levels = c("DW", "SCC", "RP"))
pines_long$SITEAREA <- factor(pines_long$SITEAREA, 
                         levels = c("DW1", "DW2A", "DW2B",
                                    "DW2C", "DW3", "SCC5", "RP1", 
                                    "RP2"))

#Converting census start date column to Date format
#(useful for plotting)
pines_long <- pines_long|>
  mutate(START_DATE = as.Date(START_DATE))

#Creating dataframe for number of cones found after fire
cones_after_fire <- plot_info|>
  dplyr::select(OLDCODE, CONES1, PLTAREA)|>
  rename(PLOTCODE = OLDCODE)

#Calculating number of new recruits per census
recruits <- pines_long|>
  group_by(PLOTCODE, AREA,
     SITEAREA, CENSUS_NUM, START_DATE, SDLAREA)|>
  summarise(NumRecruits = sum(COHORT==CENSUS_NUM))

#recruitment in first year of study 
# i.e. ~ 1 yr after fire (census 1 and 2)
recruits1 = recruits|>
  filter(CENSUS_NUM %in% c(1,2))|>
  group_by(PLOTCODE, AREA, SITEAREA, SDLAREA)|>
  #summing up the number of recruits founds in census 1 and 2
  summarise(NumRecruits = sum(NumRecruits))|>
  #joining with cones after fire dataframe
  left_join(cones_after_fire)|>
  ungroup()

#Creating dataframe for survival of first 2 cohorts 
cohorts1_2 <- pines_long|>
  filter(COHORT %in% c(1,2) & CENSUS_NUM %in% 1:3)|>
  dplyr::select(IND_ID, AREA, SITEAREA, 
    PLOTCODE, SDLAREA, COHORT, CENSUS_NUM, STAT, HT) |>
  pivot_wider(names_from = CENSUS_NUM, 
    values_from = c(STAT, HT))|>
  #creating column for status at detection
  # i.e. census 1 STAT for cohort 1, census 2 STAT for cohort 2
  mutate(stat_det = if_else(COHORT == 1, STAT_1, STAT_2)) |>
   #keeping only seedlings that were alive when first counted
  filter(stat_det == 1) |>   
  #renaming STAT3 column
  rename(alive3 = STAT_3)|>
  #removing empty height columsn (HT_1 and 2)
  dplyr::select(!c(HT_1, HT_2))|>
  #renaming HT_3 column
  rename(HT = HT_3)

#Dataframe for size (height) of first two cohorts at census 3
cohorts1_2_alive3 <- cohorts1_2|>
  #keeping only individuals that were alive at census 3
  filter(alive3==1)

#Dataframe for size of new recruits at first detection
#for cohort 3 
cohort3_ht0 <- pines_long|>
  filter(COHORT ==3  & CENSUS_NUM ==3)|>
  filter(STAT == 1)|>
  filter(!is.na(HT))

#Year 1 recruitment EDA----------------------------------
ggplot(recruits1, aes(CONES1, NumRecruits))+
  geom_point(pch = 21, size = 3)+
  theme_bw()+xlab("Number of cones found after fire")+
  ylab("Number of 1st Yr Recruits")+
  scale_fill_discrete(name ="Plot Area (m^2)")

#Plotting of number of 1st yr recruit vs number of cones
ggplot(recruits1, aes(CONES1, NumRecruits))+
  geom_point(pch = 21, size = 3, aes(fill = factor(SDLAREA)))+
  theme_bw()+xlab("Number of cones found after fire")+
  ylab("Number of 1st Yr Recruits")+
  scale_fill_discrete(name ="Plot Area (m^2)")

#Plotting of density  of 1st yr recruits vs number of cones
#with seperate colors for plot size
ggplot(recruits1, aes(CONES1, NumRecruits/SDLAREA))+
  geom_point(pch = 21, size = 3, aes(fill = factor(SDLAREA)))+
  xlab("Number of cones found after fire")+
  ylab("Number of 1st Yr Recruits per unit area")+
  scale_fill_discrete(name ="Plot Area (m^2)")+
  theme_bw() #4m2 plots seem to have higher intercept for the recruit density vs cones relationship 

#Plotting of density  of 1st yr recruits vs number of cones
#with separate colors for AREA (population type)
ggplot(recruits1, aes(CONES1, NumRecruits/SDLAREA))+
  geom_point(pch = 21, size = 3,alpha = 0.8, aes(fill = AREA))+
  xlab("Number of cones found after fire")+
  ylab("Number of 1st Yr Recruits per unit area")+
  theme_bw() #hard to say whether there are differences b/w different AREAs 

#Adding a binary column for plot size in 1st yr recruits dataframe
#since the 4m^2 plots seem to have a average recruit density
#it might be important to include plot size as a predictor
recruits1$SmallPlot <- as.numeric(recruits1$SDLAREA==4)

#Year 1 recruitment model fitting----------------------------------------------
#Recruitment in 1st year of study (~1 yr after fire):
#modelled by regressing the number of cones found 
#against number of recruits in the 1st year (census 1+2)
#after controlling for plotsize (SDLAREA)

#Model with plot area and number of cones only
mRecruit1_none <- glm.nb(NumRecruits ~
    0 + SDLAREA + I(SDLAREA * CONES1), 
data = recruits1, link = identity)

#model with interaction with AREA
mRecruit1_AREA <- glm.nb(NumRecruits ~
    0 + SDLAREA + I(SDLAREA * CONES1):AREA, 
data = recruits1, link = identity)

#model with additional additive term for the smaller plots 
mRecruit1_SmlPlt <- glm.nb(NumRecruits ~ 
  0 + SDLAREA + I(SDLAREA * SmallPlot) +
      I(SDLAREA * CONES1),
                 data = recruits1, link = identity)

#Year 1 recruitment model comparison----------------------------------------------
#list of models
rec1_models <- list(mRecruit1_none=mRecruit1_none,
        mRecruit1_AREA = mRecruit1_AREA,
        mRecruit1_SmlPlt = mRecruit1_SmlPlt)

#checking for convergence issues
sapply(rec1_models, \(x) x$converged) #no issues

#checking residuals
residuals <- lapply(rec1_models, simulateResiduals)
for (m in names(residuals)) {
  plot(residuals[[m]], title = m, quantreg = T)
  dev.new()
} #1st two models show trends in the quantile plots, while the 3rd model has no issues

#comparing AICs
rec1_AICs <- sapply(rec1_models, AIC) 
print(rec1_AICs) #model with plot size intercept has lowest AIC 
            #(and is also the one with no issues in residuals)

#selecting final model (lowest AIC)
mRecruit1 <- rec1_models[[which.min(rec1_AICs)]]

#Year 1 survival model------------------------------------------
#Survival in year 1 of study, i.e. ~1 yr after fire
#modelled by looking at proportion of cohort 1 & 2 
# that survived upto census 3

#Fitting model
mSurv1 <- glmmTMB(alive3 ~ AREA + (1|PLOTCODE), 
  family = binomial,data = cohorts1_2)
summary(mSurv1)

#checking convergence
mSurv1$fit$convergence #no issues

#checking residuals
simulateResiduals(mSurv1, plot = T) #no major issues

# Plotting predictions
cohorts1_2$p_hat <- predict(mSurv1, type = "response")  #model predictions
survival_1_summary <- cohorts1_2|>
  group_by(AREA,PLOTCODE, p_hat)|>
  #calculating observed survival
  summarise(Survived = sum(alive3), n=n())|>
  mutate(p_obs = Survived/n)
plot(survival_1_summary$p_obs, survival_1_summary$p_hat,
    xlab ="Observed Survival", ylab = "Predicted Survival")
abline(a = 0, b = 1)

#Year 1 size distribution EDA--------------------------------------
#Exploratory analysis onthe size distribution of 
# cohorts 1 and 2 at census 3

#Proportion of seedlings in each SITEAREA
#  with height measurements
prop_meas <- cohorts1_2_alive3|>
  group_by(AREA, SITEAREA)|>
  summarise(NumAlive =n(), NumMeasured = sum(!is.na(HT)),
  FracMeasured = NumMeasured/NumAlive, 
  MeanHeight = mean(HT, na.rm=T),
  SDHeight =sd(HT, na.rm =T))
prop_meas 
#uneven coverage of height measurements within DW:
#DW3 has the highest number & fraction of seedlings measured
#and has a lower mean height than rest of the DW SITEAREAs
#Year 1 size distribution Model--------------------------------------
#Model fitting for size distribution of seedlings at ~1 year age
#i.e. size distribution of cohort 1 & 2 at census 3

#Filtering dataframe to keep only seedlings with height measurements
cohorts1_2_alive3_htmeas <- cohorts1_2_alive3|>
  filter(!is.na(HT))|>
  dplyr::select(IND_ID, AREA, SITEAREA, HT)

#Fitting model - 
# Using SITEAREA as the predictor, rather than AREA 
#to account for even data collection

#Model 1 - log normal distribution
#with different variances for each SITEAREA
mSize1_1 <- glmmTMB(log(HT)~SITEAREA, 
  dispformula =  ~SITEAREA, family = "gaussian", 
  data = cohorts1_2_alive3_htmeas)

#checking residuals
sim1 <-simulateResiduals(mSize1_1)
plot(sim1) #KS test is significant

# checking QQ plot
z <- residuals(sim1,
   quantileFunction = qnorm, outlierValues = c(-4, 4))
par(mfrow = c(1, 3))
for (a in levels(cohorts1_2_alive3_htmeas$AREA)) {
  zz <- z[cohorts1_2_alive3_htmeas$AREA == a]
  qqnorm(zz, main = a); qqline(zz)
} #QQ plots indicate left skew

#Model 2 - skewed log normal distribution
#with different variances for each SITEAREA
mSize1_2 <- gamlss(log(HT) ~ SITEAREA, 
  sigma.formula = ~ SITEAREA,
   nu.formula = ~ AREA,
             family = SN1,
            data = cohorts1_2_alive3_htmeas,
             control = 
              gamlss.control(n.cyc = 500, trace = FALSE))

mSize1_2$converged

#checking QQ plots
z2 <- resid(mSize1_2)
par(mfrow = c(1, 3))
for (a in levels(cohorts1_2_alive3_htmeas$AREA)) {
  zz <- z2[cohorts1_2_alive3_htmeas$AREA == a]
  qqnorm(zz, main = a); qqline(zz)
} 
#using a skew model did not seem to make a difference
#perhaps the model was not able to find the right parameter values

#checking model estimates of the skew parameter (nu)
unique(data.frame(AREA = cohorts1_2_alive3_htmeas$AREA, 
  nu = fitted(mSize1_2, "nu"))) 
#nu values are small and postive (should be negative for left skew)

#Model 3 -  re-trying the skew model with better starting values for nu

#starting values for nu
nu0_area <- c(DW = -2, SCC = -2, RP = -2)

#calculating starting values for the rest of the parameters
st <- cohorts1_2_alive3_htmeas |>
  group_by(SITEAREA) |>
  mutate(m = mean(log(HT)), s = sd(log(HT))) |>
  ungroup() |>
  mutate(nu0   = nu0_area[as.character(AREA)],
         delta = nu0 / sqrt(1 + nu0^2),
         sig0  = s / sqrt(1 - 2 * delta^2 / pi),
         mu0   = m - sig0 * delta * sqrt(2 / pi))

#fitting model
mSize1_3 <- gamlss(log(HT) ~ SITEAREA, 
            sigma.formula = ~ SITEAREA,
             nu.formula = ~ AREA,
             family = SN1, data = cohorts1_2_alive3_htmeas,
             mu.start = st$mu0, sigma.start = st$sig0, 
             nu.start = st$nu0,
             method = mixed(5, 200),
             control =
               gamlss.control(n.cyc = 500, trace = FALSE))

#checking for convergence
mSize1_3$converged

#checking nu estimates
coef(mSize1_3, what = "nu") #they are now negative

z3 <- resid(mSize1_3)
par(mfrow = c(1, 3))
for (a in levels(cohorts1_2_alive3_htmeas$AREA)) {
  zz <- z3[cohorts1_2_alive3_htmeas$AREA == a]
  qqnorm(zz, main = a); qqline(zz)
} #residuals quantiles now mostly align with the expected line


#Model comparison - AICs
AIC(mSize1_1)
GAIC(mSize1_2)
GAIC(mSize1_3) #lowest AIC

#final model
mSize1 <- mSize1_3

