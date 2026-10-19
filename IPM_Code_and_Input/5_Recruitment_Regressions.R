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
  select(OLDCODE, CONES1, PLTAREA)|>
  rename(PLOTCODE = OLDCODE)

#Calculating number of new recruits per census
recruits <- pines_long|>
  group_by(PLOTCODE, AREA,
     SITEAREA, CENSUS_NUM, START_DATE, SDLAREA)|>
  summarise(NumRecruits = sum(COHORT==CENSUS_NUM))

#recruitment in first year (census 1 and 2)
recruits1 = recruits|>
  filter(CENSUS_NUM %in% c(1,2))|>
  group_by(PLOTCODE, AREA, SITEAREA, SDLAREA)|>
  #summing up the number of recruits in census 1 and 2
  summarise(NumRecruits = sum(NumRecruits))|>
  #joining with cones after fire dataframe
  left_join(cones_after_fire)|>
  ungroup()
recruits1$NumRecruits
#Fitting first year recruitment model----------------------------------------------
#First year recruitment:
#modelled by regressing the number of cones found 
#and number of seedlings
#after controlling for plotsize (SDLAREA)
mRecruit1 <- glm.nb(NumRecruits ~
    0 + SDLAREA + I(SDLAREA * CONES1), 
data = recruits1, link = identity)
mRecruit1

mRecruit1_2 <- glm.nb(NumRecruits ~
    0 + SDLAREA + I(SDLAREA * CONES1):AREA, 
data = recruits1, link = identity)

AIC(mRecruit1, mRecruit1_2)
mRecruit1_2 <- update(mRecruit1,
   data = subset(recruits1, SDLAREA == 25))
res <- simulateResiduals(mRecruit1_2, plot=T)
hist(res)
summary(mRecruit1_2)
