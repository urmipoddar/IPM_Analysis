## ************************************************************************** ##
## Author: Matthew Aiello-Lammens
## Modified by: Urmi Poddar

## Purpose:
## Read in and clean pine demography data
##
## ************************************************************************** ##

#Loading data and packages-----------------------------------------

## Load packages 
library(tidyverse)
library(readxl)

## Read in the Excel document of all pitch pine census data.
## warnings removed by setting 'guess_max = 7000'
pines <- read_xlsx("data/pine-demography-data-all.xlsx",
                   sheet = "CensusData", na = "-", guess_max = 7000)

#Cleaning/formatting identifying info--------------------------------------------------------------------------
## Make a data set that is all of the pines before any further filtering
pines_all <- pines

## Remove entries that have values of NA for `PLOTCODE`
print(paste("Number of rows removed because PLOTCODE is NA:", length(which(is.na(pines$PLOTCODE)))))
pines <- pines %>% filter(!is.na(PLOTCODE))

## Make SITEAREA and AREA columns. This strips off the plot number from each `PLOTCODE`
pines$SITEAREA <- sub(pattern = "\\-.*", replacement = "", pines$PLOTCODE, perl = TRUE)
pines$AREA <- sub(pattern = "[0-9].*", replacement = "", pines$SITEAREA, perl = TRUE)

## Change the order of AREA factor
pines$AREA <- factor(pines$AREA, levels = c("DW", "SCC", "RP"))
pines$SITEAREA <- factor(pines$SITEAREA, 
                         levels = c("DW1", "DW2A", "DW2B",
                                    "DW2C", "DW3", "SCC5", "RP1", 
                                    "RP2"))

## Remove Pines with no documented COHORT number
print(paste("Number of rows removed because COHORT is NA:", length(which(is.na(pines$COHORT))) ))
pines <- pines %>% filter(!is.na(COHORT))

## Read in plot information

## NOTE: In order to get this plot information I had to export data from the 
## MS Access database to a csv file. This was done using the mdb-tools package in bash.
## ```
## mdb-export ./Pine_project/census\ 2003/pinedatams.mdb plotdatams > data/plot_information.csv
## ```

## Meta-data for these data can be found in the file:
## ./Pine_project/Documentation/Seedling files/plotdatahistory.doc

plot_information <- read.csv("data/plot_information.csv")

pines <-
  left_join(pines, dplyr::select(plot_information,
                                 OLDCODE, PLTAREA, SDLAREA, TMT),
            by = c("PLOTCODE" = "OLDCODE"))

## ADDED BY URMI (everything till the end of this chunk)

#Adding a unique id for each individual
pines <- pines|>
  mutate(IND_ID = 1:nrow(pines))

#reordering columns
pines <- pines|>
  relocate(IND_ID)

#removing columns that don't contain useful data
#i.e., individual cohort columns, notes columns & 
#columns for names of observer/data recorder
pines <- pines|>
  select(-matches("Cohort_|NOTES|REC|MEAS|DATE"))

#Cleaning census status (alive/dead)---------------------------------
##  Cleaning of 'Census Status' Notations

# Before cleaning, the columns labeled STAT[census number] represent the survival 
# status of each individual plant, with values of 0 = alive, 1 = dead, and 2 = missing. 
# These values require cleaning to address the following issues:
#   
# * When a plant is marked with a 1, subsequent status may be left blank, 
#   resulting in an NA value.
# * When a plant is marked as 2, but then refound in a later census, the value 
#   of 2 should be changed to 0 to reflect that the plant was in fact alive during 
#   that census.
# * When a plant is marked as 2 in two censuses in a row, we assume it is dead, 
#   and thus the value should be marked as 1 in the first census it was missing in.
# 
# **General STAT cleaning work flow**
#   
# 1. Separate out all STAT columns from the pines data.frame
# 2. Use a for loop to test multiple scenarios across the STAT data.frame
#   a. The order in which you test the scenarios matters
#   b. In most cases, the test condition is applied to a whole STAT vector, then 
#      the new STAT vector is created using min/max values associated with the 
#      test condition
# 3. Store the results of these scenario tests in a new STAT data.frame, this 
#   way you can overwrite those results
## Loop through the `pine_stat` data.frame line by line, and edit line as needed

pines_stat <- dplyr::select(pines, contains("STAT"))


#first removing unusual values from the STAT15 column
#and converting it to numeric
##ADDED BY URMI
pines_stat <- pines_stat|>
  mutate(STAT15 = case_when(
    STAT15 == "3"| STAT15=="?" ~ "2",
    .default = STAT15))|>
  mutate(STAT15 = as.numeric(STAT15))

pines_stat_new <- pines_stat

for(ll in 1:nrow(pines_stat)){
  
  # Get the stat line
  stat <- pines_stat[ll, ]
  
  # If a pine was observed dead at any point, then mark all
  # subsequent times as  This prevents non-fire resprouting,
  # but based on data exploration, only one pine showed
  # a reasonable possibility of being a resprout
  if(any(stat == 1, na.rm = TRUE)){
    death_census <- min(which(stat==1))
    if(death_census<15){
      stat[(death_census+1):length(stat)] <- 1
    }
  }
  
  # If a pine was observed as missing, then found alive, then 
  # mark the missing periods as alive
  if(any(stat == 2, na.rm = TRUE)){
    while(min(which(stat == 2)) < max(which(stat == 0))){
      stat[min(which(stat == 2))] <- 0
    }
  }
  
  # If after the above steps, there are still missing periods
  if(any(stat == 2, na.rm = TRUE)){
    
    # Missing then later observed as alive -
    # change all prior 2s to 0s
    if(max(which(stat == 2)) < max(which(stat == 0))){
      stat[min(which(stat == 2)):max(which(stat == 0))] <- 0
    }
  }
  
  if(any(stat == 2, na.rm = TRUE)){
    # Missing then never observed as alive after -
    # change everything after max 2 to 1
    if(max(which(stat == 2)) > max(which(stat == 0))){
      stat[max(which(stat == 2)):length(stat)] <- 1
    } else {
      stat <- rep(NA, length(stat))
    }
  }
  
  # At this point, we can assume any missings are actually mortalitys
  if(any(stat == 2, na.rm = TRUE)){
    stat[stat == 2] <- 1
  }
  
  # Redo the very first condition test, now that all 2s have been
  # accounted for
  if(any(stat == 1, na.rm = TRUE)){
    death_census <- min(which(stat==1))
    if(death_census<15){
      stat[(death_census+1):length(stat)] <- 1
    }
  }
  
  
  # Add line to pines_stat_new
  pines_stat_new[ll, ] <- stat
  print(ll) ##ADDED BY URMI
}

# Clean up some unnecessary variables
rm(stat,ll)

## ADDED BY URMI
#changing alive to 1 and dead to 0
pines_stat_new_recoded <- pines_stat_new|>
  mutate(across(STAT1:STAT15,
                ~case_when(.x==0 ~1,
                           .x == 1 ~0, 
                           .default = .x)))|>
  #renaming with the format STAT__census
  rename_with(
    ~ sub("^STAT([0-9]+)$",
          "STAT__\\1", .x),
    starts_with("STAT"))

## Update pines STAT columns
pines <- dplyr::select(pines, -contains("STAT"))
pines <- cbind(pines, pines_stat_new_recoded)

#Cleaning cone data-------------------------------------------------
## STANDARDIZING SEROTINY OBSERVATIONS

pines$SEROTINY15[which(pines$SEROTINY15 == "6")] <- NA
# The changes below are based on comparing the observations of CENSUS 14 and 15,
# and assuming that "n", "non-serot", and "ns" are equivelant
pines$SEROTINY14[which(pines$SEROTINY14 == "non-serot")] <- "n"
pines$SEROTINY14[which(pines$SEROTINY14 == "ns")] <- "n"

## Now convert CENSUS14 "n" and "s" to "O" and "S" to match CENSUS15
pines$SEROTINY14[which(pines$SEROTINY14 == "s")] <- "S"
pines$SEROTINY14[which(pines$SEROTINY14 == "n")] <- "O"

#ADDED BY URMI - checking if correct
pines$SEROTINY14|>unique()
pines$SEROTINY15|>unique()

#ADDED BY URMI - #renaming with the format SEROTINY__census
pines <- pines|>
  rename_with(
    ~ sub("^SEROTINY([0-9]+)$",
          "SEROTINY__\\1", .x),
    starts_with("SEROTINY"))

## STANDARDIZING CONE COLUMNS

#standardizing names of mature cone count columns
#MODIFIED BY URMI (more efficient code, slightly different col names)
pines <-
  pines|>
  rename(CONE_Mature__12 =`CONE12_>1YR`,
         CONE_Mature__13 = `CONE13_>1YR`)|>
  rowwise()|>
  mutate(CONE_Mature__14 = sum(CONE14_GOLD, CONE14_GRAY, na.rm = T), 
         CONE_Mature__15 = sum(CONE15_2YR, `CONE15_>2YR`, na.rm = T),
         .keep = "unused")|>
  ungroup()

#standardizing names of current year & 1 yr cone count columns
## ADDED BY URMI
pines <- pines|>
  rename_with(
    ~ sub("^CONE([0-9]+)(?:_BABY|_NEW)?$", "CONE_New__\\1", .x),
    starts_with("CONE"))|>
  rename_with(
    ~ sub("^CONE([0-9]+)(?:_1YR|_GREEN)$","CONE_1Yr__\\1", .x),
    starts_with("CONE"))

#ensuring all cone columns are numeric
pines|>
  select(contains("CONE"))|>
  summarise(across(everything(), ~is.numeric(.x)))|>
  as.data.frame() #CONE_New__12 is not numeric
pines$CONE_New__12|>unique() #checking why its not numeric
pines$CONE_New__12 <- as.numeric(pines$CONE_New__12)
pines|>
  select(contains("CONE"))|>
  summarise(across(everything(), ~is.numeric(.x)))|>
  as.data.frame()
  
# Set values to NA if the tree in the census year is not alive
## MODIFIED BY URMI (Now 1 = alive & 0 = dead, renamed columns)
pines$CONE_Mature__13[which(pines$STAT__13 != 1 & 
                            !is.na(pines$CONE_Mature__13))] <- NA
pines$CONE_Mature__14[which(pines$STAT__14 != 1 & 
                            !is.na(pines$CONE_Mature__14))] <- NA
pines$CONE_Mature__15[which(pines$STAT__15 != 1 &
                            !is.na(pines$CONE_Mature__15))] <- NA
all(is.na(pines$CONE_Mature__12[which(pines$STAT__12 != 1)]))
all(is.na(pines$CONE_Mature__13[which(pines$STAT__13 != 1)]))
all(is.na(pines$CONE_Mature__14[which(pines$STAT__14 != 1)]))
all(is.na(pines$CONE_Mature__15[which(pines$STAT__15 != 1)]))

# Set values to 0 if the status of the tree is alive in the census year, but
# the current cone count is NA
## MODIFIED BY URMI ( Now 1 = alive & 0 = dead)
pines$CONE_Mature__12[which(pines$STAT__12 == 1 &
                            is.na(pines$CONE_Mature__12))] <- 0
pines$CONE_Mature__13[which(pines$STAT__13 == 1 &
                            is.na(pines$CONE_Mature__13))] <- 0
pines$CONE_Mature__14[which(pines$STAT__14 == 1 &
                            is.na(pines$CONE_Mature__14))] <- 0
pines$CONE_Mature__15[which(pines$STAT__15 == 1 & 
                            is.na(pines$CONE_Mature__15))] <- 0

#Cleaning stem diameter data--------------------------------------------------
## THIS CODE CHUNK ADDED BY URMI

#Ensuring all diameter at ankle ht columns are numeric
pines|>
  select(matches("DIAM|DAH"))|>
  summarise(across(everything(),
                   ~is.numeric(.x)))#DAH15 is not numeric
pines$DAH15|>unique() #checking why it is not numeric
#there is a value entered as "1.36/ 1.72"
#Fixing this value and converting to numeric
pines$DAH15[which(pines$DAH15=="1.36/ 1.72" )] <- 
  as.character(round(1.36/1.72, digits = 2))
pines$DAH15 <- as.numeric(pines$DAH15)
pines|>
  select(matches("DIAM|DAH"))|>
  summarise(across(everything(),
                   ~is.numeric(.x)))#everything is numeric now

#standardizing names and units of diameter at ankle ht columns
#DIAM columns are same as DAH but are in mm instead of cm
#converting everything to cm and renaming DIAM to DAH
#also ensuring column names are in the format DAH__census
pines <- pines|>
  mutate(across(contains("DIAM"), ~.x/10))|>
  rename_with( ~ sub("^DIAM([0-9]+)$", "DAH__\\1", .x),
    starts_with("DIAM"))|>
  rename_with( ~ sub("^DAH([0-9]+)$", "DAH__\\1", .x),
               starts_with("DAH"))

#ensuring that if a tree is not alive, its DAH is NA
dah_cols  <- grep("^DAH__[0-9]+$", names(pines), value = TRUE)
stat_cols <- sub("^DAH__", "STAT__", dah_cols)

pines[dah_cols] <- Map(
  function(dah, stat) {ifelse(stat != 1, NA, dah)},
  pines[dah_cols],
  pines[stat_cols])

#Looking for rows where DAH at census t is less than at census t-1
#which represent errors
dah <- pines|>
  select(dah_cols[order(as.numeric(sub("^DAH__", "", dah_cols)))])
#dah <- dah|>select(!DAH9)
dah_flag <- apply(dah, 1, function(x) {
  any(diff(x) < -0.05, na.rm = TRUE)})
dah[dah_flag,]|>View()

#ensuring there are no zero values 
any(pines[dah_cols]==0, na.rm =T) #there is at least 1 zero
#changing zeros to NAs
pines[dah_cols] <- pines[dah_cols]|>
  mutate(across(everything(), 
                ~if_else(.x==0, NA, .x)))
any(pines[dah_cols]==0, na.rm =T) #all zeros removed

## Repeating the same steps for diameter at breast ht (DBH)

#first checking if all columns are numeric
pines|>
  select(matches("DBH"))|>
  summarise(across(everything(),
        ~is.numeric(.x)))
#DBH14 & 15 are not numeric
#checking why they are not numeric
pines$DBH14|>unique() 
pines$DBH15|>unique()
#in DBH15, there is a value entered as "6.2, 2.9, 2.1"
#based on this tree's DAH15, I assume that this is supposed to be 6.2
#in DBH14, there doesn't seem to be any obvious errors
#so I will simply coerce that column to numeric
pines$DBH14 <- as.numeric(pines$DBH14)
pines$DBH15[which(pines$DBH15=="6.2, 2.9, 2.1")] <- "6.2"
pines$DBH15 <- as.numeric(pines$DBH15)
pines|>
  select(matches("DBH"))|>
  summarise(across(everything(),
           ~is.numeric(.x))) #everything is numeric now

#checking units
pines|>
  select(matches("DBH"))|>
  summarise(across(everything(),
     ~median(.x, na.rm = T))) #everything seems to be in the same units

#renaming to follow the format DBH__census
pines <- pines|>
  rename_with( ~ sub("^DBH([0-9]+)$", "DBH__\\1", .x),
               starts_with("DBH"))

#ensuring that if a tree is not alive, its DBH is NA
dbh_cols  <- grep("^DBH__[0-9]+$", names(pines), value = TRUE)
stat_cols <- sub("^DBH__", "STAT__", dbh_cols)

pines[dbh_cols] <- Map(
  function(dbh, stat) {ifelse(stat != 1, NA, dbh)},
  pines[dbh_cols],
  pines[stat_cols])

#ensuring there are no zero values 
any(pines[dbh_cols]==0, na.rm =T) #there's at least 1 zero
#changing zeros to NAs
pines[dbh_cols] <- pines[dbh_cols]|>
  mutate(across(everything(), 
                ~if_else(.x==0, NA, .x)))
any(pines[dbh_cols]==0, na.rm =T) #all zeros are removed now

#Cleaning crown diameter data-----------------------------------
#THIS CODE CHUNK ADDED BY URMI

#Ensuring all crown diameter columns are numeric
pines|>
  select(matches("WD"))|>
  summarise(across(everything(),
                   ~is.numeric(.x)))#everything is numeric

#ensuring all columns are in the same units
pines|>
  select(matches("WD"))|>
  summarise(across(everything(),
            ~median(.x, na.rm = T)))#WD_12s seem to be in mm
#converting WD1_12 and WD2_12 to cm
pines <- pines|>
  mutate(across(c(WD1_12, WD2_12), 
                ~.x/10))

#renaming columns based on the format WD1__census 
pines <- pines|>
  rename_with(~ sub("^WD([1|2])_([0-9]+)$", "WD\\1__\\2", .x),
              starts_with("WD"))

#ensuring that if a tree is not alive, its WD is NA
wd_cols  <- grep("^WD[1|2]__[0-9]+$", names(pines), value = TRUE)
stat_cols <- sub("^WD[1|2]__", "STAT__",wd_cols)

pines[wd_cols] <- Map(
  function(wd, stat) {ifelse(stat != 1, NA, wd)},
  pines[wd_cols],
  pines[stat_cols])

#checking for zeros
any(pines[wd_cols]==0, na.rm = T) #there are at least 1 zero
#changing zeros to NA
pines[wd_cols] <- pines[wd_cols]|>
  mutate(across(everything(), 
                ~if_else(.x ==0, NA, .x)))
any(pines[wd_cols]==0, na.rm = T) #all zeros removed now

#Cleaning height data--------------------------------------------------
## THIS CODE CHUNK ADDED BY URMI

#First I replace the HTB15 column with 15W/subs
#The HTB15 column seems incorrect and 15W/subs seems to be the corrected version
pines <- pines|>
  select(!HTB15)|>
  rename(HTB15 = `15W/subs`)

#comparing HTN vs HTB (both were measured in census 9)
ggplot(pines, aes(HTB9, HTN9))+
  geom_abline(slope=1, intercept = 0, lwd = 1, col="red")+
  geom_point(pch=21, fill ="grey", alpha= 0.5)+
  theme_bw()

lm(HTN9~0+HTB9, data = pines)|>summary()

#removing the HTN9 column (since there is also a HTB9 column)
pines <- pines|>
  select(!HTN9)

#Ensuring ht columns are numeric
pines|>
  select(matches("HT"))|>
  summarise(across(everything(),
                   ~is.numeric(.x)))#HTB12  is not numeric
pines$HTB12|>unique() #checking why it is not numeric
#there is a value entered as " "28..4""
#Fixing this value and converting to numeric
pines$HTB12[which(pines$HTB12== "28..4" )] <- 
  as.character(28.4)
pines$HTB12 <- as.numeric(pines$HTB12)
pines|>
  select(matches("HT"))|>
  summarise(across(everything(),
                   ~is.numeric(.x)))#everything is numeric now

#checking units
pines|>
  select(matches("HT"))|>
  summarise(across(everything(),
     ~median(.x, na.rm=T)))#looks like they're in the same units

#standardizing names to HT__census format
pines <- pines|>
  rename_with(~ sub("^HTB([0-9]+)$", "HT__\\1", .x),
    starts_with("HT"))|>
  rename_with(~ sub("^HT([0-9]+)$", "HT__\\1", .x),
              starts_with("HT"))

#ensuring that if a tree is not alive, its height is NA
ht_cols  <- grep("^HT__[0-9]+$", names(pines), value = TRUE)
stat_cols <- sub("^HT__", "STAT__", ht_cols)

pines[ht_cols] <- Map(
  function(ht, stat) {ifelse(stat != 1, NA, ht)},
  pines[ht_cols],
  pines[stat_cols])

#Looking for rows where height at census t is less than at census t-1
#which likely represent errors
ht <- pines|>
  select(ht_cols[order(as.numeric(sub("^HT__", "", ht_cols)))])

ht_flagged <- apply(ht, 1, function(x) {
  any(diff(x) < -5, na.rm = TRUE)})
ht[ht_flagged,]|>View()

#ensuring there are no zero values 
any(pines[ht_cols]==0, na.rm =T) #there are no zeros

#Converting to long format & Adding census dates----------------------------------------------
##THIS CODE CHUNK ADDED BY URMI

#loading census dates
dates <- read_xlsx("data/pine-demography-data-all.xlsx",
                   sheet = "CensusDates", skip = 2)

#tranforming to long format
#and only keeping the start date of each census
census_dates <- dates|>
  filter(`...1` == "Start on")|>
  pivot_longer( cols = -`...1`,
    names_to = "census",
    values_to = "START_DATE")|>
  mutate(CENSUS_NUM = row_number(),
    CENSUS_MONTH = month(START_DATE),
    CENSUS_YEAR = year(START_DATE))|>
  select(CENSUS_NUM, START_DATE, CENSUS_MONTH, CENSUS_YEAR)

#converting pine data to long form 
#and only retaining site info, size & demographic data columns
pines_long <- pines|>
  select(IND_ID:TAG,SITEAREA:TMT,
         matches("STAT"),
         matches("HT|DAH|DBH|WD|CONE|SEROTINY"))|>
  pivot_longer(cols = 
                 matches("STAT|HT|DAH|DBH|WD|CONE|SEROTINY"),
               names_to = c(".value", "CENSUS_NUM"),
               names_sep="__")|>
  mutate(CENSUS_NUM = as.numeric(CENSUS_NUM))

#adding census dates to pine data
pines_long <- pines_long|>
  left_join(census_dates)

#Exporting to csv-------------------------------------------
## THIS CODE CHUNK ADDED BY URMI

#exporting in wide format
write.csv(pines,"Data/pine-demography-data-all_cleaned.csv",
          row.names = F, na = "")

#exporting in long form - site/census info, size & demographic data only
write.csv(pines_long,"Data/pine_demography_cleaned_long.csv",
          row.names = F, na ="")
