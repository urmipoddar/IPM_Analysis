## ************************************************************************** ##
## Author: Urmi Poddar

## Purpose:
## Exploratory data analysis on pine demographic data
##
## ************************************************************************** ##
#Loading data and packages-----------------------------------------
library(tidyverse)
library(ggExtra)
library(ggrepel)
library(ggpubr)

pines_long <- read.csv("Data/pine_demography_cleaned_long.csv")
plot_info <- read.csv("Data/plot_information.csv")

#Data formatting---------------------------------------------------
#adding plot information
plot_info <- plot_info|>
  select(!c(SUBPLOT, AREA, PLTAREA, SDLAREA, TMT))

#Converting census start date column to Date format
#(useful for plotting)
pines_long <- pines_long|>
  mutate(START_DATE = as.Date(START_DATE))

#Converting area and site area to factors
#(useful for plotting)
pines_long$AREA <- factor(pines_long$AREA, 
                          levels = c("DW", "SCC", "RP"))
pines_long$SITEAREA <- factor(pines_long$SITEAREA, 
                         levels = c("DW1", "DW2A", "DW2B",
                                    "DW2C", "DW3", "SCC5", "RP1", 
                                    "RP2"))
#Calculating sample size and data availability------------------------------------
#sample size = number of tagged individuals alive at each census

#Calculating the number of individuals alive
#in each census x site area
num_alive <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1)|>
  summarise(NumAlive = n_distinct(IND_ID))

#repeating the same calculation but for cohort 1 only
num_alive_c1 <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1 & COHORT==1)|>
  summarise(NumAlive = n_distinct(IND_ID), n=n())

#counting number of individuals for which size  data was collected
#in each census year and site
data_avail <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1 )|>
  summarise(NumAlive = n_distinct(IND_ID), 
            NumHT = sum(!is.na(HT)),
            NumDAH = sum(!is.na(DAH)),
            NumDBH =sum(!is.na(DBH)),
            NumWD1 = sum(!is.na(WD1)),
            NumWD2 = sum(!is.na(WD2)))|>
  mutate(across(NumHT:NumWD2, 
                ~round(.x/NumAlive, digits = 2),
                .names = "{str_replace(.col, 'Num', 'Frac')}"))

#repeating the same calculation but for cohort 1 only
data_avail_c1 <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1 & COHORT==1 )|>
  summarise(NumAlive = n_distinct(IND_ID), 
            NumHT = sum(!is.na(HT)),
            NumDAH = sum(!is.na(DAH)),
            NumDBH =sum(!is.na(DBH)),
            NumWD1 = sum(!is.na(WD1)),
            NumWD2 = sum(!is.na(WD2)))|>
  mutate(across(NumHT:NumWD2, 
                ~round(.x/NumAlive, digits = 2),
                .names = "{str_replace(.col, 'Num', 'Frac')}"))


#counting number of individuals for which cone data was collected
#in each census year and site
conedata_avail <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1 )|>
  summarise(NumAlive = n_distinct(IND_ID), 
            NumCone_1Yr = sum(!is.na(CONE_1Yr)),
            NumCone_Mature = sum(!is.na(CONE_Mature)),
            NumCone_New =sum(!is.na(CONE_New)))|>
  mutate(across(NumCone_1Yr:NumCone_New, 
                ~round(.x/NumAlive, digits = 2),
                .names = "{str_replace(.col, 'Num', 'Frac')}"))

#Repeating same calculation for cohort 1 only
conedata_avail_c1 <- pines_long|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  filter(STAT == 1 & COHORT==1)|>
  summarise(NumAlive = n_distinct(IND_ID), 
            NumCone_1Yr = sum(!is.na(CONE_1Yr)),
            NumCone_Mature = sum(!is.na(CONE_Mature)),
            NumCone_New =sum(!is.na(CONE_New)))|>
  mutate(across(NumCone_1Yr:NumCone_New, 
                ~round(.x/NumAlive, digits = 2),
                .names = "{str_replace(.col, 'Num', 'Frac')}"))

#Calculating survival and recruitment-------------------------------------
#Caclulating the number of trees dying b/w two censuses
num_dead = pines_long|>
  arrange(CENSUS_NUM, IND_ID)|>
  group_by(IND_ID)|>
  #creating a column for status in previous census
  mutate(STAT_prev = lag(STAT))|>
  ungroup()|>
  filter(CENSUS_NUM>1)|>
  group_by(CENSUS_NUM, START_DATE, AREA, SITEAREA)|>
  #counting number of trees alive in previous census and dead in current census
  filter(STAT_prev==1 & STAT == 0)|>
  summarise(NumDead = n_distinct(IND_ID))

#Calculating the proportion of population dying b/w two censuses
num_dead <- num_dead|>
  #joining with number alive in each census
  right_join(num_alive)|> 
  #creating column for number alive in previous census
  arrange(AREA, SITEAREA, CENSUS_NUM)|>
  group_by(SITEAREA)|>
  mutate(NumAlive_prev = lag(NumAlive))|>
  ungroup()|>
  #proportion of population that died
  mutate(PropDead = NumDead/NumAlive_prev)|>
  filter(CENSUS_NUM>1)|>
  #changing NAs to zero
  mutate(PropDead = if_else(is.na(PropDead), 0, PropDead))

#Calculating number of new recruits per census
recruits <- pines_long|>
  #filter(TMT=="NONE")|>
  group_by(PLOTCODE, AREA, SITEAREA, CENSUS_NUM, START_DATE)|>
  summarise(NumRecruits = sum(COHORT==CENSUS_NUM))|>
  left_join(plot_info,
            by = c("PLOTCODE" = "OLDCODE"))

#Graphs of sample size and data availability--------------------------------------------
#Number of individuals alive at each census
p_alive <- ggplot(num_alive, aes(START_DATE,NumAlive, 
                      group=SITEAREA, fill = AREA))+
  geom_point(pch = 21, size = 2, alpha = 0.8)+
  geom_line(aes(col = AREA))+
  theme_classic()+
  xlab("Census date")+ylab("Number alive")
p_alive

#Number of individuals alive at each census
#for cohort 1 only
p_alive_c1 <- ggplot(num_alive_c1, aes(START_DATE,NumAlive, 
            group=SITEAREA, fill = AREA))+
  geom_point(pch = 21, size = 2, alpha= 0.8)+
  geom_line(aes(col = AREA))+ylim(0, 1500)+
  theme_classic()+ggtitle("Cohort 1 only")+
  xlab("Census start date")+ylab("Number alive")+
  theme(plot.title = element_text(hjust = 0.5))
p_alive_c1

p_alive_multiplot <- 
  num_alive|>
  ggplot( aes(START_DATE,NumAlive, 
            group=SITEAREA, fill =AREA))+
  geom_point(pch = 21, size = 2, alpha= 0.8)+
  geom_line(aes(col = AREA))+
  theme_bw()+
  xlab("Census start date")+ylab("Number alive")+
  facet_wrap(vars(AREA))+ylim(0, 1500)+
  theme(plot.title = element_text(hjust = 0.5))
p_alive_multiplot

#% of individuals whose size was measured at each census
# for each size metric
p_size_meas <- data_avail_c1|>
  group_by(CENSUS_NUM, START_DATE)|>
  #calculating mean availability for each census
  summarise(across(starts_with("Frac"), ~mean(.x,na.rm=T)))|>
  pivot_longer(cols = starts_with("Frac"))|>
  #plotting
  ggplot(aes(factor(CENSUS_NUM), value*100,
             fill = name))+
  geom_bar(stat="identity", position="dodge", color="black")+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_discrete(name = "Size variable", 
              labels = c("DAH", "DBH","Height",
          "Canopy Maj. Axis", "Canopy Min. Axis"))+
  xlab("Census Number")+ylab("% of Individuals Measured")+
  ggtitle("Cohort 1 only")+
  theme(plot.title = element_text(hjust = 0.5))
p_size_meas

#similar to previous but with seperate panels for each site
p_size_meas_site <- data_avail_c1|>
  group_by(CENSUS_NUM, START_DATE, AREA)|>
  #mean availability for each census
  summarise(across(starts_with("Frac"), ~mean(.x,na.rm=T)))|>
  pivot_longer(cols = starts_with("Frac"))|>
  #plotting
  ggplot(aes(factor(CENSUS_NUM), value*100,
             fill = name))+
  geom_bar(stat="identity", position="dodge", color="black")+
  theme_bw()+theme(legend.position = "top")+
   scale_fill_discrete(name = "Size variable", 
              labels = c("DAH", "DBH","Height",
          "Canopy Maj. Axis", "Canopy Min. Axis"))+
  xlab("Census Number")+ylab("% of Individuals Measured")+
  facet_grid(vars(AREA))
p_size_meas_site


#% of individuals whose cone production was measured at each census
# for each cone metric
p_cone_meas <- conedata_avail|>
  group_by(CENSUS_NUM, START_DATE)|>
  #calculating mean availability for each census
  summarise(across(starts_with("Frac"), ~mean(.x,na.rm=T)))|>
  pivot_longer(cols = starts_with("Frac"))|>
  #plotting
  ggplot(aes(factor(CENSUS_NUM), value*100,
             fill = name))+
  geom_bar(stat="identity", position="dodge", color="black")+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_discrete(name = "",
                      labels = c("1 Yr old cones",
                                 "New cones", 
                                 "Mature cones"))+
  xlab("Census Number")+ylab("% of Individuals Measured")+
  theme(plot.title = element_text(hjust = 0.5))
p_cone_meas

#similar to previous but with seperate panels for each site
p_cone_meas_site <- conedata_avail|>
  group_by(CENSUS_NUM, START_DATE, AREA)|>
  #mean availability for each census
  summarise(across(starts_with("Frac"), ~mean(.x,na.rm=T)))|>
  pivot_longer(cols = starts_with("Frac"))|>
  #plotting
  ggplot(aes(factor(CENSUS_NUM), value*100,
             fill = name))+
  geom_bar(stat="identity", position="dodge", color="black")+
  theme_bw()+theme(legend.position = "top")+
   scale_fill_discrete(name = "",
                      labels = c("1 Yr old cones",
                       "New cones","Mature cones"))+
  xlab("Census Number")+ylab("% of Individuals Measured")+
  facet_grid(vars(AREA))
p_cone_meas_site

#Graphs of size distributions-----------------------------
##HEIGHT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#histograms of height over time
pines_long|>
  filter(COHORT==1)|>
  filter(!is.na(HT))|>
  ggplot(aes(HT, fill=AREA))+
  geom_histogram(color = "black", alpha= 0.8,
                 aes(y = after_stat(density)))+
  facet_grid(vars(CENSUS_NUM))+theme_bw()

#boxplots of height over time
#by census number
p_ht_census <- pines_long|>
  filter(COHORT==1)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), HT, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+
  xlab("Census number")+ylab("Height (cm)")
p_ht_census

#boxplots of height over time
#by census date
p_ht_date <- pines_long|>
  filter(COHORT==1)|>
  group_by(CENSUS_NUM, AREA, START_DATE)|>
  summarise(MeanHt = mean(HT, na.rm=T),SDHt = sd(HT, na.rm=T))|>
  ggplot(aes(START_DATE, MeanHt, group=AREA, col=AREA))+
  geom_point(size = 2)+geom_line(lwd = 1)+
  geom_errorbar(aes(ymin=MeanHt-SDHt, 
                    ymax=MeanHt+SDHt),
                alpha=0.5, lwd=0.7)+
  theme_bw()+xlab("Census date")+ylab("Height (cm)")
p_ht_date

#DAH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#histograms of DAH over time
pines_long|>
  filter(COHORT==1)|>
  filter(!is.na(DAH))|>
  ggplot(aes(DAH, fill=AREA))+
  geom_histogram(color = "black", alpha= 0.8,
                 aes(y = after_stat(density)))+
  facet_grid(vars(CENSUS_NUM))+theme_bw()

#boxplot of DAH over time
#by census number
p_dah_census <- pines_long|>
  filter(COHORT==1 & CENSUS_NUM>8)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", 
        remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), DAH,
             group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+
  xlab("Census date")+ylab("DAH (cm)")
p_dah_census

#boxplots of DAH over time
#by census number & only keeping censuses where DAH was measured
p_dah_census_v2 <- pines_long|>
  filter(COHORT==1)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), DAH, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+
  xlab("Census date")+ylab("DAH (cm)")
p_dah_census_v2

#Canopy width~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#boxplot of canopy major axis (WD1) over time
#by census number
p_wd1_census <- pines_long|>
  filter(COHORT==1 )|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), WD1, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+
  xlab("Census date")+ylab("Canopy Major Axis Length (cm)")
p_wd1_census

#boxplot of canopy minor axis (WD1) over time
#by census number
p_wd2_census <- pines_long|>
  filter(COHORT==1 )|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), WD2, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+
  xlab("Census date")+ylab("Canopy Minor Axis Length (cm)")
p_wd2_census

#boxplot of canopy area over time
#by census number
p_CA_census <- pines_long|>
  filter(COHORT==1 )|>
  unite("CENSUS_AREA",CENSUS_NUM,
        AREA, sep="_", remove = F)|>
  mutate(CANOPY_AREA = pi*(WD1/2)*(WD2/2))|>
  ggplot(aes(factor(CENSUS_NUM), CANOPY_AREA,
             group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, 
               outlier.size = 1)+
  theme_bw()+xlab("Census date")+
    ylab("Canopy Area (cm^2)")+
  ylim(0,20000 )
p_CA_census

#histogram of canopy area over time
pines_long|>
  filter(COHORT==1)|>
  mutate(CANOPY_AREA = pi*(WD1/2)*(WD2/2))|>
  filter(!is.na(CANOPY_AREA))|>
  ggplot(aes(CANOPY_AREA, fill=AREA))+
  geom_histogram(color = "black", alpha= 0.8,
                 aes(y = after_stat(density)))+
  facet_grid(vars(CENSUS_NUM))+theme_bw()

#Graphs of relationship b/w size variables---------------------
#Height vs DAH, separated by area
p_ht_dah <- pines_long|>
  filter(COHORT==1)|>
  ggplot(aes(HT, DAH, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("Height (cm)")+ylab("DAH (cm)")+
  theme(legend.position =  "top")
p_ht_dah

ggMarginal(p_ht_dah, type = "histogram")

cor.test(pines_long$HT, pines_long$DAH)
summary(lm(DAH~0+HT, data = pines_long))
summary(lm(DAH~0+HT*AREA, data = pines_long))

#Height vs DAH, separated by census
p_ht_dah_census <- pines_long|>
  filter(COHORT==1 & !is.na(DAH))|>
  ggplot(aes(HT, DAH, fill=AREA))+facet_wrap(vars(CENSUS_NUM))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_bw()+xlab("Height (cm)")+ylab("DAH (cm)")
p_ht_dah_census
summary(lm(DAH~0+HT*START_DATE, data = pines_long))
anova(lm(DAH~0+HT*START_DATE, data = pines_long))

#Height vs canopy major axis, separated by area
p_ht_wd1 <- pines_long|>
  filter(COHORT==1)|>
  ggplot(aes(HT, WD1, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("Height (cm)")+
  ylab("Canopy Major Axis (cm)")+
  theme(legend.position =  "top")+xlim(0, 400)
p_ht_wd1

#Height vs canopy diameter, separated by area
#canopy diameter = geometric mean of major and minor axis
p_ht_CD <- pines_long|>
  filter(COHORT==1)|>
  mutate(CANOPY_DIAM = sqrt(WD1*WD2))|>
  ggplot(aes(HT, CANOPY_DIAM, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("Height (cm)")+
  ylab("Canopy Diamter (cm)")+
  theme(legend.position =  "top")+xlim(0, 400)
p_ht_CD

ggMarginal(p_ht_CD, type = "histogram")
pines_long|>
  filter(COHORT==1)|>  
  mutate(CANOPY_DIAM = sqrt(WD1*WD2))%>%
  lm(CANOPY_DIAM~0+HT, data=.)|>
  summary()

#Height vs canopy area, separated by area
p_ht_CA <- pines_long|>
  filter(COHORT==1)|>
  mutate(CANOPY_AREA = pi*(WD1/2)*(WD2/2))|>
  ggplot(aes((HT), CANOPY_AREA, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("Height (cm)")+
  ylab("Canopy Area (cm^2)")+
  theme(legend.position =  "top")+xlim(0, 400)
p_ht_CA

ggMarginal(p_ht_CA, type = "histogram")

#DAH vs canopy major axis, separated by area
p_dah_wd1 <- pines_long|>
  filter(COHORT==1)|>
  ggplot(aes(DAH, WD1, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("DAH (cm)")+
  ylab("Canopy Major Axis (cm^2)")+
  theme(legend.position =  "top")+
  xlim(0, 11)
p_dah_wd1

#DAH vs canopy minor axis, separated by area
p_dah_wd2 <- pines_long|>
  filter(COHORT==1)|>
  ggplot(aes(DAH, WD2, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("DAH (cm)")+
  ylab("Canopy Minor Axis (cm^2)")+
  theme(legend.position =  "top")+
  xlim(0, 11)
p_dah_wd2

#DAH vs canopy diameter, separated by area
#canopy diameter = geometric mean of major and minor axis
p_dah_CD <- pines_long|>
  filter(COHORT==1)|>
  mutate(CANOPY_DIAM = sqrt(WD1*WD2))|>
  ggplot(aes(DAH, CANOPY_DIAM, fill = AREA))+
  geom_point(pch = 21, alpha= 0.5, size = 2)+
  theme_classic()+xlab("DAH (cm)")+
  ylab("Canopy Diameter (cm)")+
  theme(legend.position =  "top")+
  xlim(0, 10)
p_dah_CD
ggMarginal(p_dah_CD, type = "histogram")
pines_long|>
  filter(COHORT==1)|>
  mutate(CANOPY_DIAM = sqrt(WD1*WD2))%>%
  lm(CANOPY_DIAM~0+DAH, data = .)|>
  summary()

#PCA of size variables--------------------------------------
size_vars <- pines_long|>
  filter(COHORT==1)|>
  #adding rownames - individual ID + census number
  unite("ID_CENSUS", IND_ID, CENSUS_NUM, sep="_", remove = F)|>
  column_to_rownames(var ="ID_CENSUS")|>
  #selecting size columns
  select(HT, DAH, WD1, WD2)|>
  na.omit()

#PCA
pca_size_vars = prcomp(size_vars, scale. = T)
summary(pca_size_vars)
(pca_size_vars$sdev)^2 / sum((pca_size_vars$sdev)^2) * 100

#Extracting scores
pca_scores <- pca_size_vars$x[, 1:2]|>
  as.data.frame()|>
  rownames_to_column(var = "ID_CENSUS")|>
  separate_wider_delim(ID_CENSUS, delim ="_",
      names = c("IND_ID", "CENSUS_NUM"))|>
  mutate(IND_ID = as.numeric(IND_ID),
         CENSUS_NUM = as.numeric(CENSUS_NUM))
pca_rotations <- pca_size_vars$rotation[,1:2]|>
  as.data.frame()|>
  rownames_to_column(var = "Variable")

#Plotting - biplot
size_pca_biplot <- ggplot(pca_scores, aes(PC1, PC2))+
  geom_hline(yintercept = 0, lty="dashed", col = "grey40")+
  geom_vline(xintercept =  0,lty ="dashed", col = "grey40")+
  geom_point(pch=21, alpha=0.3, fill="skyblue")+
  geom_segment(data = pca_rotations,
               col="darkred", lwd=0.7,
               arrow = arrow(length = unit(0.2, "cm")),
               aes(x=0, y=0, xend=PC1*5, yend = PC2*5))+
  geom_label_repel(data = pca_rotations, col="darkred",
               aes(x=PC1*5, y = PC2*5.5, 
                   label=Variable), alpha = 0.4)+
  theme_bw()+coord_fixed()+
  xlab("PC1 (88.2%)")+ylab("PC2 (6.4%)")
size_pca_biplot

#PCA scores vs site
pca_site <- pca_scores|>
  left_join(pines_long)|>
  ggplot( aes(PC1, PC2, fill = AREA))+
  geom_hline(yintercept = 0, lty="dashed", col = "grey40")+
  geom_vline(xintercept =  0,lty ="dashed", col = "grey40")+
  geom_point(pch=21, alpha=0.5)+
  theme_bw()+coord_fixed()+
  xlab("PC1 (88.2%)")+ylab("PC2 (6.4%)")
pca_site

#PCA scores vs census number
pca_scores|>
  left_join(pines_long)|>
  ggplot( aes(factor(CENSUS_NUM),PC1, 
              fill = factor(AREA)))+
  geom_boxplot()+theme_bw()
pca_scores|>
  left_join(pines_long)|>
  ggplot( aes(factor(CENSUS_NUM),PC2, fill = factor(AREA)))+
  geom_boxplot()+theme_bw()

#Graphs of cone production-------------------------
#cone production over time - current yr cones
p_ConeNew_census <- pines_long|>
  #filter(COHORT==1)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), CONE_New, group=CENSUS_AREA,
             fill = AREA))+ylim(0,30)+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+xlab("Census number")+
  ylab("Num. of New Cones")
p_ConeNew_census

#cone production over time - 1 yr old cones
p_Cone1yr_census <- pines_long|>
  #filter(COHORT==1)|>
  #filter(CENSUS_NUM>10)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), CONE_1Yr, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+xlab("Census number")+ylab("Num. of 1 Yr Old Cones")
p_Cone1yr_census

#cone production over time - mature cones
p_ConeMature_census <- pines_long|>
  #filter(COHORT==1)|>
  #filter(CENSUS_NUM>10)|>
  unite("CENSUS_AREA",CENSUS_NUM, AREA, sep="_", remove = F)|>
  ggplot(aes(factor(CENSUS_NUM), CONE_Mature, group=CENSUS_AREA,
             fill = AREA))+
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1)+
  theme_bw()+xlab("Census number")+ylab("Num. of Mature Cones")
p_ConeMature_census

#Cone production vs height - current year cones
ggplot(pines_long, aes(HT, CONE_New,
                       fill = AREA))+
  geom_point(pch = 21)+theme_bw()+
  xlab("Height")+ylab("New Cone Count")+
  facet_wrap(vars(AREA))
pines_long|>
  mutate(CONE_New = as.integer(CONE_New))%>%
glm(CONE_New~HT*AREA, data = ., 
    family ="poisson")|>summary()

#Cone production vs height - 1 year cones
ggplot(pines_long, aes(HT, CONE_1Yr,
                       fill = AREA))+
  geom_point(pch = 21)+theme_bw()+
  xlab("Height")+ylab("1 Year Old Cone Count")+
  facet_wrap(vars(AREA))

#Cone production vs height - mature cones
ggplot(pines_long, aes(HT, CONE_Mature,
                       fill = AREA))+
  geom_point(pch = 21)+theme_bw()+
  xlab("Height")+ylab("Mature Cone Count")+
  facet_wrap(vars(AREA))
glm(CONE_Mature ~ HT*AREA, data = pines_long,
    family = "poisson")|>summary()

#Cone production vs DAH - mature cones
ggplot(pines_long, aes(DAH, CONE_Mature,
                       fill = AREA))+
  geom_point(pch = 21)+theme_bw()+
  xlab("DAH")+ylab("Mature Cone Count")+
  facet_wrap(vars(AREA))

#relationship b/w Num. of curr. yr cones & Num. of mature cones in next census
pines_long|>
  filter(COHORT==1 & CENSUS_NUM>10)|>
  arrange(IND_ID, CENSUS_NUM)|>
  group_by(IND_ID)|>
  mutate(CONE_Mature_next = lead(CONE_Mature))|>
  ggplot(aes( CONE_New, CONE_Mature_next))+
  geom_abline(slope = 1, intercept =0, col = "darkblue", alpha=0.5)+
  geom_point(pch = 21, fill ="grey", size =2, alpha=0.5)+
  xlab("Number of Current Yr Cones")+
  ylab("Number of Mature Cones in next census")+
  theme_bw()

# pines_long|>
#   filter(COHORT==1 & CENSUS_NUM>10)|>
#   arrange(IND_ID, CENSUS_NUM)|>
#   group_by(IND_ID)|>
#   mutate(HT_next = lead(HT))|>
#   ggplot(aes(HT,HT_next))+
#   geom_abline(slope = 1, intercept =0, col = "darkblue", alpha=0.5)+
#   geom_point(pch = 21, fill ="grey", size =2, alpha=0.5)+
#   theme_bw()

#Serotiny
pines_long|>
  filter(SEROTINY !="")|>
  group_by(AREA, CENSUS_NUM, SEROTINY)|>
  summarise(NumInds = n(), NumCones = sum(CONE_Mature))|>
  ggplot(aes(factor(CENSUS_NUM), NumInds, fill = SEROTINY))+
  geom_bar(stat = "identity")+
  facet_grid(vars(AREA))+theme_bw()+
  xlab("Census Number")+ylab("Number of Individuals")+
  scale_fill_discrete(labels =c("Mixed", "Open", "Serotinous"))

#Graphs of survival and recruitment------------------------------------
#Proportion dying between census
p_mort <-num_dead|>
  arrange(AREA, SITEAREA, CENSUS_NUM)|>
  ggplot(aes(CENSUS_NUM, PropDead, fill=AREA, 
        group=SITEAREA))+
  geom_point(pch = 21, alpha = 0.5)+theme_bw()+
  geom_line(alpha=0.5, aes(col=AREA))+
  facet_wrap(~AREA,nrow=3)+
  xlab("Census Number")+
  ylab("Proportion of population that died b/w censuses")
p_mort

#Proportion of cohort 1 surviving over time
num_alive_c1|>
  ungroup()|>
  arrange(AREA, SITEAREA, CENSUS_NUM)|>
  group_by(AREA, SITEAREA)|>
  mutate(PropSurv = NumAlive/first(NumAlive))|>
  arrange(AREA, SITEAREA, CENSUS_NUM)|>
  ggplot(aes(CENSUS_NUM, PropSurv, fill=SITEAREA, 
             group=SITEAREA))+
  geom_point(pch = 21, size = 3,alpha = 0.5)+theme_bw()+
  geom_line(alpha=0.5, aes(col=SITEAREA))+
  #facet_wrap(~AREA,nrow=3)+
  xlab("Census Number")+
  ylab("Proportion of cohort 1 surviving")


#Number of new recruits per census
p_recruits <- ggplot(recruits, 
          aes(CENSUS_NUM, NumRecruits,
          group=PLOTCODE, fill = AREA ))+
  geom_line(alpha = 0.2)+
  geom_point(pch=21, alpha = 0.5)+
  theme_bw()+facet_wrap(vars(AREA))+
  xlab("Census Number")+ylab("Number of Recruits")
p_recruits

#Plotting number of new recruits per census
#for census 2 onwards
p_recruits_2plus_censusnum <- recruits|>
  filter(CENSUS_NUM>1)|>
  ggplot(aes(CENSUS_NUM, NumRecruits,
             group=PLOTCODE,fill = AREA ))+
  geom_point(pch=21, alpha = 0.4, size =2)+
  theme_bw()+facet_wrap(vars(AREA))+
  ylim(0,25)+xlab("Census Number")+
  ylab("Number of Recruits")
p_recruits_2plus_censusnum

p_recruits_2plus_date <- recruits|>
  filter(CENSUS_NUM>1)|>
  ggplot(aes(START_DATE, NumRecruits,
             group=PLOTCODE,fill = AREA ))+
  geom_point(pch=21, alpha = 0.4, size =2)+
  theme_bw()+facet_wrap(vars(AREA))+
  ylim(0,25)+xlab("Date")+
  ylab("Number of Recruits")
p_recruits_2plus_date

#boxplots of number of recruits per census
p_recruits_boxplot <- recruits|>
  filter(CENSUS_NUM>1)|>
  ggplot(aes(factor(CENSUS_NUM), NumRecruits,fill = AREA ))+
  geom_boxplot()+
  theme_bw()+facet_wrap(vars(AREA))
p_recruits_boxplot

#Relationship b/w number of recruits and cones left after fire
recruits|>
  filter(CENSUS_NUM ==1)|>
  ggplot(aes(CONES1, NumRecruits))+
  geom_point()+theme_bw()+
  xlab("Number of Alive Cones after Fire")+
  ylab("Number of Recruits in Census 1")+
  ggtitle("Census 1")#+xlim(0,1)+ylim(0,100)

recruits|>
  filter(CENSUS_NUM==1)%>%
  glm(NumRecruits~0+CONES1, family="poisson", data = .)|>
  summary()


#Relationship b/w number of recruits and top alive trees after fire
recruits|>
  filter(CENSUS_NUM ==1)|>
  ggplot(aes( TOPALIV1, NumRecruits))+
  geom_point()+theme_bw()+
  xlab("Number of Top Alive Trees Found after Fire")+
  ylab("Number of Recruits in Census 1")+
  ggtitle("Census 1")

recruits|>
  filter(CENSUS_NUM==1)%>%
  glm(NumRecruits~TOPALIV1, family="poisson", data = .)|>
  summary()


#Vital rates trial: survival----------------------------------------------
#Plotting height of surviving vs dead trees
pines_long|>
  #filter(COHORT==1)|>
  group_by(IND_ID)|>
  mutate(STAT_next = lead(STAT))|>
  filter(STAT ==1 & !is.na(STAT_next))|>
  ggplot(aes(factor(STAT_next), (HT), group=STAT_next))+
  geom_boxplot()+
  theme_bw()+facet_wrap(vars(AREA))

#Plotting height of surviving vs dead trees
pines_long|>
  filter(COHORT==1)|>
  group_by(IND_ID)|>
  mutate(STAT_next = lead(STAT))|>
  filter(STAT == 1 & !is.na(STAT_next))|>
  ggplot(aes( factor(STAT_next), log(DAH), group=STAT_next))+
  geom_boxplot()+
  theme_bw()+facet_wrap(vars(AREA))

#Creating dataframe for survival vs height plotting
fire_date <- as.Date("06/01/95",format= "%m/%d/%y")
survival <- pines_long|>
  group_by(IND_ID)|>
  #pairing survival in census t with survival in census t+1
  mutate(STAT_next = lead(STAT),
         CENSUS_next = lead(CENSUS_NUM),
         Date_next = lead(START_DATE))|>
  #only keeping trees that were alivecensus t
  #and had height measurements in census t
  #and were had their survival recorded in census t+1
  filter(STAT==1 & !is.na(STAT_next) & !is.na(HT))|>
  #calculating time interval b/w censuses, in years
  mutate(CensusInterval = 
           time_length(difftime(Date_next, START_DATE),"years"))|>
  mutate(CensusInterval = round(CensusInterval, digits = 1))|>
  #calculating time since fire, in years
  mutate(TimeSinceFire = 
           time_length(difftime(START_DATE, fire_date),"years"))|>
  mutate(TimeSinceFire = round(TimeSinceFire, digits = 1))

ggplot(survival, aes(HT,STAT_next))+
  geom_point()+theme_bw()+
  facet_wrap(~CensusInterval+AREA, scales = "free_x")
ggplot(survival, aes(HT,STAT_next))+
  geom_point()+theme_bw()+
  facet_wrap(~CensusInterval+TimeSinceFire, scales = "free_x")


#Vital rates trial: growth----------------------------------------------
#Calculating growth increment
growth = pines_long|>
  #filter(COHORT==1)|>
  group_by(IND_ID)|>
  #keeping only those censuses where height was measured
  filter(CENSUS_NUM %in% c(3, 5, 7, 9, 11:14))|>
  #pairing each census with values in the next census
  mutate(STAT_next = lead(STAT),
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
  mutate(HT_next_annual = HT+deltaHT_annual)

#Regression of height in year t vs height in year t+1
m_ht_annual = lm(HT_next_annual~0+HT, data = growth)
#Plotting height in year t vs height in year t+1
plot(growth$HT, growth$HT_next_annual, 
     xlab = "Height in year t (cm)",
     ylab ="Height in year t+1 (cm)", pch = 21, 
     bg =alpha("grey", 0.2),
     col = alpha("black", 0.5))
lines(growth$HT, predict(m_ht_annual), col ="red", lwd = 1)

#Regression of height in year t vs annual growth increment
m_incr_annual = lm(deltaHT_annual~0+HT, data = growth)
#Plotting height in year t vs annual growth increment
plot(growth$HT, growth$deltaHT_annual, 
     xlab = "Height in year t (cm)",
     ylab ="Annual height increment (cm)", pch = 21, 
     bg =alpha("grey", 0.2),
     col = alpha("black", 0.5))
lines(growth$HT, predict(m_incr_annual), col ="blue", lwd = 1)


#Regression of log height in year t vs log ht in year t+1
m_ht_log = lm(log(HT_next_annual)~log(HT), data = growth)
#Plotting height in year t vs annual growth increment
plot(log(growth$HT), log(growth$HT_next_annual), 
     xlab = "Log height in year t (cm)",
     ylab ="Log height in year t+1 (cm)", pch = 21, 
     bg =alpha("grey", 0.2),
     col = alpha("black", 0.5))
lines(log(growth$HT), 
    predict(m_ht_log), col ="blue", lwd = 1)

#Plots separated by area (population type)
growth|>
  ggplot(aes(HT, HT_next_annual, fill = AREA))+
  geom_point(pch = 21, alpha = 0.5)+
  facet_wrap(vars(AREA))+
  theme_bw()+xlab("Height in year t (cm)")+
  ylab("Height in year t+1 (cm)")


growth|>
  ggplot(aes(HT, deltaHT_annual/HT, fill = AREA))+
  geom_point(pch = 21, alpha = 0.2)+
  facet_wrap(vars(AREA))+
  theme_bw()+ylim(0,5)+xlab("Height (cm)")+
  ylab("Annual relative growth rate")

growth|>
  ggplot(aes(HT, deltaHT_annual, fill = AREA))+
  geom_point(pch = 21, alpha = 0.2)+
  facet_wrap(vars(AREA), scales = "free_x")+
  theme_bw()+xlab("Height (cm)")+
  ylim(-50, 100)+
  ylab("Annual growth increment")

