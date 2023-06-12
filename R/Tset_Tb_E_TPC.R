# Tb estimates from lab and GAMM Analysis
# Tset, Tb, Tpredict, TPC, 
rm(list=ls())
pacman::p_load("dplyr", "tidyverse", "lubridate", "ggplot2", "plotly", "pbapply","car", "readr", "lmerTest", "emmeans", "AICcmodavg", "nls.multstart", "broom", "purrr", "plotrix", "AICcmodavg", "MuMIn", "stargazer", "jtools", "huxtable", "officer", "grid", "png", "ggimage", "rTPC", "nls.multstart", "ggrepel", "mgcv", "gratia", "gammit", "lmerTest","performance")

######################################################################################  
############################  (1) TSET by sex - LAB ################################### 
######################################################################################  
# import data
Tsel <- read.csv(file = "R/Raw_data/Final_TselTb_file.csv" ) %>% 
  mutate(temp = mean_temp) %>% 
  filter(Sex != "ZZf")

# differences in IQR by sex? - no differences for upper or lower 
IQR.Summary <- Tsel %>%
  group_by(id, Sex) %>%
  summarise(IQR = IQR(temp), 
            Mean.Tb = mean(temp)) %>% 
  mutate(range = IQR/2,
         upper.IQR = Mean.Tb + range,
         lower.IQR = Mean.Tb - range) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()
# save this file 
saveRDS(IQR.Summary, file = "R/Final.Analysis/Final.Figure.data/IQR.Final.RDS")

# Mean differences? - yes, between females and males. Females select warmer temps
Mean.Tel <- lm(Mean.Tb ~Sex, data = IQR.Summary)
summary(Mean.Tel)
anova(Mean.Tel)
plot(residuals(Mean.Tel))
emmeans(Mean.Tel, pairwise ~ Sex)

# differences in upper? - no differences
IQR.upper <- lm(upper.IQR ~Sex,  data = IQR.Summary)
anova(IQR.upper)
emmeans(IQR.upper, pairwise ~ Sex)
# differences in lower? - yes there are differences between females and males
IQR.lower <- lm(lower.IQR ~Sex,  data = IQR.Summary)
anova(IQR.lower)
emmeans(IQR.lower, pairwise ~ Sex,)



######################################################################################  
############################  (2) Tb Predict - LAB ################################### 
######################################################################################  
# estimating Tb correction
lab.tb <- read.csv("R/Raw_data/TselvsTsurf.csv", header = T) %>% 
  drop_na() %>% 
  select(c(Lizard, Date, Time, Date.Time, Tsurf.Temp, Tsel.temp, Time.of.Day)) 
lab.tb$Date.Time <- dmy_hm(lab.tb$Date.Time)
lab.tb$hr <- hour(lab.tb$Date.Time)

# plot of relationship from lab
ggplot(lab.tb, aes(x=Tsurf.Temp, y=Tsel.temp))+
  geom_point()+
  xlab(expression(Skin~Temp~T[surf])) +
  ylab(expression(Body~Temp~T[b]))+
  stat_smooth(method = "lm" , color="red") +
  theme_bw()

# paired t test for Tb and tsurf
t.test <- t.test(lab.tb$Tsurf.Temp, lab.tb$Tsel.temp, paired = TRUE, alternative = "two.sided")
# model for predicting Tb = a + Tsurf*b
# where a is the intercept and b is the slope
mod.time <- lm(Tsurf.Temp ~Tsel.temp, data = lab.tb)
summary(mod.time)
r.squaredLR(mod.time)
confint(mod.time)
#making final table anova table for regression 
table <- summ(mod.time, digits = 4)
table
anova(mod.time)

######## 
# Tb Tset and Tsurf Tset before correction
######## 
# Tb - Tset
range(lab.tb$Tsel.temp)
Tset.tb.sum <- lab.tb %>% 
  summarise(meantb = mean(Tsel.temp),
            se = std.error(Tsel.temp))
Tset.Tb <- quantile(lab.tb$Tsel.temp)

# Tsurf - Tset
mean(lab.tb$Tsurf.Temp)
Tset.Tsurf.sum <- lab.tb %>% 
  summarise(meantb = mean(Tsurf.Temp),
            se = std.error(Tsurf.Temp))
Tset.Tsurf <-  quantile(lab.tb$Tsurf.Temp)

######## 
### apply regression Tsurf to Tb correction for tsel estimation from tsurf
######## 
Tsel.final.dat <- lab.tb %>% 
  mutate(Tb.lab.predict = -1.770131 + (Tsurf.Temp*1.057887))

# Tsel final with correction that will be used in the field
# Tset Tsurf
Tset.Tsurf.sum <- Tsel.final.dat %>% 
  summarise(meantb = mean(Tb.lab.predict),
            se = std.error(Tb.lab.predict))
Tset.Tsurf.Final <-  quantile(lab.tb$Tsurf.Temp)



######################################################################################  
############################  (3) Field accelerometer data  ##########################
############################  temperatures Spring 2018 & Winter 2019 #################
############################  prepping data for activity, TB, and GAMM's #############
###################################################################################### 
###########
# 1) arranging data - check for date and tracking days
other.accel.data <- readRDS(file = "R/Raw_data/Tsurf_field.RDS") 
str(other.accel.data)
### DATES
max(other.accel.data$Date)
# fitler out dates past 1 Spetember 2019 and times between 5am - 2100
accel.all <- other.accel.data %>% 
  filter(Date <= "2019-09-01") %>% 
  filter(HR <= 21 & HR >= 5) 

# TRACKING DAYS
# filter out animals that had less than 7 days of data
accel.all.filtered <- accel.all %>% 
  group_by(POVI, Sex, Date) %>% 
  summarise(Tsurf = mean(Tsurf))
# now get the freq of individuals 
freq <- count(accel.all.filtered, 'POVI') %>% 
  arrange(n)
# drop: POVI_23, POVI_25, POVI_33, POVI_2, POVI_40, POVI_6 
accel.all.filtered.drop <- accel.all %>% 
  filter(POVI != "POVI_23" & POVI != "POVI_25" & 
           POVI != "POVI_33" & POVI != "POVI_2" &
           POVI != "POVI_40" & POVI != "POVI_6")


########
# 2) apply Tb predicition from lab regression: Tb = a + Tsurf*b
# where a = -1.770131; b = 1.057887
# also calculating vector ms2 by taking sqrt of (xms2 + Yms2)
# note this is for every 2 mins!
# NEEDS TO BE SUMMARISED BY HOUR FOR TB data and this will be used for general activity plots!
########
# Tb predicition and activity
pref.other <- accel.all.filtered.drop %>% 
  mutate(Tb = -1.770131 + (Tsurf*1.057887),
         vect.ms2 = sqrt(X.min.ms2+Y.min.ms2))
# save for Activity analysis by min & we will use this df for summarising Tb's
saveRDS(pref.other, "R/Final.Analysis/Final.Figure.data/Final.Activity.mins.data.RDS")


######## 
# 3) DATA USED FOR GAMM ANALYSIS : TSET AND Performance in the field - use later
######## 
# extracting 95% of vect.ms2 by hr of each day for each individual - GAMM
pref.other.95 <- pref.other %>% 
  group_by(POVI, HR, Date, Month, Sex, Season) %>% 
  summarise(percent95.ms2 = quantile(vect.ms2, probs = .95),
            Tb = mean(Tb))
pref.other.95$temp <- as.factor(round(pref.other.95$Tb, digits = 0))
# summarise by ID, temp, season, Sex - GAMM
pref.other.95.final <- pref.other.95 %>% 
  group_by(POVI, temp, Season, Sex) %>% 
  summarise(percent95.ms2 = quantile(percent95.ms2, probs = .95),
            Tb = mean(Tb))
# filter out activity values that have less than 5 Tb measurements for TPC curves which are values between -4 and 45
pref.other.95.final <- pref.other.95.final %>% filter(Tb >= -4) %>% filter(Tb<=49)
# save for GAM(M) data
saveRDS(pref.other.95.final, "R/Final.Analysis/Final.Figure.data/GAMM.data.RDS")



######################################################################################  
############################  (4) Tsel Corrected   ################################### 
############################ individual db, E - FIELD  ############################### 
######################################################################################  
Tb.predict.dat <- readRDS(file = "R/Final.Analysis/Final.Figure.data/Final.Activity.mins.data.RDS") %>% filter(Sex != "ZZf")
Tb.predict.dat$Hour <- as.factor(Tb.predict.dat$HR)

###############
# 1) summarise by hour to meet data
###############
Tb.predict.sum <- Tb.predict.dat %>% 
  group_by(POVI, Date) %>% 
  summarise(Tb = mean(Tb),
            HR = unique(HR),
            Season = unique(Season),
            Sex = unique(Sex))
  
################
# 2) Ex function: using tset for db values for individuals by sex and season
################
# set up Ex function
Ex <-function(x)
{
  for (i in 1:length(x)) {
    if(x[i] >= 0.5) 
    {
      x [i]<- 0
    }
    else
    {
      x[i] <- 1
    }
  }
  return(x)
}
# Use Tset to  to field data because there are differences in sex
# Sex      IQR  Mean.Tb    range upper.IQR lower.IQR
# ZWf 6.833333 30.38976 3.416667  33.80642  26.97309
# ZZm 6.700000 28.68080 3.350000  32.03080  25.33080

###############
# 3) Tset to  to field data for males and females
###############
# females
db.ZWf<- Tb.predict.sum %>%
  filter(Sex == "ZWf") %>% 
  dplyr::mutate(resid = ifelse(between(Tb, 26.97, 33.80), 0, NA),
                resid = ifelse(Tb < 26.97, 26.97 - Tb, resid),
                resid = ifelse(Tb > 33.80, Tb - 33.80, resid))
# indicate individual EX
db.ZWf$Ex_db <- Ex(db.ZWf$resid)

# males
db.ZZm<- Tb.predict.sum %>%
  filter(Sex == "ZZm") %>% 
  dplyr::mutate(resid = ifelse(between(Tb, 25.33, 32.03), 0, NA),
                resid = ifelse(Tb < 25.33, 26.97309 - Tb, resid),
                resid = ifelse(Tb > 32.03, Tb - 32.03, resid))
# indicate individual EX
db.ZZm$Ex_db <- Ex(db.ZZm$resid)


###############
# 4) combing df and adding month col
###############
db.sex.all <- rbind(db.ZWf, db.ZZm) %>% 
  rename(individual_db = resid)
# final analysis where we get the average individual db by date and hr to meet with Te data
db.sex.all.final <- db.sex.all %>% 
  group_by(POVI, Date, HR) %>% 
  summarise(individual_db = mean(individual_db),
            Tb = mean(Tb), 
            Season = unique(Season),
            Sex = unique(Sex))
# add month for later
db.sex.all.final$Month <- format(db.sex.all.final$Date,"%B")
# now save the final data file 
saveRDS(db.sex.all.final, "R/Final.Analysis/Final.Figure.data/db.final.RDS")



######################################################################################  
####################### (5)Te and de values by season and sex   ######################
######################################################################################  
###############
# 1) Set up data for analysis, set up dates, set up times to match accelerometers
###############
Te <-read.csv(file = "R/Raw_data/Te.All.Final.csv") %>% 
  rename(Te=Value)
Te$Season <- plyr::revalue(Te$Month, c("June" = "Winter",
                                       "July" = "Winter",
                                       "August" = "Winter",
                                       "September" = "Spring", 
                                       "October" = "Spring", 
                                       "November"= "Spring",
                                       "December" = "Summer", 
                                       "January" =  "Summer",
                                       "February" = "Summer",
                                       "March" = "Autumn", 
                                       "April" = "Autumn", 
                                       "May" = "Autumn"))
#format date and filter to fit accelerometer date "2019-09-01" and times between 0500 - 2100
Te$Date <- dmy(Te$Date)
min(Te$Date)
max(Te$Date)
Te <- Te %>% 
  filter(Date < "2019/09/01") %>% 
  filter(Time >= 5) %>% 
  filter(Time <= 21) %>% 
  mutate(HR = Time)
min(Te$HR)
max(Te$HR)
  
# pre calibration fig
Fig <- Te %>%
  group_by(Date, Season)%>%
  summarise(Max = max(Te),
            Min = min(Te),
            Mean = mean(Te)) 
Fig.tb <- db.sex.all.final %>% 
  group_by(Date, Season) %>% 
  summarise(Max = max(Tb),
            Min = min(Tb),
            Mean = mean(Te))

###############
# Te - Mean, min, max fig
# Te ribbon and d
ggplot(Fig, aes( x = Date, y = Mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.2)  +
  geom_line(data = Fig.tb, aes(x = Date, y = Mean), linetype = "longdash", color = "green")+
  geom_line(data = Fig.tb, aes(x = Date, y = Min), linetype = "dotted", color = "blue")+
  geom_line(data = Fig.tb, aes(x = Date, y = Max), linetype = "dotted", color = "red")+
  theme(aspect.ratio = 1, legend.position = c(0.8, 0.2))+
  theme_bw()

###############
# 2) TE calibration (Tb vs Tb dead)
###############
Te.Calibration <- read.csv(file = "R/Final.Analysis/Final.Figure.data/Model.vs.TB.Clean.2min.csv")
Te.Calibration.Results <- lm(Tb_dead ~ Model_temp, Te.Calibration)
Te.Calibration.Results.anova <- anova(Te.Calibration.Results)
# TE model adjustment Te = a + TbCarcass*b
# where a is the intercept= 8.159055 and b is the slope 0.711953
summary(Te.Calibration.Results)
r.squaredLR(Te.Calibration.Results)
#making final table 
Te.cal.table <- summ(Te.Calibration.Results, digits = 4)
Te.cal.table
# plot of relationship between carcass and Te
ggplot(Te.Calibration, aes(x=Tb_dead, y=Model_temp))+
  geom_point()+
  xlab(expression(Carcass~Temp~T[surf])) +
  ylab(expression(Model~Temp~T[b]))+
  stat_smooth(method = "lm" , color="red") +
  theme_bw()

###############
# 3) TE calibration from carcass 
# model adjustment Te = a + TbCarcass*b
# where a is the intercept= 8.159055 and b is the slope 0.711953
###############
Te <- Te %>% 
  mutate(Te.corrected = round(8.159055 + (Te*0.711953), digits =2))
Te$Te = NULL # drop old name
# rename Te corrected to Te to fit analysis below
Te$Te = Te$Te.corrected
Te$hr = as.character(Te$Time)
# calibration figure
Fig <- Te %>% 
  group_by(Date, Season)%>%
  summarise(Max = max(Te),
            Min = min(Te),
            Mean = mean(Te)) 
# Te figure after calibrations figure
ggplot(Fig, aes( x = Date, y = Mean)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.2) +
  theme(aspect.ratio = 1)

#############
# 4) figure out de/Ex and summarize after <- VERY IMPORTANT-
# de is the  Te outside window for females and males separately 
# because they have different Tsel
#############
# de - females
de.female <- Te %>% 
  dplyr::mutate(resid = ifelse(between(Te, 26.97, 33.80), 0, NA),
                resid = ifelse(Te < 26.97, 26.97 - Te, resid),
                resid = ifelse(Te > 33.80, Te - 33.80, resid)) %>% 
  mutate(Te.Analysis = "de_ZWf")
# add Ex before summary
de.female$Ex = Ex(de.female$resid)
# now summarise de by Model date hr and habitat
de.female.final <- de.female %>% 
  group_by(Date, Type, Model_., hr) %>% 
  summarise(Te = mean(Te),
            de = mean(resid),
            Ex = mean(Ex),
            Season = unique(Season),
            Te.Analysis = unique(Te.Analysis),
            HR = unique(Time))

# de- males
de.male <- Te %>% 
  dplyr::mutate(resid = ifelse(between(Te, 25.33, 32.03), 0, NA),
                resid = ifelse(Te < 25.33, 26.97309 - Te, resid),
                resid = ifelse(Te > 32.03, Te - 32.03, resid)) %>% 
  mutate(Te.Analysis = "de_ZZm")
# add Ex before summary
de.male$Ex = Ex(de.male$resid)
# now summarise de by Model date hr and habitat
de.male.final <- de.male %>% 
  group_by(Date, Type, Model_., hr) %>% 
  summarise(Te = mean(Te),
            de = mean(resid),
            Ex = mean(Ex),
            Season = unique(Season),
            Te.Analysis = unique(Te.Analysis),
            HR = unique(Time))
# rbind male and add month for later
de.sex.final <- rbind(de.female.final, de.male.final)
de.sex.final$Month =  format(de.sex.final$Date,"%B")
# save file
saveRDS(de.sex.final, file = "R/Final.Analysis/Final.Figure.data/Te.sex.final.RDS")



######################################################################################  
####################### (6) E and Eindex E = 1-(db/de) ##############################
######################################################################################  
##########
# E and Eindex
# E = 1-(db/de)
##########
# this db to be calculated by individual lizard and season
# then apply seasonal de average to the df- this will give you se errors
# IMPORTANT THIS WILL BE DIFFERENT BY SEX BECAUSE OF TSEL
# 1) calculate Te summary across seasons: - data is already summarised by Model date hr and habitat
Te.summary <- de.sex.final %>% 
  group_by(Te.Analysis, Season) %>% 
  summarise(grandmean_de_season=mean(de))
##########
# Te.Analysis Season grandmean_de_season - used fo E
# de_ZWf      Autumn                3.78
# de_ZWf      Spring                2.33
# de_ZWf      Summer                2.99
# de_ZWf      Winter                8.48
# de_ZZm      Autumn                3.99
# de_ZZm      Spring                2.81
# de_ZZm      Summer                3.97
# de_ZZm      Winter                8.47

##########
#2) calculate seasonal mean db for lizards db is already arranged by POVI, Date, HR
individual_season <- db.sex.all %>% 
  group_by(POVI, Season, Sex) %>% 
  summarise(mean_individual_mean_db = mean(individual_db))
##########

##########
# 3) apply the de to each season for each sex
# FEMALES
ZWf_individual_season <- individual_season %>%
  filter(Sex=="ZWf")
# females for each season 
# Autumn                3.78
ZWf_individual_Autumn<- ZWf_individual_season %>%
  filter(Season=="Autumn") %>% 
  mutate(grandmean_de_season = 3.78,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# Spring                2.33
ZWf_individual_Spring<- ZWf_individual_season %>%
  filter(Season=="Spring") %>% 
  mutate(grandmean_de_season = 2.33,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# Summer                2.99
ZWf_individual_Summer <- ZWf_individual_season %>%
  filter(Season=="Summer") %>%
  mutate(grandmean_de_season = 2.99,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# de_ZWf      Winter                8.48
ZWf_individual_Winter <- ZWf_individual_season %>%
  filter(Season=="Winter") %>% 
  mutate(grandmean_de_season = 8.48,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))

# MALES
ZZm_individual_season <- individual_season %>%filter(Sex=="ZZm") 
# de_ZZm      Autumn                3.99
ZZm_individual_Autumn<- ZZm_individual_season %>%
  filter(Season=="Autumn") %>% 
  mutate(grandmean_de_season = 3.99,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# de_ZZm      Spring                2.81
ZZm_individual_Spring<- ZZm_individual_season %>%
  filter(Season=="Spring") %>% 
  mutate(grandmean_de_season = 2.81,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# de_ZZm      Summer                3.97
ZZm_individual_Summer <- ZZm_individual_season %>%
  filter(Season=="Summer") %>% 
  mutate(grandmean_de_season = 3.97,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
# de_ZZm      Winter                8.47
ZZm_individual_Winter <- ZZm_individual_season %>%
  filter(Season=="Winter") %>% 
  mutate(grandmean_de_season = 8.47,
         indvidual_season_E = 1-(mean_individual_mean_db/grandmean_de_season))
##########

##########
# 3) individual E combind and save for final analysis
E_final_season = rbind(ZWf_individual_Autumn, ZWf_individual_Spring, ZWf_individual_Summer, ZWf_individual_Winter, 
                       ZZm_individual_Autumn, ZZm_individual_Spring, ZZm_individual_Summer, ZZm_individual_Winter)
##########

##########
# 4) Check how they are different across season and sex accounting for individuals
E.mod <- lmer(indvidual_season_E  ~ Sex + Season + Season*Sex + (1|POVI), E_final_season)
E.mod.Sum <- anova(E.mod) # significant all around
E.anova.pairwise <- plot(emmeans(E.mod, pairwise ~ Sex |Season))
##########

##########
# 5) Save data with individuals for MARK and for rmd file analysis 
saveRDS(E_final_season, file = "R/Final.Analysis/Final.Figure.data/E_season_sex_individual.RDS")
##########



############################################################################### 
############## (7) Field accelerometer data for avg mins moved ################
############################################################################### 
mydata <- readRDS(file = "R/Final.Analysis/Final.Figure.data/avg.min.final.clean.RDS")
mydata$Hour <- as.factor(mydata$HR)
mydata <-mutate(mydata, Season = factor(Season, levels = c("Spring","Summer", "Autumn", "Winter"))) %>%
  group_by(Season)

#####################
## Activity function - mins moved
#####################
# activity (mins) function
act <-function(x)
{for (i in 1:length(x)) {
  if(x[i] <= 0) 
  { x [i]<- 0
  }
  else
  {x[i] <- 1}}
  return(x)
}
mydata$active.mins <- act(mydata$acc.min.sum)
mydata<-mydata %>% filter(Sex != "ZZf")

###############
# 1) Activity data for and plot analysis
# min moved - check data
max(mydata$active.mins)
mydata.final <- mydata %>%  
  group_by(Date, HR, POVI) %>% 
  summarise(freq = sum(active.mins),
            temp = mean(Tsurf),
            Sex = unique(Sex),
            Season = unique(Season))
mean(mydata.final$freq)
# by season, sex hr
mydata.final.mm <- mydata.final %>% 
  group_by(Season, Sex, HR) %>% 
  summarise(Tsurf = mean(temp),
            Movement.HR = mean(freq),
            SE = std.error(freq))
# save data for figure
saveRDS(mydata.final.mm, file = "R/Final.Analysis/Final.Figure.data/activity.fig.data.rds")

# movement plots by season
ggplot(mydata.final.mm, aes(x=HR, y=Movement.HR, group=Sex, color=Sex)) + 
  geom_pointrange(aes(ymin=Movement.HR-SE, ymax=Movement.HR+SE))+
  facet_grid(~Season) +
  ylab("Avg number of active mins") +
  xlab("Hour of day") 

# summarise activity by POVI
final.POVI.Activity <- mydata.final %>% 
  group_by(POVI, Sex) %>% 
  summarise(Avg.mins.moved = mean(freq))
final.POVI.Activity

##### Activity analysis summary/final data 
mydata.final <- mydata %>%  
group_by(Date, HR, POVI) %>% 
  summarise(freq = sum(active.mins),
            temp = mean(Tsurf),
            Sex = unique(Sex),
            Season = unique(Season))
mean(mydata.final$freq)
# by season, sex hr
act.analysis.dat <- mydata.final %>% 
  group_by(Season, Sex, POVI) %>% 
  summarise(Tsurf = mean(temp),
            Movement.HR = mean(freq),
            SE = std.error(freq))
saveRDS(act.analysis.dat, "R/Final.Analysis/Final.Figure.data/activity.season.sex.analysis.data.RDS")
###############


###############
# 2) Activity (mins moved) analysis
act.analysis.dat <- readRDS(file = "R/Final.Analysis/Final.Figure.data/activity.season.sex.analysis.data.RDS")
# Analysis - season effects once interaction is removed
activity.analysis.mod <- lmer(Movement.HR ~ Sex + Season + Season*Sex + (1|POVI), act.analysis.dat)
anova(activity.analysis.mod) # no interaction effect
activity.analysis.mod.final <- lmer(Movement.HR ~ Sex + Season + (1|POVI), act.analysis.dat)
anova(activity.analysis.mod.final) 
###############



############################################################################ 
########################## (8) GAM(M) ANLAYSIS #############################
############################################################################
# arranging datra for model
pref.95 <- readRDS(file = "R/Final.Analysis/Final.Figure.data/GAMM.data.RDS")
pref.95$temp <- as.numeric(as.character(pref.95$temp))
pref.95$percent95.ms2 <- as.numeric(as.character(pref.95$percent95.ms2))
pref.95$POVI <- as.factor(pref.95$POVI)
pref.95$Sex <- as.factor(pref.95$Sex)
pref.95$Season <- as.factor(pref.95$Season)
data.gam <- pref.95 %>%  filter(Sex != "ZZf")
write_rds(data.gam, file = "R/Final.Analysis/Final.Figure.data/data.temp.gam.analysis.rds")
data.gam <- readRDS(file = "R/Final.Analysis/Final.Figure.data/data.temp.gam.analysis.rds")

############################
## 1)  Building GAMM Models
# smooth on temp
# Temp <- gam(percent95.ms2 ~ s(temp), method = 'REML', data = data.gam)

# smooth for individual and their temp
# Temp.Povi <- gam(percent95.ms2 ~ s(temp, by = POVI) + s(POVI, bs = 're'),  method = 'REML', data = data.gam)

# m2 smooth by temp for sex - ID as random
# Sex <- gam(percent95.ms2 ~ Sex + s(temp, by = Sex) + s(POVI, bs = 're'), method = 'REML', data = data.gam)

# smooth for individual and their temp - ID as random
# Sex.Povi <-gam(percent95.ms2 ~ Sex + s(temp, by = POVI) + s(POVI, bs = 're'),  method = 'REML', data = data.gam)

# smooth for temp by season - ID as random
# Season <- gam(percent95.ms2 ~ Season + s(temp, by = Season) + s(POVI, bs = 're'), method = 'REML', data = data.gam)

# smooth for individual by season - ID as random
# Season.Povi <-gam(percent95.ms2 ~ Season + s(temp, by = POVI) + s(POVI, bs = 're'),  method = 'REML', data = data.gam)

# smooth by temp for season and sex - ID as random
# Season.Sex <- gam(percent95.ms2 ~ Season + Sex + s(temp) + s(POVI, bs = 're'), method = 'REML', data = data.gam)

# smooth for individual by season and sex - ID as random
# Season.Sex.Povi <-gam(percent95.ms2 ~ Season + Sex + s(temp, by = POVI) + s(POVI, bs = 're'),  method = 'REML', data = data.gam)

# smooth  by temp for season, sex, and their interaction - ID as random
# Season.Sex.Interaction <- gam(percent95.ms2 ~ Season + Sex + Season*Sex + s(temp) + s(POVI, bs = 're'), method = 'REML', data = data.gam)

# smooth for individual by season, sex, and their interaction - ID as random
# Season.Sex.Interaction.Povi <- gam(percent95.ms2 ~ Season + Sex + Season*Sex + s(temp, by = POVI) + s(POVI, bs = 're'), method = 'REML', data = data.gam)
############################


############################
# 2)  Saving & Read in GAMM Models
# save models
saveRDS(Temp, file = "R/Mod/Temp.rds")
saveRDS(TR/emp.Povi, file = "Models/Temp.Povi.rds")
saveRDS(Sex, file = "R/Models/Sex.rds")
saveRDS(Sex.Povi, file = "R/Models/Sex.Povi.rds")
saveRDS(Season, file = "R/Models/Season.rds")
saveRDS(Season.Povi, file = "R/Models/Season.Povi.rds")
saveRDS(Season.Sex, file = "R/Models/Season.Sex.rds")
saveRDS(Season.Sex.Povi, file = "R/Models/Season.Sex.Povi.rds")
saveRDS(Season.Sex.Interaction, file = "R/Models/Season.Sex.Interaction.rds")
saveRDS(Season.Sex.Interaction.Povi, file = "R/Models/Season.Sex.Interaction.Povi.rds")
# read models
Temp <- readRDS(file = "R/Models/Temp.rds")
Temp.Povi <- readRDS(file = "R/Models/Temp.Povi.rds")
Sex <- readRDS(file = "R/Models/Sex.rds")
Sex.Povi <- readRDS(file = "R/Models/Sex.Povi.rds")
Season <- readRDS(file = "R/Models/Season.rds")
Season.Povi <- readRDS(file = "R/Models/Season.Povi.rds")
Season.Sex <- readRDS(file = "R/Models/Season.Sex.rds")
Season.Sex.Povi <- readRDS(file = "R/Models/Season.Sex.Povi.rds")
Season.Sex.Interaction <- readRDS(file = "R/Models/Season.Sex.Interaction.rds")
Season.Sex.Interaction.Povi <- readRDS(file = "R/Models/Season.Sex.Interaction.Povi.rds")
############################


############################
## 3)  GAMM Model - CHECKS
summary(Temp)  #  % of deviance -- POVI as random effect
summary(Temp.Povi) # % of deviance -- with POVI as random, and fitted curve by POVI, no fixed effects
summary(Sex) # % of deviance -- with POVI as random, sex as fixed factors and Sex fitted
summary(Sex.Povi) # % of deviance -- with POVI as random, sex as fixed factors and fitted curve by POVI
summary(Season) # % of deviance -- with POVI as random, season as fixed factors and season fitted
summary(Season.Povi) #  % with POVI as random, season as fixed factors and fitted curve by POVI
summary(Season.Sex)  #  % of deviance --   with POVI as random, Sex and fitted curve by temp
summary(Season.Sex.Povi) #  % of deviance --  with POVI as random, Sexand fitted curve by POVI
summary(Season.Sex.Interaction) #  % of deviance --  with POVI as random, Sex, season, and interaction and fitted curve by temp
summary(Season.Sex.Interaction.Povi) #  % of deviance --  with POVI as random, Sex, season, and interaction and fitted curve by ID 
############################


############################
# 4)  ANOVA TABLE GAMM Models
# # make a table of comparing all models to model 1. 
# i.e. does adding terms significantly change our understanding of the data set?
# ANOVA Table
anova.table <- rbind(anova(Temp,Temp.Povi, test="Chisq"),
                     anova(Temp,Sex, test="Chisq")[2,],
                     anova(Temp,Sex.Povi, test="Chisq")[2,],
                     anova(Temp,Season, test="Chisq")[2,],
                     anova(Temp,Season.Povi, test="Chisq")[2,],
                     anova(Temp,Season.Sex, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Povi, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Interaction, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Interaction.Povi, test="Chisq")[2,])

anova.table <- cbind(anova.table, AIC(Temp, Temp.Povi, Sex, Sex.Povi, Season, Season.Povi, Season.Sex, Season.Sex.Povi, Season.Sex.Interaction, Season.Sex.Interaction.Povi)[,2])
anova.table <- cbind(c("Temp","Temp.Povi","Sex","Sex.Povi","Season","Season.Povi","Season.Sex","Season.Sex.Povi","Season.Sex.Interaction","Season.Sex.Interaction.Povi"), anova.table)
anova.table[,c(2:5,7)] <- round(anova.table[,c(2:5,7)],2)
anova.table[,6] <- signif(anova.table[,6],2)
rownames(anova.table) <- NULL
anova.table$DevEx <- c(summary(Temp)$dev.expl, summary(Temp.Povi)$dev.expl,
                       summary(Sex)$dev.expl, summary(Sex.Povi)$dev.expl,
                       summary(Season)$dev.expl, summary(Season.Povi)$dev.expl,
                       summary(Season.Sex)$dev.expl, summary(Season.Sex.Povi)$dev.expl,
                       summary(Season.Sex.Interaction)$dev.expl, summary(Season.Sex.Interaction.Povi)$dev.expl)
anova.table$"Dev.Expl(%)" <-  round(anova.table$DevEx*100, 2)
names(anova.table)[c(1:4,7, 8)] <- c("Model", "Residual DF", "Residual Deviance", "DF", "AICc", "Deviance Explained" )
anova.table$`Deviance Explained` <- round(anova.table$`Deviance Explained`, 2)
############################


############################
# 5)  ID TOP MOD AND SAVE
anova.table %>% arrange(AICc)
view(anova.table) # looks like the m7 that accounts for  season and sex and the interaction, with random intercept per individual, and smooth per individual
# save table -  Go in and change names of cols for manuscript
write.csv(anova.table, "R/Final.Analysis/Final.Figure.data/GAMM.anova.table.csv") 
############################



############################################################################ 
########################## (9) TPC From TOP Mod #############################
############################################################################
#########
# 1) TPC parameters; To, Tb from top model
Season.Sex.Interaction.Povi <- readRDS(file = "R/Models/Season.Sex.Interaction.Povi.rds")
anova(Season.Sex.Interaction.Povi)
summary.gam(Season.Sex.Interaction.Povi)
gam.check(Season.Sex.Interaction.Povi)
data.gam$fit <- predict(Season.Sex.Interaction.Povi, type = "response")
data.gam$se <- predict(Season.Sex.Interaction.Povi, type = "se")
#########

#########
# 2) summarizing To by id, Sex, Season
Topt <- data.gam %>% 
  group_by(POVI, Sex, Season) %>% 
  slice(Topt = which.max(fit)) 

Topt.final <- Topt %>% 
  group_by(POVI, Sex, Season) %>% 
  summarise(Topt.temp = mean(Tb),
            Topt = mean(fit))
#########


#########
# 3) Topt temp analysis
To.mod <- lmer(Topt.temp  ~ Season + Sex + Season*Sex  + (1|POVI), Topt.final)
To.mod.sum <- anova(To.mod) # no differences 
Topt.Sex.Season <- emmeans(To.mod, pairwise ~ Sex |Season) 
# extracting season sex contrasts
Topt.Sex.Season.df <- as.data.frame(Topt.Sex.Season$emmeans) %>% 
  select(c(Season, emmean, SE, Sex, lower.CL, upper.CL))
# extracting seasonal contrasts
Topt.Season <- emmeans(To.mod, pairwise ~ Season)
Topt.Season.df <- as.data.frame(Topt.Season$contrasts)
#########

#########
# 4) Pmax analysis
Pmax.mod <- lmer(Topt  ~ Season + Sex + Season*Sex  + (1|POVI), Topt.final)
Pmax.sum <- anova(Pmax.mod) # no differences between sex, but significant seasona and interaction
Pmax.Sex.Season <- emmeans(Pmax.mod, pairwise ~ Sex |Season) 
# extracting season sex contrasts
Pmax.Sex.Season.df <- as.data.frame(Pmax.Sex.Season$emmeans) %>% 
  select(c(Season, emmean, SE, Sex, lower.CL, upper.CL))
# extracting seasonal contrasts
Pmax.Season <- emmeans(Pmax.mod, pairwise ~ Season)
Pmax.Season.df <- as.data.frame(Pmax.Season$contrasts)
#########



###########################################################################
################## (10) Analysis| plots | Tables ##########################
###########################################################################

#############
#Table S1
# Tsel table by sex
# differences overall Tb by sex in the lab? 
IQR.Summary <- readRDS(file = "R/Final.Analysis/Final.Figure.data/IQR.Final.RDS")

# Mean - significant females are higher
Mean.Tel <- lm(Mean.Tb ~Sex, data = IQR.Summary)
Mean.Tel.anova <- anova(Mean.Tel)
Mean.Tel.emmeans <- emmeans(Mean.Tel, pairwise ~ Sex)
Mean.Tsel.anova.pairwise <- emmeans(Mean.Tel, pairwise ~ Sex)
Mean.Tsel.dat <- as.data.frame(Mean.Tsel.anova.pairwise$emmeans) %>% 
  select(-c(df, lower.CL, upper.CL)) %>% 
  pivot_wider(names_from = Sex, values_from = c(emmean, SE)) %>% 
  rename("ZZf Mean" = emmean_ZZf, "ZWf Mean" = emmean_ZWf,"ZZm Mean" = emmean_ZZm, 
         "ZZf SE" = SE_ZZf, "ZWf SE" = SE_ZWf, "ZZm SE" = SE_ZZm) %>% 
  mutate("Gradient Index" = "Mean", "p value" = "p < 0.05") %>% 
  select("Gradient Index", "ZZf Mean", "ZZf SE","ZWf Mean", "ZWf SE","ZZm Mean", "ZZm SE", "p value")

# upper- NS: same across the board
IQR.upper <- lm(upper.IQR ~Sex,  data = IQR.Summary)
IQR.upper.anova <- anova(IQR.upper)
IQR.upper.emmeans <- emmeans(IQR.upper, pairwise ~ Sex)
IQR.upper.dat <- as.data.frame(IQR.upper.emmeans$emmeans) %>% 
  select(-c(df, lower.CL, upper.CL)) %>% 
  pivot_wider(names_from = Sex, values_from = c(emmean, SE)) %>% 
  rename("ZZf Mean" = emmean_ZZf, "ZWf Mean" = emmean_ZWf,"ZZm Mean" = emmean_ZZm, 
         "ZZf SE" = SE_ZZf, "ZWf SE" = SE_ZWf, "ZZm SE" = SE_ZZm) %>% 
  mutate("Gradient Index" = "75% quantile", "p value" = "p = 0.05") %>% 
  select("Gradient Index", "ZZf Mean", "ZZf SE","ZWf Mean", "ZWf SE","ZZm Mean", "ZZm SE", "p value")

# lower- significant females are higher
IQR.lower <- lm(lower.IQR ~Sex,  data = IQR.Summary)
IQR.lower.anova <- anova(IQR.lower)
IQR.lower.emmeans <- emmeans(IQR.lower, pairwise ~ Sex,)
IQR.lower.dat <- as.data.frame(IQR.lower.emmeans$emmeans) %>% 
  select(-c(df, lower.CL, upper.CL)) %>% 
  pivot_wider(names_from = Sex, values_from = c(emmean, SE)) %>% 
  rename("ZZf Mean" = emmean_ZZf, "ZWf Mean" = emmean_ZWf,"ZZm Mean" = emmean_ZZm, 
         "ZZf SE" = SE_ZZf, "ZWf SE" = SE_ZWf, "ZZm SE" = SE_ZZm) %>% 
  mutate("Gradient Index" = "25% quantile", "p value" = "p < 0.05") %>% 
  select("Gradient Index", "ZZf Mean", "ZZf SE","ZWf Mean", "ZWf SE","ZZm Mean", "ZZm SE", "p value")
# final table
Tsel.Table <- rbind( IQR.lower.dat, Mean.Tsel.dat, IQR.upper.dat)
Tsel.Table
saveRDS(Tsel.Table, file = "Final.Analysis/Final.Figure.data/Tsel.table.RDS")
#############


#############
# Fig 1 - Tb within Tset
# Tb vs Te within Tsel across seasons
# Tb by season and sex
Tb.mod <- lmer(Tb  ~ Sex + Season + Season*Sex + (1|POVI), db.sex.all.final)
Tb.mod.sum <- anova(Tb.mod) # NS Sex; Significant interaction and season
Tb.anova.pairwise <- emmeans(Tb.mod, pairwise ~ Sex |Season) # Summer differences
Tb.Sex.Season.df <- as.data.frame(Tb.anova.pairwise$emmeans) %>% 
  select(c(Season, emmean, SE, Sex, asymp.UCL , asymp.LCL))
Tb.Sex.Season.df <- Tb.Sex.Season.df %>% 
  rename(upper.CL = asymp.UCL, lower.CL = asymp.LCL)
  
# Te by season and sex
Te.dat <- Te %>% group_by(Season, Model_.) %>% summarise(Te = mean(Te))
Te.mod <- lmer(Te  ~ Season  + (1|Model_.), Te.dat)
Te.mod.sum <- anova(Te.mod) # All different
Te.anova.pairwise <- emmeans(Te.mod, pairwise ~ Season) 
Te.Sex.Season.df <- as.data.frame(Te.anova.pairwise$emmeans) %>% 
  mutate(Sex = "Te") %>% 
  select(c(Season, emmean, SE, Sex, lower.CL, upper.CL))
# final df for figure
Tb.vs.Te.fig.dat <- rbind(Te.Sex.Season.df, Tb.Sex.Season.df)
Tb.vs.Te.fig.dat <- mutate(Tb.vs.Te.fig.dat, Season = factor(Season,
                              levels = c("Spring",
                                         "Summer",
                                         "Autumn",
                                         "Winter"))) %>%
  group_by(Season)

# Final figure
# Sex      IQR  Mean.Tb    range upper.IQR lower.IQR
# ZZm 6.700000 28.68080 3.350000  32.03080  25.33080
# bring in activity data
Male.Tb.Te <- Tb.vs.Te.fig.dat %>% filter(Sex != "ZWf")
Male.Tb.vs.Te.fig <- ggplot(Male.Tb.Te, aes(x=Season, y=emmean, colour=Sex, group=Sex, shape = Sex)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), colour="black", width=.1) +
  geom_line(aes(linetype=Sex)) +
  geom_point(aes(shape=Sex, color=Sex), size=2)+
  scale_shape_manual(values=c(8, 21))+
  scale_color_manual(values=c("black","black")) +
  geom_hline(yintercept=32.03080, linetype="dashed", color = "grey34")+
  geom_hline(yintercept=25.33080, linetype="dashed", color = "grey34")+ 
  scale_y_continuous(breaks=c(16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38))+
  coord_cartesian(ylim = c(16, 38))+
  labs(x = "Season", y = "Temperature (°C)")+
  theme_bw() +
  theme(legend.position = c(.85, .9))+ 
  theme(legend.title = element_blank())  

# ZWf 6.833333 30.38976 3.416667  33.80642  26.97309
Female.Tb.Te <- Tb.vs.Te.fig.dat %>% filter(Sex != "ZZm")
Female.Tb.vs.Te.fig <- ggplot(Female.Tb.Te, aes(x=Season, y=emmean, colour=Sex, group=Sex, shape = Sex)) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), colour="black", width=.1) +
  geom_line(aes(linetype=Sex)) +
  geom_point(aes(shape=Sex, color=Sex), size=2)+
  scale_shape_manual(values=c(8, 19))+
  scale_color_manual(values=c("black", "#999999")) +
  geom_hline(yintercept=33.80642, linetype="dashed", color = "grey34")+
  geom_hline(yintercept=26.97309, linetype="dashed", color = "grey34")+ 
  scale_y_continuous(breaks=c(16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38))+
  coord_cartesian(ylim = c(16, 38))+
  labs(x = "Season", y = "Temperature (°C)")+
  theme_bw() +
  theme(legend.position = c(.85, .9))+ 
  theme(legend.title = element_blank())

# Fig1 Final: Tb.fig within tset
Fig.1.Tb.Te.figure <- ggpubr::ggarrange(Male.Tb.vs.Te.fig, Female.Tb.vs.Te.fig,
                    labels = c("A - Male", "B - Female"),
                    font.label = (size = 12), 
                    hjust = - 0.95,
                    vjust = 2.0 ,
                    ncol = 2, nrow = 1)
Fig.1.Tb.Te.figure


#############
# Fig 2 - Accuracy of thermoregulation 
# db by season and sex by individual
E_final_season <- readRDS(file = "R/Final.Analysis/Final.Figure.data/E_season_sex_individual.RDS")
db.mod <- lmer(mean_individual_mean_db  ~ Sex + Season + Season*Sex + (1|POVI), E_final_season)
db.mod.Sum <- anova(db.mod) # season and season*sex sig 
db.mod.anova.pairwise <- emmeans(db.mod, pairwise ~ Sex |Season) 
season.db <- as.data.frame(db.mod.anova.pairwise$emmeans)# significant difference spring

Fig.2.Accuracy.Thermoregulation <- ggplot(season.db, aes(x=Season, y=emmean, colour=Sex, group=Sex, shape = Sex)) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), colour="black", width=.1) +
  geom_line(aes(linetype=Sex)) +
  geom_point(aes(shape=Sex, color=Sex), size=2)+
  scale_shape_manual(values=c(19, 21))+
  scale_color_manual(values=c("#999999", "black")) +
  labs(x = "Season", y = "db")+
  theme_bw() +
  theme(legend.position = c(.95, .93))+ 
  theme(legend.title = element_blank())
Fig.2.Accuracy.Thermoregulation


#############
##### Fig 3 - E index by season and se
E_final_season <- readRDS(file = "R/Final.Analysis/Final.Figure.data/E_season_sex_individual.RDS")
E.mod <- lmer(indvidual_season_E  ~ Sex + Season + Season*Sex + (1|POVI), E_final_season)
E.mod.Sum <- anova(E.mod) # significant all around
E.anova.pairwise <- emmeans(E.mod, pairwise ~ Sex |Season)
Season <- as.data.frame(E.anova.pairwise$emmeans)
# final Fig 4
Fig.4.E.index.fig <- ggplot(Season, aes(x=Season, y=emmean, fill=Sex)) + 
  scale_y_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+
  coord_cartesian(ylim = c(-0.3, 1.0))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=emmean -SE, ymax=emmean +SE), width=.2,
                position=position_dodge(.9))+
  
  labs(x="Season", y = "E index")+
  theme_classic() +
  scale_fill_manual(values=c("#999999", "white"))
Fig.4.E.index.fig



#############
##### Fig 4 - Tb, Te, Activity - grand means -----
# males
# Tb data
Tb.predict.dat <- readRDS(file = "R/Final.Analysis/Final.Figure.data/Final.Activity.mins.data.RDS") %>% filter(Sex != "ZZf")
Tb.predict.dat$Hour <- as.factor(Tb.predict.dat$HR)
Male.Tb <- Tb.predict.dat%>% filter(Sex != "ZWf") %>% 
  group_by(Date, Hour) %>% 
  ungroup() %>% 
  group_by(Hour, Season, Sex) %>% 
  summarise(Tb = mean(Tb),
            HR = unique(HR))
Male.Tb$Season <- factor(Male.Tb$Season, c("Spring", "Summer", "Autumn", "Winter"))
# activity data
Male.Activity <- mydata.final.mm %>% filter(Sex != "ZWf")
Male.Activity$Season <- factor(Male.Activity$Season, c("Spring", "Summer", "Autumn", "Winter"))
# Te data
Te.data.fig <- Te %>% 
  group_by(hr, Season) %>% 
  summarise(Te = mean(Te),
            HR = unique(HR))
Te.data.fig$Season <- factor(Te.data.fig$Season, c("Spring", "Summer", "Autumn", "Winter"))

# male plot
male.tb.te.act.fig <- ggplot() +  
  geom_point(data = Male.Tb, aes(HR,Tb), 
             fill = "white", color = "black",
             size = 2.2, shape = 21)+
  geom_point(data = Male.Activity, aes(HR,Movement.HR), 
             fill = "black", color = "black", 
             size = 2, shape = 17)+
  geom_point(data = Te.data.fig, aes(HR,Te), 
             fill = "black", color = "black", 
             size = 2, shape = 8) +
  geom_line(data = Male.Tb, aes(HR,Tb), size = 0.3)+
  geom_line(data = Male.Activity, aes(HR,Movement.HR), size = 0.3)+
  geom_line(data = Te.data.fig, aes(HR,Te), size = 0.3)+
  geom_hline(yintercept=32.03080, linetype="dashed", color = "grey34")+
  geom_hline(yintercept=25.33080, linetype="dashed", color = "grey34")+ 
  labs(x = NULL, y = "Temperature (°C) or Activity (min)") +
  facet_grid(~Season) +
  theme_bw()

# female
# data
Female.Tb <- Tb.predict.dat%>% filter(Sex != "ZZm") %>% 
  group_by(Date, Hour) %>% 
  ungroup() %>% 
  group_by(Hour, Season, Sex) %>% 
  summarise(Tb = mean(Tb),
            HR = unique(HR))
Female.Tb$Season <- factor(Female.Tb$Season, c("Spring", "Summer", "Autumn", "Winter"))
# activity data
Female.Activity <- mydata.final.mm %>% filter(Sex != "ZZm")
Female.Activity$Season <- factor(Female.Activity$Season, c("Spring", "Summer", "Autumn", "Winter"))
female.tb.te.act.fig <- ggplot() +  
  geom_point(data = Female.Tb, aes(HR,Tb), 
             fill = "#999999", color = "black",
             size = 2.2, shape = 21)+
  geom_point(data = Female.Activity, aes(HR,Movement.HR), 
             fill = "black", color = "black", 
             size = 2, shape = 17)+
  geom_point(data = Te.data.fig, aes(HR,Te), 
             fill = "black", color = "black", 
             size = 2, shape = 8) +
  geom_line(data = Female.Tb, aes(HR,Tb), size = 0.3)+
  geom_line(data = Female.Activity, aes(HR,Movement.HR), size = 0.3)+
  geom_line(data = Te.data.fig, aes(HR,Te), size = 0.3)+
  geom_hline(yintercept=33.80642, linetype="dashed", color = "grey34")+
  geom_hline(yintercept=26.97309, linetype="dashed", color = "grey34")+ 
  labs(x = "Time of day (hr)", y = "Temperature (°C) or Activity (min)") +
  facet_grid(~Season) +
  theme_bw()


#############
# Fig4 Final: Tb, Te, Activity
Fig.5.Tb.Te.Activity.fig<- ggpubr::ggarrange(male.tb.te.act.fig, female.tb.te.act.fig,
                                  font.label = (size = 12), 
                                  hjust = - 0.65,
                                  vjust = 2.0 ,
                                  labels = c( "A - Male", "B - Female"),
                                  ncol = 1, nrow = 2)

Fig.5.Tb.Te.Activity.fig



#############
##### Table 2- GAMM table
# db by season and sex by individual
Temp <- readRDS(file = "R/Models/Temp.rds")
Temp.Povi <- readRDS(file = "R/Models/Temp.Povi.rds")
Sex <- readRDS(file = "R/Models/Sex.rds")
Sex.Povi <- readRDS(file = "R/Models/Sex.Povi.rds")
Season <- readRDS(file = "R/Models/Season.rds")
Season.Povi <- readRDS(file = "R/Models/Season.Povi.rds")
Season.Sex <- readRDS(file = "R/Models/Season.Sex.rds")
Season.Sex.Povi <- readRDS(file = "R/Models/Season.Sex.Povi.rds")
Season.Sex.Interaction <- readRDS(file = "R/Models/Season.Sex.Interaction.rds")
Season.Sex.Interaction.Povi <- readRDS(file = "R/Models/Season.Sex.Interaction.Povi.rds")

# # make a table of comparing all models to model 1. 
# i.e. does adding terms significantly change our understanding of the data set?
# ANOVA Table
anova.table <- rbind(anova(Temp,Temp.Povi, test="Chisq"),
                     anova(Temp,Sex, test="Chisq")[2,],
                     anova(Temp,Sex.Povi, test="Chisq")[2,],
                     anova(Temp,Season, test="Chisq")[2,],
                     anova(Temp,Season.Povi, test="Chisq")[2,],
                     anova(Temp,Season.Sex, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Povi, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Interaction, test="Chisq")[2,],
                     anova(Temp,Season.Sex.Interaction.Povi, test="Chisq")[2,])

anova.table <- cbind(anova.table, AIC(Temp, Temp.Povi, Sex, Sex.Povi, Season, Season.Povi, Season.Sex, Season.Sex.Povi, Season.Sex.Interaction, Season.Sex.Interaction.Povi)[,2])
anova.table <- cbind(c("Temp","Temp.Povi","Sex","Sex.Povi","Season","Season.Povi","Season.Sex","Season.Sex.Povi","Season.Sex.Interaction","Season.Sex.Interaction.Povi"), anova.table)
anova.table[,c(2:5,7)] <- round(anova.table[,c(2:5,7)],2)
anova.table[,6] <- signif(anova.table[,6],2)
rownames(anova.table) <- NULL
anova.table$DevEx <- c(summary(Temp)$dev.expl, summary(Temp.Povi)$dev.expl,
                       summary(Sex)$dev.expl, summary(Sex.Povi)$dev.expl,
                       summary(Season)$dev.expl, summary(Season.Povi)$dev.expl,
                       summary(Season.Sex)$dev.expl, summary(Season.Sex.Povi)$dev.expl,
                       summary(Season.Sex.Interaction)$dev.expl, summary(Season.Sex.Interaction.Povi)$dev.expl)
anova.table$"Dev.Expl(%)" <-  round(anova.table$DevEx*100, 2)
names(anova.table)[c(1:4,7, 8)] <- c("Model", "Residual DF", "Residual Deviance", "DF", "AICc", "Deviance Explained" )
anova.table$`Deviance Explained` <- round(anova.table$`Deviance Explained`, 2)

anova.table %>% arrange(AICc)
view(anova.table) 


#############
##### Fig6 - GAM with interaction
anova(Season.Sex.Interaction.Povi, test="Chi")
summary.gam(Season.Sex.Interaction.Povi)
gam.check(Season.Sex.Interaction.Povi)
data.gam$fit <- predict(Season.Sex.Interaction.Povi, type = "response")
data.gam$se <- predict(Season.Sex.Interaction.Povi, type = "se")
data.gam$se
# set up data
data.gam$Season <- factor(data.gam$Season, c("Spring", "Summer", "Autumn", "Winter"))
Spring <- data.gam %>% filter(Season == "Spring")
Spring.plot <- ggplot(Spring, aes(x=temp, y=percent95.ms2, color=Sex))+  
  geom_point(aes(shape=Sex, color=Sex), size=2, alpha =.6)+
  scale_shape_manual(values=c(21, 19))+
  scale_color_manual(values=c("#999999", "grey11")) +
  ylab(bquote("95th"~"percentile"~"of"~"acceleration ms"^-2)) +
  xlab(expression(T[b~Predict]))+
  coord_cartesian(xlim =c(5, 55), ylim = c(0, 2.8))+
  scale_x_continuous(breaks=seq(0,55,5))+
  scale_y_continuous(breaks = seq(0, 2.8, .2))+
  geom_smooth(method = "gam")+ 
  theme_bw()+
  ggtitle("Spring")+
  theme(legend.position="none",
        plot.title = element_text(size=14, face="bold", hjust = 0.5))
####### Summer 
Summer <- data.gam %>% filter(Season == "Summer")
Summer.plot <- ggplot(Summer, aes(x=temp, y=percent95.ms2, color=Sex))+
  geom_point(aes(shape=Sex, color=Sex), size=2, alpha =.6)+
  scale_shape_manual(values=c(21, 19))+
  scale_color_manual(values=c("#999999", "grey11")) +
  xlab(expression(T[b~Predict]))+
  coord_cartesian(xlim =c(5, 55), ylim = c(0, 2.8))+
  scale_x_continuous(breaks=seq(0,55,5))+
  scale_y_continuous(breaks = seq(0, 2.8, .2))+
  geom_smooth(method = "gam")+ 
  theme_bw()+
  ggtitle("Summer") +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))
Summer.plot
####### Autumn 
Autumn <- data.gam %>% filter(Season == "Autumn")
Autumn.plot <- ggplot(Autumn, aes(x=temp, y=percent95.ms2, color=Sex))+
  geom_point(aes(shape=Sex, color=Sex), size=2, alpha =.6)+
  scale_shape_manual(values=c(21, 19))+
  scale_color_manual(values=c("#999999", "grey11")) +
  xlab(expression(T[b~Predict])) +
  coord_cartesian(xlim =c(5, 55), ylim = c(0, 2.8))+
  scale_x_continuous(breaks=seq(0,55,5))+
  scale_y_continuous(breaks = seq(0, 2.8, .2))+
  geom_smooth(method = "gam")+ 
  theme_bw()+
  ggtitle("Autumn") +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))
Autumn.plot
####### Winter 
Winter <- data.gam %>% filter(Season == "Winter")
Winter.plot <- ggplot(Winter, aes(x=temp, y=percent95.ms2, color=Sex))+
  geom_point(aes(shape=Sex, color=Sex), size=2, alpha =.6)+
  scale_shape_manual(values=c(21, 19))+
  scale_color_manual(values=c("#999999", "grey11")) +
  xlab(expression(T[b~Predict])) +
  coord_cartesian(xlim =c(5, 55), ylim = c(0, 3.0))+
  scale_x_continuous(breaks=seq(0,55,5))+
  scale_y_continuous(breaks = seq(0, 3.0, .2))+
  geom_smooth(method = "gam")+ 
  theme_bw()+
  ggtitle("Winter") +
  theme(axis.title.y = element_blank(),
        legend.position = c(.80, .80), 
        plot.title = element_text(size=12, face="bold", hjust = 0.5))
Winter.plot

# all plots
p<- cowplot::plot_grid(Spring.plot, Summer.plot, Autumn.plot, Winter.plot, ncol=4)
p


library(survminer)
# cairo_ps(filename = "Performance-curves.pdf",
         #width = 14, height = 8, pointsize = 12,
         #fallback_resolution = 300)
print(p)
dev.off()






# individuals by season
ggplot(data = data.gam) + 
  geom_line(aes(x = temp, y =fit , group = POVI, color= Sex))+
  scale_x_continuous(breaks=seq(0,55,10))+
  scale_y_continuous(breaks = seq(0, 3.0, .2))+
  facet_wrap(~Season)+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))



















##########################################################
##########################################################
##########################################################
############# OLD CODE CHUNCKS  ##########################
##########################################################
##########################################################
##########################################################
# filter out dates for MARK data for later
# we want spring
Mark.Spring <- accel.all %>% filter(Season)


# CTMin and Max - eyeball, will fix this later
fig <- plot_ly(data = Summer, x = ~temp, y = ~percent95.ms2, color = ~Sex)


# de.male.sum <- de.male %>% 
#  group_by(Season) %>% 
# summarise(de = round(mean(resid), 3), 
# Mean.Te = round(mean(Te), 3))
# de.male.winter <- de.male %>% filter(Season == "Winter") %>% mutate(de_season_value = 13.9)
# de.male.spring <- de.male %>% filter(Season == "Spring") %>% mutate(de_season_value = 6.54)
# de.male.summer <- de.male %>% filter(Season == "Summer") %>% mutate(de_season_value = 4.89)
# de.male.autumn <- de.male %>% filter(Season == "Autumn") %>% mutate(de_season_value = 6.9)
# de.male.final <- rbind(de.male.winter, de.male.spring, de.male.summer, de.male.autumn)

# use this for season averages for later
# de.female.sum <- de.female %>% 
#  group_by(Season) %>% 
#  summarise(de = round(mean(resid), 3), 
# Mean.Te = round(mean(Te), 3))
# de.female.winter <- de.female %>% filter(Season == "Winter") %>% mutate(de_season_value = 13.4)
# de.female.spring <- de.female %>% filter(Season == "Spring") %>% mutate(de_season_value = 6.25)
# de.female.summer <- de.female %>% filter(Season == "Summer") %>% mutate(de_season_value = 4.03)
# de.female.autumn <- de.female %>% filter(Season == "Autumn") %>% mutate(de_season_value = 6.68)


####################  
# Ex by (individual or Models) by individual and season
####################  
# Calculating Total Hr's of the day of Tb's within Tset and averaging across individuals and season
Tb_AVG.HRs.Within.Tsel <- db.sex.all.final %>%
  group_by(Date, POVI, Sex, Season) %>% 
  summarise(Ind_proportion_within = mean(Ex)*24) %>% 
  ungroup() %>% 
  group_by(Sex, Season) %>% 
  summarise(Mean_Hr_within = mean(Ind_proportion_within),
            SE_Hr_within = std.error(Ind_proportion_within))


# Calculating Total Hr's of the day of Tb's within Tset and averaging across individuals and season
Te_AVG.HRs.Within.Tsel <- Te %>%
  group_by(Model_., Type, Season) %>% 
  summarise(Mod_proportion_within = mean(Ex)*24) %>% 
  ungroup() %>% 
  group_by(Season) %>% 
  summarise(Mean_Hr_within = mean(Mod_proportion_within),
            SE_Hr_within = std.error(Mod_proportion_within))









# 1) individual: Topt for individuals for mark rapture data for just POVI and temp
# name: Temp.Povi
POVI <- data.gam
gam.check(Temp.Povi)
POVI$fit <- predict(Temp.Povi, type = "response")
POVI$se <- predict(Temp.Povi, type = "se")
POVI$se

# add in death data column
Dead <- read.csv(file = "MARK_Data/Dead.Animals.csv") %>% 
  dplyr::select(c(Dead, POVI))
POVI <- left_join(POVI, Dead, by = "POVI")


# quick plot for individuals
ggplot(data = POVI) + 
  geom_line(aes(x = temp, y =fit , group = POVI, color= Dead))+
  scale_x_continuous(breaks=seq(0,55,10))+
  scale_y_continuous(breaks = seq(0, 3.0, .2))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))

# Topt for POVI 
POVI.Topt <- POVI %>% 
  group_by(POVI) %>% 
  slice(To = which.max(fit))

# Tb of opt plot from summary stats
ggplot(POVI.Topt, aes(x = Dead, y = Tb, color = Dead)) +
  geom_jitter(shape = 16, size = 2, alpha = 0.4, width = 0.2) +
  geom_boxplot(alpha = 0, size = 0.75, width = 0.4) +
  facet_grid(~Sex)+
  theme_bw() +
  theme(panel.grid = element_blank())

# Ms2 plot of opt temp plot from summary stats
ggplot(POVI.Topt, aes(x = Dead, y = fit)) +
  geom_jitter(shape = 16, size = 2, alpha = 0.4, width = 0.2) +
  geom_boxplot(alpha = 0, size = 0.75, width = 0.4) +
  theme_bw() +
  facet_grid(~Sex)+
  theme(panel.grid = element_blank())

  






summary(lm(Pref ~ Sex + Season + Season*Sex, data = To.final))


# TbPredict average
Temp.season.avg <- pref.other.95.final %>% 
  group_by(Season, Sex) %>% 
  summarise(SE = std.error(Tb), 
            Tb = mean(Tb))

ggplot(Temp.season.avg, aes(x=Sex, y=Tb, group=Sex, color=Sex)) + 
  geom_pointrange(aes(ymin=Tb-SE, ymax=Tb+SE))+
  facet_grid(~Season) +
  ylab("Avg body temperature") +
  xlab("Sex") 


####################  
#####  Field accelerometer data for avg mins moved
####################
# set wd
setwd("~/Dropbox/PhD Pogona/Dissertation/4.Thermal_chapter/R")
mydata <- readRDS(file = "Ta/Final.Activity.mins.data.RDS") %>% 
  filter(Sex !="ZZf")

# arrange season for plots  
mydata <-mutate(mydata, Season = factor(Season, levels = c("Spring","Summer", "Autumn", "Winter"))) %>%
  group_by(Season)

# activity (mins) function
act <-function(x)
{for (i in 1:length(x)) {
  if(x[i] <= 0) 
  { x [i]<- 0
  }
  else
  {x[i] <- 1}}
  return(x)
}

# min moved
mydata$active.mins <- act(mydata$acc.min.sum)
max(mydata$active.mins)
mydata.final <- mydata %>%  
  group_by(Date, HR, POVI) %>% 
  summarise(freq = sum(active.mins),
            temp = mean(Tsurf),
            Sex = unique(Sex),
            Season = unique(Season))
mean(mydata.final$freq)
# by season, sex hr
mydata.final.mm <- mydata.final %>% 
  group_by(Season, Sex, HR) %>% 
  summarise(Tsurf = mean(temp),
            Movement.HR = mean(freq),
            SE = std.error(freq))

# movement plots by season
ggplot(mydata.final.mm, aes(x=HR, y=Movement.HR, group=Sex, color=Sex)) + 
  geom_pointrange(aes(ymin=Movement.HR-SE, ymax=Movement.HR+SE))+
  facet_grid(~Season) +
  ylab("Avg number of active mins") +
  xlab("Hour of day") 

# with gam smooth
ggplot(mydata.final.mm, aes(x=HR, y=Movement.HR, color=Sex))+
  geom_point()+
  stat_smooth(geom = "smooth")+
  xlim(5, 22)+
  facet_grid(~Season) +
  ylab("Avg number of active mins") +
  xlab("Hour of day") 

# summarise by POVI
final.POVI.Activity <- mydata.final %>% 
  group_by(POVI, Sex) %>% 
  summarise(Avg.mins.moved = mean(freq))
final.POVI.Activity

##########
# merge with MARK data
##########
MARK <- read.csv("MARK_Data/Mark.Season.csv")
MARK
Thermal.Data
Thermal.Mark.Final <- left_join(MARK, Thermal.Data, by = "POVI") %>% 
  drop_na() %>% 
  format(digits = 2)
# add activity line
Thermal.Mark.Activity.Final<- left_join(Thermal.Mark.Final, final.POVI.Activity, by = "POVI") %>% format(digits = 2)


# save data and export out for mark analysis
write.csv(Thermal.Mark.Activity.Final, file = "MARK_Data/Thermal.Mark.Final.csv")

##########
# Season.Sex.Interaction.Povi for performance across seasons, not used in mark analysis
##########
gam.check(Season.Sex.Interaction.Povi)
data.gam$fit <- predict(Season.Sex.Interaction.Povi, type = "response")
data.gam$se <- predict(Season.Sex.Interaction.Povi, type = "se")
data.gam$se
# set up data
data.gam$Season <- factor(data.gam$Season, c("Spring", "Summer", "Autumn", "Winter"))

ggplot(data = data.gam) + 
  geom_line(aes(x = temp, y =fit , group = POVI, color= Sex))+
  scale_x_continuous(breaks=seq(0,55,10))+
  scale_y_continuous(breaks = seq(0, 3.0, .2))+
  facet_wrap(~Season)+
  theme_bw()+
  theme(legend.position="none",
        axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))
