#######################
##### Field Tb vs Tb predict and Biophysical model comparison for:
##### 1) Tb predict vs Tb core field relationship
##### 2) Te carcass adjustment vs non-thermoregulating lizard NicheMapR
pacman::p_load(dplyr, ggplot2,lubridate, NicheMapR)

######################
### Bring in data of animals that had backpacks and inplanted ibuttons:
Tsurf.Tb <- read.csv('R/Final.Analysis/Final.Figure.data/Field_Tb_vs_Tb_predict_SI.csv')


###############
# Analysis
###############
mod <- lm(Tb_obs ~Tb_predict, data = Tsurf.Tb)
summary(mod) 
residualPlot(mod)
anova(mod)


# regression plot
# all individuals
RegressionALL <- ggplot(Tsurf.Tb, aes(Tb_obs, Tb_predict)) +
  geom_point() +
  geom_smooth(method = "lm")
RegressionALL
# accounting for individual
RegressionPOVI <-ggplot(Tsurf.Tb, aes(Tb_obs, Tb_predict, color = POVI)) +
  geom_point() +
  geom_smooth(method = "lm")
RegressionPOVI



######################
##### Question 1) Tb predict vs Tb estimation of animals w
# bring in hourly Tb predict data

######
# model and figure
data_field_tb <- read.csv('R/Final.Analysis/Final.Figure.data/Field_Tb_vs_Tb_predict_SI.csv')
mod <- lm(Tb_obs ~Tb_predict, data = data_field_tb)
summary_mod<-summary(mod) 
anova(mod)
# Extract 
slope <- summary_mod$coefficients[2, 1]
intercept <- summary_mod$coefficients[1, 1]
r_squared <- summary_mod$r.squared
p_value <- summary_mod$coefficients[2, 4]

# Plot
plot(data_field_tb$Tb_predict, data_field_tb$Tb_obs,
     xlim = c(15, 45), ylim = c(15, 45), 
     xlab = expression(T[b * ",predict"] ~ " temperature (°C)"), 
     ylab = expression(T[b * ",observed"] ~ " temperature (°C)"),
     pch = 16, col = rgb(0, 0, 0, 0.3))

# Add the one-to-one line
abline(0, 1, col = "red", lwd = 1)
# Add text for the linear regression equation, R-squared, and p-value
text(20, 40, paste0("Linear regression: y = ", round(slope, 2), "x + ", round(intercept, 2)), pos = 1)
text(20, 38, paste0(expression(R^2), " = ", round(r_squared, 2)), pos = 1)
text(20, 36, paste0("P-value < ", format.pval(p_value, digits = 3)), pos = 1)

######################
##### Question 2) Te carcass adjustment in open habitat vs non-thermoregulating lizard NicheMapR
micro<-readRDS('R/Final.Analysis/Final.Figure.data/micro_Te_Cunnamulla.RDS')

# Question 2: ecto_open_Cunnamulla - simulating predicted Te in open - dead animal
ecto2 <- ectotherm(live = 0,
                   Ww_g = 315, # mass of cylinder if filled with water
                   alpha_max = 0.7,
                   alpha_min = 0.7,
                   shape = 1,
                   postur = 0,
                   pct_cond = 25,
                   shape_b = 6.25,
                   fatosk = 0.5,
                   fatosb = 0.25,
                   pct_wet = 0)


Te <-read.csv(file = "R/Final.Analysis/Final.Figure.data/Te_biophsical_mod_dat.csv") 
Te$date_hr = ymd_hms(Te$Date_Time)
Te$hr = hour(Te$date_hr)
Te$Date = ymd(Te$Date)
Te <- Te %>% 
  filter(hr >= 5) %>% 
  filter(hr <= 21) 
  
########
### Open models and summarize hourly temps across all 10
# need to filter to first full day of measurement 
Te_open <- Te %>% 
  filter(Type == "Open") %>% 
  group_by(date_hr, Date, hr) %>% 
  filter(Date >= "2018-10-10") %>% # first full measurement of Te open
  summarise(Te = mean(Te),
            HR = unique(Time)) 
# final df for merge
te_open_final <- Te_open %>% 
  ungroup() %>% 
  dplyr::select(Te, date_hr, Date, hr) %>% 
  rename(temp = Te,
         date = Date) 

########
### Ecotherm dead  (i.e Te of animal dead in open habitat)
ecto_dead <- as.data.frame(ecto2$environ)
# micro habitat
dates <- as.data.frame(micro$dates)
# make df with dates and temps and arrange cols for filter
ecto_dead['dates'] <- dates
ecto_dead$hr <- lubridate::hour(ecto_dead$dates)
ecto_dead$date <- lubridate::date(ecto_dead$dates)
# filter dates and time that has all models
ecot_dead_final <- ecto_dead %>% 
  filter(date >= "2018-10-10" & date <= "2019-08-31")
ecot_dead_final$date = ymd(ecot_dead_final$date)
# final df before merge
ecto_open_final <- ecot_dead_final %>% 
  dplyr::select(TC, dates, hr, date) %>% 
  dplyr::rename(temp = TC,
                date_hr = dates) %>% 
  unique(hr = hr)

########
### Te predicted vs Te observed final merge 
open_observed_predicted <- dplyr::inner_join(x = te_open_final, y = ecto_open_final, by = c("date", "hr"))
# rename cols
open_observed_predicted <- open_observed_predicted %>% 
  rename(Te_observed = temp.x, 
         Te_predicted = temp.y,
         date_hr = date_hr.x)
######
# Figure
ggplot(open_observed_predicted, aes(x = Te_observed, y = Te_predicted)) +
  geom_point(alpha = 0.4)+
  geom_abline(intercept = 0, slope = 1, linewidth = 0.4)+
  xlim(-10, 80) +
  ylim(-10, 80) +
  theme_bw()


