#    
# Objective : estimation of the unknown epidemiological parameters
#                        
# Created by: Touria Jdid
#
# R version 4.3.1 
# 
library(tidyverse) 
library(dplyr)
library(tidybayes)
library(gridExtra)
library(ggplot2);
theme_set(theme_bw()) # Set the ggplot2 theme globally
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#### Incidence data
newcases_TN <- read.csv(file="C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/data/newcases_TN.csv", 
                        header=TRUE, sep=",")
cases <- newcases_TN$daily_confirmed_cases[450:625]

date <- seq(1,length(cases))
df_TN <- data.frame(date=date, cases = cases)

# Plotting COVID-19 data
df_TN %>% 
  ggplot() + 
  geom_point(mapping = aes(x = date, y = cases, color = "COVID-19 data"), size = 1.5) +
  scale_color_manual(values="black")+
  labs(y="Number of reported cases")+
  #theme(legend.position="bottom", legend.title = element_blank())+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"),
        legend.position=c(0.12, 0.8), legend.title = element_blank(),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))

# times
n_days = length(cases)
t0 = 0
#t = seq(1, n_days, by = 1)
t = 1:n_days

# Initial conditions
pop = 6910840 ## Total population in Tennessee
E0 = 35870
IA0 = 2701
IS0 = 15890
IQ0 = 18905
IH0 = 28686
R0 = 30417
D0 = 11372
S0 = pop - (E0 + IA0 + IS0 + IQ0 + IH0 + R0 + D0);
y0 = c(S0, E0, IA0, IS0, IQ0, IH0, R0, D0)

################
##### Scenario_1
pct1_train <- 0.30
n1_train <- floor(n_days * pct1_train)
n1_pred <- n_days - n1_train
times1_pred <- t[(n1_train + 1):n_days]
y1_train <- cases[1:n1_train]

# Put data into list
data1 <- list(n_days = n1_train, n_pred = n1_pred, 
              y_init = y0, cases = y1_train, 
              t0 = 0, ts = t[1:n1_train], ts_pred = times1_pred, n_difeq=8)

## Compile the model 
m1 <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/model_pred.stan")


## sampling 
time.start_fit <- Sys.time()
fit_model1 <- sampling(m1, 
                       data = data1,
                       chains = 4, 
                       iter = 2000, 
                       seed=2022)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit
 
# Save the model
save(fit_model1, file='C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred1')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred1')

#### Check for divergent transitions
rstan::check_divergences(fit_model1)
check_hmc_diagnostics(fit_model1)

##
plot_1 <- ggplot(NULL, aes(t, cases)) + 
  geom_ribbon(data = smr_y1, aes(ymin = X2.5., ymax = X97.5., fill = "Estimated 97.5% CrI")) +
  geom_line(data = smr_y1, mapping = aes(x = t, y = X50., color = "Model estimation"), linetype = 1, size = 0.5) +
  
  geom_ribbon(data = smr_y11_pred, aes(ymin = X2.5., ymax = X97.5., fill = "Predicted 97.5% CrI")) +
  geom_line(data = smr_y11_pred, mapping = aes(x = t, y = X50., color = "Model prediction"), linetype = 1, size = 0.5) +
  
  geom_point(mapping = aes(shape = "COVID-19 data"), size = 0.5)+
  geom_vline(xintercept = 52.5, linetype = 2, color= "red", size = 0.5) +
  #geom_vline(xintercept = 89, linetype = 2, color= "orange", size = 0.5) +
  
  labs(x = "Time (days)", y = "Daily infections")+
  #scale_x_date(date_labels = "%m-%Y")+
  #scale_x_date(date_labels = "%m-%d-%Y")+
  
  scale_color_manual(values=c("blue", "green"))+
  scale_fill_manual(breaks = c("Estimated 97.5% CrI", "Predicted 97.5% CrI"),
                    values = c("#0000FF19", "#00FF0019"),
                    labels=c("Estimated 97.5% CrI", "Predicted 97.5% CrI")) +
  
  #theme(legend.position="bottom", legend.title = element_blank())
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=8),
        legend.position="none", legend.title = element_blank(),
        legend.text=element_text(size = 5),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_1)

################
##### Scenario_2
pct2_train <- 0.50
n2_train <- floor(n_days * pct2_train)
n2_pred <- n_days - n2_train
times2_pred <- t[(n2_train + 1):n_days]
y2_train <- cases[1:n2_train]

# Put data into list
data2 <- list(n_days = n2_train, n_pred = n2_pred, 
              y_init = y0, cases = y2_train, 
              t0 = 0, ts = t[1:n2_train], ts_pred = times2_pred, n_difeq=8)

## Compile the model 
m2 <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/model_pred.stan")

## sampling 
time.start_fit <- Sys.time()
fit_model2 <- sampling(m2, 
                       data = data2, 
                       chains = 4, 
                       iter = 2000, 
                       seed=2022)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit

# Save the model
save(fit_model2, file='C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred2')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred2')

#### Check for divergent transitions
rstan::check_divergences(fit_model1)
check_hmc_diagnostics(fit_model1)

##
plot_2 <- ggplot(NULL, aes(t, cases)) + 
  geom_ribbon(data = smr_y2, aes(ymin = X2.5., ymax = X97.5., fill = "Estimated 97.5% CrI")) +
  geom_line(data = smr_y2, mapping = aes(x = t, y = X50., color = "Model estimation"), linetype = 1, size = 0.5) +
  
  geom_ribbon(data = smr_y22_pred, aes(ymin = X2.5., ymax = X97.5., fill = "Predicted 97.5% CrI")) +
  geom_line(data = smr_y22_pred, mapping = aes(x = t, y = X50., color = "Model prediction"), linetype = 1, size = 0.5) +
  
  geom_point(mapping = aes(shape = "COVID-19 data"), size = 0.5)+
  geom_vline(xintercept = 88.5, linetype = 2, color= "red", size = 0.5) +
  
  labs(x = "Time (days)", y = "Daily infections")+
  #scale_x_date(date_labels = "%m-%Y")+
  #scale_x_date(date_labels = "%m-%d-%Y")+
  
  scale_color_manual(values=c("blue", "green"))+
  scale_fill_manual(breaks = c("Estimated 97.5% CrI", "Predicted 97.5% CrI"),
                    values = c("#0000FF19", "#00FF0019"),
                    labels=c("Estimated 97.5% CrI", "Predicted 97.5% CrI")) +
  #theme(legend.position="bottom", legend.title = element_blank())
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=8),
        legend.position="none", legend.title = element_blank(),
        legend.text=element_text(size = 5),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_2)

################
##### Scenario_3
pct3_train <- 0.70
n3_train <- floor(n_days * pct3_train)
n3_pred <- n_days - n3_train
times3_pred <- t[(n3_train + 1):n_days]
y3_train <- cases[1:n3_train]

# Put data into list
data3 <- list(n_days = n3_train, n_pred = n3_pred, 
              y_init = y0, cases = y3_train, 
              t0 = 0, ts = t[1:n3_train], ts_pred = times3_pred, n_difeq=8)

## Compile the model 
m3 <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/model_pred.stan")

## sampling 
time.start_fit <- Sys.time()
fit_model3 <- sampling(m3, 
                       data = data3, 
                       chains = 4, 
                       iter = 2000, 
                       seed=2022)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit

# Save the model
save(fit_model3, file='C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred3')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred3')

#### Check for divergent transitions
rstan::check_divergences(fit_model1)
check_hmc_diagnostics(fit_model1)

##
plot_3 <- ggplot(NULL, aes(t, cases)) + 
  geom_ribbon(data = smr_y3, aes(ymin = X2.5., ymax = X97.5., fill = "Estimated 97.5% CrI")) +
  geom_line(data = smr_y3, mapping = aes(x = t, y = X50., color = "Model estimation"), linetype = 1, size = 0.5) +
  
  geom_ribbon(data = smr_y33_pred, aes(ymin = X2.5., ymax = X97.5., fill = "Predicted 97.5% CrI")) +
  geom_line(data = smr_y33_pred, mapping = aes(x = t, y = X50., color = "Model prediction"), linetype = 1, size = 0.5) +
  
  geom_point(mapping = aes(shape = "COVID-19 data"), size = 0.5)+
  geom_vline(xintercept = 123.5, linetype = 2, color= "red", size = 0.5) +
  
  labs(x = "Time (days)", y = "Daily infections")+
  #scale_x_date(date_labels = "%m-%Y")+
  #scale_x_date(date_labels = "%m-%d-%Y")+
  
  scale_color_manual(values=c("blue", "green"))+
  scale_fill_manual(breaks = c("Estimated 97.5% CrI", "Predicted 97.5% CrI"),
                    values = c("#0000FF19", "#00FF0019"),
                    labels=c("Estimated 97.5% CrI", "Predicted 97.5% CrI")) +
  #theme(legend.position="bottom", legend.title = element_blank())
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=8),
        legend.position="none", legend.title = element_blank(),
        legend.text=element_text(size = 5),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_3)

################
##### Scenario_4
pct4_train <- 0.90
n4_train <- floor(n_days * pct4_train)
n4_pred <- n_days - n4_train
times4_pred <- t[(n4_train + 1):n_days]
y4_train <- cases[1:n4_train]

# Put data into list
data4 <- list(n_days = n4_train, n_pred = n4_pred, 
              y_init = y0, cases = y4_train, 
              t0 = 0, ts = t[1:n4_train], ts_pred = times4_pred, n_difeq=8)

## Compile the model 
m4 <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/model_pred.stan")

## sampling 
time.start_fit <- Sys.time()
fit_model4 <- sampling(m4, 
                       data = data4, 
                       chains = 4, 
                       iter = 2000, 
                       seed=2022)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit

# Save the model
save(fit_model4, file='C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred4')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/fit_pred4')

#### Check for divergent transitions
rstan::check_divergences(fit_model1)
check_hmc_diagnostics(fit_model1)

##
plot_4 <- ggplot(NULL, aes(t, cases)) + 
  geom_ribbon(data = smr_y4, aes(ymin = X2.5., ymax = X97.5., fill = "Estimated 97.5% CrI")) +
  geom_line(data = smr_y4, mapping = aes(x = t, y = X50., color = "Model estimation"), linetype = 1, size = 0.5) +
  
  geom_ribbon(data = smr_y44_pred, aes(ymin = X2.5., ymax = X97.5., fill = "Predicted 97.5% CrI")) +
  geom_line(data = smr_y44_pred, mapping = aes(x = t, y = X50., color = "Model prediction"), linetype = 1, size = 0.5) +
  
  geom_point(mapping = aes(shape = "COVID-19 data"), size = 0.5)+
  geom_vline(xintercept = 158.5, linetype = 2, color= "red", size = 0.5) +
  
  labs(x = "Time (days)", y = "Daily infections")+
  #scale_x_date(date_labels = "%m-%Y")+
  #scale_x_date(date_labels = "%m-%d-%Y")+
  
  scale_color_manual(values=c("blue", "green"))+
  scale_fill_manual(breaks = c("Estimated 97.5% CrI", "Predicted 97.5% CrI"),
                    values = c("#0000FF19", "#00FF0019"),
                    labels=c("Estimated 97.5% CrI", "Predicted 97.5% CrI")) +
  #theme(legend.position="bottom", legend.title = element_blank())
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=8),
        legend.position="none", legend.title = element_blank(),
        legend.text=element_text(size = 5),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_4)

## Subplots
library(cowplot)
plot_grid(plot_1, plot_2, plot_3, plot_4, labels = c('(a)', '(b)', '(c)', '(d)'), 
          label_size = 10, hjust = -0.1, vjust = 1.5)
ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/Figure 14.pdf", 
       dpi = 600, width = 15, height = 10, units = "cm")


library(patchwork)
plot_1 + plot_2 + plot_3 + plot_4 + plot_layout(ncol = 2)
ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/prediction/Figure 14.pdf", 
       dpi = 600, width = 14, height = 9, units = "cm")


