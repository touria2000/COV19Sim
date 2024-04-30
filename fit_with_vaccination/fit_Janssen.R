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
V0 = 25460
S0 = pop - (E0 + IA0 + IS0 + IQ0 + IH0 + R0 + D0 + V0)
y0 = c(S0, E0, IA0, IS0, IQ0, IH0, R0, D0, V0)

# Put data into list
data = list(n_days = n_days, n_difeq=9, y_init = y0, t0 = t0, ts = t, cases = cases)

## Compile the model
m_j <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/models/model_Janssen.stan")

## sampling 
time.start_fit <- Sys.time()
fit_Janssen <- sampling(m_j, 
                      data = data, 
                      chains = 4, 
                      iter = 2000, 
                      seed=2022)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit

#### Check for divergent transitions
rstan::check_divergences(fit_Janssen)
check_hmc_diagnostics(fit_Janssen)

# Save the model
save(fit_Janssen, file = 'C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/saved_models/fit_Janssen')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/saved_models/fit_Janssen')

#### Print the parameter estimates
pars1 = c("beta_S", "beta_A", "gamma_S", "gamma_A", "gamma_H", "gamma_Q", "delta_S", "delta_Q", "delta_H",
          "sigma_S", "sigma_Q","eta" , "theta_V", "phi_inv", "R_0","lp__")
print(fit_Janssen, pars = pars1)
print(fit_Janssen, pars = pars1, digits = 4, format = "f")

pars2 = c("recovery_time_S", "recovery_time_A", "recovery_time_H", "recovery_time_Q", "phi")
print(fit_model, pars = pars2, digits = 4, format = "f")

pars3 = c("E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0", "V0")
print(fit_Janssen, pars = pars3)
print(fit_Janssen, pars = pars3, digits = 4, format = "f")

library(bayesplot)
color_scheme_set("mix-brightblue-red")
mcmc_trace(fit_Janssen, pars = "theta_V") 

#### Post analysis
# Extract posterior distributions for parameters and initial conditions
posts_fit <- rstan::extract(fit_Janssen)
posts_fit <- extract(fit_Janssen)


########## Checking the model
#### A posterior predictive check indicates if the model matches the data or not

#date format:yyyy/m/d
di = as.Date("2021-06-04")
df = as.Date("2021-11-26")
Date = seq(di, df, 1)

smr_pred <- cbind(as.data.frame(summary(fit_Janssen, 
                                        pars = "pred_cases", 
                                        probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975))$summary), 
                  Date=Date, cases = cases)

colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

blue75 <- rgb(0, 0, 255, max = 255, alpha = 50, names = "blue75")
blue97 <- rgb(0, 0, 255, max = 255, alpha = 25, names = "blue97")

plot_1 <- ggplot(smr_pred, mapping = aes(x = Date)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = "97.5% CrI")) +
  geom_ribbon(aes(ymin = X25., ymax = X75., fill = "50% CrI")) +
  geom_line(mapping = aes(x = Date, y = X50., color = "Model prediction"), linetype = 1, size = 0.8) +
  geom_point(mapping = aes(y = cases, shape = "COVID-19 data"), size = 0.8) +
  labs(x = "Date", y = "Daily confirmed cases")+
  scale_x_date(date_labels = "%b %d", date_breaks = "40 days")+
  scale_color_manual(values="blue")+
  scale_fill_manual(breaks = c("50% CrI", "97.5% CrI"),
                    values = c("#0000FF32", "#0000FF19"),
                    labels=c("50% CrI", "97.5% CrI")) +
  theme(axis.text=element_text(size=7), 
        axis.title=element_text(size=10,face="bold"),
        #legend.position = c(0.15, 0.65), legend.title = element_blank(),
        legend.position = "none", legend.title = element_blank(),
        legend.text = element_text(size =6),
        legend.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_1)

########## results of fitting the model
di = as.Date("2021-06-04")
df = as.Date("2021-11-26")
Date = seq(di, df, 1)

par <- lapply(t, function(i){sprintf("y_hat[%s,4]", i)}) #number of infected for each day
smr_y <-  cbind(as.data.frame(summary(fit_Janssen,
                                      pars = par, 
                                      probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975))$summary),
                Date=Date, cases = cases)
colnames(smr_y) <- make.names(colnames(smr_y)) # to remove % in the col names

plot_2 <- ggplot(smr_y, mapping = aes(x = Date)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = "97.5% CrI")) +
  geom_line(mapping = aes(x = Date, y = X50., color = "Model estimation"), linetype = 1, size = 0.8) + 
  geom_point(mapping = aes(y = cases, shape = "COVID-19 data"), size = 0.8) +
  labs(x = "Date", y = "Daily confirmed cases")+
  scale_x_date(date_labels = "%b %d", date_breaks = "40 days")+
  scale_color_manual(values="blue")+
  scale_fill_manual(breaks = "97.5% CrI",
                    values = "#0000FF19") +
  theme(axis.text=element_text(size=7), 
        axis.title=element_text(size=10,face="bold"),
        #legend.position=c(0.15, 0.70), legend.title = element_blank(),
        legend.position="none", legend.title = element_blank(),
        legend.text=element_text(size =6),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot_2)

## Export ggplot2 plots
## Subplots
library(cowplot)
plot_grid(plot_1, plot_2, labels = c('(a)', '(b)'), 
          label_size = 10, hjust = -0.1, vjust = 1.5)
ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/plots/Figure 12.pdf", 
       dpi = 600, width = 20, height = 7, units = "cm")

# Plot vaccination rate
#####
posterior1 <- extract(fit_Janssen, inc_warmup = TRUE, permuted = FALSE)

library(bayesplot)
color_scheme_set("mix-brightblue-red")
p <- mcmc_trace(posterior1,  pars = "theta_V", n_warmup = 1000,
                facet_args = list(ncol = 1, 
                                  strip.position = "left", 
                                  labeller = label_parsed))
p + facet_text(size = 10)

##
posterior2 <- as.matrix(fit_Janssen)

color_scheme_set("mix-brightblue-red")
pp <- mcmc_areas(posterior2, pars = "theta_V", prob = 0.9)
pp + facet_text(size = 12)

## Subplots
library(cowplot)
plot_row <- plot_grid(p, pp, labels = c('(a)', '(b)'), label_size = 10, hjust = -0.1, vjust = 1.5)
# now add the title
title <- ggdraw() + 
  draw_label("Janssen COVID-19 vaccine", size = 8, fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/plots/Figure 10.pdf", 
       dpi = 600, width = 16, height = 5, units = "cm")



