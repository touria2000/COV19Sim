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
S0 = pop - (E0 + IA0 + IS0 + IQ0 + IH0 + R0 + D0)
y0 = c(S0, E0, IA0, IS0, IQ0, IH0, R0, D0)

# Put data into list
data = list(n_days = n_days, n_difeq=8, y_init = y0, t0 = t0, ts = t, cases = cases)

## Compile the model 
m <- stan_model("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/models/model_v1.stan")

## sampling 
time.start_fit <- Sys.time()
fit_model <- sampling(m, 
                      data = data, 
                      chains = 4, 
                      iter = 2000, 
                      seed=2023)
time.end_fit <- Sys.time()
duration_fit <- time.end_fit - time.start_fit

#### Check for divergent transitions
rstan::check_divergences(fit_model)
check_hmc_diagnostics(fit_model)

# Save the model
save(fit_model, file='C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/saved_models/fit1')
# Load the model
load('C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/saved_models/fit1')

########## Checking the model
#### A posterior predictive check indicates if the model matches the data or not

#date format:yyyy/m/d
di = as.Date("2021-06-04")
df = as.Date("2021-11-26")
Date = seq(di, df, 1)

smr_pred <- cbind(as.data.frame(summary(fit_model, 
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

########## Plot a posterior predictive check for each chain
fit_model %>% 
  spread_draws(pred_cases[n_days]) %>% 
  left_join(tibble(cases = cases, n_days = 1:length(cases))) %>% 
  group_by(n_days, .chain) %>% 
  summarise(cases = mean(cases), pred_mean = mean(pred_cases), pred_9 = quantile(pred_cases, 0.95), pred_1 = quantile(pred_cases, 0.05)) %>% 
  ggplot(aes(x = n_days)) +
  geom_ribbon(aes(ymin = pred_1, ymax = pred_9), fill = "#0000FF19", alpha=0.2)+
  geom_line(mapping = aes(y=pred_mean), size=0.8, color = "blue")+
  geom_point(mapping = aes(y=cases), size=0.5)+
  labs(x = "Date (days)", y = "Predicted incidence")+
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~.chain)

ggsave("C:/Users/Benbrahim/Documents/Touria/R_Project/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/fit_chain.pdf", 
       dpi = 600, width = 18, height = 10, units = "cm")

########## Print the parameter estimates
pars1 = c("beta_S", "beta_A", "gamma_S", "gamma_A", "gamma_H", "gamma_Q", 
          "delta_S", "delta_Q", "delta_H","sigma_S", "sigma_Q","eta" , 
          "phi_inv", "phi", "R_0","lp__")
print(fit_model, pars = pars1)
print(fit_model, pars = pars1, digits = 4, format = "f")

pars2 = c("recovery_time_S", "recovery_time_A", "recovery_time_H", "recovery_time_Q")
print(fit_model, pars = pars2, digits = 4, format = "f")

pars3 = c("E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0")
print(fit_model, pars = pars2)
print(fit_model, pars = pars2, digits = 4, format = "f")

#
fit_model_summary <- summary(fit_model, pars = pars1)$summary
print(fit_model_summary, scientific = FALSE, digits = 5)

########## Densities
pars4 = c("beta_S", "beta_A", "gamma_S", "gamma_A", "gamma_H", "gamma_Q", 
          "delta_S", "delta_Q", "delta_H","sigma_S", "sigma_Q","eta" , 
          "phi_inv", "phi","lp__")

stan_dens(fit_model, pars = pars4, separate_chains = TRUE)

library(bayesplot)
posterior_1 <- as.array(fit_model)
#
color_scheme_set("brightblue")
mcmc_dens(posterior_1, pars = pars4)
#
color_scheme_set("mix-brightblue-red")
mcmc_dens_overlay(posterior_1, pars = pars4)

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/Figure 4.pdf", 
       dpi = 600, width = 30, height = 18, units = "cm")

########## Trace plot
traceplot(fit_model, pars = pars4)

library(bayesplot)
color_scheme_set("mix-brightblue-red")
mcmc_trace(fit_model, pars =pars4) 

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/Figure 5.pdf", 
       dpi = 600, width = 30, height = 18, units = "cm")

########## plotting CI for each parameter
pars5 = c("beta_S", "beta_A", "gamma_S", "gamma_A")
pars6 = c("gamma_H", "gamma_Q", "delta_S", "delta_Q", "delta_H",
          "sigma_S", "sigma_Q","eta", "phi_inv")

plot(fit_model, pars = pars5)
plot(fit_model, pars = pars6)

#### Central posterior uncertainty intervals
library(bayesplot)
posterior_1 <- as.array(fit_model)

# MCMC intervals
color_scheme_set("brightblue")
mcmc_intervals(posterior_1, pars = pars5)

color_scheme_set("brightblue")
mcmc_intervals(posterior_1, pars = pars6)

# R_hat
rhat(fit_model, pars = pars5)

# color by rhat value
color_scheme_set("brightblue")
fake_rhat_values <- c(1.001, 1.001, 0.999, 0.999)
mcmc_intervals(posterior_1, pars = pars5,
               rhat = fake_rhat_values)

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/Figure 2.pdf", 
       dpi = 600, width = 16, height = 6, units = "cm")

rhat(fit_model, pars = pars6)

# color by rhat value
color_scheme_set("brightblue")
fake_rhat_values <- c(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999)
mcmc_intervals(posterior_1, pars = pars6,
               rhat = fake_rhat_values)

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/Figure 3.pdf", 
       dpi = 600, width = 16, height = 10, units = "cm")

########## Basic reproduction number
library(bayesplot)

posterior2 <- extract(fit_model, inc_warmup = FALSE, permuted = FALSE)
color_scheme_set("mix-brightblue-red")
p <- mcmc_trace(posterior2,  pars = "R_0",
                facet_args = list(ncol = 1, 
                                  #nrow = 1,
                                  strip.position = "left", 
                                  labeller = label_parsed))
p + facet_text(size = 10)


posterior3 <- as.matrix(fit_model)
color_scheme_set("mix-brightblue-red")
pp <- mcmc_areas(posterior3, pars = "R_0", prob = 0.9)
pp + facet_text(size = 12)

## Subplots
library(cowplot)
plot_row <- plot_grid(p, pp, labels = c('(a)', '(b)'), label_size = 10, hjust = -0.1, vjust = 1.5)
# now add the title
title <- ggdraw() + 
  draw_label("Basic reproduction number", size = 8, fontface = 'bold', x = 0, hjust = 0) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))

ggsave("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/plots/Figure 7.pdf", 
       dpi = 600, width = 16, height = 5, units = "cm")


########## Post analysis
# Extract posterior distributions for parameters and initial conditions
posts_fit <- rstan::extract(fit_model)
posts_fit <- extract(fit_model)

# Save Estimate parameters and initial conditions
## Method_1
# Parameters
df_theta <- data.frame(beta_S = posts_fit$beta_S, 
                       beta_A = posts_fit$beta_A,
                       gamma_S = posts_fit$gamma_S, 
                       gamma_A = posts_fit$gamma_A, 
                       gamma_H = posts_fit$gamma_H, 
                       gamma_Q = posts_fit$gamma_Q,
                       delta_S = posts_fit$delta_S,
                       delta_Q = posts_fit$delta_Q,
                       delta_H = posts_fit$delta_H, 
                       sigma_S = posts_fit$sigma_S, 
                       sigma_Q = posts_fit$sigma_Q, 
                       eta = posts_fit$d_eta,
                       R_0 = posts_fit$R_0)

write.table(df_theta, 
            file = "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/estimation/theta_values.csv", 
            sep = ",", row.names = FALSE)

# Initial conditions
Mat_IC <- posts_fit$y_init
write.table(Mat_IC, 
            file = "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/estimation/IC_no.csv", 
            sep = ",", 
            col.names = c("S0", "E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0"),
            row.names = FALSE)

# Fit
df_fit <- data.frame(beta_S = posts_fit$beta_S, 
                     beta_A = posts_fit$beta_A,
                     gamma_S = posts_fit$gamma_S, 
                     gamma_A = posts_fit$gamma_A, 
                     gamma_H = posts_fit$gamma_H, 
                     gamma_Q = posts_fit$gamma_Q,
                     delta_S = posts_fit$delta_S,
                     delta_Q = posts_fit$delta_Q,
                     delta_H = posts_fit$delta_H, 
                     sigma_S = posts_fit$sigma_S, 
                     sigma_Q = posts_fit$sigma_Q, 
                     eta = posts_fit$d_eta,
                     phi_inv = posts_fit$phi_inv,
                     phi = posts_fit$phi,
                     log_post = posts_fit$lp__,
                     E0 = posts_fit$E0,
                     IA0 = posts_fit$IA0,
                     IS0 = posts_fit$IS0,
                     IQ0 = posts_fit$Q0,
                     IH0 = posts_fit$H0,
                     R0 = posts_fit$R0,
                     D0 = posts_fit$D0)

write.table(df_fit, file = "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/estimation/fit_estimates.csv", 
            sep = ",", row.names = FALSE)

## Method_2 
source("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/pars_estimates.R")

summary <- summary_statistics(fit_model)
df_theta <- df_theta(posts_fit)
Mat_IC <- Mat_IC(posts_fit)
df_fit <- df_fit(posts_fit)


