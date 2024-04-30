#    
# Objective : estimation of the unknown epidemiological parameters
#                        
# Created by: Touria Jdid
#
# R version 4.3.1 
# 
library(deSolve)
library(ggplot2);
theme_set(theme_bw()) # Set the ggplot2 theme globally

##### Define ode system
SEIQHRD <- function(t, y, par){ 
  ## State variables
  S = y[1]
  E = y[2]
  IA = y[3]
  IS = y[4]
  IQ = y[5]
  IH = y[6]
  R = y[7]
  D = y[8]
 
  ## Model parameters
  alpha = par[1]
  rho = par[2]
  omega_R = par[3]
  m_H = par[4]
  m_Q = par[5]
  beta_S = par[6]
  beta_A = par[7]
  gamma_S = par[8]
  gamma_A = par[9]
  gamma_H = par[10]
  gamma_Q = par[11]
  delta_S = par[12]
  delta_Q = par[13]
  delta_H = par[14]
  sigma_S = par[15]
  sigma_Q = par[16]
  eta = par[17]
  
  # Forces of infection
  N = S + E + IA + IS + IQ + IH + R + D; 
  lmbda  = (beta_A * IA + beta_S * IS + m_Q*beta_S*IQ + m_H*beta_S*IH) / N; 
  
  # ODE system
  dS = omega_R*R - lmbda*S;
  dE = lmbda*S - alpha*E;
  dIA = rho*alpha*E - gamma_A*IA;
  dIS = (1 - rho)*alpha*E - (delta_S + sigma_S + gamma_S + eta)*IS;
  dIQ = eta*IS - (delta_Q + sigma_Q + gamma_Q)*IQ;
  dIH = sigma_S*IS + sigma_Q*IQ - (delta_H + gamma_H)*IH;
  dR = gamma_A*IA + gamma_Q*IQ + gamma_S*IS + gamma_H*IH - omega_R*R;
  dD = delta_S*IS + delta_Q*IQ + delta_H*IH;
  
  return(list(c(dS, dE, dIA, dIS, dQ, dH, dR, dD)))
}

##### Initial conditions
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

#date format:yyyy/m/d
di = as.Date("2021-06-04")
df = as.Date("2021-11-26")
Date1 = seq(di, df, 1)
time = seq(0, length(Date1)-1, 1)

## Load parameters
theta_val <- read.csv("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/theta_values.csv", 
                      header=TRUE, sep=",")

par = c(alpha = 1.0 / 5.2,
        rho = 0.179,
        omega_R = 1.0 / 183,
        m_H = 0.01,
        m_Q = 0.1,
        beta_S = mean(theta_val$beta_S),
        beta_A = mean(theta_val$beta_A),
        gamma_S = mean(theta_val$gamma_S),
        gamma_A = mean(theta_val$gamma_A),
        gamma_H = mean(theta_val$gamma_H),
        gamma_Q= mean(theta_val$gamma_Q),
        delta_S = mean(theta_val$delta_S),
        delta_Q = mean(theta_val$delta_Q),
        delta_H = mean(theta_val$delta_H), 
        sigma_S = mean(theta_val$sigma_S),
        sigma_Q = mean(theta_val$sigma_Q),
        eta = mean(theta_val$eta)
        )

##### Solve ODE system
out = as.matrix(ode(y = y0, times = time, func = SEIQHRD, parms = par, method = "rk4"))
out.df <- as.data.frame(out)

