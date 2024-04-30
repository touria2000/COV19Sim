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
SEIQHRDV <- function(t, y, par){ 
  ## State variables
  S = y[1]
  E = y[2]
  IA = y[3]
  IS = y[4]
  IQ = y[5]
  IH = y[6]
  R = y[7]
  D = y[8]
  V = y[9]
 
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
  omega_V =  par[18]
  xi_V = par[19]
  theta_V = par[20]
  
  # Forces of infection
  N = S + E + IA + IS + IQ + IH + R + D + V
  lmbda  = (beta_A * IA + beta_S * IS + m_Q*beta_S*IQ + m_H*beta_S*IH) / N
  
  # ODE system
  dS = omega_R*R + omega_V*V - (lmbda + theta_V)*S
  dE = lmbda*S - (1 - xi_V)*lmbda*V - alpha*E
  dIA = rho*alpha*E - gamma_A*IA
  dIS = (1 - rho)*alpha*E - (delta_S + sigma_S + gamma_S + eta)*IS
  dIQ = eta*IS - (delta_Q + sigma_Q + gamma_Q)*IQ
  dIH = sigma_S*IS + sigma_Q*IQ - (delta_H + gamma_H)*IH
  dR = gamma_A*IA + gamma_Q*IQ + gamma_S*IS + gamma_H*IH - omega_R*R
  dD = delta_S*IS + delta_Q*IQ + delta_H*IH
  dV = theta_V*S - (omega_V + (1 - xi_V)*lmbda)*V
  
  return(list(c(dS, dE, dIA, dIS, dQ, dH, dR, dD, dV)))
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
V0 = 25460
S0 = pop - (E0 + IA0 + IS0 + IQ0 + IH0 + R0 + D0 + V0)
y0 = c(S0, E0, IA0, IS0, IQ0, IH0, R0, D0, V0)

#date format:yyyy/m/d
di = as.Date("2021-06-04")
df = as.Date("2021-11-26")
Date1 = seq(di, df, 1)
time = seq(0, length(Date1)-1, 1)

## Load parameters
theta_val <- read.csv("C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_with_vaccination/estimation/theta_values.csv", 
                      header=TRUE, sep=",")

#### Scenario_1: Impact of vaccination rate
# Define a vector of values for theta_V
theta_vec<- c(0.00, 0.006, 0.008, 0.01)

# Define a list to hold the results
res1 <- vector(length(theta_vec),mode="list")

for (k in seq_along(theta_vec)){ #range of values for theta_v
  res1[[k]] <- ode(y = y0, times = time, func = ode1, 
                   parms = c(alpha = 1.0 / 5.2,
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
                             eta = mean(theta_val$eta),
                             omega_V = 1.0 / 183,
                             xi_V = 0.95,
                             theta_V = theta_vec[k]), 
                   method = "rk4")
}

# To get a matrix of the last states from each run
t(sapply(res1,tail,1))

# Make the results in one long data frame
names(res1) <- theta_vec  ## to get theta_v value incorporated in results
df_ode1 <- dplyr::bind_rows(lapply(res1,as.data.frame),.id="theta_V")

# Converting theta_V from a character vector to a numeric vector
df_ode1$theta_V <- as.numeric(df_ode1$theta_V)



#### Scenario_2: Impact of vaccine efficacy
# Define a vector of values for EV_vec
EV_vec<- c(0.00, 0.67, 0.80, 0.95)

# Define a list to hold the results
res2 <- vector(length(EV_vec),mode="list")

for (k in seq_along(EV_vec)){ #range of values for theta_v
  res2[[k]] <- ode(y = y0, times = time, func = ode1, 
                   parms = c(alpha = 1.0 / 5.2,
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
                             eta = mean(theta_val$eta),
                             omega_V = 1.0 / 183,
                             xi_V = EV_vec[k],
                             theta_V = 0.0051), 
                   method = "rk4")
}

# To get a matrix of the last states from each run
t(sapply(res2,tail,1))

# Make the results in one long data frame
names(res2) <- theta_vec  ## to get theta_v value incorporated in results
df_ode2 <- dplyr::bind_rows(lapply(res2,as.data.frame),.id="E_V")

# Converting theta_V from a character vector to a numeric vector
df_ode2$E_V <- as.numeric(df_ode2$E_V)

