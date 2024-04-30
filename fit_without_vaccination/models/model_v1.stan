functions {
  real[] SEIQHRD(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
      // State variables
      real S = y[1];
      real E = y[2];
      real IA = y[3];
      real IS = y[4];
      real IQ = y[5];
      real IH = y[6];
      real R = y[7];
      real D = y[8];

      // Fixed parameters
      real alpha = 1.0 / 5.2;
      real rho = 0.179;
      real omega_R = 1.0 / 183;
      real m_H = 0.01;
      real m_Q =0.1;
      
      // Parameters to monitor
      real beta_S = theta[1];
      real beta_A = theta[2];
      real gamma_S = theta[3];
      real gamma_A = theta[4];
      real gamma_H = theta[5];
      real gamma_Q = theta[6];
      real delta_S = theta[7];
      real delta_Q = theta[8];
      real delta_H = theta[9];
      real sigma_S = theta[10];
      real sigma_Q = theta[11];
      real eta = theta[12];

      // Force of infection
      real N = S + E + IA + IS + IQ + IH + R + D; 
      real lmbda  = (beta_A * IA + beta_S * IS + m_Q*beta_S*IQ + m_H*beta_S*IH) / N; 
      
      // ODE system
      real dS_dt = omega_R*R - lmbda*S;
      real dE_dt = lmbda*S - alpha*E;
      real dIA_dt = rho*alpha*E - gamma_A*IA;
      real dIS_dt = (1 - rho)*alpha*E - (delta_S + sigma_S + gamma_S + eta)*IS;
      real dIQ_dt = eta*IS - (delta_Q + sigma_Q + gamma_Q)*IQ;
      real dIH_dt = sigma_S*IS + sigma_Q*IQ - (delta_H + gamma_H)*IH;
      real dR_dt = gamma_A*IA + gamma_Q*IQ + gamma_S*IS + gamma_H*IH - omega_R*R;
      real dD_dt = delta_S*IS + delta_Q*IQ + delta_H*IH;
      
      return {dS_dt, dE_dt, dIA_dt, dIS_dt, dIQ_dt, dIH_dt, dR_dt, dD_dt};
  }
}

data {
  int<lower = 1> n_days;    // number of observed days
  int<lower = 1> n_difeq;   // number of differential equations
  real y_init[n_difeq];     // initial conditions for state variables
  real t0;                  // initial time point (zero)
  real ts[n_days];          // observed time points 
  int cases[n_days];        // data, daily number of infected individuals
  //int<lower = 1> pop;       // population size
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> beta_S;
  real<lower=0> beta_A;
  real<lower=0> gamma_S;
  real<lower=0> gamma_A;  
  real<lower=0> gamma_H;
  real<lower=0> gamma_Q;
  real<lower=0.0, upper=0.004> delta_S;
  real<lower=0.0, upper=0.004> delta_Q;
  real<lower=0.0, upper=0.04> delta_H;
  real<lower=1.0 / 15, upper=1.0 / 7> sigma_S;
  real<lower=1.0 / 8, upper=1.0 / 6> sigma_Q;
  real<lower=1.0 / 10, upper=1.0 / 4> eta;
  
  real<lower=0> phi_inv;
}

transformed parameters{
  real y_hat[n_days, n_difeq];  // solution from the ODE solver
  real incidence[n_days];
  real phi = 1.0 / phi_inv;
  real alpha = 1.0 / 5.2;
  real rho = 0.179;

  {
    real theta[12];
    theta[1] = beta_S;
    theta[2] = beta_A;
    theta[3] = gamma_S;
    theta[4] = gamma_A;
    theta[5] = gamma_H;
    theta[6] = gamma_Q;
    theta[7] = delta_S;
    theta[8] = delta_Q;
    theta[9] = delta_H;
    theta[10] = sigma_S;
    theta[11] = sigma_Q;
    theta[12] = eta;

    y_hat = integrate_ode_rk45(SEIQHRD, y_init, t0, ts, theta, x_r, x_i);
  }
  
  // Calculate reported incidence from the model for the time interval we want to fit to
  for (i in 1:n_days){
    incidence[i] = (1-rho)*alpha* y_hat[i,2];
  } 
}  

model {
  // Priors
  beta_S ~ normal(1.5, 0.2);
  beta_A ~ normal(1.5, 0.2);
  gamma_S ~ normal(0.14 , 0.1);   
  gamma_A ~ normal(0.14 , 0.1);   
  gamma_H ~ normal(0.071 , 0.01); 
  gamma_Q ~ normal(0.071 , 0.01); 
  delta_S ~ uniform(0.0, 0.004);
  delta_Q ~ uniform(0.0, 0.004);
  delta_H ~ uniform(0.0, 0.04);
  sigma_S ~ uniform(1.0 / 15, 1.0 / 7);
  sigma_Q ~ uniform(1.0 / 8 , 1.0 / 6);
  eta ~ uniform(1.0 / 10, 1.0 / 4);

  phi_inv ~ exponential(5);
  
  // Observation model
  cases[1:(n_days)] ~ neg_binomial_2(incidence, phi);
}

generated quantities {
  // recovery rates
  real recovery_time_S = 1 / gamma_S;
  real recovery_time_A = 1 / gamma_A; 
  real recovery_time_H = 1 / gamma_H;
  real recovery_time_Q = 1 / gamma_Q;
  
  // Basic reproduction number
  real R_0;     
  real m_H = 0.01;
  real m_Q =0.1;
  real psi_1 = delta_S + sigma_S + gamma_S + eta;
  real psi_2 = delta_Q + sigma_Q + gamma_Q;
  R_0 = (rho*beta_A)/gamma_A + (((1-rho)*beta_S)/psi_1)*(1+(eta*m_Q)/psi_2 + (sigma_S + (eta*sigma_Q)/psi_2)*(m_H/(delta_H + gamma_H)));

  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(incidence, phi);
}



