
summary_statistics <- function(fit_model){
  parameters = c("beta_S", "beta_A", "gamma_S", "gamma_A", "gamma_H", "gamma_Q", 
                 "delta_S", "delta_Q", "delta_H","sigma_S", "sigma_Q","eta", "R_0", 
                 "phi_inv", "phi", "R_0","lp__",
                 "E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0")
  fit_summary <- summary(fit_model, pars = parameters)$summary
  df_pars <- data.frame(fit_summary)
  
  names(df_pars)[4] <-"2.5%"
  names(df_pars)[5] <-"25%"
  names(df_pars)[6] <-"50%"
  names(df_pars)[7] <-"75%"
  names(df_pars)[8] <-"97.5%"
  
  colnames(df_pars) <- c('mean', 'se_mean','sd',
                         '2.5%', '25%','50%',
                         '75%' ,'97.5%', 'n_eff', 'Rhat')
  
  row.names(df_pars) <- c("beta_S", "beta_A", "gamma_S", "gamma_A", "gamma_H", "gamma_Q", 
                          "delta_S", "delta_Q", "delta_H","sigma_S", "sigma_Q","eta", "R_0", 
                          "phi_inv", "phi", "R_0","lp__",
                          "E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0")
  
  write.csv(df_pars, "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/pars_summary.csv")
}

df_theta <- function(posts_fit){
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
  
  write.table(df_theta, file = "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/theta_values.csv", 
              sep = ",", row.names = FALSE)
}

Mat_IC <- function(posts_fit){
  Mat_IC <- posts_fit$y_init
  write.table(Mat_IC, file = "C:/Users/Touria/Documents/R_Project/NovelCovid19_VaccinationModel/fit_without_vaccination/IC.csv", 
              sep = ",", 
              col.names = c("S0", "E0", "IA0", "IS0", "IQ0", "IH0", "R0", "D0"),
              row.names = FALSE)
}

df_fit <- function(posts_fit){
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
                       R_0 = posts_fit$R_0,
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
  
  write.table(df_fit, file = "C:/Users/Touria/Documents/NovelCovid19_VaccinationModel/fit_without_vaccination/fit_estimates.csv", 
              sep = ",", row.names = FALSE)
}


