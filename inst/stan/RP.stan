// Royston-Parmar splines model 

functions {
  real rps_lpdf(vector t, vector d, vector gamma, matrix B, matrix DB, vector linpred) {
    // t = vector of observed times
    // d = event indicator (=1 if event happened and 0 if censored)
    // gamma = M+2 vector of coefficients for the flexible part
    // B = matrix of basis
    // DB = matrix of derivatives for the basis
    // linpred = fixed effect part
    vector[num_elements(t)] eta;
    vector[num_elements(t)] eta_prime;
    vector[num_elements(t)] log_lik;
    real lprob;
    
    eta = B*gamma + linpred;
    eta_prime = DB*gamma;
    log_lik = d .* (-log(t) + log(eta_prime) + eta) - exp(eta);
    lprob = sum(log_lik);
    return lprob;
  }
}

data {
  int<lower=1> n;                   // number of observations
  int<lower=0> M;                   // number of internal knots for the splines model
  int<lower=1> H;                   // number of covariates in the (time-independent) linear predictor
  vector<lower=0>[n] t;             // observed times (including censored values)
  vector<lower=0,upper=1>[n] d;     // censoring indicator: 1 if fully observed, 0 if censored
  matrix[n,H] X;                    // matrix of covariates for the (time-independent) linear predictor
  matrix[n,M+2] B;                  // matrix with basis
  matrix[n,M+2] DB;                 // matrix with derivatives of the basis
  vector[H] mu_beta;                // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  vector[M+2] mu_gamma;             // mean of the splines coefficients
  vector<lower=0>[M+2] sigma_gamma; // sd of the splines coefficients
}

parameters {
  vector[M+2] gamma;
  vector[H] beta;
}

model {
  // Priors
  gamma ~ normal(mu_gamma,sigma_gamma);
  beta ~ normal(mu_beta,sigma_beta);
  
  // Data model
  t ~ rps(d,gamma,B,DB,X*beta);
}
