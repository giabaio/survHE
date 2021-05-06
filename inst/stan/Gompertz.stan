// Gompertz survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, real shape, vector rate) {
    vector[num_elements(t)] log_h;
    log_h = log(rate) + (shape * t);
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, real shape, vector rate) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -rate[i]/shape * (exp(shape * t[i]) - 1);
    }
    return log_S;
  }
  
  // Defines the sampling distribution
  real surv_gompertz_lpdf (vector t, vector d, real shape, vector rate) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,shape,rate) + log_S(t,shape,rate);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                          // number of observations
  vector[n] t;                    // observed times
  vector[n] d;                    // censoring indicator (1=observed, 0=censored)
  int H;                          // number of covariates
  matrix[n,H] X;                  // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	            // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta; // sd of the covariates coefficients
  real<lower=0> a_alpha;          // mu_alpha
  real<lower=0> b_alpha;          // sigma_alpha
}

parameters {
  vector[H] beta;                 // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;            // shape parameter
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  linpred = X*beta;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }
}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_gompertz(d,alpha,mu);
}

generated quantities {
  real rate;                        // rate parameter
  rate = exp(beta[1]);
}
