// log-Normal survival model

functions {
  // Defines the log survival
  vector log_S (vector t, vector mean, real sd) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = log(1-Phi((log(t[i])-mean[i])/sd));
    }
    return log_S;
  }
  
  // Defines the log hazard
  vector log_h (vector t, vector mean, real sd) {
    vector[num_elements(t)] log_h;
    vector[num_elements(t)] ls;
    ls = log_S(t,mean,sd);
    for (i in 1:num_elements(t)) {
      log_h[i] = lognormal_lpdf(t[i]|mean[i],sd) - ls[i];
    }
    return log_h;
  }

  // Defines the sampling distribution
  real surv_lognormal_lpdf (vector t, vector d, vector mean, real sd) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,mean,sd) + log_S(t,mean,sd);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  int H;                  // number of covariates
  matrix[n,H] X;          // matrix of covariates (with n rows and H columns)
  vector[H] mu_beta;	    // mean of the covariates coefficients
  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients
  real a_alpha;			      // lower bound for the sd of the data			  
  real b_alpha;			      // upper bound for the sd of the data
}

parameters {
  vector[H] beta;         // Coefficients in the linear predictor (including intercept)
  real<lower=0> alpha;    // log-sd parameter
}

model {
  alpha ~ uniform(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  t ~ surv_lognormal(d,X*beta,alpha);
}

generated quantities {
  real meanlog;
  meanlog = beta[1];
}
