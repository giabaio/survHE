// Gamma survival model
data {
  int<lower=1> n_obs;                     // number of observed cases
  int<lower=0> n_cens;                    // number of censored cases
  vector<lower=0>[n_obs] t;               // fully observed times
  vector<lower=0>[n_cens] d;              // observed censoring times
  int<lower=1> H;                         // number of covariates (including intercept)
  matrix[n_obs,H] X_obs;                  // matrix of categorical covariates for the valid cases (0/1 => dummy variables)
  matrix[n_cens,H] X_cens;                // matrix of categorical covariates for the censored cases (0/1 => dummy variables)
  vector[H] mu_beta;                      // vector of means for the covariates
  vector<lower=0>[H] sigma_beta;          // vector of sd for the covariates
  real<lower=0> a_alpha;		  		        // first parameter for the shape distribution
  real<lower=0> b_alpha;		  		        // second parameter for the shape distribution
}

parameters {
  real<lower=0> alpha;                    // shape of the Gamma distribution
  vector[H] beta;                         // coefficients for the covariates
  vector<lower=1>[n_cens] cens;           // censoring variable (latent)
}

transformed parameters {
  vector[n_obs] loglambda_obs;            // loglinear predictor for the observed cases
  vector[n_cens] loglambda_cens;          // loglinear predictor for the censored cases
  vector[n_obs] lambda_obs;               // rescaled predictor (rate) for the observed cases
  vector[n_cens] lambda_cens;             // rescaled predictor (rate) for the censored cases
  // Only run this if n_cens>0, in which case needs to account for censoring
  if(n_cens>0) {
    loglambda_cens = X_cens*beta + log(d);
    for (i in 1:n_cens) {
      lambda_cens[i] = exp(loglambda_cens[i]);
    }
  }
  loglambda_obs = X_obs*beta;
  for (i in 1:n_obs) {
    lambda_obs[i] = exp(loglambda_obs[i]);
  }
}

model {
  // Prior distributions
  alpha ~ gamma(a_alpha,b_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  // Only does the censoring part if n_cens>0
  if(n_cens>0) {
    cens ~ gamma(alpha,lambda_cens);
  }
  t ~ gamma(alpha,lambda_obs);
}

generated quantities {
  real rate;
  rate = exp(beta[1]);
}
