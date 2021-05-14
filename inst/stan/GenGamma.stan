// Generalised Gamma model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised Gamma model for fully observed and censored cases

functions {
  real gen_gamma_lpdf(vector x, vector mu, real sigma, real Q) {
    // Uses the same parameterisation as flexsurv
    // mu = location
    // sigma = scale
    // Q = shape
    vector[num_elements(x)] prob;
    real lprob;
    vector[num_elements(x)] w;
    // Constructs the log-density for each observation
    w = ((log(x)-mu))/sigma;
    for (i in 1:num_elements(x)) {
      prob[i] = -log(sigma*x[i])+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w[i]-exp(Q*w[i]))-lgamma(pow(Q,-2));
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = sum((prob));
    return lprob;
  }
  
  real gen_gamma_cens_lpdf(vector x, vector mu, real sigma, real Q, vector u) {
    // Rescales the distribution accounting for right censoring
    vector[num_elements(x)] prob;
    real lprob;
    vector[num_elements(x)] w;
    vector[num_elements(x)] tr;
    // Constructs the log-density for each observation
    tr = x .* u;
    w = ((log(tr)-mu))/sigma;
    for (i in 1:num_elements(x)) {
      prob[i] = log(u[i])-log(sigma*tr[i])+log(fabs(Q))+pow(Q,-2)*log(pow(Q,-2))+pow(Q,-2)*(Q*w[i]-exp(Q*w[i]))-lgamma(pow(Q,-2));
    }
    // And the total log-density (as a sum of the individual terms)
    lprob = sum((prob));
    return lprob;
  }
}

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
  real mu_Q;                              // mean for the parameter Q
  real<lower=0> sigma_Q;                  // sd for the parameter Q
  real<lower=0> a_sigma;                  // first parameter for the scale distribution
  real<lower=0> b_sigma;                  // second parameter for the scale distribution
}

parameters {
  real Q;                                 // shape of the Generalised Gamma distribution
  real<lower=0> sigma;                    // scale of the Generalised Gamma distribution
  vector[H] beta;                         // coefficients for the covariates
  vector<lower=1>[n_cens] cens;           // censoring variable (latent)
}

model {
  // Prior distributions
  Q ~ normal(mu_Q,sigma_Q);
  sigma ~ gamma(a_sigma,b_sigma);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  // Only does the censoring part if n_cens>0
  if(n_cens>0) {
     cens ~ gen_gamma_cens(X_cens*beta,sigma,Q,d);
  }
  t ~ gen_gamma(X_obs*beta,sigma,Q);
}

generated quantities {
  real mu;
  mu = beta[1];
}
