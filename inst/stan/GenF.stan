// Generalised F model for survival analysis - uses the parameterisation of flexsurv 

// First defines the Generalised F model for fully observed and censored cases
functions {
    real genf_lpdf(vector x, vector mu, real sigma, real Q, real P) {
      // mu = location (\in R)
      // sigma = scale (>0)
      // Q = first shape (\in R)
      // P = second shape (>0)
      vector[num_elements(x)] prob;
      real lprob;
      real tmp;
      real delta;
      real s1; 
      real s2;
      vector[num_elements(x)] expw;
      tmp = pow(Q, 2) + 2*P;
      delta = pow(tmp, 0.5);
      s1 = 2/(tmp + Q*delta);
      s2 = 2/(tmp - Q*delta);
      for (i in 1:num_elements(x)) {
        expw[i] = pow(x[i],delta/sigma)*exp(-mu[i]*delta/sigma);
        prob[i] = log(delta) + s1/sigma * delta * (log(x[i]) - mu[i]) + s1 * (log(s1) - log(s2)) - log(sigma * x[i]) - 
                  (s1 + s2) * log(1 + s1 * expw[i]/s2) - lbeta(s1, s2);
      }
      lprob = sum(prob);
      return lprob;
    }
    
    real genf_cens_lpdf(vector x, vector mu, real sigma, real Q, real P, vector u) {
      vector[num_elements(x)] prob;
      real lprob;
      real tmp;
      real delta;
      real s1; 
      real s2;
      vector[num_elements(x)] expw;
      vector[num_elements(x)] tr;
      tmp = pow(Q,2)+2*P;
      delta = pow(tmp,0.5);
      s1 = 2/(tmp + Q*delta);
      s2 = 2/(tmp - Q*delta);
      tr = x .* u;
      for (i in 1:num_elements(x)) {
        expw[i] = pow(tr[i],delta/sigma)*exp(-mu[i]*delta/sigma);
        prob[i] = log(u[i]) + log(delta) + s1/sigma*delta*(log(tr[i]) - mu[i]) + s1 * (log(s1) - log(s2)) - log(sigma*tr[i]) - 
                  (s1 + s2)*log(1 + s1 * expw[i]/s2) - lbeta(s1, s2);
      }
      lprob=sum(prob);
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
  real<lower=0> a_sigma;		 		          // first parameter for the shape parameter
  real<lower=0> b_sigma;		 		          // second parameter for the shape parameter
  real mu_P;						                  // mean for the parameter logP
  real<lower=0> sigma_P;		  		        // sd for the parameter logP
  real mu_Q; 				  			              // mean for the parameter Q
  real<lower=0> sigma_Q;		 		          // sd for the parameter Q
}

parameters {
  real<lower=0> sigma;                    // scale of the Generalised F distribution
  real Q;                                 // first degree of freedom
  real logP;                              // second degree of freedom (log-scale)
  vector[H] beta;                         // coefficients for the covariates
  vector<lower=1>[n_cens] cens;           // censoring variable (latent)
}

transformed parameters {
  real<lower=0> P;                        // second degree of freedom (natural scale)
  P = exp(logP);
}

model {
  // Prior distributions
  sigma ~ gamma(a_sigma,b_sigma);
  logP ~ normal(mu_P,sigma_P);
  Q ~ normal(mu_Q,sigma_Q);
  beta ~ normal(mu_beta,sigma_beta);
  // Data model
  cens ~ genf_cens(X_cens*beta,sigma,Q,P,d);
  t ~ genf(X_obs*beta,sigma,Q,P);
}

generated quantities {
  real mu;
  mu = beta[1];
}
