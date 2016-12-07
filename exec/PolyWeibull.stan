// Poly-Weibull model for survival analysis - based on Demiris et al (2015), SMMR 24(2), 287-301

functions {
    real polyweibull_lpdf(vector t, vector d, vector shape, matrix rate, int M) {
      matrix[num_elements(t),M] h;
      matrix[num_elements(t),M] log_S;
  	  vector[num_elements(t)] log_lik;
      real prob;
      for (j in 1:M) {
          for (i in 1:num_elements(t)) {
              h[i,j] = shape[j]*rate[i,j]*pow(t[i],(shape[j]-1));
              log_S[i,j] = rate[i,j]*pow(t[i],shape[j]);
          }
      }
      for (i in 1:num_elements(t)) {
        log_lik[i] = d[i]*log(sum(h[i,])) - sum(log_S[i,]);
      }
      prob = sum(log_lik);
      return prob;
    }
}

// Poly-Weibull survival model
data {
    int n;                            // number of observations
    vector[n] t;                      // observed times
    vector[n] d;                      // censoring indicator (1=observed, 0=censored)
    int H;                            // number of covariates 
    int M;                            // number of components in the PW model
    matrix[n,H] X[M];                 // array of covariates matrix of covariates (each with n rows and H columns, the max number of covariates)
    matrix[H,M] mu_beta;              // mean of the covariates coefficients
    matrix<lower=0>[H,M] sigma_beta;  // sd of the covariates coefficients
}

parameters {
    vector[H] beta[M];              // array of vectors of the coefficients for the linear predictors
    positive_ordered[M] alpha;      // ordered vector of shapes
}

transformed parameters {
    matrix[n,M] loglambda;
    matrix[n,M] lambda;
    for (m in 1:M) {
        loglambda[,m] = X[m,,]*beta[m,];
        for (i in 1:n) {
            lambda[i,m] = exp(loglambda[i,m]);
        }
    }
}

model {
  // Prior distributions
  for (m in 1:M) {
      beta[m,] ~ normal(mu_beta[m,],sigma_beta[m,]);
  }
  t ~ polyweibull(d,alpha,lambda,M);
}
