// Poly-Weibull model for survival analysis - based on Demiris et al (2015), SMMR 24(2), 287-301

functions {
  /*
  /* This is a more complex version of the definition of the log-density, using two external functions to compute
     the log sum (by columns) of the hazards and the sum (by columns) of the log survival. Then the actual sampling
     code can all be vectorised. However, it is only marginally faster than the plain, simpler code that just defines
     all the relevant quantities (including loops) in the sampling code definition.
  /*
    vector sum_hazard(vector t, vector shape, matrix rate, int M) {
      matrix[num_elements(t),M] h;
      vector[num_elements(t)] sh;
      for (j in 1:M) {
        for (i in 1:num_elements(t)) {
          h[i,j] = shape[j]*rate[i,j]*pow(t[i],(shape[j]-1));
        }
      }
      for (i in 1:num_elements(t)) {
        sh[i] = log(sum(h[i]));
      }
      return sh;
    }
    
    vector sum_log_survival(vector t, vector shape, matrix rate, int M) {
      matrix[num_elements(t),M] log_S;
      vector[num_elements(t)] sls;
      for (j in 1:M) {
        for (i in 1:num_elements(t)) {
          log_S[i,j] = rate[i,j]*pow(t[i],shape[j]);
        }
      }
      for (i in 1:num_elements(t)) {
        sls[i] = sum(log_S[i]);
      }
      return sls;
    }
    
    real polyweibull_lpdf(vector t, vector d, vector shape, matrix rate, int M) {
      vector[num_elements(t)] log_h;
      vector[num_elements(t)] log_S;
  	  vector[num_elements(t)] log_lik;
      real prob;
      log_h = sum_hazard(t,shape,rate,M);
      log_S = sum_log_survival(t,shape,rate,M);
      log_lik = d .* log_h - log_S;
      prob = sum(log_lik);
      return prob;
    }
    */
    
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
    int<lower=2> H;                   // number of covariates 
    int<lower=2> M;                   // number of components in the PW model
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
      beta[m,] ~ normal(mu_beta[,m],sigma_beta[,m]);
  }
  t ~ polyweibull(d,alpha,lambda,M);
}
