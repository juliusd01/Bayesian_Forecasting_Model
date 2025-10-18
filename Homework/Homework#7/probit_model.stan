//
// trying to fit a probit model
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;  // Number of observations
  int<lower=0,upper=1> y[N];  // Binary outcome
  int<lower=1> K;  // Number of predictors
  matrix[N, K] X;  // Predictor matrix
}

parameters {
  vector[K] beta;  // Regression coefficients
}

model {
  beta ~ normal(0, 5);  // Weakly informative priors
  y ~ bernoulli(Phi(X * beta));  // Probit likelihood using normal CDF
}