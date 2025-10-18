rm(list=ls())
library(bayesm)
source('varselv3.R')

#set.seed(66)

hilfseed=.Random.seed

simbprobit=function(X,beta) {
##  function to simulate from binary probit including x variable
z=(X%*%beta+rnorm(nrow(X)))
y=ifelse(z<0,0,1)
list(X=X,y=y,z=z,beta=beta)
}

nobs=100
X=cbind(rep(1,nobs),runif(nobs),runif(nobs))
beta=c(-3,3,6)
nvar=ncol(X)
simout=simbprobit(X,beta)

y=simout$y

R=200000 ## 
# different number of 'phantom' explanatory vars, to check whether the correct ones can be chosen
tk=10

###standard analysis

Xt=cbind(simout$X,matrix(runif(nobs*tk),ncol=tk))
Data1=list(X=Xt,y=simout$y)
Mcmc1=list(R=R,keep=5)


#standard weakly infomative prior
A1= .01  # #.01 #.05 #.1 #1
Prior1=list(A=A1*diag(nvar+tk),mu=rep(0,nvar+tk),aBeta=2,bBeta=2)#aBeta=1000,bBeta=1)#
outbinom2=rbprobitGibbs(Data=Data1,Mcmc=Mcmc1,Prior=Prior1)
summary(outbinom2$betadraw)
#plot(outbinom$betadraw)


####model with variable selection
####model with variable selection
####model with variable selection
####model with variable selection
A1= .01 #.05 #.1 #1
# and different startings for tau
#Init1=list(tau=c(1,rep(1,nvar+tk-1)),z=double(length(simout$z)))  # start with everything in
Init1=list(tau=c(1,rep(0,nvar+tk-1)),z=double(length(simout$z)))  # start with everything out
Prior1=list(A=A1*diag(nvar+tk),mu=rep(0,nvar+tk),aBeta=2,bBeta=2)#aBeta=1000,bBeta=1)#
out=probvs(Data=Data1,Mcmc=Mcmc1,Prior=Prior1,Init=Init1)

summary(out$betadraw)

windows()
matplot(out$betadraw[,1:3],type="l",ylim=c(-20,25))
abline(beta[1],0,lwd=2, lty=2, col='grey')
abline(beta[2],0,lwd=2, lty=2, col='dark red')
abline(beta[3],0,lwd=2, lty=2, col='dark green')


matplot(outbinom2$betadraw[,1:3],type="l",ylim=c(-20,25))
abline(beta[1],0,lwd=2, lty=2, col='grey')
abline(beta[2],0,lwd=2, lty=2, col='dark red')
abline(beta[3],0,lwd=2, lty=2, col='dark green')



#--------------------------------
# SECOND TASK
#i) Netwon-Raftery
#  Define the log-likelihood function for the probit model
probit_log_likelihood <- function(beta, y, X) {
  linear_predictor <- X %*% beta  # X is the design matrix, beta is the parameter vector
  prob <- pnorm(linear_predictor)  # Probit link function
  log_lik <- sum(dbinom(y, size = 1, prob = prob, log = TRUE))  # Bernoulli log-likelihood
  return(log_lik)
}

# Compute log-likelihoods for each posterior sample
log_likelihoods <- apply(outbinom2$betadraw, 1, function(beta) {
  probit_log_likelihood(beta, y = Data1$y, X = Data1$X)
})

# Compute the log marginal density using logMargDenNR
log_marginal_density <- logMargDenNR(log_likelihoods)




# ii) Gelfand-Dey
# Define IndMWnvec first
IndMWnvec <- function(x, mu, rooti) {
  z <- t(rooti) %*% (x - mu)
  log_density <- -(nrow(x)/2) * log(2 * pi) - 0.5 * colSums(z * z) + sum(log(diag(rooti)))
  return(log_density)
}

# Compute posterior mean and covariance
posterior_beta <- outbinom2$betadraw  # Ensure this is a matrix (iterations x parameters)
posterior_mean <- colMeans(posterior_beta)
posterior_cov <- cov(posterior_beta)

# Compute inverse Cholesky factor of posterior covariance
chol_posterior <- chol(posterior_cov)
rooti_posterior <- solve(chol_posterior)  # Inverse of Cholesky

# Log-density of posterior samples under q(theta) ~ MVN(posterior_mean, posterior_cov)
log_q <- IndMWnvec(t(posterior_beta), mu = posterior_mean, rooti = rooti_posterior)

# Log-prior density (Prior: beta ~ N(0, solve(A)))
prior_A <- Prior1$A  # Prior precision matrix (e.g., 0.01 * diag(nvar + tk))
prior_cov <- solve(prior_A)  # Prior covariance
chol_prior <- chol(prior_cov)
rooti_prior <- solve(chol_prior)
log_prior <- IndMWnvec(t(posterior_beta), mu = Prior1$mu, rooti = rooti_prior)

# Compute lqp' = log(q) - (log_likelihood + log_prior)
lqp_prime <- log_q - (log_likelihoods + log_prior)
max_lqp <- max(lqp_prime)

# Final log marginal density
log_marginal_GD <- -max_lqp - log(mean(exp(lqp_prime - max_lqp)))




# iii) Importance Sampling
# Generate samples from q(theta)
num_samples <- nrow(posterior_beta)
importance_samples <- MASS::mvrnorm(n = num_samples, mu = posterior_mean, Sigma = posterior_cov)

# Compute log-likelihood, log-prior, and log-q for importance samples
log_lik_importance <- apply(importance_samples, 1, function(beta) {
  probit_log_likelihood(beta, y = Data1$y, X = Data1$X)
})

log_q_importance <- IndMWnvec(t(importance_samples), mu = posterior_mean, rooti = rooti_posterior)
log_prior_importance <- IndMWnvec(t(importance_samples), mu = prior_mean, rooti = rooti_prior)

# Compute weights and log marginal density
weights <- exp(log_lik_importance + log_prior_importance - log_q_importance)
log_marginal_IS <- log(mean(weights))



# iv) Bridge Sampling
library(rstan)
library(mvtnorm)
library(bridgesampling)

# Stan model (save as "probit_model.stan")
model_code <- "
data {
  int<lower=1> N;
  int<lower=0,upper=1> y[N];
  int<lower=1> K;
  matrix[N, K] X;
}
parameters {
  vector[K] beta;
}
model {
  beta ~ normal(0, sqrt(1/0.01));  // Prior: beta ~ N(0, 100)
  y ~ bernoulli(Phi(X * beta));
}
"

# Fit model in Stan
stan_data <- list(
  N = length(Data1$y),
  y = as.vector(Data1$y),
  K = ncol(Data1$X),
  X = Data1$X
)

stan_fit <- stan(model_code = model_code, data = stan_data, iter = 2000, chains = 4)
bridge_result <- bridge_sampler(stan_fit)
log_marginal_bridge <- bridge_result$logml




