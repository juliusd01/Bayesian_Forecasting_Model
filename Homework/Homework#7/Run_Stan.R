rm(list=ls())
# Load required libraries
library(rstan)

# Set seed for reproducibility
set.seed(66)

# Data Generating Process (Probit Model)
simbprobit <- function(X, beta) {
  z <- (X %*% beta + rnorm(nrow(X)))  # Latent variable
  y <- ifelse(z < 0, 0, 1)  # Binary outcome
  return(list(X = X, y = y, z = z, beta = beta))
}

# Generate Data
nobs <- 100  # Small dataset
X <- cbind(rep(1, nobs), runif(nobs), runif(nobs))  # Design matrix with intercept
beta <- c(-3, 3, 6)  # True coefficients
simout <- simbprobit(X, beta)

y <- as.vector(simout$y)

# Prepare data for Stan
Data_Stan <- list(
  N = nobs,  
  K = ncol(X),  
  X = X,  
  y = y
)

# Run Stan Model
fit <- stan(file = "probit_model.stan", data = Data_Stan, iter = 20000, chains = 4)

# Print results
print(fit, pars = "beta")

# Extract posterior samples
posterior_samples <- extract(fit)$beta


windows()
matplot(posterior_samples[,1:3],type="l",ylim=c(-20,25))
abline(beta[1],0,lwd=2, lty=2, col='grey')
abline(beta[2],0,lwd=2, lty=2, col='dark red')
abline(beta[3],0,lwd=2, lty=2, col='dark green')
