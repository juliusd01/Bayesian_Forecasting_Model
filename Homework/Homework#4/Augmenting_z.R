# We need to install truncnorm
# install.packages("truncnorm")

simulate_mnprobit_data_augmented <- function(n, p, beta) {
  # Generate predictors for each latent variable
  x <- array(rnorm(n * 3 * p), dim = c(n, 3, p))
  
  # Generate epsilon (multivariate normal errors with identity covariance)
  epsilon <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  # Compute latent variables z
  z <- matrix(NA, nrow = n, ncol = 3)  # Placeholder for latent variables
  for (i in 1:n) {
    z[i, ] <- x[i, , ] %*% beta + epsilon[i, ]
  }
  
  # Compute observed outcome y_i = argmax(z_1i, z_2i, z_3i)
  y <- apply(z, 1, which.max)
  
  # Data augmentation: Resample z's based on the observed y
  z_augmented <- matrix(NA, nrow = n, ncol = 3)  # Placeholder for augmented z's
  
  for (i in 1:n) {
    j1 <- y[i]  # Chosen alternative
    
    # Resample z for the chosen alternative j1 (unrestricted)
    z_augmented[i, j1] <- rnorm(1, mean = x[i, j1, ] %*% beta, sd = 1)
    
    # Resample z for non-chosen alternatives (truncated normal)
    for (k in setdiff(1:3, j1)) {
      mean_k <- x[i, k, ] %*% beta
      z_augmented[i, k] <- truncnorm::rtruncnorm(1, a = -Inf, b = z_augmented[i, j1], mean = mean_k, sd = 1)
    }
  }
  
  list(x = x, y = y, z = z, z_augmented = z_augmented)
}

# Parameters for simulation
set.seed(123)
n <- 100  # Number of observations
p <- 5    # Number of predictors
beta <- runif(p, -1, 1)  # Random beta coefficients

# Simulate data with augmentation
sim_data <- simulate_mnprobit_data_augmented(n, p, beta)

# Display results
cat("First 10 observed outcomes:\n")
print(sim_data$y[1:10])
cat("\nLatent variables (first 10 rows):\n")
print(sim_data$z[1:10, ])
cat("\nAugmented latent variables (first 10 rows):\n")
print(sim_data$z_augmented[1:10, ])
