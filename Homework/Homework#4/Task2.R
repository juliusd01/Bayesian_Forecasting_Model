simulate_mnprobit_data <- function(n, p, beta) {
  # Generate predictors for each latent variable
  x <- array(rnorm(n * 3 * p), dim = c(n, 3, p))
  
  # Generate epsilon (multivariate normal errors)
  epsilon <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  # Compute latent variables z
  z <- matrix(NA, nrow = n, ncol = 3)  # Placeholder for latent variables
  for (i in 1:n) {
    # Compute z_j,i for each observation
    z[i, ] <- x[i, , ] %*% beta + epsilon[i, ]
  }
  
  # Compute observed outcome y_i = argmax(z_1i, z_2i, z_3i)
  y <- apply(z, 1, which.max)
  
  list(x = x, y = y, z = z)
}

# Parameters for simulation
set.seed(123)
n <- 100  # Number of observations
p <- 5    # Number of predictors
beta <- runif(p, -1, 1)  # Random beta coefficients

# Simulate data
sim_data <- simulate_mnprobit_data(n, p, beta)

# Display first few results
str(sim_data$x[1,,])  # First observation's predictors
sim_data$y[1:10]      # First 10 observed outcomes
sim_data$z[1:10, ]    # First 10 latent variable rows