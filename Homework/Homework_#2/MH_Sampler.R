rm(list=ls(all=TRUE))

# Load necessary libraries
library(stats)
library(coda)  # For autocorrelation and diagnostics

# Parameters for the beta prior (alpha, beta) and data
alpha_prior <- 2   # Prior parameter alpha
beta_prior <- 5    # Prior parameter beta
N <- 400            # Number of binomial trials
p_true <- 0.1      # True probability of success for simulation
y <- rbinom(1, N, p_true)  # Simulated observed data: number of successes

# Function for log-posterior (up to constant)
log_posterior <- function(p, y, N, alpha_prior, beta_prior) {
  if (p < 0 || p > 1) return(-Inf)  # Reject invalid p
  log_likelihood <- dbinom(y, N, p, log = TRUE)
  log_prior <- dbeta(p, alpha_prior, beta_prior, log = TRUE)
  return(log_likelihood + log_prior)  # Posterior up to a constant
}

# Random walk MH sampler function with console output
mh_sampler <- function(y, N, alpha_prior, beta_prior, n_iter, init_p, step_size) {
  # Initialize
  samples <- numeric(n_iter)
  samples[1] <- init_p
  current_log_post <- log_posterior(init_p, y, N, alpha_prior, beta_prior)
  accepted <- 0  # Track acceptance
  
  cat("Starting Metropolis-Hastings sampling...\n")
  cat("=========================================\n")
  cat(sprintf("Initial probability p: %.4f\n", init_p))
  cat(sprintf("Total iterations: %d\n", n_iter))
  cat("Step size:", step_size, "\n\n")
  
  for (i in 2:n_iter) {
    # Propose a new value
    proposed_p <- samples[i - 1] + rnorm(1, mean = 0, sd = step_size)
    
    # Calculate the posterior for the proposed value
    proposed_log_post <- log_posterior(proposed_p, y, N, alpha_prior, beta_prior)
    
    # Acceptance ratio --> Adding and subtracting the normalizing
    # constant obviously gives the same result
    log_accept_ratio <- proposed_log_post - current_log_post
    if (log(runif(1)) < log_accept_ratio) {
      samples[i] <- proposed_p
      current_log_post <- proposed_log_post
      accepted <- accepted + 1  # Increment acceptance count
    } else {
      samples[i] <- samples[i - 1]
    }
    
    # Progress updates at intervals
    if (i %% (n_iter / 10) == 0) {
      cat(sprintf("Iteration %d: Acceptance rate = %.2f%%, Running mean of p = %.4f\n",
                  i, (accepted / i) * 100, mean(samples[1:i])))
    }
  }
  
  cat("\nSampling complete!\n")
  cat(sprintf("Final acceptance rate: %.2f%%\n", (accepted / n_iter) * 100))
  
  return(samples)
}

# Run the sampler
set.seed(123)
n_iter <- 10000
init_p <- 0.5  # Initial value for p
step_size <- 0.05  # Experiment with different step sizes
burnin <- n_iter*0.2

samples <- mh_sampler(y, N, alpha_prior, beta_prior, n_iter, init_p, step_size)

# Diagnostic Plots
par(mfrow = c(2, 2))

plot(samples, type = 'l', main = 'Trace Plot', xlab = 'Iteration', ylab = 'p')
abline(v = burnin, col = 'red', lty = 2)  # Add vertical line for burn-in
legend("topright", legend = "Burn-in End", col = "red", lty = 2)

# Histogram of the sampled p values (Posterior Distribution)
hist(samples, breaks = 30, probability = TRUE, main = 'Posterior Distribution', xlab = 'p')
lines(density(samples), col = 'blue', lwd = 2)

# Autocorrelation plot
acf(samples, main = "Autocorrelation of Samples")

# Effective sample size calculation
ess <- effectiveSize(samples)
# Add effective sample size to the plot
text(10, 0.8 * max(acf(samples)$acf), paste("ESS:", round(ess)), col = 'red')

# Generate summary statistics
posterior_summary <- function(samples) {
  # Moments
  mean_p <- mean(samples)
  sd_p <- sd(samples)
  n_samples <- length(samples)
  se <- sd_p / sqrt(n_samples)
  ess <- effectiveSize(samples)  # Effective sample size
  
  # Relative efficiency
  rel_eff <- ess / n_samples
  
  # Quantiles
  quantiles <- quantile(samples, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  
  # Create a summary data frame for moments
  moments_summary <- data.frame(
    Metric = c("Mean", "Standard Deviation", "Standard Error", "Effective Sample Size", "Relative Efficiency"),
    Value = c(mean_p, sd_p, se, ess, rel_eff)
  )
  
  # Create a summary data frame for quantiles
  quantiles_summary <- data.frame(
    Quantile = c("2.5%", "5%", "50%", "95%", "97.5%"),
    Value = quantiles
  )
  
  # Output moments summary table
  cat("Summary of Posterior Marginal Distributions (Moments)\n")
  print(kable(moments_summary, format = "pipe", caption = "Moments Summary"))
  
  # Output quantiles summary table
  cat("\nQuantiles Summary\n")
  print(kable(quantiles_summary, format = "pipe", caption = "Quantiles Summary"))
  
}

posterior_summary(samples = samples[(burnin + 1):n_iter])