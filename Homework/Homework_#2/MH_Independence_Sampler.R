# Clear the global environment
rm(list = ls(all = TRUE))

# Load necessary libraries
library(stats)
library(coda)  # For autocorrelation and diagnostics
library(knitr)

# Set the random seed for reproducibility
set.seed(12)

# Parameters for the beta prior (alpha, beta) and data
alpha_prior <- 2   # Prior parameter alpha
beta_prior <- 5    # Prior parameter beta
N <- 50            # Number of binomial trials
p_true <- 0.1      # True probability of success for simulation
y <- rbinom(1, N, p_true)  # Simulated observed data: number of successes

# Function for log-posterior (up to constant)
log_posterior <- function(p, y, N, alpha_prior, beta_prior) {
  if (p < 0 || p > 1) return(-Inf)  # Reject invalid p
  log_likelihood <- dbinom(y, N, p, log = TRUE)
  log_prior <- dbeta(p, alpha_prior, beta_prior, log = TRUE)
  return(log_likelihood + log_prior)  # Posterior up to a constant
}

# Independence Metropolis-Hastings sampler function
mh_sampler_independence <- function(y, N, alpha_prior, beta_prior, n_iter, init_p, proposal_dist, proposal_params) {
  # Initialize
  samples <- numeric(n_iter)
  samples[1] <- init_p
  current_log_post <- log_posterior(init_p, y, N, alpha_prior, beta_prior)
  accepted <- 0
  
  for (i in 2:n_iter) {
    # Generate proposal based on specified distribution
    proposed_p <- switch(proposal_dist,
                         "uniform" = runif(1, proposal_params[1], proposal_params[2]),
                         "normal" = rnorm(1, mean = proposal_params[1], sd = proposal_params[2]),
                         "beta" = rbeta(1, proposal_params[1], proposal_params[2]))
    
    # Calculate the posterior for the proposed value
    proposed_log_post <- log_posterior(proposed_p, y, N, alpha_prior, beta_prior)
    
    # Calculate proposal density for MH ratio
    if (proposal_dist == "uniform") {
      log_q_current <- dunif(samples[i - 1], proposal_params[1], proposal_params[2], log = TRUE)
      log_q_proposed <- dunif(proposed_p, proposal_params[1], proposal_params[2], log = TRUE)
    } else if (proposal_dist == "normal") {
      log_q_current <- dnorm(samples[i - 1], mean = proposal_params[1], sd = proposal_params[2], log = TRUE)
      log_q_proposed <- dnorm(proposed_p, mean = proposal_params[1], sd = proposal_params[2], log = TRUE)
    } else if (proposal_dist == "beta") {
      log_q_current <- dbeta(samples[i - 1], proposal_params[1], proposal_params[2], log = TRUE)
      log_q_proposed <- dbeta(proposed_p, proposal_params[1], proposal_params[2], log = TRUE)
    }
    
    # Acceptance ratio for independence sampler
    log_accept_ratio <- (proposed_log_post + log_q_current) - (current_log_post + log_q_proposed)
    
    # Accept/reject step
    if (log(runif(1)) < log_accept_ratio) {
      samples[i] <- proposed_p
      current_log_post <- proposed_log_post
      accepted <- accepted + 1
    } else {
      samples[i] <- samples[i - 1]
    }
  }
  
  cat(sprintf("Acceptance rate for %s proposal: %.2f%%\n", proposal_dist, (accepted / n_iter) * 100))
  return(samples)
}

# Number of iterations and burn-in period
n_iter <- 10000
init_p <- 0.5
burnin <- n_iter * 0.2

# Run the sampler with each proposal distribution
samples_uniform <- mh_sampler_independence(y, N, alpha_prior, beta_prior, n_iter, init_p, "uniform", c(0, 1))
samples_normal <- mh_sampler_independence(y, N, alpha_prior, beta_prior, n_iter, init_p, "normal", c(0.4, 0.1))
samples_beta <- mh_sampler_independence(y, N, alpha_prior, beta_prior, n_iter, init_p, "beta", c(alpha_prior, beta_prior))

# Run the sampler using posterior distribution as the proposal (should yield acceptance rate ~ 1)
alpha_post <- alpha_prior + y
beta_post <- beta_prior + N - y
samples_posterior <- mh_sampler_independence(y, N, alpha_prior, beta_prior, n_iter, init_p, "beta", c(alpha_post, beta_post))

# Diagnostic Plots
par(mfrow = c(4, 2))

# Define named list for samples for each proposal
sample_list <- list(
  "Uniform Proposal" = samples_uniform,
  "Normal Proposal" = samples_normal,
  "Beta Proposal" = samples_beta,
  "Posterior-based Proposal" = samples_posterior
)

# Plot trace and posterior distribution for each proposal type with proper titles
for (name in names(sample_list)) {
  samples <- sample_list[[name]]
  
  # Trace Plot
  plot(samples, type = 'l', main = paste("Trace Plot:", name), xlab = 'Iteration', ylab = 'p')
  abline(v = burnin, col = 'red', lty = 2)
  
  # Posterior Distribution Histogram
  hist(samples, breaks = 30, probability = TRUE, main = paste("Posterior Distribution:", name), xlab = 'p')
  lines(density(samples), col = 'blue', lwd = 2)
}

# Autocorrelation and Effective Sample Size (ESS) for each proposal
for (name in names(sample_list)) {
  samples <- sample_list[[name]]
  acf(samples, main = paste("Autocorrelation:", name))
  ess <- effectiveSize(samples[(burnin + 1):n_iter])
  cat(sprintf("Effective Sample Size for %s: %f\n", name, ess))
}

# Summary of Posterior Marginal Distributions (Moments and Quantiles)
posterior_summary <- function(samples) {
  mean_p <- mean(samples)
  sd_p <- sd(samples)
  n_samples <- length(samples)
  se <- sd_p / sqrt(n_samples)
  ess <- effectiveSize(samples)
  rel_eff <- ess / n_samples
  quantiles <- quantile(samples, probs = c(0.025, 0.5, 0.975))
  
  moments_summary <- data.frame(
    Metric = c("Mean", "Standard Deviation", "Standard Error", "Effective Sample Size", "Relative Efficiency"),
    Value = c(mean_p, sd_p, se, ess, rel_eff)
  )
  
  quantiles_summary <- data.frame(
    Quantile = c("2.5%", "50%", "97.5%"),
    Value = quantiles
  )
  
  cat("Moments Summary\n")
  print(kable(moments_summary, format = "pipe"))
  
  cat("\nQuantiles Summary\n")
  print(kable(quantiles_summary, format = "pipe"))
}

# Generate summary for each proposal's samples
for (name in names(sample_list)) {
  cat(sprintf("\nSummary for %s:\n", name))
  posterior_summary(sample_list[[name]][(burnin + 1):n_iter])
}
