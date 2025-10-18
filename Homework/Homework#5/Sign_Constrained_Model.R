##################################################################################################
#####Simulation & Estimation - Constrained vs. Unconstrained (log-normal)#########################
##################################################################################################

rm(list=ls())

# LOAD LIBRARIES REQUIRED TO CREATE THE SIMULATED DATA. YOU MAY NEED TO INSTALL THESE PACKAGES.
library(devtools)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(ggplot2)
#library(tikzDevice)
#library(plyr)
library(latex2exp)
library(pracma)


set.seed(77)

###Setting
#Sample size
nunits = 1000

#Variance-covariance matrix of betastars
sigma = diag(c(0.3,.5, 1),nrow=3,ncol=3) #Specify variance covariance matrix of normally distributed RV 
tsigma = t(sigma)
sigma[lower.tri(sigma)] = tsigma[lower.tri(tsigma)]
#Expected values 
avgbeta = c(-0.2,-2, 3) ###price, log delta brand A, brand B,
#Draw
betastar = mvrnorm(n=nunits, avgbeta, sigma) #draw true betas from assumed distribution of heterogeneity
beta=betastar
#Price --> Price sensitivity is negative
beta[,1] = -exp(betastar[,1])

#Brand --> combine additivate and multiplicative effects
beta[,2] = betastar[,3] + exp(betastar[,2])

###Compare estimates to population data 
#Draw population
betastar_pop = mvrnorm(n=200000, avgbeta, sigma) #draw true betas from assumed distribution of heterogeneity
beta_pop=betastar_pop
#Price
beta_pop[,1] =-exp(betastar_pop[,1])
#Brand
beta_pop[,2] = betastar_pop[,3] + exp(betastar_pop[,2])

###Summary Population distribution
summary_population_distribution <- array(0,dim=c(7,3))
rownames(summary_population_distribution) = c("1%","25%","50%","75%","99%","Mean","Variance")
colnames(summary_population_distribution) = c("Price","brand A", "brand B")
summary_population_distribution[,1] = c(quantile(beta_pop[,1],probs=c(0.01,0.25,0.5,0.75,0.99)),mean(beta_pop[,1]),var(beta_pop[,1]))
summary_population_distribution[,2] = c(quantile(beta_pop[,2],probs=c(0.01,0.25,0.5,0.75,0.99)),mean(beta_pop[,2]),var(beta_pop[,2]))
summary_population_distribution[,3] = c(quantile(beta_pop[,3],probs=c(0.01,0.25,0.5,0.75,0.99)),mean(beta_pop[,3]),var(beta_pop[,3]))
summary_population_distribution

windows()
par(mfrow=c(2,2))
plot(density(beta_pop[,1]), main = "price");grid()
plot(density(beta_pop[,2]), main = "brand A");grid()
plot(density(beta_pop[,3]), main = "brand B");grid()
smoothScatter(beta_pop[,2],beta_pop[,3], main = "brand A versus B");grid();abline(0,1)



###Function to create artificial data set
simmnlv2 = function(n,beta){
  
  #
  # p. rossi 2004
  # Modified by Max Pachali & Thomas Otter 2014/2015/2016
  #
  # Purpose: simulate from MNL (including X values)
  #
  # Arguments:
  # p is number of alternatives per choice task in a alternative complete setting
  # n is number of choice tasks
  # beta is true parm value
  #
  # Output:
  # list of X 
  # y (indicator of choice --> 1, ...,p)
  # prob is a n x p matrix of choice probs
  #
  
  k = length(beta)
  ###MATRIX OF ALL POSSIBLE ATTRIBUTE LEVEL COMBINATIONS (ONLY ONE BRAND RELATIVE TO OUTSIDE)
  X_full_original = matrix(c(1,0,0,0,1,0),ncol=2)
  ###DEDUCE p
  p = dim(X_full_original)[1]
  ###REPLICATE ORIGINAL MATRIX n TIMES --> n*p rows, the third option is outside option (no brand is chosen)
  for(i in 1:n){
    if(i==1){
      X_full = X_full_original
    }else{
      X_full = rbind(X_full,X_full_original)
    }
  }
  ###ADD PRICE --> random prices between 0.5 and 3 for all n*p rows
  X_full = cbind(runif(n*p,min=0.5,max=3),X_full)
  ###Make prices for brand A larger --> Vertical differentiation
  for(i in 1:n){
    X_full[i*p-2,1] = X_full[i*p-1,1] + runif(1,min=0.6,max=3)
  }
  ###DELETE PRICES FOR OUTSIDE OPTION --> Prices are set to zero for outside option
  for(i in 1:n){
    X_full[i*p,1] = 0
  }
  ###CONSTRUCT PROBABILITIES
  Xbeta=X_full%*%beta
  p=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=p) #reshape Xbeta matrix to be of (nxp)-dimension
  Prob=exp(Xbeta)
  iota=rep(1,p)
  denom=Prob%*%iota
  # Softmax function for choice probabilities
  Prob=Prob/as.vector(denom)
  ###DRAW CHOICES Y
  y=vector("double",n)
  ind=1:p
  for (i in 1:n){
    yvec=rmultinom(1,1,Prob[i,])
    y[i]=ind%*%yvec
  }
  
  return(list(y=y,X=X_full,beta=beta,prob=Prob))
}

###CREATE MULTINOMIAL LOGIT y AND X FOR EACH UNIT USING BETA
cmax_1 = 7
cmax_2 = 60    
cmax <- cmax_1    # set cmax <- cmax_2 if you want to switch to model with 60 choices per unit
lgtdata_sim=NULL
hdata=NULL
for (i in 1:nunits) {
  hdata[[i]] = simmnlv2(cmax,beta[i,])
}
###Create list 
for(i in 1:nunits){
  lgtdata_sim[[i]]=list(X=hdata[[i]]$X,y=hdata[[i]]$y) 
}

#View(lgtdata_sim[[1]]$X)
#View(lgtdata_sim[[1]]$y)

N = length(lgtdata_sim)
t = cmax
p = dim(lgtdata_sim[[1]]$X)[1]/cmax

###Prepare estimation data
E_Data=list(p=p,lgtdata=lgtdata_sim) #construct a list with all relevant data for estimation sample

Mcmc = list(R=12000,keep=2)
Prior = list(ncomp=1, SignRes = c(-1, 0, 0))
out_constrained = rhierMnlRwMixture(Data=E_Data,Mcmc=Mcmc,Prior=Prior)

# posterior draws of mu -> population level
bbarmc=matrix(double(6000*3),ncol=3)

for (r in 1:6000){
  bbarmc[r,] = out_constrained$nmix$compdraw[[r]][[1]]$mu
}

betadrawconverged = out_constrained$betadraw[,,3001:6000]

# individual level 
betaexchange <- array(aperm(betadrawconverged, perm=c( 1 , 3 , 2 )),
                      dim=c(dim(betadrawconverged)[1] * dim(betadrawconverged)[3],
                            dim (betadrawconverged)[2]))


# PLOTS for INTERPRETATION

# assess convergence with log likelihood
windows()
plot(out_constrained$loglike, type = 'l')

# Density plots for price, brand A, and brand B coefficients
windows()
par(mfrow = c(2, 2))  # 2x2 layout for plots
plot(density(betaexchange[, 1]), main = "Price Coefficient (β1)", col = "blue", lwd = 2); grid()
plot(density(betaexchange[, 2]), main = "Brand A Coefficient (β2)", col = "red", lwd = 2); grid()
plot(density(betaexchange[, 3]), main = "Brand B Coefficient (β3)", col = "green", lwd = 2); grid()
smoothScatter(betaexchange[, 2], betaexchange[, 3],
              main = "Brand A vs Brand B (β2 vs β3)", xlab = "β2 (Brand A)", ylab = "β3 (Brand B)"); 
abline(0, 1, col = "black", lwd = 2)  # Line where β2 = β3


# Summary statistics for the individual-level parameters
beta_summary <- apply(betaexchange, 2, function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x), min = min(x), max = max(x))
})
rownames(beta_summary) <- c("Price (β1)", "Brand A (β2)", "Brand B (β3)")
beta_summary

# Plot posterior means of population-level preferences
windows()
matplot(bbarmc, type = "l", col = 1:3, lty = 1, lwd = 2,
        main = "Population-Level Posterior Means", xlab = "Iteration", ylab = "Posterior Means")
legend("center", legend = c("Price (β1)", "Brand A (β2)", "Brand B (β3)"),
       col = 1:3, lty = 1, lwd = 2)
windows()
plot(bbarmc[,2],bbarmc[,3]);abline(0,1)

# Calculate summary statistics (mean, median, sd, quantiles) for each parameter in bbarmc
summary_stats <- t(apply(bbarmc, 2, function(x) c(quantile(x, probs = c(0.01, 0.25, 0.5, 0.75, 0.99)), mean = mean(x), median = median(x), sd = sd(x))))

# Convert the matrix to a data frame and add row names for clarity
summary_stats_df <- t(data.frame(summary_stats))
rownames(summary_stats_df) <- c("1%", "25%", "50%", "75%", "99%", "Mean", "Median", "SD")
colnames(summary_stats_df) <- c("β1 (Price)", "β2 (Brand A)", "β3 (Brand B)")

# Print the summary table
print(summary_stats_df)





# PROFIT MAXIMIZATION

# Function to compute the demand for Brand A based on prices and model coefficients
calculate_demand_A <- function(p_A, p_B, beta) {
  # Utility for Brand A
  utility_A = beta[1] * p_A + beta[2] * p_B
  
  # Choice probability for Brand A (logit choice probability)
  prob_A = exp(utility_A) / (1 + exp(utility_A))
  return(prob_A)
}

# Function to compute the profit for Brand A
calculate_profit_A <- function(p_A, p_B, beta, marginal_cost_A) {
  # Calculate the demand for Brand A
  demand_A = calculate_demand_A(p_A, p_B, beta)
  
  # Profit = (Price - Cost) * Demand
  profit_A = (p_A - marginal_cost_A) * demand_A
  return(profit_A)
}

# Optimization function to find the optimal price for Brand A
optimize_price_A <- function(beta_samples, p_B, marginal_cost_A) {
  optimal_prices <- numeric(nrow(beta_samples))  # Store optimal prices for each posterior draw
  
  # Iterate over each posterior sample (each row corresponds to a set of parameters)
  for (i in 1:nrow(beta_samples)) {
    # Define a function to maximize for this posterior sample
    profit_function <- function(p_A) {
      calculate_profit_A(p_A, p_B, beta_samples[i, ], marginal_cost_A)
    }
    
    # Use optim() to maximize profit with respect to p_A (initial guess: p_A = 1)
    result <- optim(1, profit_function, lower=0, upper=5, method="L-BFGS-B", control = list(fnscale = -1))
    
    # Store the optimal price for this posterior sample
    optimal_prices[i] <- result$par
  }
  
  # Return the optimal prices for all posterior samples
  return(optimal_prices)
}

# Set fixed price for Brand B
p_B <- 0.5

# Marginal cost for Brand A
marginal_cost_A <- 1

# Compute optimal prices for Brand A using the population-level posterior (bbarmc)
optimal_prices_A_pop <- optimize_price_A(bbarmc, p_B, marginal_cost_A)

# Compute optimal prices for Brand A using individual-level posterior (betaexchange)
#optimal_prices_A_ind <- optimize_price_A(betaexchange, p_B, marginal_cost_A)

# Calculate the average optimal price for both population-level and individual-level prices
mean_optimal_price_pop <- mean(optimal_prices_A_pop)
#mean_optimal_price_ind <- mean(optimal_prices_A_ind)

# Print out the results
cat("Mean Optimal Price for Brand A (Population-Level):", mean_optimal_price_pop, "\n")
#cat("Mean Optimal Price for Brand A (Individual-Level):", mean_optimal_price_ind, "\n")

# Visualizing the profit curve
profit_curve_posterior <- function(price_range, p_B, beta_samples, marginal_cost_A) {
  n_samples <- nrow(beta_samples)
  n_prices <- length(price_range)
  profits <- matrix(0, nrow = n_samples, ncol = n_prices)
  
  for (i in 1:n_samples) {
    for (j in 1:n_prices) {
      profits[i, j] <- calculate_profit_A(price_range[j], p_B, beta_samples[i, ], marginal_cost_A)
    }
  }
  return(profits)
}

# Define price range for Brand A
price_range <- seq(0.5, 20, length.out = 100)

# Compute profits for all posterior samples (using population-level posterior bbarmc)
profits_posterior <- profit_curve_posterior(price_range, p_B, bbarmc, marginal_cost_A)

# Calculate mean and 95% credible intervals for profits at each price level
mean_profits <- apply(profits_posterior, 2, mean)
lower_bound <- apply(profits_posterior, 2, quantile, probs = 0.025)
upper_bound <- apply(profits_posterior, 2, quantile, probs = 0.975)

# Plot the mean profit curve with credible intervals
plot(price_range, mean_profits, type = "l", col = "blue", lwd = 2,
     xlab = "Price of Brand A", ylab = "Profit", main = "Profit Curve with Posterior Uncertainty")
polygon(c(price_range, rev(price_range)), c(upper_bound, rev(lower_bound)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(price_range, mean_profits, col = "blue", lwd = 2)
grid()