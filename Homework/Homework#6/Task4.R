

# Load the bayesm package
library(bayesm)


E_Data <- load("G:/Dokumente/UNI/Goethe/1. Semester/Bayesian Modelling/Homework6/Test_2.RData")
p <- 3

# Combine y vectors and X matrices from all respondents
y_combined <- unlist(lapply(lgtdata, function(x) x$y))            # Combine all y vectors
X_combined <- do.call(rbind, lapply(lgtdata, function(x) x$X))   # Combine all X matrices

# Verify dimensions
n <- length(y_combined)  # Total number of choice tasks
cat("Total number of choice tasks (n):", n, "\n")
cat("Dimensions of combined X (should be n*p x k):", dim(X_combined), "\n")

# Create the Data object
Data <- list(y = y_combined, X = X_combined, p = p)

# Inspect the structure of the Data object
str(Data)

R=240000
keep=5


# Set the prior parameters
Prior <- list(
  betabar = rep(0, ncol(Data$X)), # Prior mean for beta coefficients
  A = diag(0.01, ncol(Data$X))   # Prior precision matrix (weakly informative)
)

# Set MCMC parameters
Mcmc <- list(
  R = R,    # Total number of MCMC iterations
  keep = keep, # Number of iterations to keep for posterior analysis
  nprint = 5000  # Print progress every 5000 iterations
)

# Call the rmnlIndepMetrop function
set.seed(123) # Set seed for reproducibility
output <- rmnlIndepMetrop(Data = Data, Prior = Prior, Mcmc = Mcmc)

# Summarize the results
cat("Summary of MCMC output:\n")
summary(output$betadraw)

# Analyze graphically
plot(output$betadraw)