rm(list=ls(all=TRUE))
library(bayesm)
source("SingleMoveMH.R")
source("rmnlRWMetrop_noTuning.r")
# The above two lines did not work for me, so I just copied the code in R in full awareness,
# that it looks a bit messy this way.


SingleMoveMH <- function(Data, Prior, Mcmc) {
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }
  # Extract data and check dimensions
  if (missing(Data)) pandterm("Requires Data argument -- list of p, y, X")
  if (is.null(Data$X)) pandterm("Requires Data element X")
  if (is.null(Data$y)) pandterm("Requires Data element y")
  if (is.null(Data$p)) pandterm("Requires Data element p")
  
  X <- Data$X
  y <- Data$y
  p <- Data$p
  nvar <- ncol(X)
  nobs <- length(y)
  
  if (length(y) != (nrow(X) / p)) pandterm("length(y) ne nrow(X)/p")
  if (sum(y %in% (1:p)) < nobs) pandterm("invalid values in y vector -- must be integers in 1:p")
  
  # Prior and MCMC settings
  if (missing(Prior)) {
    betabar <- rep(0, nvar)
    A <- 0.01 * diag(nvar)
  } else {
    betabar <- ifelse(is.null(Prior$betabar), rep(0, nvar), Prior$betabar)
    A <- ifelse(is.null(Prior$A), 0.01 * diag(nvar), Prior$A)
  }
  
  if (missing(Mcmc)) pandterm("requires Mcmc argument")
  R <- Mcmc$R
  keep <- ifelse(is.null(Mcmc$keep), 1, Mcmc$keep)
  sbeta <- ifelse(is.null(Mcmc$sbeta), 1.0, Mcmc$sbeta)
  nprint <- ifelse(is.null(Mcmc$nprint), 100, Mcmc$nprint)
  
  cat("Starting Single-Move Metropolis Sampler for Multinomial Logit Model", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  print(betabar)
  print(A)
  
  # Initialization
  betadraw <- matrix(0, nrow = floor(R / keep), ncol = nvar)
  loglike <- double(floor(R / keep))
  beta <- rep(0, nvar)
  naccept <- rep(0, nvar)  # Track acceptance rate for each beta element
  
  # Cholesky decomposition for prior covariance
  priorcov <- chol2inv(chol(A))
  rootpi <- backsolve(chol(priorcov), diag(nvar))
  
  # Initial log-likelihood and posterior
  oldloglike <- llmnl(beta, y, X)
  oldlpost <- -0.5 * sum((rootpi %*% (beta - betabar))^2) + oldloglike
  
  itime <- proc.time()[3]
  cat("MCMC Iteration (est time to end - min)", fill = TRUE)
  
  # Main MCMC loop
  for (rep in 1:R) {
    for (j in 1:nvar) {  # Loop over each element of beta !MAIN DIFFERENCE!!!
      betac <- beta  # Copy current beta
      repeat {
        betac[j] <- beta[j] + sbeta * rnorm(1)  # Propose new value for jth element
        
        # Apply truncation restrictions:
        if (j == 1 && betac[j] > 0) break  # First element must be negative
        if (j == 2 && betac[j] < 0) break  # Second element must be positive
        if (j == 3 && betac[j] < 0) break  # Third element must be positive
      }
      
      # Compute log-likelihood and posterior for proposed beta
      cloglike <- llmnl(betac, y, X)
      clpost <- -0.5 * sum((rootpi %*% (betac - betabar))^2) + cloglike
      
      # Compute acceptance probability
      ldiff <- clpost - oldlpost
      alpha <- min(1, exp(ldiff))
      
      if (runif(1) < alpha) {
        beta[j] <- betac[j]  # Accept
        oldloglike <- cloglike
        oldlpost <- clpost
        naccept[j] <- naccept[j] + 1
      }
    }
    
    # Save draws
    if (rep %% keep == 0) {
      betadraw[rep / keep, ] <- beta
      loglike[rep / keep] <- oldloglike
    }
    
    # Print progress
    if (rep %% nprint == 0) {
      ctime <- proc.time()[3]
      timetoend <- ((ctime - itime) / rep) * (R - rep)
      cat(" ", rep, " (", round(timetoend / 60, 1), ")", fill = TRUE)
    }
  }
  
  ctime <- proc.time()[3]
  cat("  Total Time Elapsed: ", round((ctime - itime) / 60, 2), "\n")
  
  attributes(betadraw)$class <- c("bayesm.mat", "mcmc")
  attributes(betadraw)$mcpar <- c(1, R, keep)
  return(list(betadraw = betadraw, loglike = loglike, acceptr = naccept / R))
}







simmnl= function(p,n,beta) {
  #   note: create X array with 1 alt.spec vars
  k=length(beta)
  X1=matrix(runif(n*p,min=-1,max=1),ncol=p)
  #X2=matrix(runif(n*p,min=-1,max=1),ncol=p)
  #X1=matrix(runif(n*p*(length(beta)-p+1),min=-1,max=1),ncol=p*(length(beta)-p+1))  2015-12-04 generalize the above two lines
  #X=createX(p,na=2,nd=NULL,Xd=NULL,Xa=cbind(X1,X2),base=1)
  X=createX(p,na=1,nd=NULL,Xd=NULL,Xa=cbind(X1),base=1)
  Xbeta=X%*%beta # now do probs
  p=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=p)
  Prob=exp(Xbeta)
  iota=c(rep(1,p))
  denom=Prob%*%iota
  Prob=Prob/as.vector(denom)
  # draw y
  y=vector("double",n)
  ind=1:p
  for (i in 1:n) 
  { yvec=rmultinom(1,1,Prob[i,]); y[i]=ind%*%yvec }
  return(list(y=y,X=X,beta=beta,prob=Prob))
}



set.seed(66)
n=500; p=3 # number of alternatives
#beta=c(-0.5,0.5,0.2,0.3)
beta=c(-0.5,0.5,0.2)
simout=simmnl(p,n,beta)
A=diag(c(rep(.01,length(beta)))); betabar=rep(0,length(beta))

Data=list(y=simout$y,X=simout$X,p=p)
Mcmc=list(R=20000,nprint=10000,keep=1)
Prior=list(A=A,betabar=betabar)


# Single-Move MH Sampler
Mcmc = list(R = 20000, nprint = 1000, keep = 1, sbeta = 0.1)  
outSM = SingleMoveMH(Data = Data, Prior = Prior, Mcmc = Mcmc)
#windows()
plot(outSM$betadraw, tvalues = beta)

Mcmc = list(R = 20000, nprint = 1000, keep = 1, sbeta = 0.005)  
outSM_small = SingleMoveMH(Data = Data, Prior = Prior, Mcmc = Mcmc)
plot(outSM_small$betadraw, tvalues=beta)

Mcmc = list(R = 20000, nprint = 1000, keep = 1, sbeta = 1)  
outSM_sbig = SingleMoveMH(Data = Data, Prior = Prior, Mcmc = Mcmc)
plot(outSM_sbig$betadraw, tvalues=beta)




cat("Acceptance Rate: ", round(outSM$acceptr * 100, 2), "%\n")
summary(outSM_small$betadraw,tvalues=beta)
cat("Acceptance Rate: ", round(outSM_small$acceptr * 100, 2), "%\n")
summary(outSM_sbig$betadraw,tvalues=beta)
cat("Acceptance Rate: ", round(outSM_sbig$acceptr * 100, 2), "%\n")