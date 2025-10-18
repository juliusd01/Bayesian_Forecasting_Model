###Clear working directory
rm(list=ls())

###Load packages


###set plot flag
plotflag = F  # if F, no plots are are displayed

###set seed flag
seedflag = T # if F, random seeds are used



###########################################################
######## Simulating data from a binomial likelihood #######
###########################################################


###Set seed (only for purposes of *exact* replication.  Beware of seeds in real work!!)
if (seedflag){set.seed(66)}

p=.3        ## the probability of success    !!!!!!!!!!!!!!!!!
N= 500       ## the number of binomial trials !!!!!!!!!!!!!!!!!

## 1) draw N independent random unform numbers on the interval [0,1]
## 2) check for each number if smaller than p.  --> p "successes" in expectation
y=runif(N)<p


##############################
# subjective prior parameters
# subjective prior parameters
# subjective prior parameters
# subjective prior parameters
##############################
# values close to one define a neutral (non-informative) prior
# larger values --> more informative priors
# larger (smaller) ratios alpha/beta --> higher (lower) prior success expectations
# values close to zero produce a bathtub-shaped prior favoring parameter values at the boundary

alpha=70 #!!!!!!!!!!!!!!!
beta=30  #!!!!!!!!!!!!!!!


###########################################################
######## Preparing for analysis ###########################

## define (non-normalized) likelihood function
binloglik <- function(p,y){
  N=length(y)
  Y=sum(y)
  ll=Y*log(p) + (N-Y)*log(1-p)
}

## define (non-normalized) likelihood function
binloglikopt <- function(pstar,y){
  p=exp(pstar)/(1+exp(pstar))
  N=length(y)
  Y=sum(y)
  ll=Y*log(p) + (N-Y)*log(1-p)
}


## define (non-normalized) log prior
logprior <- function(alpha,beta,theta){
  #if ((theta = 0) | (theta = 1)){
  #  lprior=-.Machine$double.xmax
  #}else{
  lprior=(alpha-1)*log(theta) + (beta-1)*log(1-theta)
}

## define log-posterior
logpost <- function(theta,y,alpha, beta){
  lpost = binloglik(theta,y) + logprior(alpha,beta,theta)
}


## define log-posterior with transformed parameter
## (direct optimization on the interval [0,1], with non-finite functional values at the bounds crashes)
logpostopt <- function(thetastar,y,alpha, beta){
  theta=exp(thetastar)/(1+exp(thetastar))
  lpost = binloglik(theta,y) + logprior(alpha,beta,theta)
}


###################################################
## classical and Bayesian analysis by optimization
## classical and Bayesian analysis by optimization
## classical and Bayesian analysis by optimization
## plot log-likelihood-profile, log-prior and log-posterior
###################################################
thetagrid=seq(0,1,.001) # a grid of parameter values

outlik=double(length(thetagrid)) 
outlprior=double(length(thetagrid))
outlpost=double(length(thetagrid))

# evaluate likelihood, prior and posterior for each parameter value (on the log-scale)
for (i in 1:length(thetagrid)){
  outlik[i]=binloglik(thetagrid[i],y)
  outlprior[i]=logprior(alpha,beta,thetagrid[i])
  outlpost[i]=logpost(thetagrid[i],y,alpha,beta)
}


# find maxima
mlik=which(outlik==max(outlik, na.rm=T))
mprior=which(outlprior==max(outlprior, na.rm=T))
mpost=which(outlpost==max(outlpost, na.rm=T))

# use optimizer to maximize log-likelihood
# use optimizer to maximize log-likelihood
# use optimizer to maximize log-likelihood
# here I use the function binloglikopt that works on thetastar where theta = exp(thetastar)/(1+exp(thetastar))
# this turns a constrained problem into an unconstrained one
outopt=optim(par=0,binloglikopt,gr=NULL,method = "BFGS", y=y,
             hessian=TRUE, control = list(fnscale=-1,trace=1))



# compute numerical standard for *transformed* parameter from the Hessian
hilf <- function(x){
  tryCatch(
    # This is what I want to do...
    {
      solve(-x)^.5
    },error=function(error_message){
      # ... but if an error occurs, tell me what happened: 
      .Machine$double.xmax
    }
  )
}
sellopt=hilf(outopt$hessian)
# computing ML-estimate
thetahatll=exp(outopt$par)/(1+exp(outopt$par))
# use optimHess and the function logpost to get the standard error for theta
hilf2 <- function(x,y){
  tryCatch(
    # This is what I want to do...
    {
      out=optimHess(par=x,binloglik,gr=NULL, y=y, 
                    control = list(fnscale=-1,trace=1))
      return(out)
    },error=function(error_message){
      # ... but if an error occurs, tell me what happened: 
      -.Machine$double.xmax
    }
  )
}
out=hilf2(x=thetahatll,y=y)
# numerical standard error
sell=hilf(out)


hilf3 <- function(x,y,alpha,beta){
  tryCatch(
    # This is what I want to do...
    {
      out=optimHess(par=x,logpost,gr=NULL, y=y, alpha=alpha, beta=beta, 
                    control = list(fnscale=-1,trace=1))
      return(out)
    },error=function(error_message){
      # ... but if an error occurs, tell me what happened: 
      -.Machine$double.xmax
    }
  )
}




# use optimizer to maximize log-posterior
# use optimizer to maximize log-posterior
# use optimizer to maximize log-posterior
# here I use the function logpostopt that works on thetastar where theta = exp(thetastar)/(1+exp(thetastar))
# this turns a constrained problem into an unconstrained one
outopt=optim(par=0,logpostopt,gr=NULL,method = "BFGS", y=y, alpha=alpha, beta=beta, 
          hessian=TRUE, control = list(fnscale=-1,trace=1))
# compute numerical standard for *transformed* parameter from the Hessian
seopt=solve(-outopt$hessian)^.5
# computing ML-estimate
thetahat=exp(outopt$par)/(1+exp(outopt$par))
# use optimHess and the function logpost to get the standard error for theta
out=hilf3(x=thetahat,y=y,alpha=alpha,beta=beta)
# numerical standard error
se=hilf(out)


####################################################
# Draw from prior
####################################################
R=100000
draws=rbeta(R,alpha,beta)
priormean=mean(draws) # computing the posterior mean
# a 95% posterior credible interval
upperpr=sort(draws)[ceiling(R*.975)]
lowerpr=sort(draws)[ceiling(R*.025)]


####################################################
# Bayesian analysis by iid simulation from posterior
###################################################
R=100000
draws=rbeta(R,sum(y)+alpha,N-sum(y)+beta)
postmean=mean(draws) # computing the posterior mean
# a 95% posterior credible interval
upper=sort(draws)[ceiling(R*.975)]
lower=sort(draws)[ceiling(R*.025)]


####################################################
# plot log-likelihood, log-prior, and log-posterior
####################################################

if (plotflag){
  windows()
  par(mfrow=c(2,2))
  plot(thetagrid,outlik,main='Loglik-profile',xlab = 'theta', ylab = 'log-likelihood');grid();
  abline(v=thetagrid[mlik], col='red', lwd=2)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  
  plot(thetagrid,outlprior,main='Logprior-profile',xlab = 'theta, greendash::priorCI', ylab = 'log-prior');grid();
  abline(v=thetagrid[mprior], col='green', lwd=2)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)

  plot(thetagrid,outlpost,main='Logposterior-profile',xlab = 'theta', ylab = 'log-posterior');grid();
  abline(v=thetagrid[mpost], col='blue', lwd=2)
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  
  plot(thetagrid,outlpost,main='Logposterior-profile',xlab = 'theta', ylab = 'log-posterior');grid()
  abline(v=thetagrid[mlik], col='red', lwd=3)
  abline(v=thetagrid[mprior], col='green', lwd=3)
  abline(v=thetagrid[mpost], col='blue', lwd=3)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  
  
  ##########################################################################
  # overlay with results from optimization and compare to simulation results
  ##########################################################################
  windows()
  par(mfrow=c(2,1))
  plot(thetagrid,outlpost,main='Logposterior-profile',
       xlab = 'theta, g::prior, r::lik, b::post, p::optim', ylab = 'log-posterior');grid()
  abline(v=thetagrid[mlik], col='red', lwd=2)
  abline(v=thetagrid[mprior], col='green', lwd=2)
  abline(v=thetagrid[mpost], col='blue', lwd=2)
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
    if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  
  hist(draws, main = 'posterior draws', xlim=c(0,1), xlab= 'theta, r::lik, m::postmean, mdash::posteriorCI')
  abline(v=postmean, col='magenta', lwd=4, lty=1)
  abline(v=upper, col='magenta', lwd=3, lty=2)
  abline(v=lower, col='magenta', lwd=3, lty=3)
  abline(v=thetagrid[mlik], col='red', lwd=2)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  
  
  # this can be adjusted to zoom into where the maxima occur
  # this can be adjusted to zoom into where the maxima occur
  # this can be adjusted to zoom into where the maxima occur
  indcenter=which((thetagrid-postmean)^2==min((thetagrid-postmean)^2))
  indlower=which((thetagrid-c(.75*lower))^2==min((thetagrid-c(.75*lower))^2))
  indhigher=which((thetagrid-c(1.5*upper))^2==min((thetagrid-c(1.5*upper))^2))
  
  index=(indlower:indhigher)  ## the indexing is into thetagrid
  
  windows()
  par(mfrow=c(2,2))
  plot(thetagrid[index],outlik[index],main='Loglik-profile',xlab = 'theta', ylab = 'log-likelihood');grid();
  abline(v=thetagrid[mlik], col='red', lwd=2)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  
  
  plot(thetagrid[index],outlprior[index],main='Logprior-profile',xlab = 'theta', ylab = 'log-prior');grid();
  abline(v=thetagrid[mprior], col='green', lwd=2)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)

  plot(thetagrid[index],outlpost[index],main='Logposterior-profile',xlab = 'theta', ylab = 'log-posterior');grid();
  abline(v=thetagrid[mpost], col='blue', lwd=2)
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  
  
  plot(thetagrid[index],outlpost[index],main='Logposterior-profile',xlab = 'theta', ylab = 'log-posterior');grid()
  abline(v=thetagrid[mlik], col='red', lwd=2)
  abline(v=thetagrid[mprior], col='green', lwd=2)
  abline(v=thetagrid[mpost], col='blue', lwd=2)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  
  
  
  windows()
  par(mfrow=c(2,1))
  plot(thetagrid[index],outlpost[index],main='Logposterior-profile',xlab = 'theta, g::prior, r::lik, b::post, p::optim', ylab = 'log-posterior');grid()
  abline(v=thetagrid[mlik], col='red', lwd=2)
  abline(v=thetagrid[mprior], col='green', lwd=2)
  abline(v=thetagrid[mpost], col='blue', lwd=2)
  #abline(v=thetahat, col='purple3', lwd=4, lty=2)
  abline(v=thetahat+2*se, col='purple3', lwd=3, lty=2)
  abline(v=thetahat-2*se, col='purple3', lwd=3, lty=3)
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  
  hist(draws, main = 'posterior draws', xlim=c(thetagrid[index[1]],thetagrid[index[length(index)]]),
       xlab= 'theta, r::lik, m::postmean, mdash::posteriorCI')
  abline(v=postmean, col='magenta', lwd=4, lty=1)
  abline(v=upper, col='magenta', lwd=3, lty=2)
  abline(v=lower, col='magenta', lwd=3, lty=3)
  abline(v=thetagrid[mlik], col='red', lwd=2)
  if (is.finite(sell)){
    abline(v=thetahatll+2*sell, col='red', lwd=3, lty=2)
    abline(v=thetahatll-2*sell, col='red', lwd=3, lty=3)
  }
  abline(v=upperpr, col='green', lwd=3, lty=2)
  abline(v=lowerpr, col='green', lwd=3, lty=3)
  
  
  ##########################################################################
  # is the posterior normal?
  ##########################################################################
  windows()
  qqnorm(draws, pch = 1, frame = FALSE)
  qqline(draws, col = "steelblue", lwd = 2)
}


##########################################################################
# Point summary of inference
##########################################################################
summarye=rbind(cbind(thetahatll-2*sell,thetahatll,thetahatll+2*sell,4*sell),
               cbind(thetahat-2*se,thetahat,thetahat+2*se,4*se),
               cbind(lower,postmean,upper,upper-lower))

colnames(summarye) <- c('2.5%', 'point-estimate', '97.5%', 'CI-length')
row.names(summarye) <- c('ML', 'Bayes-normal-approx', 'Bayes')



##########################################################################
# Display numeric results
##########################################################################

print(round(summarye,digits = 3))


