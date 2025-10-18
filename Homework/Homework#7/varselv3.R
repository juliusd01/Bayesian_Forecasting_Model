probvs=function(Data, Prior, Mcmc, Init){
  
  
	# variable selection a la Wachtel and Otter 2013 built into Rossi's rbprobitGibbs (see the bayesm package)
	# this version does the following
	# -	use p(tau_k) ~ beta(2,2) as a prior for all k = 1, ..., p
	# -	update this prior to p(tau_k=1) = beta(2+sum(tau), 2 + length(tau) - sum(tau)) 
	# -	sample a new p(tau_k=1) from this distribution
	# -	use this p(tau_k=1) in the metropolis ratio
	
  breg1 = function(root, X, y, Abetabar) {
    cov = crossprod(root, root)
    betatilde = cov %*% (crossprod(X, y) + Abetabar)
    return(list(btilde=betatilde,beta=betatilde + t(root) %*% rnorm(length(betatilde))))
  }	
  
	X=Data$X
	y=Data$y
	A=Prior$A
	aBeta=Prior$aBeta
	bBeta=Prior$bBeta
	betabar=Prior$mu
	R=Mcmc$R
	keep=Mcmc$keep
	z=Init$z
	tau=Init$tau
	# Getting the dimensions
	nobs=length(y)
	nvarmax=dim(X)[2]
		
	# Allocating memory for the draws
	betadraw=matrix(rep(0,nvarmax*(floor(R/keep))),ncol=nvarmax)
	taudraw=matrix(rep(0,(nvarmax-1)*(floor(R/keep))),ncol=nvarmax-1)
	loglike=rep(0,floor(R/keep))
	alphadraw=matrix(rep(0,(floor(R/keep))),ncol=1)
	naccept=0
	
	# Initializing the chain
	beta=rep(0,nvarmax)
    a = ifelse(y == 0, -100, 0)
    b = ifelse(y == 0, 0, 100)
	Atau=A[as.logical(tau),as.logical(tau)]	# indexed A by the corresponding tau
	sigma=c(rep(1,nobs))	# fix sigma to one
	Xtau = X[,as.logical(tau)]
	alph=rep(rbeta(1,aBeta,bBeta),nvarmax-1)
	
	
	itime = proc.time()[3]
	cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
	for (rep in 1:R){
		# updating augmented latent "utilities" z | beta, y
		betatau=beta[as.logical(tau)]
		mu=Xtau%*%as.matrix(betatau)
		z = rtrun(mu, sigma, a, b)
		# computing moments for tau sampling | z, theta (alph)
		root = chol(chol2inv(chol((crossprod(Xtau, Xtau) + Atau))))
		betabartau=betabar[as.logical(tau)]
		Abetabar = crossprod(Atau, betabartau)
		bregout=breg1(root, Xtau, z, Abetabar)
		res=z-Xtau%*%as.matrix(bregout$btilde)
		s2=crossprod(res,res)
		lold=.5*determinant(as.matrix(Atau),log=T)$modulus[1]-s2/2-.5*determinant(as.matrix(crossprod(Xtau, Xtau) + Atau),log=T)$modulus[1]
		
		# tau proposal and moments at proposal
		ctau=rbinom(nvarmax-1,1,alph) # generate proposal from hiearchical prior
		Atau=A[as.logical(c(1,ctau)),as.logical(c(1,ctau))]	# indexed A by the corresponding ctau
		sigma=c(rep(1,nobs))	# fix sigma to one
		Xtau = X[,as.logical(c(1,ctau))]
		betabartau=betabar[as.logical(c(1,ctau))]
		root = chol(chol2inv(chol((crossprod(Xtau, Xtau) + Atau))))
		Abetabar = crossprod(Atau, betabartau)
		bregout=breg1(root, Xtau, z, Abetabar)
		res=z-Xtau%*%bregout$btilde
		s2=crossprod(res,res)
		clpost=.5*determinant(as.matrix(Atau),log=T)$modulus[1]-s2/2-.5*determinant(as.matrix(crossprod(Xtau, Xtau) + Atau),log=T)$modulus[1]
		ldiff=clpost-lold
		alpha = min(1, exp(ldiff))
		if (alpha < 1) {
		  unif = runif(1)
		}
		else {
		  unif = 0
		}
		if (unif <= alpha) {
		  lold=clpost
		  tau=c(1,ctau)
		  naccept = naccept + 1
		  beta=rep(0,nvarmax)  # reset beta to zero (only active betas will be updated below)
		}
		# updating prior selection probabilities 
		alph=rep(rbeta(1,aBeta+sum(tau[-1]),bBeta+nvarmax-1-sum(tau[-1])),nvarmax-1)
				Atau=A[as.logical(tau),as.logical(tau)]	# indexed A by the corresponding tau
		# updating beta | tau, z
				sigma=c(rep(1,nobs))	# fix sigma to one
		Xtau = X[,as.logical(tau)]
		betabartau=betabar[as.logical(tau)]
		root = chol(chol2inv(chol((crossprod(Xtau, Xtau) + Atau))))
		Abetabar = crossprod(Atau, betabartau)
		bregout=breg1(root, Xtau, z, Abetabar)
		beta[as.logical(tau)] = bregout$beta
		if (rep%%100 == 0) {
		  ctime = proc.time()[3]
		  timetoend = ((ctime - itime)/rep) * (R - rep)
		  cat(" ", rep, " (", round(timetoend/60, 1), ")", 
		      fill = TRUE)
		}
        if (rep%%keep == 0) {
             mkeep = rep/keep
             betadraw[mkeep, ] = beta
			 taudraw[mkeep, ] = tau[-1]
			 loglike[mkeep] = lold
			 alphadraw[mkeep] = alph[1]
         }
	}
    ctime = proc.time()[3]
    cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
        "\n")
    attributes(betadraw)$class = c("bayesm.mat", "mcmc")
    attributes(betadraw)$mcpar = c(1, R, keep)

    return(list(betadraw = betadraw, loglike = loglike, acceptr = naccept/R, taudraw=taudraw, alphadraw=alphadraw))
	}
