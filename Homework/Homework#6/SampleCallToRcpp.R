

#####################################
####sign and order constrained############
########sign and order constrained############
############sign and order constrained############
################sign and order constrained############
####################sign and order constrained############
E_Data <- load("G:/Dokumente/UNI/Goethe/1. Semester/Bayesian Modelling/Homework6/Test_2.RData")

R=240000
keep=5

###Source cpp function for estimation
Rcpp::sourceCpp("rhierMnlRwMixture_rcpp_loop_Sim_modHP_solution.cpp",showOutput = FALSE)
source('rhierMnlRwMixture_main.R')


# MODEL 5i)
#number of constrained coefficients (here price, L1 and L2)
nvar_c = 0

Prior = list(ncomp=1)
Mcmc = list(R=R,keep=keep)

Data <- list(
  p = 3,          
  lgtdata = E_Data,
  Z = NULL
)

# Hierarchical Model
out_order = rhierMnlRwMixture_SR(Data=Data,
                                 Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,flag="approx")



# MODEL 5ii)
nvar_c = 1

###Prior setting
Amu = diag(1/10, nrow = nvar_c, ncol = nvar_c)
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
nu = 15 + nvar_c
V = nu * diag(nvar_c)*0.5

Prior = list(ncomp=1, Amu = Amu, mustarbarc = mustarbarc, nu = nu, V = V)
Mcmc = list(R=R,keep=keep)

# Hierarchical Model
out_order = rhierMnlRwMixture_SR(Data=E_Data,
                                 Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,flag="approx")





  
