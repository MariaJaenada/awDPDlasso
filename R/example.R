rm(list=ls())

library(MASS)
library(glmnet)

source("dpd_estimation.R")
source("logit.R")
source("hgic_dpd_logit.R")
source("hgic_dpd_linear.R")

###################################
# example for linear regression #
###################################
p=500
n=100
alpha = 0.3 #DPD loss tunning parameter

num_signal = 1 #number of nonzero coefficients (x3)
epsilon = .1 #outliers percent


## paramter beta and covariance matrix of covariates
tB = rep(0,p)
for(b in 1:num_signal){tB[(20*(b-1)+1):(20*(b-1)+5)] = c(3,1.5,0,0,2)}

Sigma = diag(p)
for(i in 2:p){
  for(j in 1:(i-1)){
    Sigma[i,j] = .5^abs(i-j)
    Sigma[j,i] = Sigma[i,j]
  }
}

#Generate data
X = mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
Y = drop(X %*% tB + rnorm(n)*.5)

## add outliers on Y
nout = ceiling(100*epsilon)
if (nout > 0){Y[1:nout] = Y[1:nout] + 20 +rnorm(nout)}

z <- awDPDlasso.lm(X, as.matrix(Y), intercept=T, init.mode="RLARS", init.model=NULL, penalty = "adaptive",
              weight.mode = "lasso",lambda.mode="lambda0", lmin=0.01, lmax=0.1, gam=alpha, nlambda = 50)
z$best$beta[z$best$beta !=0]

z <- awDPDlasso.lm(X, as.matrix(Y), intercept=T, init.mode="RLARS", init.model=NULL, penalty = "adaptive",
                   weight.mode = "adaptive",lambda.mode="lambda0", lmin=0.01, lmax=0.1, gam=alpha, nlambda = 50)
z$best$beta[z$best$beta !=0]

z <- awDPDlasso.lm(X, as.matrix(Y), intercept=T, init.mode="RLARS", init.model=NULL, penalty = "adaptive",
                   weight.mode = "scad",lambda.mode="lambda0", lmin=0.01, lmax=0.1, gam=alpha, nlambda = 50)
z$best$beta[z$best$beta !=0]

###################################
# example for logistic regression #
###################################

k=500 #number of variables
n=100 #sample size
alpha = 0.3 #DPD loss tunning parameter

num_signal = 1 #number of nonzero coefficients
epsilon = 0.05 #outliers percent


## True regression parameters
tB = rep(0,k)
for(b in 1:num_signal){tB[(20*(b-1)+1):(20*(b-1)+5)] = c(5,5,0,0,5)}

## Variance-covariance matrix of explanatory variables
Sigma = diag(k)
  for(i in 2:k){
    for(j in 1:(i-1)){
      Sigma[i,j] = 0.1^abs(i-j)
      Sigma[j,i] = Sigma[i,j]
    }
  }
  
#simulate data
X = mvrnorm(n, mu=rep(0,k), Sigma=Sigma)

pi = sapply(drop(X %*% tB), invlogit)
pi[pi<1e-4]=0; pi[pi>(1-1e-4)]=1 ###avoiding large estimates 
    
if(epsilon >0){
  raffle = rbinom(size=1,n=n,p=(1-epsilon))
  Y1 =sapply(pi,function(x) rbinom(size=1,n=1,p=x)) # pure data
  Y2 = sapply(pi,function(x) rbinom(size=1,n=1,p=(1-x))) # contamination
  Y  = raffle*Y1+(1-raffle)*Y2
}else{Y = sapply(pi,function(x) rbinom(size=1,n=1,p=x))}

X=cbind(1,X)

# initial estimator
LS_LASSO<-cv.glmnet(X[,2:(k+1)], Y, family = "binomial", intercept=T)
beta.init = coef(LS_LASSO)[1:(k+1)]
beta.init[beta.init!=0]

#different adaptive estimators based on DPD loss 
lmin = 0.005;lmax=0.5
#lasso
z <- awDPDlasso.logistic(X,Y, beta.init, weight.mode = "lasso",LR.loss = "LS", lambda.mode="given",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]

#adaptive
z <- awDPDlasso.logistic(X,Y, beta.init, weight.mode = "adaptive", LR.loss = "LS", lambda.mode="lambda0",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]
      
#scad
z <- awDPDlasso.logistic(X,Y, beta.init, weight.mode = "scad", LR.loss = "LS", lambda.mode="lambda0",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]