# Function to choose the value of the regularization parameter lambda using HGIC criterion
# X : n x k matrix containing observations
# Y : n binary vector containing response ; only 0 and 1 are allowed
# beta.init: initial estimator
# weight.mode: adaptive weight mode (lasso, adaptive, scad)
# LR. loss: loss used to solve the auxiliar least squares problem. Defeault LS (least squares)
# alpha : DPD tunning parameter controlling trade-off between efficiency and robustness
# lambda: penalty parameter controlling the shrinkage
# lambda.mode: grid mode (given or an improvement)
# nlambda : length of the grid


awDPDlasso.logistic <- function(X,Y, beta.init,
                             weight.mode = "adaptive",
                             LR.loss = "LS", 
                             lambda.mode=c("lambda0","given"),
                             lmax=.2,lmin=0.01, nlambda=50,
                             alpha =0.1){
  
  
  ####################################################
  # detect the data type and tuning parameters value #
  ####################################################
  if(is.matrix(X)=="FALSE" || is.vector(Y)=="FALSE" ){ stop(message="The data type of X is only Matrix and Y is only a vector") }
  if(alpha <0 ){ stop(message="Robust tuning parameters alpha is only positive.") }
  
  nlambda <- as.integer(nlambda)
  if(nlambda < 1 ){stop(message="The number of grids for sparse tuning parameter is more than 1.")}
  if(lmin < 0 || lmax <0){ stop(message="lmin and lmax are positive.") }
  
  N = nrow(X); k = ncol(X)
  init.time <- proc.time()[3]
  
  ########################
  # regularization grids #
  ########################
  lmode <- match.arg(lambda.mode)
  
  if(lmode=="lambda0"){
    
    LS_LASSO<-cv.glmnet(X, Y, family = "binomial", intercept=F)
    lambdamax0<-LS_LASSO$lambda.min 
    lambdamax1<-try(optim.lam.logit(X,Y,beta.init, weight.mode,LR.loss,lambdamax0, alpha)) #Improvement
    
    if(class(lambdamax1)=='try-error'){lambdamax1<-lambdamax0; warning('Using approximate lambdamax')}
    
    lambda<-seq(1e-4,lambdamax1,length.out=nlambda)
    lambda<-lambda[2:nlambda]
    
  }else{
    lambda <- exp(seq( log(lmin)  ,log(lmax)  ,length=nlambda))
  }
  
  ###################
  # caluculate HGIC #
  ###################
  
  Cn = log(log(N))
  mult = Cn*log(k-1)
  # mult = (log(N)/N+0.5*log(k-1)/N) for EBIC criterion
  
  res.reggam <- list()
  HGIC.vec = rep(NA, length(lambda))
  
  for(ind in 1:length(lambda)){
    res <- logistic_estimation(X,Y,beta.init,weight.mode,LR.loss, alpha,lambda[ind])
    if(sum(res == "NaN") ==0){
      res.reggam[[ind]] <- list(beta=res)
      HGIC.vec[ind] = -2*loglikelihood(X,Y,res) +sum(res[2:k] != 0)*mult
    }
  }
  
  # for large gamma, some values in HBIC.vec are NaN so control for that
  best=NA # if all are NaN then no best model
  which.good = which(!is.na(HGIC.vec))
  if(length(which.good)>0){
    best_HGIC = min(HGIC.vec[which.good])
    len = length(which(HGIC.vec == best_HGIC))
    best = res.reggam[[which(HGIC.vec == best_HGIC)[len]]] # get good entries and best model among them
  }
  return(list(init.time=proc.time()[3]-init.time, lambda=lambda, error=HGIC.vec, 
              lambda.min = lambda[[which(HGIC.vec == best_HGIC)[len]]],
              fits=res.reggam, best=best))
  
}


optim.lam.logit<-function(X,Y,beta.init, weight.mode, LR.loss,lambdamax, alpha){
  #Finds (approximately) smallest lambda such that the slope estimate is zero
  #INPUT
  #X, y: data set
  #lambdamax: initial guess for lambda
  #alpha : parameter of DPD loss
  #weight.mode : lasso, adaptive or scad
  #OUTPUT
  #lambda: smallest lambda such that the slope estimate is zero
  
  N<-nrow(X)
  k<-ncol(X)
  lambda2<-lambdamax
  
  beta.n <- logistic_estimation(X,Y, beta.init,weight.mode,LR.loss, alpha,lambda2)[2:(k)]
  zeros.n<-sum(beta.n==0)
  
  while (zeros.n<(k-1) & lambda2<max(N,1e2)){
    lambda2<-2*lambda2
    beta.n <-logistic_estimation(X,Y,beta.init,weight.mode, LR.loss,alpha,lambda2)[2:(k)]
    zeros.n<-sum(abs(beta.n)==0)
  }
  lambda1<-lambdamax/2
  beta.o<-logistic_estimation(X,Y,beta.init,weight.mode,LR.loss,alpha,lambda1)[2:(k)]
  zeros.o<-sum(abs(beta.o)==0)
  while(zeros.o==(k-1)){
    lambda1<-lambda1/2
    beta.o<-logistic_estimation(X,Y,beta.init,weight.mode,LR.loss,alpha,lambda1)[2:(k)]
    zeros.o<-sum(abs(beta.o)==0)
  }
  for (j in 1:5){
    lambda<-0.5*(lambda1+lambda2)
    beta.n<-logistic_estimation(X,Y,beta.init,weight.mode,LR.loss,alpha,lambda)[2:(k)]
    zeros.n<-sum(beta.n==0)
    if (zeros.n<(k-1)){
      lambda1<-lambda
    }else{
      lambda2<-lambda
    }
  }
  return(lambda2)
}


