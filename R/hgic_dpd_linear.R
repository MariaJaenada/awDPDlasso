awDPDlasso.lm <- function(X,Y, init.mode=c("RLARS","DPD-lasso"), init.model=NULL,
                               penalty = "lasso", weight.mode = "lasso",
                               lambda.mode=c("lambda0","given"),
                               lmax=1,lmin=0.05, nlambda=50,
                               seed = 1024, intercept=T, alpha=1){
  
  
  #########################################
  # detect the data type and tuning parameters value
  #########################################
  if(is.matrix(X)=="FALSE" || is.matrix(Y)=="FALSE" ){
    stop(message="The data type of X and Y is only Matrix.")
  }
  
  nlambda <- as.integer(nlambda)
  if(nlambda < 1 ){
    stop(message="The number of grids for sparse tuning parameter is more than 1.")
  }
 
  if(lmin < 0 || lmax <0){
    stop(message="lmin and lmax are positive.")
  }
  
  #########################################
  # choose a initial point
  #########################################
  init <- match.arg(init.mode)
  
  if(init == "RLARS"){
    RLARS = init.model
    init.time = NULL
    if(is.null(RLARS)){
      init.time = system.time(RLARS <- rlars(Y~X, seed=seed))}
    
    beta0.init <- intercept*(coef(RLARS)[1])
    beta.init <- as.matrix(coef(RLARS)[-1])
    sigma.init <- getScale(RLARS)
    
    weight.init = weight.mode
  }else if(init == "DPD-lasso"){
    
    RLARS = init.model
    if(is.null(RLARS)){init.time = system.time(RLARS <- rlars(Y~X,seed = seed))}
    beta0.init.lasso <- intercept*(coef(RLARS)[1])
    beta.init.lasso <- as.matrix(coef(RLARS)[-1])
    sigma.init.lasso <- getScale(RLARS)
    weight.init = "lasso"
    
   }else{
    X1 = cbind(1,X)
    b = ginv(crossprod(X1)) %*% t(X1) %*% Y
    beta0.init = intercept*b[1]
    beta.init = as.matrix(b[-1])
    sigma.init = 1.4826*mad(Y - X1 %*% b)
  }
  
  ########################################
  # regularization grids
  ########################################
  lmode <- match.arg(lambda.mode)
  
  if(lmode=="lambda0"){
    lambdamax0<-lambda0(X,Y) #Initial candidate
    lambdamax1<-try(optim.lam(X,Y,lambdamax0, 
                              beta0 = intercept*(coef(RLARS)[1]), beta= as.matrix(coef(RLARS)[-1]),sigma = getScale(RLARS),
                              gam, penalty, weight.mode, intercept))#Improvement
    if(class(lambdamax1)=='try-error'){
      lambdamax1<-lambdamax0
      warning('Using approximate lambdamax')
    }
    lambda<-seq(1e-8,lambdamax1,length.out=nlambda) #return sequence which annull all coeff when init is RLARS
    if(dim(X)[2]>=dim(X)[1]){lambda<-lambda[2:nlambda]}
    
  }else{
    lambda <- exp(seq( log(lmin)  ,log(lmax)  ,length=nlambda))
  }
  
  #######################################
  # caluculate HBIC and best beta
  #######################################
  n = nrow(X); p = ncol(X)
  Cn = log(log(n))
  mult = Cn*log(p)/n
  
  # calculate solutions
  res.reggam <- list()
  HBIC.vec = rep(NA, length(lambda))
  
  for(k in 1:length(lambda)){
    
    #initialize with the same lambda
    if(init == "DPD-lasso"){
      
      DPD <- dpd_estimator(X, as.matrix(Y), beta.init.lasso, beta0.init.lasso, sigma.init.lasso, 
                           penalty = "lasso", weight.mode = "lasso", lambda = lambda[k], gamma = gam, intercept = intercept)
  
      beta.init <- DPD$beta
      beta0.init <- DPD$beta0
      sigma.init <- DPD$sigma
      
      if(beta0.init != "NaN"){
        
        if(sum(beta.init) !=0 ){ res <- dpd_estimator(X, as.matrix(Y), beta.init, beta0.init, sigma.init, penalty = penalty, weight.mode = weight.mode, lambda = lambda[k], gamma = gam, intercept = intercept)
        
        }else{res<- list(beta0 = beta0.init, beta = beta.init, sigma = sigma.init, converge = F)}
        
      }else{res <- list(beta0 = "NaN", beta = "NaN", sigma = "NaN", converge = F)}
      
    }else{ 
      #beta init RLARS
      res <- dpd_estimator(X, as.matrix(Y), beta.init, beta0.init, sigma.init, penalty = penalty, weight.mode = weight.mode, lambda = lambda[k], gamma = gam, intercept = intercept)
    }
    if( !is.na(res$beta0) & res$converge){
      res.reggam[[k]] <- list(beta0=res$beta0, beta=res$beta, sigma=res$sigma)
      #print(res$beta[res$beta!=0])
      M = sum(res$beta != 0)
      SSE = res$sigma^2 # robust measure of variance
      #SSE2 = sd(Y-cbind(1,X)%*%c(res$beta0,res$beta))
      HBIC.vec[k] = log(SSE) + M*mult
      if(all(res$beta==0)){break}
      
    }
  }
  
  
  # for large gamma, some values in HBIC.vec are NaN so control for that
  best=NA # if all are NaN then no best model
  which.good = which(!is.na(HBIC.vec))
  if(length(which.good)>0){
    
    best_HBIC = min(HBIC.vec[which.good])
    best = res.reggam[[which(HBIC.vec == best_HBIC)[1]]] # get good entries and best model among them
  }
  print(lambda[[which(HBIC.vec == best_HBIC)[1]]])
  
  return(list(init.time=init.time[3], lambda=lambda, error=HBIC.vec,
              fits=res.reggam, best=best))
  
}


optim.lam<-function(X,Y,beta0=NULL, beta=NULL, sigma=NULL, lambdamax, gam, penalty, weight.mode,intercept, seed){
  #Finds (approximately) smallest lambda such that the slope estimate is zero
  #INPUT
  #X, y: data set
  #lambdamax: initial guess for lambda
  #gam : parameter of DPD loss
  #weight.mode : scad or adapt
  #OUTPUT
  #lambda: smallest lambda such that the slope estimate is zero
  
  n<-nrow(X)
  p<-ncol(X)

  lambda2<-lambdamax
  
  if(is.null(beta0)){
    RLARS <- rlars(Y~X)
    beta0 <- intercept*(coef(RLARS)[1])
    beta <- as.matrix(coef(RLARS)[-1])
    sigma <- getScale(RLARS)
  }
  
  beta.n<- dpd_estimator(X, as.matrix(Y), beta, beta0, sigma, penalty = penalty, weight.mode = weight.mode, lambda = lambda2, gamma = gam, intercept = intercept)$beta
  zeros.n<-sum(beta.n==0)
  
  while (zeros.n<p & lambda2<min(n,10)){
    lambda2<-2*lambda2
    beta.n<-dpd_estimator(X, as.matrix(Y), beta, beta0, sigma, penalty = penalty, weight.mode = weight.mode, lambda = lambda2, gamma = gam, intercept = intercept)$beta
    zeros.n<-sum(beta.n==0)
  }
  
  if(lambda2>=min(n,10)){return(lambda2)}
  
  lambda1<-lambdamax/2
  beta.o<- dpd_estimator(X, as.matrix(Y),  beta, beta0, sigma, penalty = penalty, weight.mode = weight.mode, lambda = lambda1, gamma = gam, intercept = intercept)$beta
  zeros.o<-sum(beta.o==0)
  
  while(zeros.o==p){
    lambda1<-lambda1/2
    beta.o<-dpd_estimator(X, as.matrix(Y),  beta, beta0, sigma, penalty = penalty, weight.mode = weight.mode, lambda = lambda1, gamma = gam, intercept = intercept)$beta
    zeros.o<-sum(beta.o==0)
  }
  
  for (j in 1:5){
    lambda<-0.5*(lambda1+lambda2)
    beta.n<-dpd_estimator(X, as.matrix(Y),  beta, beta0, sigma, penalty = penalty, weight.mode = weight.mode, lambda = lambda, gamma = gam, intercept = intercept)$beta
    zeros.n<-sum(beta.n==0)
    if (zeros.n<p){
      lambda1<-lambda
    }else{
      lambda2<-lambda
    }
  }
  return(lambda2)
}


