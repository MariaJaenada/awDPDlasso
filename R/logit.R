#main function using IRLS algorithm for logistic regression with DPD loss
# X : n x k matrix containing observations
# Y : n binary vector containing response ; only 0 and 1 are allowed
# beta.init: initial estimator
# weight.mode: adaptive weight mode (lasso, adaptive, scad)
# LR. loss: loss used to solve the auxiliar least squares problem
# alpha : DPD tunning parameter controlling trade-off between efficiency and robustness
# lambda: penalty parameter controlling the shrinkage

#invlogit function
invlogit <- function(x){return(exp(x)/(1+exp(x)))}

# DPD loss function in LR 
LR_DPD_loss <- function(X,Y, beta, alpha){
  n <-dim(X)[1]
  pi <- sapply(X%*%beta,invlogit) #dim x : n x p
  res <- (pi^(alpha+1)+(1-pi)^(alpha+1)) - (1+1/alpha)*(Y*pi^alpha+(1-Y)*(1-pi)^alpha)
  return((1/n^alpha)*mean(res))
}

#loglikelihood loss
loglikelihood <- function(X,Y,beta){
  pi <- sapply(X%*%beta,invlogit) #dim x : n x p
  nonzero <- (pi != 1 & pi !=0)
  Y<- Y[nonzero]
  pi <-pi[nonzero]
  ll <-Y*log(pi)+(1-Y)*log(1-pi)
  return(sum(ll))
}

weight.function <- function(beta, weight.mode = c("lasso", "adaptive", "scad"),lambda, N = 100){
  
  k=length(beta)
  
  # calculate weight
  w = rep(0,k)
  weight <- match.arg(weight.mode)
  
  if(weight =="lasso"){ 
    w = w +1/N
  }else if(weight == "adaptive"){
    
    for(i in 1:k){if(beta[i] != 0){w[i] = 1/abs(beta[i])} } #calculate weights of adaptitive LASSO
    w[w == 0] = 10*max(w) #for zero coefficients the greatest penalization
    w = w/max(w)
    
  }else if(weight == "scad"){
    
    a=3.7 #empirical
    for(i in 1:k){
      if(abs(beta[i]) <= lambda){
        w[i] = 1
      }else if(lambda < abs(beta[i]) & abs(beta[i]) <= a*lambda){
        w[i] = (a*lambda - abs(beta[i]))/((a-1)*lambda)
      }
    }
    
  }else{ stop("you need to specify a valid weight mode") }
  
  return(w)
}

LR_penalized_DPD <- function(X,Y, beta, beta.init, alpha,lambda, weight.mode){
  n <- dim(X)[1]; k<-dim(X)[2]  ### X contains 1 columns to intercept
  w <- weight.function(beta.init[2:k], weight.mode,lambda,n)
  return(LR_DPD_loss(X,Y, beta, alpha)+ lambda*sum(abs(w*beta[2:k])))
}


# IRLS iterate
IRLS_direction <- function(X,Y,  beta.init, weight.mode = c("lasso", "adaptive", "scad"),LR.loss=c("LS","DPD"), alpha, lambda){
  
  ## define regression parameters and derivatives ##
  n<- dim(X)[1] ; k <- dim(X)[2]
  
  pi <- sapply(drop(X%*%beta.init), invlogit)
  pi[pi<1e-4]=0; pi[pi> (1-1e-4)]=1 #avoiding too large coefficients
  weight <- match.arg(weight.mode)
  LR <- match.arg(LR.loss)
  
  # gradient and hessian
  W <- ((1+alpha)/n)*(pi^(alpha+1)*(1-pi)+(1-pi)^alpha*pi^2 -Y*(pi^alpha*(1-pi)+pi*(1-pi)^alpha))
  Lambda <- ((1+alpha)/n)*(pi^(alpha+1)*(1-pi)*((1+alpha)-(2+alpha)*pi) +
               (1-pi)^alpha*pi^2*(2-(2+alpha)*pi) -Y*(pi^alpha*(1-pi)*(alpha-(alpha+1)*pi) +
                pi*(1-pi)^alpha*(1-(alpha+1)*pi)))
  
  Lambda <- sapply(Lambda, function(x){max(0.0001,x)}) #avoiding zero division 

  if(sum(Lambda == 0.0001)>ceiling(3/4*n)){
    
    warning(message="Not convex objective function or zero hessian matrix")
    # If the Hessian matrix is not definite positive, we cant apply Newthon method 
    beta = beta.init
    direction = beta.init
    loss = LR_penalized_DPD(X,Y, beta,beta.init, alpha,lambda, weight)
    
  }else{
  
  ## weighted matrices for IRLS algorithm ##
  Lambda_sq <- diag(sqrt(Lambda))
  #We do not use the intercept to redefine the problem bc the regressor beta0 shouldnt be penalized
  X_w <- Lambda_sq%*%X
  Y_w <- diag(1/Lambda)%*%W - X%*%beta.init 
  Y_w <- Lambda_sq%*%Y_w
  
  if(LR == "LS"){
    ## use Least Squares loss for the linear regression problem ##
    w = weight.function(beta.init[2:k], weight ,lambda, n)
    
    if(sum(w)==0){ w = w+1/n; warning("penalty factor equal zero, changed to LASSO penalty")}
    
    # #solve the least square problem
    # estimate <- ncvreg(X_w, Y_w, family="gaussian", penalty="lasso", lambda = lambda, penalty.factor = c(min(w),w), intercept = F)
    # direction = -drop(estimate$beta)[2:(k+1)]
    
    z <- glmnet(X_w, Y_w, standardize = T, alpha = 1, intercept = F,
                family = "gaussian", lambda = lambda,
                penalty.factor = c(min(w),w))
    direction = - drop(coef(z))[2:(k+1)]
    
  }else if(LR == "DPD"){
    ## use DPD loss for the linear regression problem ##
    
    # dpd_adaptLASSO algorithm needs standarized matrices
    X_ww = scale(X_w,scale=T)
    Y_ww = scale(as.matrix(Y_w),scale=T)
  
    # special function to LR (incertecpt is a regressor but not penalized)
    z <- dpd_estimator(X = X_ww, Y = Y_ww, beta = beta.init[2:k], beta0 = beta.init[1],
                       sigma = 1, penalty = "adaptive", weight.mode = weight,
                       lambda = lambda, gamma = alpha, intercept = T)
    if(z$beta == "NaN"){
      z <- glmnet(X_w, Y_w, standardize = T, alpha = 1, intercept = F,
                family = "gaussian", lambda = lambda)
      direction = -drop(coef(z))[2:(k+1)]
      warning("DPD algorithm fails")
    }else{
      direction = - z$beta*attributes(Y_ww)["scaled:scale"][[1]]/attributes(X_ww)["scaled:scale"][[1]] #unstandarize
    }
    
    }else{ stop("Loss function for the Linear regression problem missespecified")}
  

  ## grid search for the update ##
  
  search<- function(t){return(LR_penalized_DPD(X,Y, t*beta.init +(1-t)*direction, beta.init,alpha,lambda, weight))}
  t0 <- optimize(f=search,interval=c(0,1))$minimum #grid search for the step size
  # evaluate at extremes
  if(search(0)<search(t0)){t0=0}
  if(search(1)<search(t0)){t0=1}
  #update
  beta <- t0*beta.init+(1-t0)*direction
  loss = search(t0)
  
  }
  
  return(list("beta"=beta, "direction" = direction,"loss" = loss))
}


## logistic estimation using IRLS 
logistic_estimation <- function(X,Y,beta, weight.mode, LR.loss, alpha,lambda){
  
  N = dim(X)[1] ;  k = dim(X)[2]
  
  for(m in 2:500){
    
    beta_tmp <- beta
    beta <- IRLS_direction(X,Y,beta,weight.mode,LR.loss, alpha,lambda)$beta
   
    #stopping criteria : differents criteria tested
    stop = max( abs(beta-beta_tmp) )
    #stop = (LR_penalized_DPD(X,Y, beta, beta_tmp, alpha,lambda, weight.mode) -LR_penalized_DPD(X,Y, beta_tmp,beta_tmp, alpha,lambda, weight.mode))
    
    if(all(beta[2:k]==0) | abs(stop) <1e-4 ){break}
    
  }

  return(beta)
}

