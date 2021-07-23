###### estimate function ######
#gamma: parameter DPD function
#lambda : penalty function parameter 
#weigh.mode = Ad-DPD or AW-DPD
library(ncvreg)
library(robustHD)
library(glmnet)

dpd_estimator <- function(X, Y, beta, beta0, sigma, penalty = c("lasso", "adaptive", "SCAD", "MCP"), weight.mode = c("lasso","adaptive","scad"), lambda, gamma, intercept){
  
  #if every coeff on init beta is zero, stop
  if(all(beta==0)){
    stop("null beta init- non penalized DPD")
  }
  
  N = dim(X)[1]
  p = dim(X)[2]
  
  #create intercept term
  tmp = rep(1, N)
  tmp1 = intercept*beta0*tmp
  tmp2 = X%*%beta #X matrix doesnt contain ones column 
  tmp3 = drop(tmp1 + tmp2) #y estimate
  
  converge = T
  # iterative minimization 
  for (m in 1:500){
    
    penalty_it = penalty
    if(penalty == "adaptive"){w = weight.function(beta = beta, weight.mode = weight.mode , lambda = lambda, N=N); penalty_it = "lasso" }else{w = rep(1,p)}
    
    #temporary copies 
    beta0_tmp = beta0*intercept
    beta_tmp = beta
    sigma_tmp = sigma
    
    #weights of MM-algorithm
    mu = exp(-(gamma/2)*((Y-tmp3)/sigma)^2)
    mu[abs(mu)<1e-16]= 0
    if(sum(mu) == 0){ beta <- beta0 <- sigma <- "NaN"; converge = F; break}
    mu = drop(mu/sum(mu))
    
    #update beta0: doesn't depend on penalization
    beta0 = drop(t(mu)%*%(Y - tmp2)*intercept)
    tmp1 = beta0*tmp
    
    #weigthed matrices : reduce to least squares problem
    Y_w = diag(sqrt(mu))%*%((Y-tmp1)/sigma) #Y contains incercept (mean)
    X_w =  diag(sqrt(mu))%*%(X/sigma) #without interceptcept
    
    #minimize beta with adaptative lasso
    #estimate = ncvfit(X_w, Y_w, penalty=penalty, lambda = lambda, penalty.factor = w)
    estimate = ncvreg(X_w, Y_w, family="gaussian", penalty=penalty_it, lambda = lambda, penalty.factor = w)
    #estimate <- glmnet(X_w, Y_w, intercept = FALSE, standardize = TRUE, alpha = 1, family = "gaussian", lambda = lambda, penalty.factor = w)
   
    beta = estimate$beta[1:p+1]
    #beta = coef(estimate)[1:p+1]
    #update tmp3 with new beta
    tmp2 = X%*%beta 
    tmp3 = drop(tmp1 + tmp2) 
    
    #update sigma now
    wi = exp((-gamma/2)*((Y-tmp3)/sigma)^2)
    sigma = sqrt(abs(drop(mean(wi*(Y-tmp3)^2))*(mean(wi)-gamma/(gamma+1)^(3/2))^(-1)))
    
    if(sigma < 1e-5){beta <- beta0 <- sigma <- "NaN"; converge = F; break}
    
    ##stopping criteria LASSO penalization 
    # stop2 = abs(sigma^(-gamma)*(1/(1+gamma)^(3/2)-(1/(N*gamma)*sum(exp(-(gamma/2)*((Y-tmp3)/sigma)^2)) ))

    #                          + N*lambda*sum(abs(w*beta_tmp))- N*lambda*sum(abs(w*beta)) )
    
    #sparsity condition
    if(sum(beta!=0) > n/2){converge = FALSE; break}
    # if(penalty == "adaptive"){
      stop2 = abs(sigma^(-gamma)*(1/(1+gamma)^(3/2)-(1/(N*gamma)*sum(exp(-(gamma/2)*((Y-tmp3)/sigma)^2)) ))
                  - (sigma_tmp^(-gamma)*(1/(1+gamma)^(3/2)-(1/(N*gamma)*sum(exp(-(gamma/2)*((Y-beta0_tmp*tmp-X%*%beta_tmp)/sigma_tmp)^2)) )))
                  + N*penalty.function(beta_tmp, penalty, w, lambda) - N*penalty.function(beta, penalty, w, lambda) )
    # }else{
    # stop2 <- abs( (-gamma/(1+gamma))*log(sigma_tmp) - (-gamma/(1+gamma))*log(sigma) 
    #       + (-1/gamma)*log(sum(exp(-gamma*(Y-beta0_tmp*tmp-X%*%beta_tmp)^2/(2*sigma_tmp^2))*(2*pi*sigma_tmp^2)^(-gamma/2)))
    #       -(-1/gamma)*log(sum( exp(-gamma*(Y-tmp3)^2/(2*sigma^2))*(2*pi*sigma^2)^(-gamma/2))) 
    #       + N*penalty.function(beta_tmp, penalty, w, lambda) - N*penalty.function(beta, penalty, w, lambda) )
    
    #}
    
    
    if(all(beta==0) | (stop2 < 1e-9) ){ break }
  }
  if(m==500){converge = F}
  return(list("beta0"=beta0,"beta" = beta, "sigma"= sigma, "converge" = converge))
}


weight.function <- function(beta, weight.mode = c("lasso", "adaptive", "scad"),lambda, N = 100){
  
  p=length(beta)
  
  # #calculate weight
  w = rep(0,p)
  weight <- match.arg(weight.mode)
  
  if(weight =="lasso"){ 
    w = w +1 #/N
  }else if(weight == "adaptive"){
    
    for(i in 1:p){if(beta[i] != 0){w[i] = 1/abs(beta[i])} } #calculate weights of adaptitive LASSO
    w[w == 0] = 4*max(w) #for zero coefficients the greatest penalization
    
  }else if(weight == "scad"){
    
    a=3.7 #empirical
    for(i in 1:p){
      
      if(abs(beta[i]) <= lambda){
        w[i] = 1
      }else if(lambda < abs(beta[i]) & abs(beta[i]) <= a*lambda){
        w[i] = (a*lambda - abs(beta[i]))/((a-1)*lambda)
      }
    }
    
  }else{ stop("you need to specify a valid weight mode") }

  return(w)
}

penalty.function <- function(beta, penalty, w, lambda){
  
  if(penalty == "SCAD"){ pen = sum(scad(beta,lambda))
  }else if(penalty == "MCP"){ pen = sum(mcp(beta, lambda))
  }else{ pen = lambda*sum(abs(w*beta)) }
  
  return(pen)
}
