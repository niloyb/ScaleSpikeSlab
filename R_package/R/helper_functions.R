###### Helper functions to generate synthetic data and choose hyperparameters ######
# Generate synthetic dataset

#' synthetic_data
#' @description Generates synthetic linear and logistic regression data
#' @param n number of observations
#' @param p number of covariates
#' @param s0 sparsity (number of non-zero components of the true signal)
#' @param error_std Standard deviation of the Gaussian noise (linear regression only)
#' @param type dataset type ('linear' or 'logistic')
#' @param scale design matrix X has columns mean zero and standard deviation 1 (TRUE or FALSE)
#' @param signal non-zero components of the true signal ('constant' or 'deacy')
#' @return Design matrix, response and true signal vector for linear and logistic regression
#' @export
synthetic_data <- function(n,p,s0,error_std,type='linear',scale=TRUE,signal='constant'){
  true_beta <- matrix(0,p,1)
  s0 <- min(p,s0)
  if(s0>0){
    if(signal=='constant'){true_beta[1:s0] <- 2}
    if(signal=='decay'){true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))}
  }
  
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  if(scale){X <- matrix(scale(X),n,p)}
  X_truebeta <- X%*%true_beta
  
  if(type=='linear'){
    # Error terms 
    error_terms = error_std*rnorm(n, mean = 0, sd = 1)
    y = X_truebeta + error_terms
  } else if(type=='logistic'){
    true_aug_y = rlogis(n, location = X_truebeta)
    y <- ifelse(true_aug_y>0,1,0) # Logistic response
  } else if(type=='probit'){
    true_aug_y = rnorm(n, mean = X_truebeta)
    y <- ifelse(true_aug_y>0,1,0) # Probit response
  }
  return(list(X=X, y=y, true_beta=true_beta))
}

#' spike_slab_params
#' @description Generates hyperparameters for spike-and-slab
#' @param n number of observations
#' @param p number of covariates
#' @return spike-and-slab hyperparameters q, tau0, tau1, a0, b0
#' @export
spike_slab_params <- function(n, p){
  # Choice of q, tau0, tau1: following skinny gibbs paper
  K <- max(10,log(n))
  q_seq <- seq(0.0001,(1-0.0001),0.0001)
  probs <- abs(pbinom(K,p,q_seq)-0.9)
  q_index <- which(probs==min(probs))
  if(length(q_index)>1){
    q <- 1/p
  } else {
    q <- q_seq[q_index]
  }
  tau0 <- 1/sqrt(n)
  tau1 <- 1
  a0 <- 1
  b0 <- 1
  return(list(q=q, tau0=tau0, tau1=tau1, a0=a0, b0=b0))
}
