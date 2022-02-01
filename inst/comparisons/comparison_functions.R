# Functions for comparison with alternatives
require(Rcpp)
require(RcppEigen)
# Rcpp functions
Rcpp::cppFunction('
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y){
Eigen::VectorXd output;
output = X*y;
return output;
}', depends = 'RcppEigen')
# Rcpp::cppFunction('
# Eigen::MatrixXd cpp_vec_mat_prod(const Eigen::VectorXd y, const Eigen::MatrixXd X){
# Eigen::VectorXd output;
# output = y.adjoint()*X;
# return output;
# }', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
  return Eigen::MatrixXd(X*Y);
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X){
const int n = X.cols();
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}', depends = 'RcppEigen')


## Matrix calculation helper functions ##
matrix_full <- function(Xt, d){
  n <- dim(Xt)[2]
  # NOTE: (1) For MacOS with veclib BLAS, crossprod is fast via multit-hreading
  # return(crossprod(Xt*c(1/d)^0.5) + diag(n))
  return(fcprd(Xt*c(1/d)^0.5) + diag(n))
}
matrix_precompute <- function(Xt, d, prev_d, prev_matrix){
  if(sum(d!=prev_d)==0){return(prev_matrix)}
  
  pos_indices <- (1/d-1/prev_d)>0
  if(sum(pos_indices)==0){
    pos_part <- 0
  } else{
    pos_d_update <- (1/d-1/prev_d)[pos_indices]
    pos_Xt_update <- Xt[pos_indices,,drop=FALSE]
    pos_Xt_d_update <- pos_Xt_update*(pos_d_update^0.5)
    pos_part <- fcprd(pos_Xt_d_update)
  }
  
  neg_indices <- (1/d-1/prev_d)<0
  if(sum(neg_indices)==0){
    neg_part <- 0
  } else{
    neg_d_update <- -(1/d-1/prev_d)[neg_indices]
    neg_Xt_update <- Xt[neg_indices,,drop=FALSE]
    neg_Xt_d_update <- neg_Xt_update*(neg_d_update^0.5)
    neg_part <- fcprd(neg_Xt_d_update)
  }
  return(pos_part-neg_part+prev_matrix)
}
inverse_precompute <- function(Xt, d, prev_d, prev_inverse){
  swap_indices <- (d!=prev_d)
  no_swaps <- sum(swap_indices)
  if(no_swaps==0){return(prev_inverse)}
  
  diag_diff <- (1/d-1/prev_d)[swap_indices]
  Xt_update <- Xt[swap_indices,,drop=FALSE]
  Xt_update_prev_inverse <- cpp_prod(Xt_update,prev_inverse)
  M_matrix <- diag(1/diag_diff, no_swaps)+cpp_prod(Xt_update_prev_inverse,t(Xt_update))
  # Use (M_matrix+t(M_matrix))/2 to avoid numerical error, as M_matrix symmetric
  M_matrix_inverse <- solve((M_matrix+t(M_matrix))/2)
  inverse_update <- cpp_prod(t(Xt_update_prev_inverse),
                             cpp_prod(M_matrix_inverse,Xt_update_prev_inverse))
  return(prev_inverse-inverse_update)
}

## Spike slab linear regression functions ##
# beta update
# beta update
update_beta_comparison <- 
  function(z, sigma2, X, Xt, y, prev_z, prev_matrix, 
           prev_inverse, tau0, tau1, algo, u=NULL, delta=NULL){
  p <- length(z)
  n <- length(y)
  no_swaps <- sum(z!=prev_z)
  
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
  
  if(is.null(u)){u <- rnorm(p, 0, 1)}
  u = u/(d^0.5)
  if(is.null(delta)){delta <- c(rnorm(n,0,1))}
  # v = cpp_mat_vec_prod(X,u) + delta
  v = X%*%u + delta
  
  if(algo=='ScalableSpikeSlab'){
    M_matrix <- matrix_precompute(Xt, d, prev_d, prev_matrix)
    
    if(no_swaps<n){
      M_matrix_inverse <- inverse_precompute(Xt, d, prev_d, prev_inverse)
      weighted_cprd <- NA
    } else {
      M_matrix_inverse <- chol2inv(chol(M_matrix))
    }
  } else if(algo=='Sota'){
    M_matrix <- matrix_full(Xt, d)
    M_matrix_inverse <- chol2inv(chol(M_matrix))
  } else{
    stop("algo parameter must be ScalableSpikeSlab or Sota")
  }
  
  # v_star <- cpp_mat_vec_prod(M_matrix_inverse,(y/sqrt(sigma2) - v))
  # beta <- sqrt(sigma2)*(u + (d^(-1))*cpp_mat_vec_prod(Xt,v_star))
  v_star <- M_matrix_inverse%*%(y/sqrt(sigma2) - v)
  beta <- sqrt(sigma2)*(u + (d^(-1))*(Xt%*%v_star))
  return(list('beta'=beta, 'matrix'=M_matrix, 
              'matrix_inverse'=M_matrix_inverse))
}
# z update
update_z <- function(beta, sigma2, tau0, tau1, q, u_crn=NULL){
  p <- length(beta)
  if(is.null(u_crn)){u_crn <- runif(p)}
  log_prob1 <- log(q)+dnorm(beta, sd=(tau1*sqrt(sigma2)),log = TRUE)
  log_prob2 <- log(1-q)+dnorm(beta, sd=(tau0*sqrt(sigma2)),log = TRUE)
  probs <- 1/(1+exp(log_prob2-log_prob1))
  z <- ifelse(u_crn<probs,1,0)
  return(z)
}
# sigma2 update
update_sigma2 <- function(beta, z, tau0, tau1, a0, b0, X, y, u_crn=NULL){
  if(is.null(u_crn)){u_crn <- runif(1)}
  n <- length(y)
  p <- length(beta)
  # rss <- sum((y-cpp_mat_vec_prod(X,beta))^2)
  rss <- sum((y-X%*%beta)^2)
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  beta2_d <- sum((beta)^2*d)
  sigma2 <- 1/qgamma(u_crn,shape = (a0+n+p)/2, rate = (b0+rss+beta2_d)/2)
  return(sigma2)
}

spike_slab_linear_kernel <- 
  function(beta, z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
           tau0, tau1, q, a0, b0, algo, random_samples){
    beta_output <- 
      update_beta_comparison(z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                  tau0, tau1, algo, u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples$z_u)
    sigma2_new <- update_sigma2(beta_new, z_new, tau0, tau1, a0, b0, X, y, u_crn=random_samples$sigma2_u)
    return(list('beta'=beta_new,'z'=z_new,'sigma2'=sigma2_new,
                'prev_z'=z, 'prev_matrix'=mat, 'prev_inverse'=mat_inverse))
  }

#' comparison_spike_slab_linear
#' @description Generates Markov chain targeting the posterior corresponding to
#' Bayesian linear regression with spike and slab priors
#' @param chain_length Markov chain length
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @param y Response
#' @param tau0 prior hyperparameter (non-negative real)
#' @param tau1 prior hyperparameter (non-negative real)
#' @param q prior hyperparameter (strictly between 0 and 1)
#' @param a0 prior hyperparameter (non-negative real)
#' @param b0 prior hyperparameter (non-negative real)
#' @return Output from Markov chain targeting the posterior corresponding to
#' Bayesian linear regression with spike and slab priors
#' @export
comparison_spike_slab_linear <- 
  function(chain_length,X,Xt,y,tau0,tau1,q,a0=1,b0=1,algo,rinit=NULL,
           verbose=FALSE,burnin=0,store=TRUE){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      sigma2 <- 1/rgamma(1,shape = (a0/2), rate = (b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0*sqrt(sigma2))
      beta[z==1] <- beta[z==1]*(tau1*sqrt(sigma2))
      
      prev_z <- z
      prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      prev_matrix <- matrix_full(Xt, prev_d)
      prev_inverse <- chol2inv(chol(prev_matrix))
    }
    
    if(store){
      beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      sigma2_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = 1)
    } else {
      beta_ergodic_sum <- rep(0,p)
      z_ergodic_sum <- rep(0,p)
      sigma2_ergodic_sum <- 0
    }
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), z_u=runif(p), sigma2_u=runif(1))
      new_state <- 
        spike_slab_linear_kernel(beta, z, sigma2, X, Xt, y, 
                                 prev_z, prev_matrix, prev_inverse,
                                 tau0, tau1, q, a0, b0, algo, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      sigma2 <- new_state$sigma2
      prev_z <- new_state$prev_z
      prev_matrix <- new_state$prev_matrix
      prev_inverse <- new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          beta_samples[(t-burnin),] <- beta
          z_samples[(t-burnin),] <- z
          sigma2_samples[(t-burnin),] <- sigma2
        } else{
          beta_ergodic_sum <- beta_ergodic_sum + beta
          z_ergodic_sum <- z_ergodic_sum + z
          sigma2_ergodic_sum <- sigma2_ergodic_sum + sigma2
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('beta'=beta_samples, 'z'=z_samples, 'sigma2'=sigma2_samples))
    } else {
      return(list('beta_ergodic_avg'=beta_ergodic_sum/(chain_length-burnin), 
                  'z_ergodic_avg'=z_ergodic_sum/(chain_length-burnin), 
                  'sigma2_ergodic_avg'=sigma2_ergodic_sum/(chain_length-burnin)))
    }
  }



## Spike slab logistic regression functions ##
update_e <- function(beta, sigma2, y, X, u_crn=NULL){
  n <- length(y)
  if(is.null(u_crn)){u_crn <- runif(n)}
  # means <- cpp_mat_vec_prod(X, beta)
  means <- X%*%beta
  sds <- sqrt(sigma2)
  ubs <- rep(Inf,n)
  lbs <- rep(-Inf,n)
  lbs[y==1] <- 0
  ubs[y!=1] <- 0
  e <- TruncatedNormal::qtnorm(u_crn, mu=means, sd=sds, lb=lbs, ub=ubs)
  return(e)
}
spike_slab_logistic_kernel <- 
  function(beta, z, e, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
           tau0, tau1, q, a0, b0, algo, random_samples){
    beta_output <- 
      update_beta_comparison(z, sigma2, X, Xt, e, prev_z, prev_matrix, prev_inverse,
                  tau0, tau1, algo, u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- update_e(beta_new, sigma2, y, X, u_crn=random_samples$e_u)
    nu <- 7.3
    w2 <- (pi^2)*(nu-2)/(3*nu)
    a0 <- nu
    b0 <- w2*nu
    sigma2_new <- update_sigma2(beta_new, z_new, tau0, tau1, a0, b0, X, e_new, u_crn=random_samples$sigma2_u)
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new,'sigma2'=sigma2_new,
                'prev_z'=z, 'prev_matrix'=mat, 'prev_inverse'=mat_inverse))
  }

#' comparison_spike_slab_logistic
#' @description Generates Markov chain targeting the posterior corresponding to
#' Bayesian logistic regression with spike and slab priors
#' @param chain_length Markov chain length
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @param y Response
#' @param tau0 prior hyperparameter (non-negative real)
#' @param tau1 prior hyperparameter (non-negative real)
#' @param q prior hyperparameter (strictly between 0 and 1)
#' @param a0 prior hyperparameter (non-negative real)
#' @param b0 prior hyperparameter (non-negative real)
#' @return Output from Markov chain targeting the posterior corresponding to
#' Bayesian logistic regression with spike and slab priors
#' @export
comparison_spike_slab_logistic <- 
  function(chain_length,X,Xt,y,tau0,tau1,q,a0=1,b0=1,algo, rinit=NULL,verbose=FALSE,
           burnin=0,store=TRUE){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      sigma2 <- 1/rgamma(1,shape = (a0/2), rate = (b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0*sqrt(sigma2))
      beta[z==1] <- beta[z==1]*(tau1*sqrt(sigma2))
      # e <- cpp_mat_vec_prod(X, beta) + sqrt(sigma2)*rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + sqrt(sigma2)*rnorm(n, mean = 0, sd = 1)
      
      prev_z <- z
      prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      prev_matrix <- matrix_full(Xt, prev_d)
      prev_inverse <- chol2inv(chol(prev_matrix))
    }
    
    beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
    z_samples <- matrix(NA, nrow = chain_length, ncol = p)
    e_samples <- matrix(NA, nrow = chain_length, ncol = n)
    sigma2_samples <- matrix(NA, nrow = chain_length, ncol = 1)
    
    if(store){
      beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
      z_samples <- matrix(NA, nrow = chain_length, ncol = p)
      e_samples <- matrix(NA, nrow = chain_length, ncol = n)
      sigma2_samples <- matrix(NA, nrow = chain_length, ncol = 1)
    } else {
      beta_ergodic_sum <- rep(0,p)
      z_ergodic_sum <- rep(0,p)
      e_ergodic_sum <- rep(0,n)
      sigma2_ergodic_sum <- 0
    }
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), e_u=runif(n), z_u=runif(p), sigma2_u=runif(1))
      new_state <- 
        spike_slab_logistic_kernel(beta, z, e, sigma2, X, Xt, y, 
                                   prev_z, prev_matrix, prev_inverse,
                                   tau0, tau1, q, a0, b0, algo, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      sigma2 <- new_state$sigma2
      prev_z <- new_state$prev_z
      prev_matrix <- new_state$prev_matrix
      prev_inverse <- new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          beta_samples[(t-burnin),] <- beta
          z_samples[(t-burnin),] <- z
          e_samples[(t-burnin),] <- e
          sigma2_samples[(t-burnin),] <- sigma2
        } else{
          beta_ergodic_sum <- beta_ergodic_sum + beta
          z_ergodic_sum <- z_ergodic_sum + z
          e_ergodic_sum <- e_ergodic_sum + e
          sigma2_ergodic_sum <- sigma2_ergodic_sum + sigma2
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('beta'=beta_samples, 'z'=z_samples, 'e'=e_samples, 'sigma2'=sigma2_samples))
    } else {
      return(list('beta_ergodic_avg'=beta_ergodic_sum/(chain_length-burnin), 
                  'z_ergodic_avg'=z_ergodic_sum/(chain_length-burnin), 
                  'e_ergodic_avg'=e_ergodic_sum/(chain_length-burnin), 
                  'sigma2_ergodic_avg'=sigma2_ergodic_sum/(chain_length-burnin)))
    }
  }

#' comparison_spike_slab_mcmc
#' @description Generates Markov chain targeting the posterior corresponding to
#' Bayesian linear or logistic regression with spike and slab priors
#' @param chain_length Markov chain length
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @param y Response
#' @param tau0 prior hyperparameter (non-negative real)
#' @param tau1 prior hyperparameter (non-negative real)
#' @param q prior hyperparameter (strictly between 0 and 1)
#' @param a0 prior hyperparameter (non-negative real)
#' @param b0 prior hyperparameter (non-negative real)
#' @return Output from Markov chain targeting the posterior corresponding to
#' Bayesian linear or logistic regression with spike and slab priors
#' @export
comparison_spike_slab_mcmc <- 
  function(chain_length,X,Xt,y,tau0,tau1,q,a0=1,b0=1,algo='ScalableSpikeSlab',
           type=NULL,rinit=NULL,verbose=FALSE,burnin=0,store=TRUE){
    # Binary 0,1 labels for logistic regression
    if(all(y*(1-y)==0)){
      return(comparison_spike_slab_logistic(chain_length,X,Xt,y,tau0,tau1,q,a0,b0,algo,rinit,verbose,burnin,store))
    } else {
      return(comparison_spike_slab_linear(chain_length,X,Xt,y,tau0,tau1,q,a0,b0,algo,rinit,verbose,burnin,store))
    }
  }











