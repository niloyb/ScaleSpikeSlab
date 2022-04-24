# Functions for Scalable Spike-and-Slab

## Matrix calculation helper functions ##
matrix_full <- function(Xt, d){
  n <- dim(Xt)[2]
  # NOTE: (1) For MacOS with veclib BLAS, crossprod is fast via multit-hreading
  # return(crossprod(Xt*c(1/d)^0.5) + diag(n))
  return(fcprd(Xt*c(1/d)^0.5) + diag(n))
}
matrix_precompute <- function(Xt, d, prev_d, prev_matrix, XXt, tau0, tau1){
  p <- dim(Xt)[1]
  n <- dim(Xt)[2]
  no_swaps <- sum(d!=prev_d)
  slab_indices <- d!=1/(tau0)^2
  no_slabs <- sum(slab_indices)
  no_spikes <- p-no_slabs
  
  if(no_swaps==0){return(prev_matrix)}
  if(no_slabs==0){
    tau0_mat <- tau0^2*XXt
    diag(tau0_mat) <- diag(tau0_mat)+1
    return(tau0_mat)
  }
  if(no_spikes==0){
    tau1_mat <- tau1^2*XXt
    diag(tau1_mat) <- diag(tau1_mat)+1
    return(tau1_mat)
  }
  
  # if(no_slabs==0){return(diag(n)+tau0^2*XXt)}
  # if(no_spikes==0){return(diag(n)+tau1^2*XXt)}
  
  if(no_swaps<min(no_slabs, no_spikes)){
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
  } else if (no_slabs<=min(no_swaps, no_spikes)){
    Xt_update <- Xt[slab_indices,,drop=FALSE]
    update_part <- fcprd(Xt_update)*(tau1^2-tau0^2)
    # return(diag(n)+tau0^2*XXt+update_part)
    tau0_mat <- tau0^2*XXt
    diag(tau0_mat) <- diag(tau0_mat)+1
    return(tau0_mat+update_part)
  } else{
    Xt_update <- Xt[!slab_indices,,drop=FALSE]
    update_part <- fcprd(Xt_update)*(tau0^2-tau1^2)
    # return(diag(n)+tau1^2*XXt+update_part)
    tau1_mat <- tau1^2*XXt
    diag(tau1_mat) <- diag(tau1_mat)+1
    return(tau1_mat+update_part)
  }
}

inverse_precompute <- 
  function(Xt, d, prev_d, prev_inverse, tau0_inverse, tau1_inverse, tau0, tau1){
  p <- dim(Xt)[1]
  n <- dim(Xt)[2]
  
  swap_indices <- (d!=prev_d)
  no_swaps <- sum(swap_indices)
  slab_indices <- d!=1/(tau0)^2
  no_slabs <- sum(slab_indices)
  no_spikes <- p-no_slabs
  
  # print(c(no_swaps, no_slabs, no_spikes))
  
  if(no_swaps==0){return(prev_inverse)}
  if(no_slabs==0){return(tau0_inverse)}
  if(no_spikes==0){return(tau1_inverse)}
  
  if(no_swaps<min(no_slabs, no_spikes)){
    diag_diff <- (1/d-1/prev_d)[swap_indices]
    Xt_update <- Xt[swap_indices,,drop=FALSE]
    Xt_update_prev_inverse <- cpp_prod(Xt_update,prev_inverse)
    M_matrix <- diag(1/diag_diff, no_swaps)+cpp_prod(Xt_update_prev_inverse,t(Xt_update))
    # Use (M_matrix+t(M_matrix))/2 to avoid numerical error, as M_matrix symmetric
    M_matrix_inverse <- solve((M_matrix+t(M_matrix))/2)
    inverse_update <- cpp_prod(t(Xt_update_prev_inverse),
                               cpp_prod(M_matrix_inverse,Xt_update_prev_inverse))
    return(prev_inverse-inverse_update)
  } else if (no_slabs<=min(no_swaps, no_spikes)){
    Xt_update <- Xt[slab_indices,,drop=FALSE]
    Xt_update_prev_inverse <- cpp_prod(Xt_update,tau0_inverse)
    M_matrix <- diag(no_slabs)/(tau1^2-tau0^2) + cpp_prod(Xt_update_prev_inverse,t(Xt_update))
    M_matrix_inverse <- solve((M_matrix+t(M_matrix))/2)
    inverse_update <- cpp_prod(t(Xt_update_prev_inverse),
                               cpp_prod(M_matrix_inverse,Xt_update_prev_inverse))
    return(tau0_inverse-inverse_update)
  } else{
    Xt_update <- Xt[!slab_indices,,drop=FALSE]
    Xt_update_prev_inverse <- cpp_prod(Xt_update,tau1_inverse)
    M_matrix <- diag(no_spikes)/(tau0^2-tau1^2) + cpp_prod(Xt_update_prev_inverse,t(Xt_update))
    M_matrix_inverse <- solve((M_matrix+t(M_matrix))/2)
    inverse_update <- cpp_prod(t(Xt_update_prev_inverse),
                               cpp_prod(M_matrix_inverse,Xt_update_prev_inverse))
    return(tau1_inverse-inverse_update)
  }
}

## Spike slab linear regression functions ##
# beta update
update_beta <- function(z, sigma2, X, Xt, y, prev_z, prev_matrix, 
                        prev_inverse, XXt, tau0_inverse, tau1_inverse,
                        tau0, tau1, u=NULL, delta=NULL){
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
  
  M_matrix <- matrix_precompute(Xt, d, prev_d, prev_matrix, XXt, tau0, tau1)
  
  if(no_swaps<n){
    M_matrix_inverse <- 
      inverse_precompute(Xt, d, prev_d, prev_inverse, tau0_inverse, tau1_inverse, tau0, tau1)
  } else {
    M_matrix_inverse <- chol2inv(chol(M_matrix))
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
update_sigma2_linear <- function(beta, z, tau0, tau1, a0, b0, X, y, u_crn=NULL){
  n <- length(y)
  p <- length(beta)
  # rss <- sum((y-cpp_mat_vec_prod(X,beta))^2)
  rss <- sum((y-X%*%beta)^2)
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  beta2_d <- sum((beta)^2*d)
  
  if(is.null(u_crn)){
    sigma2 <- 1/rgamma(1,shape = (a0+n+p)/2, rate = (b0+rss+beta2_d)/2)
  } else{
    # print(c((a0+n+p),(b0+rss+beta2_d)))
    log_sigma2 <- -log(qgamma(log(u_crn),shape = (a0+n+p)/2, rate = (b0+rss+beta2_d)/2, log.p = TRUE))
    sigma2 <- exp(log_sigma2)
    # sigma2 <- 1/qgamma(log(u_crn),shape = (a0+n+p)/2, rate = (b0+rss+beta2_d)/2, log.p = TRUE)
  }
  # if(is.null(u_crn)){u_crn <- runif(1)}
  # sigma2 <- 1/qgamma(u_crn,shape = (a0+n+p)/2, rate = (b0+rss+beta2_d)/2)
  return(sigma2)
}

spike_slab_linear_kernel <- 
  function(beta, z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
           XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, a0, b0, random_samples){
    beta_output <- 
      update_beta(z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                  XXt, tau0_inverse, tau1_inverse, tau0, tau1, 
                  u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples$z_u)
    sigma2_new <- update_sigma2_linear(beta_new, z_new, tau0, tau1, a0, b0, X, y, u_crn=random_samples$sigma2_u)
    return(list('beta'=beta_new,'z'=z_new,'sigma2'=sigma2_new,
                'prev_z'=z, 'prev_matrix'=mat, 'prev_inverse'=mat_inverse))
  }

#' spike_slab_linear
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
spike_slab_linear <- 
  function(chain_length,X,y,tau0,tau1,q,a0=1,b0=1,rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,
           Xt=NULL, XXt=NULL, tau0_inverse=NULL, tau1_inverse=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if(is.null(XXt)){XXt <- fcprd(Xt)}
    if(is.null(tau0_inverse)){tau0_inverse <- chol2inv(chol(diag(n)+tau0^2*XXt))}
    if(is.null(tau1_inverse)){tau1_inverse <- chol2inv(chol(diag(n)+tau1^2*XXt))}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      sigma2 <- 1/rgamma(1,shape=(a0/2), rate=(b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0*sqrt(sigma2))
      beta[z==1] <- beta[z==1]*(tau1*sqrt(sigma2))
      
      prev_z <- z
      prev_matrix <- diag(n)+tau0^2*XXt+fcprd(Xt[prev_z==1,,drop=FALSE])*(tau1^2-tau0^2)
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_matrix <- matrix_full(Xt, prev_d)
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
                                 XXt, tau0_inverse, tau1_inverse,
                                 tau0, tau1, q, a0, b0, random_samples)
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


## Spike slab probit regression functions ##
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
spike_slab_probit_kernel <- 
  function(beta, z, e, X, Xt, y, prev_z, prev_matrix, prev_inverse,
           XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, random_samples){
    beta_output <- 
      update_beta(z, sigma2=1, X, Xt, e, prev_z, prev_matrix, prev_inverse,
                  XXt, tau0_inverse, tau1_inverse, tau0, tau1, 
                  u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- update_z(beta_new, sigma2=1, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- update_e(beta_new, sigma2=1, y, X, u_crn=random_samples$e_u)
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new,
                'prev_z'=z, 'prev_matrix'=mat, 'prev_inverse'=mat_inverse))
  }
#' spike_slab_probit
#' @description Generates Markov chain targeting the posterior corresponding to
#' Bayesian probit regression with spike and slab priors
#' @param chain_length Markov chain length
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @param y Response
#' @param tau0 prior hyperparameter (non-negative real)
#' @param tau1 prior hyperparameter (non-negative real)
#' @param q prior hyperparameter (strictly between 0 and 1)
#' @return Output from Markov chain targeting the posterior corresponding to
#' Bayesian logistic regression with spike and slab priors
#' @export
spike_slab_probit <- 
  function(chain_length,X,y,tau0,tau1,q,rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,
           Xt=NULL, XXt=NULL, tau0_inverse=NULL, tau1_inverse=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if(is.null(XXt)){XXt <- fcprd(Xt)}
    if(is.null(tau0_inverse)){tau0_inverse <- chol2inv(chol(diag(n)+tau0^2*XXt))}
    if(is.null(tau1_inverse)){tau1_inverse <- chol2inv(chol(diag(n)+tau1^2*XXt))}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      # nu <- 7.3
      # w2 <- (pi^2)*(nu-2)/(3*nu)
      # a0 <- nu
      # b0 <- w2*nu
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0)
      beta[z==1] <- beta[z==1]*(tau1)
      # e <- cpp_mat_vec_prod(X, beta) + rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + rnorm(n, mean = 0, sd = 1)
      
      prev_z <- z
      prev_matrix <- diag(n)+tau0^2*XXt+fcprd(Xt[prev_z==1,,drop=FALSE])*(tau1^2-tau0^2)
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_matrix <- matrix_full(Xt, prev_d)
      prev_inverse <- chol2inv(chol(prev_matrix))
    }
    
    if(store){
      beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
    } else {
      beta_ergodic_sum <- rep(0,p)
      z_ergodic_sum <- rep(0,p)
      e_ergodic_sum <- rep(0,n)
    }
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), e_u=runif(n), z_u=runif(p))
      new_state <- 
        spike_slab_probit_kernel(beta, z, e, X, Xt, y, 
                                 prev_z, prev_matrix, prev_inverse,
                                 XXt, tau0_inverse, tau1_inverse,
                                 tau0, tau1, q, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      prev_z <- new_state$prev_z
      prev_matrix <- new_state$prev_matrix
      prev_inverse <- new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          beta_samples[(t-burnin),] <- beta
          z_samples[(t-burnin),] <- z
          e_samples[(t-burnin),] <- e
        } else{
          beta_ergodic_sum <- beta_ergodic_sum + beta
          z_ergodic_sum <- z_ergodic_sum + z
          e_ergodic_sum <- e_ergodic_sum + e
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('beta'=beta_samples, 'z'=z_samples, 'e'=e_samples))
    } else {
      return(list('beta_ergodic_avg'=beta_ergodic_sum/(chain_length-burnin), 
                  'z_ergodic_avg'=z_ergodic_sum/(chain_length-burnin), 
                  'e_ergodic_avg'=e_ergodic_sum/(chain_length-burnin)))
    }
  }

## Spike slab logistic regression functions ##
# Output D%*%M%*%D for symmetric matrix M and a diagonal matrix D
diagonal_inner_outer_prod <- function(symetric_mat,diag_vec){
  return(diag_vec*t(symetric_mat*diag_vec))
}
matrix_full_logistic <- function(Xt, d, sigma2tilde){
  n <- dim(Xt)[2]
  return(diagonal_inner_outer_prod(crossprod(Xt*sqrt(1/d)),1/sqrt(sigma2tilde))+diag(n))
}

matrix_precompute_logistic <- 
  function(Xt, d, sigma2tilde, prev_d, prev_matrix, prev_sigma2tilde, 
           XXt, tau0, tau1){
  p <- dim(Xt)[1]
  n <- dim(Xt)[2]
  no_swaps <- sum(d!=prev_d)
  slab_indices <- d!=1/(tau0)^2
  no_slabs <- sum(slab_indices)
  no_spikes <- p-no_slabs
  
  if(no_swaps==0){
    diag(prev_matrix) <- diag(prev_matrix)-1
    prev_matrix <- diagonal_inner_outer_prod(prev_matrix,sqrt(prev_sigma2tilde))
    new_matrix <- diagonal_inner_outer_prod(prev_matrix,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  }
  if(no_slabs==0){
    new_matrix <- diagonal_inner_outer_prod(tau0^2*XXt,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  }
  if(no_spikes==0){
    new_matrix <- diagonal_inner_outer_prod(tau1^2*XXt,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  }
  
  if(no_swaps<min(no_slabs, no_spikes)){
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
    diag(prev_matrix) <- diag(prev_matrix)-1
    prev_matrix <- diagonal_inner_outer_prod(prev_matrix,sqrt(prev_sigma2tilde))
    prev_matrix <- pos_part-neg_part+prev_matrix
    new_matrix <- diagonal_inner_outer_prod(prev_matrix,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  } else if (no_slabs<=min(no_swaps, no_spikes)){
    Xt_update <- Xt[slab_indices,,drop=FALSE]
    update_part <- fcprd(Xt_update)*(tau1^2-tau0^2)
    new_matrix <- diagonal_inner_outer_prod(tau0^2*XXt+update_part,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  } else{
    Xt_update <- Xt[!slab_indices,,drop=FALSE]
    update_part <- fcprd(Xt_update)*(tau0^2-tau1^2)
    new_matrix <- diagonal_inner_outer_prod(tau1^2*XXt+update_part,1/sqrt(sigma2tilde))
    diag(new_matrix) <- diag(new_matrix)+1
    return(new_matrix)
  }
}
update_beta_logistic <- 
  function(z, sigma2tilde, X, Xt, y, prev_z, prev_sigma2tilde, prev_matrix, 
           XXt, tau0, tau1, u=NULL, delta=NULL){
  p <- length(z)
  n <- length(y)
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
  
  if(is.null(u)){u <- rnorm(p, 0, 1)}
  u = u/(d^0.5)
  if(is.null(delta)){delta <- c(rnorm(n,0,1))}
  # v = cpp_mat_vec_prod(sqrt(1/sigma2tilde)*X,u) + delta
  v = sqrt(1/sigma2tilde)*X%*%u + delta
  M_matrix <-
    matrix_precompute_logistic(Xt, d, sigma2tilde,
                               prev_d, prev_matrix, prev_sigma2tilde,
                               XXt, tau0, tau1)
  # M_matrix <- matrix_full_logistic(Xt, d, sigma2tilde)
  M_matrix_inverse <- chol2inv(chol(M_matrix))
  
  # v_star <- cpp_mat_vec_prod(M_matrix_inverse,(y/sqrt(sigma2tilde) - v))
  # beta <- (u + (1/d)*cpp_mat_vec_prod(Xt,sqrt(1/sigma2tilde)*v_star))
  v_star <- M_matrix_inverse%*%(y/sqrt(sigma2tilde) - v)
  beta <- u + (1/d)*(Xt%*%(sqrt(1/sigma2tilde)*v_star))
  return(list('beta'=beta, 'matrix'=M_matrix, 
              'matrix_inverse'=M_matrix_inverse))
}
# sigma2tilde update
update_sigma2tilde_logistic <- function(beta, z, tau0, tau1, a0, b0, X, y, u_crn=NULL){
  n <- length(y)
  # rss <- (y-cpp_mat_vec_prod(X,beta))^2
  rss <- (y-X%*%beta)^2
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  
  if(is.null(u_crn)){
    sigma2tilde <- 1/rgamma(n,shape = (a0+1)/2, rate = (b0+rss)/2)
  } else{
    log_sigma2tilde <- -log(qgamma(log(u_crn),shape = (a0+1)/2, rate = (b0+rss)/2, log.p = TRUE))
    sigma2tilde <- exp(log_sigma2tilde)
    # sigma2tilde <- 1/qgamma(log(u_crn),shape = (a0+1)/2, rate = (b0+rss)/2, log.p = TRUE)
  }
  # if(is.null(u_crn)){u_crn <- runif(n)}
  # sigma2tilde <- 1/qgamma(u_crn,shape = (a0+1)/2, rate = (b0+rss)/2)
  return(as.vector(sigma2tilde))
}
spike_slab_logistic_kernel <- 
  function(beta, z, e, sigma2tilde, X, Xt, y, 
           prev_z, prev_sigma2tilde, prev_matrix,
           XXt, tau0, tau1, q, random_samples){
    beta_output <-
      update_beta_logistic(z, sigma2tilde, X, Xt, e,
                           prev_z, prev_sigma2tilde, prev_matrix,
                           XXt, tau0, tau1,
                           u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    
    z_new <- update_z(beta_new, sigma2=1, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- update_e(beta_new, sigma2=sigma2tilde, y, X, u_crn=random_samples$e_u)
    nu <- 7.3
    w2 <- (pi^2)*(nu-2)/(3*nu)
    a0 <- nu
    b0 <- w2*nu
    sigma2tilde_new <-
      update_sigma2tilde_logistic(beta_new, z_new, tau0, tau1, a0, b0, X, e_new, u_crn=random_samples$sigma2tilde_u)
    
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new,'sigma2tilde'=sigma2tilde_new,
                'prev_z'=z,'prev_sigma2tilde'=sigma2tilde,'prev_matrix'=mat,'prev_inverse'=mat_inverse))
  }

#' spike_slab_logistic
#' @description Generates Markov chain targeting the posterior corresponding to
#' Bayesian logistic regression with spike and slab priors
#' @param chain_length Markov chain length
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @param y Response
#' @param tau0 prior hyperparameter (non-negative real)
#' @param tau1 prior hyperparameter (non-negative real)
#' @param q prior hyperparameter (strictly between 0 and 1)
#' @return Output from Markov chain targeting the posterior corresponding to
#' Bayesian logistic regression with spike and slab priors
#' @export
spike_slab_logistic <- 
  function(chain_length,X,y,tau0,tau1,q,rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,
           Xt=NULL, XXt=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if(is.null(XXt)){XXt <- fcprd(Xt)}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      nu <- 7.3
      w2 <- (pi^2)*(nu-2)/(3*nu)
      a0 <- nu
      b0 <- w2*nu
      sigma2tilde <- 1/rgamma(n,shape = (a0/2), rate = (b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0)
      beta[z==1] <- beta[z==1]*(tau1)
      # e <- cpp_mat_vec_prod(X, beta) + sqrt(sigma2tilde)*rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + sqrt(sigma2tilde)*rnorm(n, mean = 0, sd = 1)
      
      prev_z <- z
      prev_sigma2tilde <- sigma2tilde
      prev_matrix <- 
        diagonal_inner_outer_prod(tau0^2*XXt+fcprd(Xt[prev_z==1,,drop=FALSE])*(tau1^2-tau0^2),
                                  1/sqrt(sigma2tilde)) + diag(n)
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_matrix <- matrix_full_logistic(Xt, prev_d, prev_sigma2tilde)
      prev_inverse <- chol2inv(chol(prev_matrix))
    }
    
    if(store){
      beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
      sigma2tilde_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
    } else {
      beta_ergodic_sum <- rep(0,p)
      z_ergodic_sum <- rep(0,p)
      e_ergodic_sum <- rep(0,n)
      sigma2tilde_ergodic_sum <- rep(0,n)
    }
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), e_u=runif(n), z_u=runif(p), sigma2tilde_u=runif(n))
      new_state <- 
        spike_slab_logistic_kernel(beta, z, e, sigma2tilde, X, Xt, y, 
                                   prev_z, prev_sigma2tilde, prev_matrix,
                                   XXt, tau0, tau1, q, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      sigma2tilde <- new_state$sigma2tilde
      prev_z <- new_state$prev_z
      prev_sigma2tilde <- new_state$prev_sigma2tilde
      prev_matrix <- new_state$prev_matrix
      prev_inverse <- new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          beta_samples[(t-burnin),] <- beta
          z_samples[(t-burnin),] <- z
          e_samples[(t-burnin),] <- e
          sigma2tilde_samples[(t-burnin),] <- sigma2tilde
        } else{
          beta_ergodic_sum <- beta_ergodic_sum + beta
          z_ergodic_sum <- z_ergodic_sum + z
          e_ergodic_sum <- e_ergodic_sum + e
          sigma2tilde_ergodic_sum <- sigma2tilde_ergodic_sum + sigma2tilde
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('beta'=beta_samples, 'z'=z_samples, 'e'=e_samples, 'sigma2tilde'=sigma2tilde_samples))
    } else {
      return(list('beta_ergodic_avg'=beta_ergodic_sum/(chain_length-burnin), 
                  'z_ergodic_avg'=z_ergodic_sum/(chain_length-burnin), 
                  'e_ergodic_avg'=e_ergodic_sum/(chain_length-burnin), 
                  'sigma2tilde_ergodic_avg'=sigma2tilde_ergodic_sum/(chain_length-burnin)))
    }
  }