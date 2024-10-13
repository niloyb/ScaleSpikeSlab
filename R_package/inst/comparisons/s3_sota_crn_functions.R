# Functions for numerical stability plot
library(ScaleSpikeSlab)
library(doParallel)
registerDoParallel(cores = detectCores()-1)
library(foreach)

################################ Skinny Gibbs functions ################################
# # Installing Skinny Gibbs
# install.packages("/Users/niloybiswas/Downloads/UASA_A_1482754_Supplement/Skinny Gibbs/skinnybasad_0.0.1.tar.gz", repos = NULL, type ="source")
library(skinnybasad)
# Editing the skinnybasad:::skinnybasad function for implementation
# Added the line: PACKAGE = "skinnybasad"
skinnybasad <- function(X, E, pr, B0, Z0, nsplit=10, modif, nburn = 1000, niter = 5000, 
                        printitrsep = 0, maxsize = 1000, a0 = 0.01, b0 = 1){
  n = as.integer(dim(X)[1])
  p = as.integer(dim(X)[2])
  printitr = 0
  if (printitrsep != 0) 
    printitr = 1
  if (!is.vector(B0)) 
    B0 = rep(0, p)
  if (!is.vector(Z0)) 
    Z0 = rep(0, p)
  del = 0.1
  
  s1 = 1
  # s1 = max(a0 * p^{2 + del}/n, 1)
  s0 = b0/n
  
  res = .C("skinnybasad", as.double(as.vector(X)), as.double(as.vector(E)), 
           as.double(as.vector(B0)), as.double(as.vector(Z0)), as.double(pr), 
           n, p, as.double(s1), as.double(s0), as.integer(nburn), 
           as.integer(niter), as.integer(nsplit), as.integer(modif), 
           as.integer(printitr), as.integer(printitrsep), as.integer(maxsize), 
           outmarZ = double(p), outsize = integer(niter + nburn),
           PACKAGE = "skinnybasad")
  return(list(marZ = res$outmarZ, allsize = res$outsize[nburn:(nburn + niter)]))
}


################################ SOTA functions ################################
# SOTA linear
update_beta_sota <- 
  function(z, sigma2, X, Xt, y, tau0, tau1, u=NULL, delta=NULL){
    p <- length(z)
    n <- length(y)
    d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
    
    if(is.null(u)){u <- rnorm(p, 0, 1)}
    u = u/(d^0.5)
    if(is.null(delta)){delta <- c(rnorm(n,0,1))}
    # v = cpp_mat_vec_prod(X,u) + delta
    v = X%*%u + delta
    
    M_matrix <- ScaleSpikeSlab:::matrix_full(Xt, d)
    M_matrix_inverse <- chol2inv(chol(M_matrix))
    
    # v_star <- cpp_mat_vec_prod(M_matrix_inverse,(y/sqrt(sigma2) - v))
    # beta <- sqrt(sigma2)*(u + (d^(-1))*cpp_mat_vec_prod(Xt,v_star))
    v_star <- M_matrix_inverse%*%(y/sqrt(sigma2) - v)
    beta <- sqrt(sigma2)*(u + (d^(-1))*(Xt%*%v_star))
    return(list('beta'=beta, 'matrix'=M_matrix, 
                'matrix_inverse'=M_matrix_inverse))
  }
sota_spike_slab_linear_kernel <- 
  function(beta, z, sigma2, X, Xt, y, tau0, tau1, q, a0, b0, random_samples){
    beta_output <- 
      update_beta_sota(z, sigma2, X, Xt, y, tau0, tau1, 
                       u=random_samples$beta_u, delta=random_samples$beta_delta)
    beta_new <- beta_output$beta
    z_new <- ScaleSpikeSlab:::update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples$z_u)
    sigma2_new <- ScaleSpikeSlab:::update_sigma2_linear(beta_new, z_new, tau0, tau1, a0, b0, X, y, u_crn=random_samples$sigma2_u)
    return(list('beta'=beta_new,'z'=z_new,'sigma2'=sigma2_new))
  }
sota_spike_slab_linear <- 
  function(chain_length,X,y,tau0,tau1,q,a0=1,b0=1,rinit=NULL,
           verbose=FALSE,burnin=0,store=TRUE,Xt=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      sigma2 <- 1/rgamma(1,shape = (a0/2), rate = (b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0*sqrt(sigma2))
      beta[z==1] <- beta[z==1]*(tau1*sqrt(sigma2))
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
        sota_spike_slab_linear_kernel(beta, z, sigma2, X, Xt, y, tau0, tau1, q, a0, b0, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      sigma2 <- new_state$sigma2
      
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

# SOTA probit
sota_spike_slab_probit_kernel <- 
  function(beta, z, e, X, Xt, y, tau0, tau1, q, random_samples){
    beta_output <- 
      update_beta_sota(z, sigma2=1, X, Xt, e, tau0, tau1, 
                       u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- ScaleSpikeSlab:::update_z(beta_new, sigma2=1, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- ScaleSpikeSlab:::update_e(beta_new, sigma2=1, y, X, u_crn=random_samples$e_u)
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new))
  }
sota_spike_slab_probit <- 
  function(chain_length,X,y,tau0,tau1,q,rinit=NULL,
           verbose=FALSE,burnin=0,store=TRUE,Xt=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      nu <- 7.3
      w2 <- (pi^2)*(nu-2)/(3*nu)
      a0 <- nu
      b0 <- w2*nu
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0)
      beta[z==1] <- beta[z==1]*(tau1)
      # e <- cpp_mat_vec_prod(X, beta) + rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + rnorm(n, mean = 0, sd = 1)
      
      # prev_z <- z
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_matrix <- ScaleSpikeSlab:::matrix_full(Xt, prev_d)
      # prev_inverse <- chol2inv(chol(prev_matrix))
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
        sota_spike_slab_probit_kernel(beta, z, e, X, Xt, y, tau0, tau1, q, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      
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



# SOTA logistic
update_beta_logistic_sota <- 
  function(z, sigma2tilde, X, Xt, y, XXt, tau0, tau1, u=NULL, delta=NULL){
    p <- length(z)
    n <- length(y)
    d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
    # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
    
    if(is.null(u)){u <- rnorm(p, 0, 1)}
    u = u/(d^0.5)
    if(is.null(delta)){delta <- c(rnorm(n,0,1))}
    # v = cpp_mat_vec_prod(sqrt(1/sigma2tilde)*X,u) + delta
    v = sqrt(1/sigma2tilde)*X%*%u + delta
    # M_matrix <-
    #   matrix_precompute_logistic(Xt, d, sigma2tilde,
    #                              prev_d, prev_matrix, prev_sigma2tilde,
    #                              XXt, tau0, tau1)
    M_matrix <- ScaleSpikeSlab:::matrix_full_logistic(Xt, d, sigma2tilde)
    M_matrix_inverse <- chol2inv(chol(M_matrix))
    
    # v_star <- cpp_mat_vec_prod(M_matrix_inverse,(y/sqrt(sigma2tilde) - v))
    # beta <- (u + (1/d)*cpp_mat_vec_prod(Xt,sqrt(1/sigma2tilde)*v_star))
    v_star <- M_matrix_inverse%*%(y/sqrt(sigma2tilde) - v)
    beta <- u + (1/d)*(Xt%*%(sqrt(1/sigma2tilde)*v_star))
    return(list('beta'=beta, 'matrix'=M_matrix, 
                'matrix_inverse'=M_matrix_inverse))
  }
sota_spike_slab_logistic_kernel <- 
  function(beta, z, e, sigma2tilde, X, Xt, y, XXt, tau0, tau1, q, random_samples){
    beta_output <-
      update_beta_logistic_sota(z, sigma2tilde, X, Xt, e, XXt, tau0, tau1,
                                u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    
    z_new <- ScaleSpikeSlab:::update_z(beta_new, sigma2=1, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- ScaleSpikeSlab:::update_e(beta_new, sigma2=sigma2tilde, y, X, u_crn=random_samples$e_u)
    
    # sigma2tilde_new <- rep(1,length(y)) 
    nu <- 7.3
    w2 <- (pi^2)*(nu-2)/(3*nu)
    a0 <- nu
    b0 <- w2*nu
    sigma2tilde_new <-
      ScaleSpikeSlab:::update_sigma2tilde_logistic(beta_new, z_new, tau0, tau1, a0, b0, X, e_new, u_crn=random_samples$sigma2tilde_u)
    
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new,'sigma2tilde'=sigma2tilde_new,
                'prev_z'=z,'prev_sigma2tilde'=sigma2tilde,'prev_matrix'=mat,'prev_inverse'=mat_inverse))
  }
sota_spike_slab_logistic <- 
  function(chain_length,X,y,tau0,tau1,q,
           rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,Xt=NULL, XXt=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if(is.null(XXt)){XXt <- ScaleSpikeSlab:::fcprd(Xt)}
    
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
      
      # prev_z <- z
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_sigma2tilde <- sigma2tilde
      # prev_matrix <- ScaleSpikeSlab:::matrix_full_logistic(Xt, prev_d, prev_sigma2tilde)
      # prev_inverse <- chol2inv(chol(prev_matrix))
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
        sota_spike_slab_logistic_kernel(beta, z, e, sigma2tilde, X, Xt, y, 
                                        XXt, tau0, tau1, q, random_samples)
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      sigma2tilde <- new_state$sigma2tilde
      
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

########################### SOTA & S3 CRN coupling #############################
s3_sota_crn_logistic <- 
  function(chain_length,X,y,tau0,tau1,q,
           rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,Xt=NULL, XXt=NULL){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if(is.null(XXt)){XXt <- ScaleSpikeSlab:::fcprd(Xt)}
    
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
      
      # prev_z <- z
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_sigma2tilde <- sigma2tilde
      # prev_matrix <- ScaleSpikeSlab:::matrix_full_logistic(Xt, prev_d, prev_sigma2tilde)
      # prev_inverse <- chol2inv(chol(prev_matrix))
      prev_z <- z
      prev_sigma2tilde <- sigma2tilde
      prev_matrix <- ScaleSpikeSlab:::diagonal_inner_outer_prod(tau0^2*XXt + ScaleSpikeSlab:::fcprd(Xt[prev_z==1,,drop=FALSE])*(tau1^2-tau0^2), 1/sqrt(sigma2tilde)) + diag(n)
      prev_inverse <- chol2inv(chol(prev_matrix))
      
      
    }
    
    if(store){
      sota_beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      sota_z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      sota_e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
      sota_sigma2tilde_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
      
      s3_beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      s3_z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      s3_e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
      s3_sigma2tilde_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
    } else {
      sota_beta_ergodic_sum <- rep(0,p)
      sota_z_ergodic_sum <- rep(0,p)
      sota_e_ergodic_sum <- rep(0,n)
      sota_sigma2tilde_ergodic_sum <- rep(0,n)
      
      s3_beta_ergodic_sum <- rep(0,p)
      s3_z_ergodic_sum <- rep(0,p)
      s3_e_ergodic_sum <- rep(0,n)
      s3_sigma2tilde_ergodic_sum <- rep(0,n)
    }
    
    sota_beta <- beta
    sota_z <- z
    sota_e <- e
    sota_sigma2tilde <- sigma2tilde
    
    s3_beta <- beta
    s3_z <- z
    s3_e <- e
    s3_sigma2tilde <- sigma2tilde
    s3_prev_z <- prev_z
    s3_prev_sigma2tilde <- prev_sigma2tilde
    s3_prev_matrix <- prev_matrix
    s3_prev_inverse <- prev_inverse
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), e_u=runif(n), z_u=runif(p), sigma2tilde_u=runif(n))
      sota_new_state <- 
        sota_spike_slab_logistic_kernel(sota_beta, sota_z, sota_e, sota_sigma2tilde, X, Xt, y, 
                                        XXt, tau0, tau1, q, random_samples)
      sota_beta <- sota_new_state$beta
      sota_z <- sota_new_state$z
      sota_e <- sota_new_state$e
      sota_sigma2tilde <- sota_new_state$sigma2tilde
      
      
      s3_new_state <- 
        ScaleSpikeSlab:::spike_slab_logistic_kernel(s3_beta, s3_z, s3_e, s3_sigma2tilde, 
                                   X, Xt, y, s3_prev_z, s3_prev_sigma2tilde, s3_prev_matrix, 
                                   XXt, tau0, tau1, q, random_samples)
      s3_beta <- s3_new_state$beta
      s3_z <- s3_new_state$z
      s3_e <- s3_new_state$e
      s3_sigma2tilde <- s3_new_state$sigma2tilde
      s3_prev_z <- s3_new_state$prev_z
      s3_prev_sigma2tilde <- s3_new_state$prev_sigma2tilde
      s3_prev_matrix <- s3_new_state$prev_matrix
      s3_prev_inverse <- s3_new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          sota_beta_samples[(t-burnin),] <- sota_beta
          sota_z_samples[(t-burnin),] <- sota_z
          sota_e_samples[(t-burnin),] <- sota_e
          sota_sigma2tilde_samples[(t-burnin),] <- sota_sigma2tilde
          
          s3_beta_samples[(t-burnin),] <- s3_beta
          s3_z_samples[(t-burnin),] <- s3_z
          s3_e_samples[(t-burnin),] <- s3_e
          s3_sigma2tilde_samples[(t-burnin),] <- s3_sigma2tilde
        } else{
          sota_beta_ergodic_sum <- sota_beta_ergodic_sum + sota_beta
          sota_z_ergodic_sum <- sota_z_ergodic_sum + sota_z
          sota_e_ergodic_sum <- sota_e_ergodic_sum + sota_e
          sota_sigma2tilde_ergodic_sum <- sota_sigma2tilde_ergodic_sum + sota_sigma2tilde
          
          s3_beta_ergodic_sum <- s3_beta_ergodic_sum + s3_beta
          s3_z_ergodic_sum <- s3_z_ergodic_sum + s3_z
          s3_e_ergodic_sum <- s3_e_ergodic_sum + s3_e
          s3_sigma2tilde_ergodic_sum <- s3_sigma2tilde_ergodic_sum + s3_sigma2tilde
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('sota'=list('beta'=sota_beta_samples, 'z'=sota_z_samples, 'e'=sota_e_samples, 'sigma2tilde'=sota_sigma2tilde_samples),
                  's3'=list('beta'=s3_beta_samples, 'z'=s3_z_samples, 'e'=s3_e_samples, 'sigma2tilde'=s3_sigma2tilde_samples)))
    } else {
      return(list('sota'=list('beta_ergodic_avg'=sota_beta_ergodic_sum/(chain_length-burnin), 
                              'z_ergodic_avg'=sota_z_ergodic_sum/(chain_length-burnin), 
                              'e_ergodic_avg'=sota_e_ergodic_sum/(chain_length-burnin), 
                              'sigma2tilde_ergodic_avg'=sota_sigma2tilde_ergodic_sum/(chain_length-burnin)),
                  's3'=list('beta_ergodic_avg'=s3_beta_ergodic_sum/(chain_length-burnin), 
                            'z_ergodic_avg'=s3_z_ergodic_sum/(chain_length-burnin), 
                            'e_ergodic_avg'=s3_e_ergodic_sum/(chain_length-burnin), 
                            'sigma2tilde_ergodic_avg'=s3_sigma2tilde_ergodic_sum/(chain_length-burnin))))
    }
  }



s3_sota_crn_probit <- 
  function(chain_length,X,y,tau0,tau1,q,rinit=NULL,
           verbose=FALSE,burnin=0,store=TRUE,Xt=NULL,XXt=NULL,
           tau0_inverse=NULL,tau1_inverse=NULL){
    
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(Xt)){Xt <- t(X)}
    if (is.null(XXt)) {XXt <- ScaleSpikeSlab:::fcprd(Xt)}
    if (is.null(tau0_inverse)) {tau0_inverse <- chol2inv(chol(diag(n) + tau0^2 * XXt))}
    if (is.null(tau1_inverse)) {tau1_inverse <- chol2inv(chol(diag(n) + tau1^2 * XXt))}
    
    if(is.null(rinit)){
      # Initializing from the prior
      z <- rbinom(p,1,q)
      nu <- 7.3
      w2 <- (pi^2)*(nu-2)/(3*nu)
      a0 <- nu
      b0 <- w2*nu
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0)
      beta[z==1] <- beta[z==1]*(tau1)
      # e <- cpp_mat_vec_prod(X, beta) + rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + rnorm(n, mean = 0, sd = 1)
      
      # prev_z <- z
      # prev_d <- as.vector(prev_z/(tau1^2)+(1-prev_z)/(tau0^2))
      # prev_matrix <- ScaleSpikeSlab:::matrix_full(Xt, prev_d)
      # prev_inverse <- chol2inv(chol(prev_matrix))
      
      prev_z <- z
      prev_matrix <- diag(n)+tau0^2*XXt + ScaleSpikeSlab:::fcprd(Xt[prev_z==1,,drop=FALSE])*(tau1^2-tau0^2)
      prev_inverse <- chol2inv(chol(prev_matrix))
    }
    
    if(store){
      sota_beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      sota_z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      sota_e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
      
      s3_beta_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      s3_z_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = p)
      s3_e_samples <- matrix(NA, nrow = (chain_length-burnin), ncol = n)
    } else {
      sota_beta_ergodic_sum <- rep(0,p)
      sota_z_ergodic_sum <- rep(0,p)
      sota_e_ergodic_sum <- rep(0,n)
      
      s3_beta_ergodic_sum <- rep(0,p)
      s3_z_ergodic_sum <- rep(0,p)
      s3_e_ergodic_sum <- rep(0,n)
    }
    
    sota_beta <- beta
    sota_z <- z
    sota_e <- e
    
    s3_beta <- beta
    s3_z <- z
    s3_e <- e
    s3_prev_z <- prev_z
    s3_prev_matrix <- prev_matrix
    s3_prev_inverse <- prev_inverse
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), e_u=runif(n), z_u=runif(p))
      sota_new_state <- 
        sota_spike_slab_probit_kernel(sota_beta, sota_z, sota_e, X, Xt, y, tau0, tau1, q, random_samples)
      sota_beta <- sota_new_state$beta
      sota_z <- sota_new_state$z
      sota_e <- sota_new_state$e
      
      s3_new_state <- 
        ScaleSpikeSlab:::spike_slab_probit_kernel(s3_beta, s3_z, s3_e, X, Xt, y, s3_prev_z, s3_prev_matrix, s3_prev_inverse, 
                                                  XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, random_samples)
      s3_beta <- s3_new_state$beta
      s3_z <- s3_new_state$z
      s3_e <- s3_new_state$e
      s3_prev_z <- s3_new_state$prev_z
      s3_prev_matrix <- s3_new_state$prev_matrix
      s3_prev_inverse <- s3_new_state$prev_inverse
      
      if(t>burnin){
        if(store){
          sota_beta_samples[(t-burnin),] <- sota_beta
          sota_z_samples[(t-burnin),] <- sota_z
          sota_e_samples[(t-burnin),] <- sota_e
          
          s3_beta_samples[(t-burnin),] <- s3_beta
          s3_z_samples[(t-burnin),] <- s3_z
          s3_e_samples[(t-burnin),] <- s3_e
        } else{
          sota_beta_ergodic_sum <- sota_beta_ergodic_sum + sota_beta
          sota_z_ergodic_sum <- sota_z_ergodic_sum + sota_z
          sota_e_ergodic_sum <- sota_e_ergodic_sum + sota_e
          
          s3_beta_ergodic_sum <- s3_beta_ergodic_sum + s3_beta
          s3_z_ergodic_sum <- s3_z_ergodic_sum + s3_z
          s3_e_ergodic_sum <- s3_e_ergodic_sum + s3_e
        }
      }
      if(verbose){print(t)}
    }
    
    if(store){
      return(list('sota'=list('beta'=sota_beta_samples, 'z'=sota_z_samples, 'e'=sota_e_samples),
                  's3'=list('beta'=s3_beta_samples, 'z'=s3_z_samples, 'e'=s3_e_samples)))
    } else {
      return(list('sota'=list('beta_ergodic_avg'=sota_beta_ergodic_sum/(chain_length-burnin), 
                              'z_ergodic_avg'=sota_z_ergodic_sum/(chain_length-burnin), 
                              'e_ergodic_avg'=sota_e_ergodic_sum/(chain_length-burnin)),
                  's3'=list('beta_ergodic_avg'=s3_beta_ergodic_sum/(chain_length-burnin), 
                            'z_ergodic_avg'=s3_z_ergodic_sum/(chain_length-burnin), 
                            'e_ergodic_avg'=s3_e_ergodic_sum/(chain_length-burnin))))
    }
  }




