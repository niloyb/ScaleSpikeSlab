# Functions for comparison with alternatives
library(ScaleSpikeSlab)
library(doParallel)
registerDoParallel(cores = detectCores()-1)
library(foreach)

############################ Skinny Gibbs functions ############################
# # Installing Skinny Gibbs
# install.packages("/Users/niloybiswas/Downloads/UASA_A_1482754_Supplement/Skinny Gibbs/skinnybasad_0.0.1.tar.gz", repos = NULL, type ="source")
library(skinnybasad)
# Editing the skinnybasad:::skinnybasad function for implementation
# Added the line: PACKAGE = "skinnybasad"
skinnybasad <- function(X, E, pr, B0, Z0, nsplit, modif, nburn = 1000, niter = 5000, 
                        printitrsep = 0, maxsize = 1000, a0 = 1, b0 = 1){
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
  s1 = max(a0 * p^{
    2 + del
  }/n, 1)
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
    sigma2_new <- ScaleSpikeSlab:::update_sigma2(beta_new, z_new, tau0, tau1, a0, b0, X, y, u_crn=random_samples$sigma2_u)
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
        sota_spike_slab_linear(beta, z, sigma2, X, Xt, y, 
                               tau0, tau1, q, a0, b0, random_samples)
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

sota_spike_slab_logistic_kernel <- 
  function(beta, z, e, sigma2, X, Xt, y, tau0, tau1, q, random_samples){
    beta_output <- 
      update_beta_sota(z, sigma2, X, Xt, e, tau0, tau1, u=random_samples$beta_u, delta=random_samples$beta_delta)
    mat <- beta_output$matrix
    mat_inverse <- beta_output$matrix_inverse
    beta_new <- beta_output$beta
    z_new <- ScaleSpikeSlab:::update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples$z_u)
    e_new <- ScaleSpikeSlab:::update_e(beta_new, sigma2, y, X, u_crn=random_samples$e_u)
    nu <- 7.3
    w2 <- (pi^2)*(nu-2)/(3*nu)
    a0 <- nu
    b0 <- w2*nu
    sigma2_new <- ScaleSpikeSlab:::update_sigma2(beta_new, z_new, tau0, tau1, a0, b0, X, e_new, u_crn=random_samples$sigma2_u)
    return(list('beta'=beta_new,'z'=z_new,'e'=e_new,'sigma2'=sigma2_new))
  }

sota_spike_slab_logistic <- 
  function(chain_length,X,y,tau0,tau1,q, rinit=NULL,verbose=FALSE,
           burnin=0,store=TRUE,Xt=NULL){
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
      sigma2 <- 1/rgamma(1,shape = (a0/2), rate = (b0/2))
      beta <- rnorm(p)
      beta[z==0] <- beta[z==0]*(tau0*sqrt(sigma2))
      beta[z==1] <- beta[z==1]*(tau1*sqrt(sigma2))
      # e <- cpp_mat_vec_prod(X, beta) + sqrt(sigma2)*rnorm(n, mean = 0, sd = 1)
      e <- X%*%beta + sqrt(sigma2)*rnorm(n, mean = 0, sd = 1)
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
        sota_spike_slab_logistic_kernel(beta, z, e, sigma2, X, Xt, y, 
                                        tau0, tau1, q, random_samples)
      
      beta <- new_state$beta
      z <- new_state$z
      e <- new_state$e
      sigma2 <- new_state$sigma2
      
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

sota_spike_slab_mcmc <- 
  function(chain_length,X,y,tau0,tau1,q,a0=1,b0=1,
           type=NULL,rinit=NULL,verbose=FALSE,burnin=0,store=TRUE,Xt=NULL){
    if(is.null(Xt)){Xt <- t(X)}
    
    # Binary 0,1 labels for logistic regression
    if(all(y*(1-y)==0)){ 
      return(sota_spike_slab_logistic(chain_length,X,y,tau0,tau1,q,rinit,verbose,burnin,store,Xt))
    } else {
      return(sota_spike_slab_linear(chain_length,X,y,tau0,tau1,q,a0,b0,rinit,verbose,burnin,store,Xt))
    }
  }



#################### SOTA, Skinny Gibbs and S^3 comparison ####################
comparison_sims <- function(n_p_error_s0_list,chain_length=1e4,burnin=5e3,no_repeats=1,
                            algos=c('ScalableSpikeSlab', 'Sota', 'SkinnyGibbs'), 
                            signal='constant', store=TRUE){
  foreach(n_p_error_s0 = n_p_error_s0_list, .combine = rbind)%:%
    foreach(i = c(1:no_repeats), .combine = rbind)%dopar%{
      n <- n_p_error_s0$n
      p <- n_p_error_s0$p
      error_std <- n_p_error_s0$error_std
      s0 <- n_p_error_s0$s0
      
      logreg_sync_data <- synthetic_data(n, p, s0, type = 'logistic', signal=signal)
      X <- logreg_sync_data$X
      X <- matrix(scale(X), n, p)
      y <- logreg_sync_data$y
      Xt <- t(X)
      signal_indices <- logreg_sync_data$true_beta!=0
      
      params <- spike_slab_params(n, p)
      
      output <- data.frame()
      
      if('ScalableSpikeSlab' %in% algos){
        ###### Scalable spike and slab
        sss_time_taken <-
          system.time(
            sss_chain <- 
              spike_slab_mcmc(chain_length=chain_length, X=X,Xt=Xt,y=y,
                                         tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                         a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE,
                                         store=store))
        if(store){
          delta <- rowSums(sss_chain$z[c(1:(chain_length-1)),]!=sss_chain$z[c(2:chain_length),])
          no_active <- rowSums(sss_chain$z[c(1:chain_length),])
          
          sss_tpr <- mean(colMeans(sss_chain$z[c(burnin:chain_length),signal_indices,drop=FALSE])>0.5)
          sss_fdr <- mean(colMeans(sss_chain$z[c(burnin:chain_length),!signal_indices,drop=FALSE])>0.5)
          sss_mse <- mean((colMeans(sss_chain$beta[c(burnin:chain_length),])-logreg_sync_data$true_beta)^2)
          
          output <- 
            rbind(output, 
                  data.frame(algo='ScalableSpikeSlab', time=as.double(sss_time_taken[1])/chain_length, 
                             tpr=sss_tpr, fdr=sss_fdr, mse=sss_mse, delta_mean=mean(delta), delta_var=var(delta),
                             no_active_mean=mean(delta), no_active_var=var(delta), n=n, p=p, s0, iteration=i))
        } else{
          output <- 
            rbind(output, data.frame(algo='ScalableSpikeSlab', 
                                     time=as.double(sss_time_taken[1])/chain_length, n=n, p=p, s0, iteration=i))
        }
      }
      
      if('Sota' %in% algos){
        sota_time_taken <-
          system.time(
            sota_chain <- 
              sota_spike_slab_mcmc(chain_length=chain_length, X=X,Xt=Xt,y=y,
                                         tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                         a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE,
                                         store=store))
        if(store){
          delta <- rowSums(sota_chain$z[c(1:(chain_length-1)),]!=sota_chain$z[c(2:chain_length),])
          no_active <- rowSums(sota_chain$z[c(1:chain_length),])
          
          sota_tpr <- mean(colMeans(sota_chain$z[c(burnin:chain_length),signal_indices,drop=FALSE])>0.5)
          sota_fdr <- mean(colMeans(sota_chain$z[c(burnin:chain_length),!signal_indices,drop=FALSE])>0.5)
          sota_mse <- mean((colMeans(sota_chain$beta[c(burnin:chain_length),])-logreg_sync_data$true_beta)^2)
          
          output <- 
            rbind(output, 
                  data.frame(algo='Sota', time=as.double(sota_time_taken[1])/chain_length, 
                             tpr=sota_tpr, fdr=sota_fdr, mse=sota_mse, delta_mean=mean(delta), delta_var=var(delta), 
                             no_active_mean=mean(no_active), no_active_var=var(no_active), n=n, p=p, s0, iteration=i))
        } else{
          output <- 
            rbind(output, 
                  data.frame(algo='Sota', time=as.double(sota_time_taken[1])/chain_length, 
                             n=n, p=p, s0, iteration=i))
        }
      }
      
      if('SkinnyGibbs' %in% algos){
        # Skinny Gibbs: Initializing from the prior
        z <- rbinom(p,1,params$q)
        sigma2 <- 1/rgamma(1,shape = (params$a0/2), rate = (params$b0/2))
        beta <- rnorm(p)
        beta[z==0] <- beta[z==0]*(params$tau0*sqrt(sigma2))
        beta[z==1] <- beta[z==1]*(params$tau1*sqrt(sigma2))
        skinnygibbs_time_taken <-
          system.time(
            skinny_chain <-
              skinnybasad(X,y,params$q,beta,z,nsplit=10,modif=1,nburn=burnin,
                          niter=(chain_length-burnin),printitrsep=1,a0=1,b0=1)
          )
        
        if(store){
          skinny_tpr <- mean(skinny_chain$marZ[signal_indices,drop=FALSE]>0.5)
          skinny_fdr <- mean(skinny_chain$marZ[!signal_indices,drop=FALSE]>0.5)
          skinny_mse <- NA
          output <- 
            rbind(output, 
                  data.frame(algo='SkinnyGibbs', time=as.double(skinnygibbs_time_taken[1])/chain_length, 
                             tpr=skinny_tpr, fdr=skinny_fdr, mse=skinny_mse,
                             delta_mean=NA, delta_var=NA, 
                             no_active_mean=NA, no_active_var=NA, 
                             n=n, p=p, s0, iteration=i))
        } else{
          output <- 
            rbind(output, 
                  data.frame(algo='SkinnyGibbs', time=as.double(skinnygibbs_time_taken[1])/chain_length, 
                             n=n, p=p, s0, iteration=i))
        }
      }
      
      print(n_p_error_s0)
      return(output)
    }
}


