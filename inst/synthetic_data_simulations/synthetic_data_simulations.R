# Synthetic dataset simulations
rm(list=ls())
source('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/FastSpikeSlab/inst/comparisons/comparison_functions.R')
# library(FastSpikeSlab)

library(dplyr)
library(ggplot2)
library(latex2exp)


###### Funcstions to generate synthetic data and choose hyperparameters ######
# Generate synthetic dataset
synthetic_data <- function(n,p,s0, type='linear'){
  true_beta <- matrix(0,p,1)
  true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))
  
  # True beta: following a example from skinny gibbs paper
  # s0 <- 4
  # true_beta[1:s0] = c(-1.5,2,-2.5,3) 
  
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  X_truebeta <- X%*%true_beta
  
  if(type=='linear'){
    # Error terms 
    error_std <- 2
    error_terms = error_std*rnorm(n, mean = 0, sd = 1)
    y = X_truebeta + error_terms
  } else if(type=='logistic'){
    true_aug_y = rlogis(n, location = X_truebeta) # recall variance is pi^2/3
    y <- ifelse(true_aug_y>0,1,0) # Logistic response
  }
  return(list(X=X, y=y, true_beta=true_beta))
}
# Choose hyperparameters
spike_slab_params <- function(n, p){
  # Choice of q, tau0, tau1: following skinny gibbs paper
  K <- max(10,log(n))
  q_seq <- seq(0,1,0.0001)
  probs <- abs(pbinom(K,p,q_seq)-0.9)
  q_index <- which(probs==min(probs))
  if(length(q_index)>1){
    q <- 1/p
  } else {
    q <- q_seq[q_index]
  }
  tau0 <- 1/sqrt(n)
  tau1 <- 1
  # tau1 <- sqrt(max(1, p^(2.1)/(100*n)))
  a0 <- 1
  b0 <- 1
  return(list(K=K, q=q, tau0=tau0, tau1=tau1, a0=a0, b0=b0))
}

######## Fix n vary p ########
n <- 100
s0 <- 5
time_taken_df <- data.frame()
for (p in seq(100,1000,100)){
  linreg_sync_data <- synthetic_data(n, p, s0)
  X <- linreg_sync_data$X
  y <- linreg_sync_data$y
  Xt <- t(X)
  
  params <- spike_slab_params(n, p)
  
  burnin <- 5e2
  chain_length <- burnin+5e2
  ###### Fast spike and slab
  # Full Gibbs with w fixed
  fastspikeslab_time_taken <-
    system.time(
      log_reg <- 
        spike_slab_mcmc(chain_length=1e3, X=X,Xt=t(X),y=y,
                        tau0=params$tau0, tau1=params$tau1, q=params$q, 
                        a0=params$a0,b0=params$b0, rinit=NULL, verbose=FALSE))

  delta <- rowSums(log_reg$z[c(1:(chain_length-1)),]!=log_reg$z[c(2:chain_length),])
  no_active <- rowSums(log_reg$z[c(1:chain_length),])
  
  time_taken_df <- 
    rbind(time_taken_df, 
          data.frame(algo='fast_spike_slab',
                     time=as.double(fastspikeslab_time_taken[1])/chain_length, 
                     delta_mean=mean(delta), delta_var=var(delta), 
                     no_active_mean=mean(no_active), no_active_var=var(no_active), 
                     n=n, p=p))
  print(p)
}

ggplot(time_taken_df %>% filter(n==100), aes(x=p, y=time*1000, linetype=algo)) + 
geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('Time per iteration (ms)')) +
theme_classic(base_size = 18) +
theme(legend.position="bottom")

ggplot(time_taken_df %>% filter(n==100), aes(x=p, y=delta_mean, linetype=algo)) + 
geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('p_t')) +
theme_classic(base_size = 18) +
theme(legend.position="bottom")

ggplot(time_taken_df %>% filter(n==100), aes(x=p, y=no_active_mean, linetype=algo)) + 
geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('p_t')) +
theme_classic(base_size = 18) +
theme(legend.position="bottom")


p <- 100
s0 <- 20
time_taken_df <- data.frame()
for (n in seq(100,2000,100)){
  lingreg_sync_data <- synthetic_data(n, p, s0)
  X <- lingreg_sync_data$X
  y <- lingreg_sync_data$y
  Xt <- t(X)
  
  params <- spike_slab_params(n, p)
  
  burnin <- 5e2
  chain_length <- burnin+5e2
  ###### Fast spike and slab
  # Full Gibbs with w fixed
  fast_time_taken <-
    system.time(
      lin_reg_test <- 
        spike_slab_mcmc(chain_length=1e3, X=X,Xt=t(X),y=y,
                        tau0=params$tau0, tau1=params$tau1, q=params$q, 
                        a0=params$a0,b0=params$b0, rinit=NULL, verbose=FALSE))
  delta <- rowSums(lin_reg_test$z[c(1:(chain_length-1)),]!=lin_reg_test$z[c(2:chain_length),])
  
  time_taken_df <- 
    rbind(time_taken_df, 
          data.frame(algo='fast_spike_slab',
                     time=as.double(fast_time_taken[1])/chain_length, 
                     delta_mean=mean(delta), delta_var=var(delta), n=n, p=p))
  
  print(n)
}

comparison_plot3 <- 
  ggplot(time_taken_df %>% filter(p==100), aes(x=n, y=time*1000, linetype=algo)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('Time per iteration (ms)')) +
  theme_classic(base_size = 18) +
  theme(legend.position="bottom")
comparison_plot3

comparison_plot4 <- 
  ggplot(time_taken_df %>% filter(p==100), aes(x=n, y=delta_mean*n, linetype=algo)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('p_t')) +
  theme_classic(base_size = 18) +
  theme(legend.position="bottom")
comparison_plot4








# n <- 100
# p <- 1000
# s0 <- 20
# lingreg_sync_data <- synthetic_data(n,p,s0)
# params <- spike_slab_params(n,p)
# lin_reg_test <- 
#   spike_slab_mcmc(chain_length=1e3, X=lingreg_sync_data$X,Xt=t(lingreg_sync_data$X),y=lingreg_sync_data$y,
#                     tau0=params$tau0, tau1=params$tau1, q=params$q, 
#                     a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE)
# matplot(lin_reg_test$beta[,c(1:20)], type='l')
# 
# n <- 100
# p <- 1000
# s0 <- 4
# logreg_sync_data <- synthetic_data(n,p,s0,type = 'logistic')
# params <- spike_slab_params(n,p)
# X <- logreg_sync_data$X
# # X <- matrix(scale(logreg_sync_data$X), n, p)
# log_reg_test <- 
#   spike_slab_mcmc(chain_length=1e3, X=X,Xt=t(X),y=logreg_sync_data$y,
#                   tau0=params$tau0, tau1=params$tau1, q=params$q, 
#                   a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE)
# matplot(log_reg_test$beta[,c(1:20)], type='l')



