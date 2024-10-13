########### Functions for k-fold cross validation prediction errors ###########
library(dplyr)

# Predict
predict_response_over_time_cumean <- function(X,beta_samples,type='linear'){
  traj_length <- dim(beta_samples)[1]
  if(type=='linear'){
    # Output predictive mean of response
    # return(as.vector(X%*%colMeans(beta_samples)))
    return(apply(X%*%t(beta_samples),1,dplyr::cummean))
  }
  if(type=='logistic'){
    # Output predictive probability of binary outcome
    # plogis is the logistic function
    # return(rowMeans(plogis(X%*%t(beta_samples))))
    return(apply(plogis(X%*%t(beta_samples)),1,dplyr::cummean))
  }
  if(type=='probit'){
    # Output predictive probability of binary outcome
    # pnorm is the std normal cdf
    # return(rowMeans(pnorm(X%*%t(beta_samples))))
    return(apply(pnorm(X%*%t(beta_samples)),1,dplyr::cummean))
  }
}

# Performance
mse_over_time <- function(predicted_response,true_response){
  return(rowMeans(sweep(predicted_response,2,true_response)^2))
}

kfold_rmse <-
  function(X,y,k=10,type='linear',chain_length=5000,burnin=1e3, no_chains=1, verbose=FALSE){
    n <- length(y)
    p <- dim(X)[2]
    #Creating folds
    folds <- matrix(sample(seq(1,n,1), size = n)[1:(floor(n/k)*k)], nrow=k)

    kfold_df <- data.frame()
    for(i in c(1:k)){
      test_ind=folds[i,]
      n_train <- n-length(test_ind)
      params <- spike_slab_params(n_train,p)
      # Choosing same hyperpamaters as skinnybasad package for SkinnyGibbs
      params$tau0 <- 1/sqrt(n_train)
      params$tau1 <- sqrt(max(0.01 * p^{2 + 0.1}/n_train, 1))
      tau0 <- params$tau0
      tau1 <- params$tau1
      q <- params$q

      sss_beta <-
        foreach(chain = c(1:no_chains), .combine = '+')%dopar%{
          if(type=='linear'){
            sss_chain <-
              spike_slab_linear(chain_length, X[-test_ind,,drop=FALSE], y[-test_ind,drop=FALSE], tau0, tau1, q, burnin=burnin, verbose=verbose)
          }
          if(type=='logistic'){
            sss_chain <-
              spike_slab_logistic(chain_length, X[-test_ind,,drop=FALSE], y[-test_ind,drop=FALSE], tau0, tau1, q, burnin=burnin, verbose=verbose)
          }
          if(type=='probit'){
            sss_chain <-
              spike_slab_probit(chain_length, X[-test_ind,,drop=FALSE], y[-test_ind,drop=FALSE], tau0, tau1, q, burnin=burnin, verbose=verbose)
          }
          return(sss_chain$beta/no_chains)
        }

      y_pred <- predict_response_over_time_cumean(X[test_ind,,drop=FALSE],sss_beta,type=type)
      rmse <- mse_over_time(y_pred,y[test_ind,drop=FALSE])^0.5
      kfold_df <- rbind(kfold_df,rmse)
      print(i)
    }
    return(kfold_df)
  }


