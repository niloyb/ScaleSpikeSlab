# Dataset functions for Section 4

import_dataset <- function(dataset_name){
  ######################### Import dataset #########################
  # Lin reg datasets
  if(dataset_name=='Riboflavin'){
    data("riboflavin")
    y <- riboflavin$y
    X <- riboflavin$x
    colnames(X) <- NULL
    rownames(X) <- NULL
  }
  if(dataset_name=='Maize'){
    load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/processed_maize_data/design_matrix_Xnew.RData')
    load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/processed_maize_data/response_ynew.RData')
  }
  # Log reg datasets
  if(dataset_name=='Malware'){
    malware_data <- read.csv('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/uci_malware/uci_malware_detection.csv', header = TRUE)
    y <- as.matrix(rep(0, nrow(malware_data)))
    y[malware_data[,1]=='malicious',] <- 1
    X <- as.matrix(malware_data[,-1])
    rownames(X) <- NULL
    colnames(X) <- NULL
  }
  if(dataset_name=='Lymph'){
    lymph_data <- read.table('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/lymph_data/lymph.dat', header = FALSE)
    y <- lymph_data[,1]
    X <- lymph_data[,-1]
    colnames(X) <- NULL
    rownames(X) <- NULL
  }
  if(dataset_name=='PCR'){
    pcr_data <- read.csv('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/pcr/pcr.csv', header = FALSE)
    # dim(pcr_data)
    X <- t(pcr_data[,-1])
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- read.delim('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/pcr/GSE3330_B6BTBRlivergenotype60mice.txt', header = TRUE)
    y <- y$X.GPAT[(2+1):(2+60)]
    # Creating binary response from continuous GPAT (if want to do logistic)
    # y <- 1*(y<quantile(y,0.4))
  }
  # Synthetic datasets
  if(dataset_name=='Synthetic_Binary'){
    # Just to sense check things work
    syn_data <- synthetic_data(1000,50000,20,2,signal ='decay',type = 'logistic')
    X <- syn_data$X
    y <- syn_data$y
  }
  if(dataset_name=='Synthetic_Continuous'){
    # Just to sense check things work
    syn_data <- synthetic_data(1000,50000,20,2,signal ='decay',type = 'linear')
    X <- syn_data$X
    y <- syn_data$y
  }
  
  # Logistic regression datasets from datamicroarray package
  if(dataset_name=='Borovecki'){
    data('borovecki', package = 'datamicroarray')
    X <- borovecki$x
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- borovecki$y
    binary_y <- rep(NA,nrow(X))
    binary_y[y==unique(as.factor(y))[1]] <- 0
    binary_y[y==unique(as.factor(y))[2]] <- 1
    y <- binary_y
  }
  if(dataset_name=='Chin'){
    data('chin', package = 'datamicroarray')
    X <- chin$x
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- chin$y
    binary_y <- rep(NA,nrow(X))
    binary_y[y==unique(as.factor(y))[1]] <- 0
    binary_y[y==unique(as.factor(y))[2]] <- 1
    y <- binary_y
  }
  if(dataset_name=='Chowdary'){
    data('chowdary', package = 'datamicroarray')
    X <- chowdary$x
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- chowdary$y
    binary_y <- rep(NA,nrow(X))
    binary_y[y==unique(as.factor(y))[1]] <- 0
    binary_y[y==unique(as.factor(y))[2]] <- 1
    y <- binary_y
  }
  if(dataset_name=='Gordon'){
    data('gordon', package = 'datamicroarray')
    X <- gordon$x
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- gordon$y
    binary_y <- rep(NA,nrow(X))
    binary_y[y==unique(as.factor(y))[1]] <- 0
    binary_y[y==unique(as.factor(y))[2]] <- 1
    y <- binary_y
  }
  
  # Removing co-variates which do not vary across observations
  col_sds <- sapply(c(1:dim(X)[2]), function(i){sd(X[,i])})
  X <- X[,col_sds!=0]
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(scale(X),n,p)
  Xt <- t(X)
  XXt <- X%*%Xt
  
  if(length(unique(y))==2){
    dataset_type <- 'binary'
  } else{
    dataset_type <- 'continuous'
  }
  
  return(list(X=X,Xt=Xt,XXt=XXt,y=y,n=n,p=p,type=dataset_type,name=dataset_name))
}

comparison_dataset_sims <- 
  function(dataset_name,chain_length=5e3,burnin=1e3,no_chains=1,
           algos=c('S3_logistic','S3_probit', 'S3_linear', 
                   'SOTA_logistic', 'SOTA_probit','SOTA_linear',
                   'SKINNY')){
    dataset <- import_dataset(dataset_name)
    X <- dataset$X
    Xt <- dataset$Xt
    XXt <- dataset$XXt
    y <- dataset$y
    n <- dataset$n
    p <- dataset$p
    
    # Select hyperparameters
    params <- spike_slab_params(n, p)
    # Choosing same hyperpamaters as skinnybasad package for SkinnyGibbs
    params$tau0 <- 1/sqrt(n)
    params$tau1 <- sqrt(max(0.01 * p^{2 + 0.1}/n, 1))
    tau0_inverse <- chol2inv(chol(diag(n) + params$tau0^2 * XXt))
    tau1_inverse <- chol2inv(chol(diag(n) + params$tau1^2 * XXt))
    
    foreach(i = c(1:no_chains), .combine = rbind)%dopar%{
      output <- data.frame()
      if('S3_linear' %in% algos){
        ###### Scalable spike and slab
        sss_time_taken <-
          system.time(
            sss_chain <- 
              spike_slab_linear(chain_length=chain_length, X=X, y=y,
                                tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE,
                                burnin=burnin, store=FALSE, Xt=Xt, XXt=XXt,
                                tau0_inverse=tau0_inverse, tau1_inverse=tau1_inverse)
          )
        output <- 
          rbind(output, 
                data.frame(algo='S3_linear',
                           time=as.double(sss_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sss_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i, 
                           burnin=burnin, chain_length=chain_length))
      }
      if('S3_logistic' %in% algos){
        ###### Scalable spike and slab
        sss_time_taken <-
          system.time(
            sss_chain <- 
              spike_slab_logistic(chain_length=chain_length, X=X, y=y,
                                  tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                  rinit=NULL, verbose=TRUE,
                                  burnin=burnin, store=FALSE, Xt=Xt, XXt=XXt)
          )
        output <- 
          rbind(output, 
                data.frame(algo='S3_logistic',
                           time=as.double(sss_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sss_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      if('S3_probit' %in% algos){
        ###### Scalable spike and slab
        sss_time_taken <-
          system.time(
            sss_chain <- 
              spike_slab_probit(chain_length=chain_length, X=X, y=y,
                                tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                rinit=NULL, verbose=TRUE,
                                burnin=burnin, store=FALSE, Xt=Xt, XXt=XXt,
                                tau0_inverse=tau0_inverse, tau1_inverse=tau1_inverse)
          )
        output <- 
          rbind(output, 
                data.frame(algo='S3_probit',
                           time=as.double(sss_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sss_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      if('SOTA_logistic' %in% algos){
        ###### Scalable spike and slab
        sota_time_taken <-
          system.time(
            sota_chain <- 
              sota_spike_slab_logistic(chain_length=chain_length, X=X, y=y,
                                       tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                       rinit=NULL, verbose=TRUE,
                                       burnin=burnin, store=FALSE, Xt=Xt, XXt=XXt)
          )
        output <- 
          rbind(output, 
                data.frame(algo='SOTA_logistic',
                           time=as.double(sota_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sota_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      
      if('SOTA_probit' %in% algos){
        ###### Scalable spike and slab
        sota_time_taken <-
          system.time(
            sota_chain <- 
              sota_spike_slab_probit(chain_length=chain_length, X=X, y=y,
                                     tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                     rinit=NULL, verbose=TRUE,
                                     burnin=burnin, store=FALSE, Xt=Xt)
          )
        output <- 
          rbind(output, 
                data.frame(algo='SOTA_probit',
                           time=as.double(sota_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sota_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      if('SOTA_linear' %in% algos){
        ###### Scalable spike and slab
        sss_time_taken <-
          system.time(
            sss_chain <- 
              sota_spike_slab_linear(chain_length=chain_length, X=X, y=y,
                                     tau0=params$tau0, tau1=params$tau1, q=params$q, 
                                     a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE,
                                     burnin=burnin, store=FALSE, Xt=Xt)
          )
        output <- 
          rbind(output, 
                data.frame(algo='SOTA_linear',
                           time=as.double(sss_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(sss_chain$z_ergodic_avg))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      if('SKINNY' %in% algos){
        ###### Skinny Gibbs Logistic chain
        z <- rbinom(p,1,params$q)
        sigma2 <- 1/rgamma(1,shape = (params$a0/2), rate = (params$b0/2))
        beta <- rnorm(p)
        beta[z==0] <- beta[z==0]*(params$tau0*sqrt(sigma2))
        beta[z==1] <- beta[z==1]*(params$tau1*sqrt(sigma2))
        
        skinny_time_taken <-
          system.time(
            skinny_chain_logistic <-
              skinnybasad(X,y,params$q,beta,z,modif=1,nburn=burnin,
                          niter=(chain_length-burnin),printitrsep=TRUE)
          )
        output <- 
          rbind(output, 
                data.frame(algo='SKINNY',
                           time=as.double(skinny_time_taken[1])/chain_length, 
                           z_ergodic_avg=I(list(as.vector(skinny_chain_logistic$marZ))),
                           dataset=dataset_name, n=n, p=p, iteration=i,
                           burnin=burnin, chain_length=chain_length))
      }
      
      print(c(dataset_name,i))
      return(output)
    }
  }



