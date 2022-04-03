# Logistic regression dataset example
rm(list=ls())
library(ScaleSpikeSlab)
source('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/comparison_functions.R')
# source('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScalableSpikeSlab/R/helper_functions.R')

library(doParallel)
registerDoParallel(cores = detectCores()-1)
library(foreach)

# Comparison of run-time for different datasets

# library(devtools)
# install_github('ramhiser/datamicroarray')
library(datamicroarray)
# datamicroarray_datasets <- c("borovecki","chin","chowdary","gordon")
# describe_data()

library(dplyr)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(ggridges)
library(ggpubr)

comparison_dataset_sims <- 
  function(dataset,chain_length=1e4,burnin=5e3,no_chains=1,
           algos=c('S3_logistic','S3_probit', 'S3_linear', 
                   'SOTA_logistic', 'SOTA_probit','SOTA_linear')){
  ######################### Import dataset #########################
  # Lin reg datasets
  if(dataset=='Riboflavin'){
    data("riboflavin")
    y <- riboflavin$y
    X <- riboflavin$x
    colnames(X) <- NULL
    rownames(X) <- NULL
  }
  if(dataset=='Maize'){
    load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/Pierre Jacob Lab/misc/high_dim_datasets/processed_maize_data/design_matrix_Xnew.RData')
    load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/Pierre Jacob Lab/misc/high_dim_datasets/processed_maize_data/response_ynew.RData')
  }
  # Log reg datasets
  if(dataset=='Malware'){
    malware_data <- read.csv('/Users/niloybiswas/Downloads/uci_malware_detection.csv', header = TRUE)
    y <- as.matrix(rep(0, nrow(malware_data)))
    y[malware_data[,1]=='malicious',] <- 1
    X <- as.matrix(malware_data[,-1])
    rownames(X) <- NULL
    colnames(X) <- NULL
  }
  if(dataset=='Lymph'){
    lymph_data <- read.table('/Users/niloybiswas/Dropbox/Apps/Overleaf/scalable_spike_slab/datasets/lymph.dat', header = FALSE)
    y <- lymph_data[,1]
    X <- lymph_data[,-1]
    colnames(X) <- NULL
    rownames(X) <- NULL
  }
  if(dataset=='PCR'){
    pcr_data <- read.csv('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/skinny_gibbs/pcr.csv', header = FALSE)
    # dim(pcr_data)
    X <- t(pcr_data[,-1])
    rownames(X) <- NULL
    colnames(X) <- NULL
    y <- read.delim('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/skinny_gibbs/GSE3330_B6BTBRlivergenotype60mice.txt', header = TRUE)
    y <- y$X.GPAT[(2+1):(2+60)]
    # Creating binary response from continuous GPAT (if want to do logistic)
    # y <- 1*(y<quantile(y,0.4))
  }
  # Synthetic datasets
  if(dataset=='Synthetic_Binary'){
    # Just to sense check things work
    syn_data <- synthetic_data(1000,50000,20,2,signal ='decay',type = 'logistic')
    X <- syn_data$X
    y <- syn_data$y
  }
  if(dataset=='Synthetic_Continuous'){
    # Just to sense check things work
    syn_data <- synthetic_data(1000,50000,20,2,signal ='decay',type = 'linear')
    X <- syn_data$X
    y <- syn_data$y
  }
  
  # Logistic regression datasets from datamicroarray package
  if(dataset=='Borovecki'){
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
  if(dataset=='Chin'){
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
  if(dataset=='Chowdary'){
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
  if(dataset=='Gordon'){
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
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Removing co-variates which do not vary across observations
  col_sds <- sapply(c(1:p), function(i){sd(X[,i])})
  X <- X[,col_sds!=0]
  X <- as.matrix(scale(X),n,p)
  Xt <- t(X)
  XXt <- X%*%Xt
  
  # Select hyperparameters
  params <- spike_slab_params(n, p)
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
                         dataset=dataset, n=n, p=p, iteration=i))
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
                         dataset=dataset, n=n, p=p, iteration=i))
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
                         dataset=dataset, n=n, p=p, iteration=i))
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
                         dataset=dataset, n=n, p=p, iteration=i))
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
                         dataset=dataset, n=n, p=p, iteration=i))
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
                         dataset=dataset, n=n, p=p, iteration=i))
    }
    print(c(dataset,i))
    return(output)
    }
}

######################## Malware Dataset Simulations ###########################
dataset <- 'Malware'
no_chains <- 10
sss_chain_length <- 1e4
sss_burnin <- 1e3
sss_output <- 
  comparison_dataset_sims(dataset,chain_length=sss_chain_length,burnin=sss_burnin,
                          no_chains=no_chains, 
                          algos=c('S3_logistic','S3_probit'))

sota_chain_length <- 1e3
sota_burnin <- 1e2
sota_output <- 
  comparison_dataset_sims(dataset,chain_length=sota_chain_length,burnin=sota_burnin,
                          no_chains=no_chains, 
                          algos=c('SOTA_logistic','SOTA_probit'))

sss_sota_comparison <- rbind(sss_output, sota_output)

sss_sota_comparison_df <- 
  sss_sota_comparison %>% 
  group_by(algo) %>% 
  select(time, z_ergodic_avg, n, p, iteration) %>%
  summarise(time_mean = mean(time), time_sd = sd(time),
            z_avg=I(list(colMeans(do.call(rbind, z_ergodic_avg)))),
            dataset,n=n, p=p)

filename1 <- paste("/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/dataset_simulations/",dataset,'_sims.RData',sep = '')
# save(file = filename1, sss_sota_comparison_df)

############### Plot of runtime comparison of different datasets ###############
datasets_list_cont <- c('PCR',"Maize",'Synthetic_Continuous')
datasets_list_binary <- c('Malware','Lymph','Synthetic_Binary','Borovecki','Chin','Chowdary','Gordon')
sss_chain_length <- 1e2
sota_chain_length <- 1e1
no_chains <- 5

time_comparison_output_linear <- foreach(dataset = datasets_list_cont, .combine = rbind)%dopar%{
  sss_output <- 
    comparison_dataset_sims(dataset,chain_length=sss_chain_length,
                            burnin=0, no_chains=no_chains, 
                            algos=c('S3_linear'))
  sota_output <- 
    comparison_dataset_sims(dataset,chain_length=sota_chain_length,
                            burnin=0, no_chains=no_chains,
                            algos=c('SOTA_linear'))
  sss_sota_comparison <- rbind(sss_output,sota_output)
  return(sss_sota_comparison)
}
time_comparison_output_probit <- foreach(dataset = datasets_list_binary, .combine = rbind)%dopar%{
  sss_output <- 
    comparison_dataset_sims(dataset,chain_length=sss_chain_length,
                            burnin=0, no_chains=no_chains, 
                            algos=c('S3_probit'))
  sota_output <- 
    comparison_dataset_sims(dataset,chain_length=sota_chain_length,
                            burnin=0, no_chains=no_chains,
                            algos=c('SOTA_probit'))
  sss_sota_comparison <- rbind(sss_output,sota_output)
  return(sss_sota_comparison)
}
time_comparison_output <- 
  rbind(time_comparison_output_linear,time_comparison_output_probit)

sss_sota_multi_dataset_time_comparison_df <- 
  time_comparison_output %>% 
  group_by(dataset, algo) %>% 
  select(time, n, p) %>%
  summarise(time_mean = mean(time), time_sd = sd(time),n=mean(n), p=mean(p), no_chains=n()) %>% 
  arrange((n^2*p))

filename2 <- paste("/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/dataset_simulations/multiple_dataset_sims_new.RData",sep = '')
# save(file = filename2, sss_sota_multi_dataset_time_comparison_df)



time_comparison <- 
  ggplot(sss_sota_multi_dataset_time_comparison_df, aes(x=dataset, y=time_mean*1000, color=algo)) + 
  geom_point(size=4, shape=4) + xlab(TeX('Datasets')) + ylab(TeX('Time per iteration (ms)')) +
  scale_x_discrete(limits=sss_sota_multi_dataset_time_comparison_df$dataset) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),labels=unname(TeX(c('S^3', 'SOTA'))),values = c('Black', 'Gray')) +
  # geom_errorbar(aes(ymax=(time_mean+time_sd)*1000, ymin=(time_mean-time_sd)*1000),position=position_dodge(.9)) +
  scale_y_continuous(trans='log10') + theme_classic(base_size = 9) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
time_comparison















