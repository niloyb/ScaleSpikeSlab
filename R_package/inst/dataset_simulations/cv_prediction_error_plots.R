rm(list=ls())
library(ScaleSpikeSlab)
library(doParallel)
library(foreach)
# no_cores <- 2
no_cores <- detectCores()-1
registerDoParallel(cores = no_cores)

library(ScaleSpikeSlab)
# setwd('github/ScaleSpikeSlab/R_package/')
source('inst/dataset_simulations/cv_prediction_error_functions.R')
source('inst/comparisons/comparison_functions.R')
source('inst/dataset_simulations/dataset_functions.R')

############### Plot of 10-fold prediction error of different datasets ###############
# datasets_list_cont <- c('PCR',"Maize",'Synthetic_Continuous')
# datasets_list_binary <- c('Malware','Lymph','Synthetic_Binary','Borovecki','Chin','Chowdary','Gordon')




dataset_name <- 'Synthetic_Binary'
dataset <- import_dataset(dataset_name)

sss_chain_length <- 5e3
sss_burnin <- 1e3
verbose <- TRUE

if(dataset$type=='binary'){
  kfold_output_probit <- kfold_rmse(dataset$X,dataset$y,k=5,type='probit',chain_length=sss_chain_length,burnin=sss_burnin,verbose=verbose)
  sss_probit_df <- data.frame(dataset=dataset_name,algo='S3_probit',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_probit)))
  sota_probit_df <- data.frame(dataset=dataset_name,algo='SOTA_probit',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_probit)))
  
  kfold_output_logistic <- kfold_rmse(dataset$X,dataset$y,k=5,type='logistic',chain_length=sss_chain_length,burnin=sss_burnin,verbose=verbose)
  sss_logistic_df <- data.frame(dataset=dataset_name,algo='S3_logistic',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_logistic)))
  sota_logistic_df <- data.frame(dataset=dataset_name,algo='SOTA_logistic',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_logistic)))
  rmse_df <- rbind(sss_probit_df, sota_probit_df, sss_logistic_df, sota_logistic_df)
}
if(dataset$type=='continuous'){
  kfold_output_linear <- 
    kfold_rmse(dataset$X,dataset$y,k=5,type='linear',chain_length=sss_chain_length,burnin=sss_burnin,verbose=verbose)
  sss_linear_df <- data.frame(dataset=dataset_name,algo='S3_linear',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_linear)))
  sota_linear_df <- data.frame(dataset=dataset_name,algo='SOTA_linear',iteration=c((sss_burnin+1):sss_chain_length),rmse=as.vector(colMeans(kfold_output_linear)))
  rmse_df <- rbind(sss_linear_df,sota_linear_df)
}

# plot(colMeans(kfold_output_probit))
# plot(colMeans(kfold_output_logistic))
# plot(colMeans(kfold_output_linear))

kfold_cv_filename <- paste('inst/dataset_simulations/',dataset_name,'_10fold_cv_rmse.Rdata',sep = '')
# save(file = kfold_cv_filename, rmse_df)
# load(kfold_cv_filename)






