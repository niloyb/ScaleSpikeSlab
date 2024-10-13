# Simulation results for different datasets
rm(list=ls())
library(ScaleSpikeSlab)
source('inst/comparisons/comparison_functions.R')
source('inst/dataset_simulations/dataset_functions.R')

library(doParallel)
library(foreach)
# no_cores <- 2
no_cores <- detectCores()-1
registerDoParallel(cores = no_cores)


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

############### Plot of runtime comparison of different datasets ###############
# datasets_list_cont <- c('PCR',"Maize",'Synthetic_Continuous')
# datasets_list_binary <- c('Malware','Lymph','Synthetic_Binary','Borovecki','Chin','Chowdary','Gordon')

datasets_list_cont <- c()
datasets_list_binary <- c('Synthetic_Binary')

sss_chain_length <- 5e3
sss_burnin <- 1e3
skinny_chain_length <- 1e3
skinny_burnin <- 1e2
sota_chain_length <- 1e1
sota_burnin <- 0
no_chains <- 1

dataset_output_cont <- data.frame()
for(dataset_name in datasets_list_cont){
  sss_output <- 
    comparison_dataset_sims(dataset_name,chain_length=sss_chain_length,
                            burnin=sss_burnin, no_chains=no_chains, 
                            algos=c('S3_linear'))
  sota_output <- 
    comparison_dataset_sims(dataset_name,chain_length=sota_chain_length,
                            burnin=0, no_chains=no_chains,
                            algos=c('SOTA_linear'))
  dataset_output_cont <- rbind(dataset_output_cont,sss_output,sota_output)
  print(dataset_name)
}

dataset_output_binary <- data.frame()
for(dataset_name in datasets_list_binary){
  sss_output <- 
    comparison_dataset_sims(dataset_name,chain_length=sss_chain_length,
                            burnin=sss_burnin, no_chains=no_chains, 
                            algos=c('S3_probit','S3_logistic'))
  skinny_output <- 
    comparison_dataset_sims(dataset_name,chain_length=skinny_chain_length,
                            burnin=skinny_burnin, no_chains=no_chains, 
                            algos=c('SKINNY'))
  sota_output <- 
    comparison_dataset_sims(dataset_name,chain_length=sota_chain_length,
                            burnin=0, no_chains=no_chains,
                            algos=c('SOTA_probit','SOTA_logistic'))
  dataset_output_binary <- rbind(dataset_output_binary,sss_output,skinny_output,sota_output)
  print(dataset_name)
}

 datasets_output <- rbind(dataset_output_cont,dataset_output_binary)

multi_datasets_df <- 
  datasets_output %>% 
  group_by(dataset, algo) %>% 
  select(time, z_ergodic_avg, n, p) %>%
  summarise(time_mean = mean(time), time_sd = sd(time),
            z_avg=I(list(colMeans(do.call(rbind, z_ergodic_avg)))),
            n=mean(n), p=mean(p), no_chains=n()) 
# %>% arrange((n^2*p))

# filename <- paste("inst/dataset_simulations/datasets_sims.RData",sep = '')
filename <- paste("inst/dataset_simulations/",dataset_name,"_sims.RData",sep = '')
# save(file = filename, multi_datasets_df)











