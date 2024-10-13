# ESS calculations
rm(list=ls())
library(ScaleSpikeSlab)
library(doParallel)
registerDoParallel(cores = detectCores()-2)
library(foreach)
library(mcmcse)

library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)

source('inst/dataset_simulations/dataset_functions.R')
source('inst/comparisons/comparison_functions.R')

# datasets_list_cont <- c('PCR',"Maize",'Synthetic_Continuous')
# datasets_list_binary <- c('Malware','Lymph','Synthetic_Binary','Borovecki','Chin','Chowdary','Gordon')

datasets_list <- c('Synthetic_Continuous','Synthetic_Binary',"Maize")
sss_chain_length <- 5e3
sss_burnin <- 1e3
sota_chain_length <- 1e2
sota_burnin <- 1e1

ess_df <- data.frame()
for(dataset_name in datasets_list){
  dataset <- import_dataset(dataset_name)
  n <- dataset$n
  p <- dataset$p
  X <- dataset$X
  y <- dataset$y
  Xt <- dataset$Xt
  XXt <- dataset$XXt
  
  # Select hyperparameters
  params <- spike_slab_params(n,p)
  # Choosing same hyperpamaters as skinnybasad package for SkinnyGibbs
  params$tau0 <- 1/sqrt(n)
  params$tau1 <- sqrt(max(0.01 * p^{2 + 0.1}/n, 1))
  tau0_inverse <- chol2inv(chol(diag(n) + params$tau0^2 * XXt))
  tau1_inverse <- chol2inv(chol(diag(n) + params$tau1^2 * XXt))
  
  if(dataset$type=='binary'){
    sss_probit_time_taken <-
      system.time(
        sss_probit <-
          spike_slab_probit(sss_chain_length, X, y, params$tau0, params$tau1, params$q, verbose = TRUE, 
                            burnin = sss_burnin, store = TRUE, Xt = Xt, XXt = XXt, tau0_inverse = tau0_inverse, 
                            tau1_inverse = tau1_inverse)
      )
    sss_logistic_time_taken <-
      system.time(
        sss_logistic <-
          spike_slab_logistic(sss_chain_length, X, y, params$tau0, params$tau1, params$q, verbose = TRUE, 
                              burnin = sss_burnin, store = TRUE, Xt = Xt, XXt = XXt)
      )
    ess_s3_probit_beta <- mcmcse::ess(sss_probit$beta)
    ess_s3_logistic_beta <- mcmcse::ess(sss_logistic$beta)
    
    
    sota_probit_time_taken <-
      system.time(
        sota_probit <-
          sota_spike_slab_probit(sota_chain_length,X,y,params$tau0,params$tau1,params$q,burnin=sota_burnin,Xt=Xt)
      )
    sota_logistic_time_taken <-
      system.time(
        sota_logistic <-
          sota_spike_slab_logistic(sota_chain_length,X,y,params$tau0,params$tau1,params$q,
                                   burnin=sota_burnin,store=TRUE,Xt=Xt,XXt=XXt)
      )
    ess_sota_probit_beta <- mcmcse::ess(sota_probit$beta)
    ess_sota_logistic_beta <- mcmcse::ess(sota_logistic$beta)
    
    ess_s3_probit_df <- data.frame(dataset=dataset_name, algo='S3_probit', ess_beta = mean(ess_s3_probit_beta), chain_length=sss_chain_length, burnin=sss_burnin,
                                time=as.double(sss_probit_time_taken[1])/sss_chain_length)
    ess_s3_logistic_df <- data.frame(dataset=dataset_name, algo='S3_logistic', ess_beta = mean(ess_s3_logistic_beta), chain_length=sss_chain_length, burnin=sss_burnin,
                                  time=as.double(sss_logistic_time_taken[1])/sss_chain_length)
    ess_sota_probit_df <- data.frame(dataset=dataset_name, algo='SOTA_probit', ess_beta = mean(ess_sota_probit_beta), chain_length=sota_chain_length, burnin=sota_burnin,
                                   time=as.double(sota_probit_time_taken[1])/sota_chain_length)
    ess_sota_logistic_df <- data.frame(dataset=dataset_name, algo='SOTA_logistic', ess_beta = mean(ess_sota_logistic_beta), chain_length=sota_chain_length, burnin=sota_burnin,
                                     time=as.double(sota_logistic_time_taken[1])/sota_chain_length)
    ess_df <- rbind(ess_df, ess_s3_probit_df, ess_s3_logistic_df, ess_sota_probit_df, ess_sota_logistic_df)
  }
  if(dataset$type=='continuous'){
    sss_linear_time_taken <-
      system.time(
        sss_linear <-
          spike_slab_linear(sss_chain_length, X, y, params$tau0, params$tau1, params$q, verbose = TRUE, 
                            burnin = sss_burnin, store = TRUE, Xt = Xt, XXt = XXt, tau0_inverse = tau0_inverse, 
                            tau1_inverse = tau1_inverse)
      )
    ess_s3_linear_beta <- mcmcse::ess(sss_linear$beta)
    sota_linear_time_taken <-
      system.time(
        sota_linear <- 
          sota_spike_slab_linear(sota_chain_length,X,y,params$tau0,params$tau1,params$q, 
                                 verbose=TRUE,burnin=sota_burnin,store=TRUE,Xt=Xt)
      )
    ess_sota_linear_beta <- mcmcse::ess(sota_linear$beta)
    
    ess_s3_linear_df <- data.frame(dataset=dataset_name, algo='S3_linear', ess_beta = mean(ess_s3_linear_beta), chain_length=sss_chain_length, burnin=sss_burnin,
                                   time=as.double(sss_linear_time_taken[1])/sss_chain_length)
    ess_sota_linear_df <- data.frame(dataset=dataset_name, algo='SOTA_linear', ess_beta = mean(ess_sota_linear_beta), chain_length=sota_chain_length, burnin=sota_burnin,
                                   time=as.double(sota_linear_time_taken[1])/sota_chain_length)
    ess_df <- rbind(ess_df, ess_s3_linear_df, ess_sota_linear_df)
  }
  print(dataset_name)
}

filename <- paste("inst/dataset_simulations/ESS_sims.RData",sep = '')
# save(file = filename, ess_df)
# load(filename)

# ess_df[ess_df$dataset=='Synthetic_Continuous',]$dataset <- "Synthetic Continuous"
# ess_df[ess_df$dataset=='Synthetic_Binary',]$dataset <- "Synthetic Binary"

s3_ess_plot <-
  ggplot(ess_df %>% filter(algo %in% c('S3_linear','S3_logistic','S3_probit')), aes(x=dataset, y=ess_beta/(chain_length-burnin), shape=algo)) +
  geom_point(size=3) + xlab(TeX('Dataset')) + ylab(TeX('ESS per iteration of $$\\beta$$')) +
  scale_shape_manual(name=TeX('Sampler'), 
                     breaks=c('S3_linear','S3_logistic','S3_probit'), 
                     labels=unname(TeX(c('$S^3$ Linear', '$S^3$ Logistic', '$S^3$ Probit'))), values = c(3,4,1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  guides(shape=guide_legend(nrow=2,ncol=3,byrow=TRUE)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# s3_ess_plot

label_title <- TeX('')
label_breaks <- c('S3_linear','S3_logistic','S3_probit','SOTA_linear','SOTA_logistic','SOTA_probit')
label_names <- unname(TeX(c('$S^3$ Linear', '$S^3$ Logistic', '$S^3$ Probit', 'SOTA Linear', 'SOTA Logistic', 'SOTA Probit')))
s3_sota_ess_plot <-
  ggplot(ess_df, aes(x=dataset, y=ess_beta/time, shape=algo, color=algo)) +
  geom_point(size=3) + xlab(TeX('Dataset')) + ylab(TeX('ESS per second of $$\\beta$$')) +
  scale_shape_manual(name=label_title, breaks=label_breaks, labels=label_names, values = c(3,4,1,3,4,1)) +
  scale_color_manual(name=label_title, breaks=label_breaks, labels=label_names, values = c('black','black','black','gray','gray','gray')) +
  scale_y_continuous(trans='log10') +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# s3_sota_ess_plot


# Combined plot
plot_combined <- ggarrange(s3_ess_plot, s3_sota_ess_plot, legend = "bottom", nrow=1)
# plot_combined
# ggsave(filename = "inst/dataset_simulations/ess_combined_plot.pdf", plot = plot_combined, width = 10, height = 4)


