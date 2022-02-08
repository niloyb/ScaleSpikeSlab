
rm(list=ls())
load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/FastSpikeSlab/inst/dataset_simulations/riboflavin.RData')

library(doParallel)
registerDoParallel(cores = detectCores()-1)
library(foreach)
library(ggplot2)
library(latex2exp)

# uci_malware_detection.csv
y <- riboflavin$y
X <- riboflavin$x
Xt <- t(X)
X <- matrix(scale(X), n, p)
n <- nrow(X)
p <- ncol(X)
params <- spike_slab_params(n=n,p=p)

no_chains <- 20
sss_test <- 
  foreach(i = c(1:no_chains), .combine=rbind)%dopar%{
  sss_chain <- spike_slab_mcmc(chain_length=1e4,burnin=1e3,X=X,Xt=Xt,y=y,
                               tau0=params$tau0,tau1=1,q=params$q,
                               verbose=TRUE,store=FALSE)
  return(as.vector(sss_chain$z_ergodic_avg))
}


riboflavin_df <- 
  data.frame(post_prob_mean=apply(sss_test,2,mean),
             post_prob_sd=apply(sss_test,2,sd),
             cov_index=c(1:p), no_chains=no_chains) %>%
  arrange(desc(post_prob_mean)) %>%
  mutate(xaxis =1:n())


ggplot(riboflavin_df, aes(x=xaxis, y=post_prob_mean)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymax=(post_prob_mean+post_prob_sd/sqrt(no_chains)), 
                    ymin=(post_prob_mean-post_prob_sd/sqrt(no_chains))),
                position=position_dodge(.9)) +
  xlab('Malware Covariates') + 
  ylab(TeX('Marginal posterior probabilities')) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 9)







plot(sort(sss_test, decreasing = TRUE), log='x')


plot(sort(sss_chain$z_ergodic_avg, decreasing = TRUE), log='x')

plot(sort(colMeans(sss_chain$z), decreasing = TRUE), log='x')


