# Fixed dataset simulations
library(ScaleSpikeSlab)
library(doParallel)
registerDoParallel(cores = detectCores()-2)
library(foreach)

library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)

source('inst/dataset_simulations/dataset_functions.R')
source('inst/comparisons/s3_sota_crn_functions.R')

dataset_name <- 'Gordon'
dataset <- import_dataset(dataset_name)
n <- dataset$n
p <- dataset$p
X <- dataset$X
y <- dataset$y
XXt <- dataset$XXt

# Select hyperparameters
params <- spike_slab_params(n,p)
# Choosing same hyperpamaters as skinnybasad package for SkinnyGibbs
params$tau0 <- 1/sqrt(n)
params$tau1 <- sqrt(max(0.01 * p^{2 + 0.1}/n, 1))
tau0_inverse <- chol2inv(chol(diag(n) + params$tau0^2 * XXt))
tau1_inverse <- chol2inv(chol(diag(n) + params$tau1^2 * XXt))

chain_length <- 5e3
burnin <- 1e3

# chains_probit <-
#   s3_sota_crn_probit(chain_length,X,y,tau0=params$tau0,tau1=params$tau1,q=params$q,
#                      rinit=NULL,verbose=TRUE,burnin=burnin,store=TRUE,Xt=NULL, XXt=NULL,
#                      tau0_inverse=NULL,tau1_inverse=NULL)
# ess_sota_probit <- mcmcse::ess(chains_probit$sota$beta)
# ess_s3_probit <- mcmcse::ess(chains_probit$s3$beta)
# mean(abs(chains_probit$sota$z-chains_probit$s3$z))
# max(abs(chains_probit$sota$beta-chains_probit$s3$beta))
# 
# chains_logistic <-
#   s3_sota_crn_logistic(chain_length,X,y,tau0=params$tau0,tau1=params$tau1,q=params$q,
#                        rinit=NULL,verbose=TRUE,burnin=burnin,store=TRUE,Xt=NULL, XXt=NULL)
# ess_sota_logistic <- mcmcse::ess(chains_logistic$sota$beta)
# ess_s3_logistic <- mcmcse::ess(chains_logistic$s3$beta)
# max(abs(chains_logistic$sota$z-chains_logistic$s3$z))
# max(abs(chains_logistic$sota$beta-chains_logistic$s3$beta))
# 
# # Skinny Gibbs Logistic chain
# z <- rbinom(p,1,params$q)
# sigma2 <- 1/rgamma(1,shape = (params$a0/2), rate = (params$b0/2))
# beta <- rnorm(p)
# beta[z==0] <- beta[z==0]*(params$tau0*sqrt(sigma2))
# beta[z==1] <- beta[z==1]*(params$tau1*sqrt(sigma2))
# skinny_chain_logistic <-
#   skinnybasad(X,y,params$q,beta,z,modif=1,nburn=burnin,
#               niter=(chain_length-burnin),printitrsep=TRUE)
# chains_output <-
#   rbind(data.frame(algo='SOTA_probit', z_avg=I(list(colMeans(chains_probit$sota$z))),
#                    beta_avg=I(list(chains_probit$sota$beta_ergodic_avg)), ess=mean(ess_sota_probit)),
#         data.frame(algo='S3_probit', z_avg=I(list(colMeans(chains_probit$s3$z))),
#                    beta_avg=I(list(chains_probit$s3$beta_ergodic_avg)), ess=mean(ess_s3_probit)),
#         data.frame(algo='SOTA_logistic', z_avg=I(list(colMeans(chains_logistic$sota$z))),
#                    beta_avg=I(list(chains_logistic$sota$beta_ergodic_avg)), ess=mean(ess_sota_logistic)),
#         data.frame(algo='S3_logistic', z_avg=I(list(colMeans(chains_logistic$s3$z))),
#                    beta_avg=I(list(chains_logistic$s3$beta_ergodic_avg)), ess=mean(ess_s3_logistic)),
#         data.frame(algo='SKINNY', z_avg=I(list(skinny_chain_logistic$marZ)),
#                    beta_avg=NA, ess=NA))

nchains <- 5
chains_output <-
  foreach(i = c(1:nchains), .combine = rbind)%dopar%{
    chains_probit <-
      s3_sota_crn_probit(chain_length,X,y,tau0=params$tau0,tau1=params$tau1,q=params$q,
                         rinit=NULL,verbose=TRUE,burnin=burnin,store=FALSE,Xt=NULL, XXt=NULL,
                         tau0_inverse=NULL,tau1_inverse=NULL)
    # ess_sota_probit <- mcmcse::ess(chains_probit$sota$beta)
    # ess_s3_probit <- mcmcse::ess(chains_probit$s3$beta)

    chains_logistic <-
      s3_sota_crn_logistic(chain_length,X,y,tau0=params$tau0,tau1=params$tau1,q=params$q,
                           rinit=NULL,verbose=TRUE,burnin=burnin,store=FALSE,Xt=NULL, XXt=NULL)
    # ess_sota_logistic <- mcmcse::ess(chains_logistic$sota$beta)
    # ess_s3_logistic <- mcmcse::ess(chains_logistic$s3$beta)

    # Skinny Gibbs Logistic chain
    z <- rbinom(p,1,params$q)
    sigma2 <- 1/rgamma(1,shape = (params$a0/2), rate = (params$b0/2))
    beta <- rnorm(p)
    beta[z==0] <- beta[z==0]*(params$tau0*sqrt(sigma2))
    beta[z==1] <- beta[z==1]*(params$tau1*sqrt(sigma2))
    skinny_chain_logistic <-
      skinnybasad(X,y,params$q,beta,z,modif=1,nburn=burnin,
                  niter=(chain_length-burnin),printitrsep=TRUE)
    return(rbind(data.frame(algo='SOTA_probit', iteration=i, z_avg=I(list(chains_probit$sota$z_ergodic_avg)),
                            beta_avg=I(list(chains_probit$sota$beta_ergodic_avg))),
                 data.frame(algo='S3_probit', iteration=i, z_avg=I(list(chains_probit$s3$z_ergodic_avg)),
                            beta_avg=I(list(chains_probit$s3$beta_ergodic_avg))),
                 data.frame(algo='SOTA_logistic', iteration=i, z_avg=I(list(chains_logistic$sota$z_ergodic_avg)),
                            beta_avg=I(list(chains_logistic$sota$beta_ergodic_avg))),
                 data.frame(algo='S3_logistic', iteration=i, z_avg=I(list(chains_logistic$s3$z_ergodic_avg)),
                            beta_avg=I(list(chains_logistic$s3$beta_ergodic_avg))),
                 data.frame(algo='SKINNY', iteration=i, z_avg=I(list(skinny_chain_logistic$marZ)),
                            beta_avg=NA)))
  }



# filename <- paste("inst/dataset_simulations/",dataset_name,"_chains.RData",sep="")
# save(file = filename, chains_output)
# load(filename)

#### Variable selection comparison
variable_section_df <- 
  chains_output %>% 
  group_by(algo) %>%
  summarise(prob_covariate_slab = rowMeans(do.call(cbind, z_avg))) %>%
  mutate(prob_covariate_slab=prob_covariate_slab, covariate_index=1:n()) %>%
  arrange(desc(prob_covariate_slab)) %>%
  mutate(xaxis =1:n())

xlabel <- TeX(paste(dataset_name,'Covariates'))
# label_name <- TeX('Sampler')
# label_breaks <- c('S3_logistic','S3_probit','SKINNY_logistic','SOTA_logistic','SOTA_probit', 'S3_linear', 'SOTA_linear')
# label_names <- unname(TeX(c('$S^3$ Logistic', '$S^3$ Probit', 'Skinny Gibbs Logistic',
#                             'SOTA Logistic', 'SOTA Probit', '$S^3$ Linear', 'SOTA Linear')))
# 
# z_variable_section <- 
#   ggplot(variable_section_df, 
#          aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo, shape=algo)) + 
#   geom_point(size=3, stroke=0.5) + xlab(xlabel) + 
#   ylab(TeX('Marginal posterior probabilities')) + scale_x_continuous(trans='log10')+
#   scale_color_manual(name=label_name, breaks=label_breaks, labels=label_names,
#                      values = c('Black','Black','Black','Gray','Gray')) +
#   scale_shape_manual(name=label_name, breaks = label_breaks,
#                      labels=label_names, 
#                      values = c(4, 1, 3, 4, 1)) +
#   scale_x_continuous(trans='log10') + theme_classic(base_size = 12) +
#   theme(legend.position = 'bottom',legend.text=element_text(size=12))
# # z_variable_section



z_variable_section_probit <-
  ggplot(variable_section_df %>% filter(algo%in%c('S3_probit','SOTA_probit')),
         aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), shape=algo)) +
  geom_point(size=6, color='gray', stroke=1) + xlab(xlabel) +
  ylab(TeX('Marginal posterior probabilities')) +
  scale_shape_manual(name=TeX('Sampler'), breaks = c("S3_probit", "SOTA_probit"),
                     labels=unname(TeX(c('$S^3$ Probit','SOTA Probit'))), values = c(4, 1, 3)) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 22) +
  theme(legend.position = 'bottom',legend.text=element_text(size=24))
# z_variable_section_probit

z_variable_section_logistic <-
  ggplot(variable_section_df %>% filter(algo%in%c('S3_logistic','SOTA_logistic', 'SKINNY')),
         aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), shape=algo)) +
  geom_point(size=6, color='black', stroke=0.5) + xlab(xlabel) +
  ylab(TeX('Marginal posterior probabilities')) +
  scale_shape_manual(name=TeX('Sampler'), breaks = c("S3_logistic", "SOTA_logistic", 'SKINNY'),
                     labels=unname(TeX(c('$S^3$','SOTA', 'Skinny Gibbs Logistic'))), values = c(4, 1, 3)) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 22) +
  theme(legend.position = 'bottom',legend.text=element_text(size=24))
# z_variable_section_logistic


#### Time per iteration plot
load("inst/dataset_simulations/Gordon_sims.RData")
level_order <- 
  factor(multi_datasets_df$algo, 
         level = c('S3_logistic','SOTA_logistic','SKINNY','S3_probit','SOTA_probit'))
label= TeX('Sampler')
label_breaks <- c('S3_logistic','SOTA_logistic','SKINNY','S3_probit','SOTA_probit')
label_names <- unname(TeX(c('$S^3$ Logistic', 'SOTA Logistic', 'Skinny Gibbs Logistic', '$S^3$ Probit','SOTA Probit')))

# Both logistic and probit
time_comparison <- 
  ggplot(multi_datasets_df, aes(x=level_order, y=time_mean*1000, color=algo, shape=algo)) + 
  geom_point(size=6, stroke=1) + 
  xlab(TeX('Sampler')) + ylab(TeX('Time per iteration (ms)')) +
  scale_color_manual(name=label, breaks=label_breaks, labels=label_names,
                     values = c('Black','Black','Black','Gray','Gray')) +
  scale_shape_manual(name=label, breaks=label_breaks, labels=label_names,
                     values = c(4,1,3,4,1)) +
  scale_x_discrete(name=TeX('Sampler'), breaks=label_breaks, labels=label_names) +
  geom_errorbar(aes(ymax=(time_mean+time_sd/no_chains)*1000, ymin=(time_mean-time_sd/no_chains)*1000),
                position=position_dodge(.9), show.legend = FALSE) +
  # scale_y_continuous(limits =c(0,40)) +
  scale_y_continuous(trans='log10') +
  theme_classic(base_size = 22) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), 
        legend.text=element_text(size=24)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, colour = 'white'))
# time_comparison 
shared_legend <- ggpubr::get_legend(time_comparison, position = NULL)

# Combined plot
plot_combined <- 
  ggarrange(z_variable_section_logistic, z_variable_section_probit,
            time_comparison, common.legend=TRUE, legend.grob=shared_legend, 
            legend = "bottom", nrow=1)
# plot_combined
# ggsave(filename = "inst/dataset_simulations/Gordon_plot_main.pdf", plot = plot_combined, width = 12, height = 7)

