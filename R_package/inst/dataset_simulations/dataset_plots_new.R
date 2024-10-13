################################################################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)

######## Multiple datasets plot ########
datasets_list_cont <- c('PCR',"Maize",'Synthetic_Continuous')
datasets_list_binary <- c('Malware','Lymph','Synthetic_Binary','Borovecki','Chin','Chowdary','Gordon')

######## Variable selection df ########
dataset_name <- 'Synthetic_Binary'
dataset_sims <- paste("inst/dataset_simulations/",dataset_name,"_sims.RData",sep = '')
# Loads multi_datasets_df
load(dataset_sims)
variable_section_df <- 
  multi_datasets_df %>% filter(algo%in%c('S3_linear','S3_logistic','S3_probit','SKINNY')) %>% 
  group_by(algo) %>%
  summarise(prob_covariate_slab = rowMeans(do.call(cbind, z_avg))) %>%
  mutate(prob_covariate_slab=prob_covariate_slab, covariate_index=1:n()) %>%
  arrange(desc(prob_covariate_slab)) %>%
  mutate(xaxis =1:n())

######## k-fold CV RMSE df ########
dataset_10fold_rmse <- paste('inst/dataset_simulations/',dataset_name,'_10fold_cv_rmse.Rdata',sep = '')
# Loads rmse_df
load(dataset_10fold_rmse)
if(dataset_name %in% datasets_list_binary){dataset_response_type <- 'binary'}
if(dataset_name %in% datasets_list_cont){dataset_response_type <- 'continuous'}
xlabel <- TeX(paste(dataset_name,'Covariates'))

# filename2 <- paste("inst/dataset_simulations/datasets_sims.RData",sep = '')
# load(filename2)
# unique(multi_datasets_df$dataset)
# dataset_name <- "Lymph"
# dataset_response_type <- 'binary'
# dataset_response_type <- 'continuous'

if(dataset_response_type == 'binary'){
  z_variable_section <- 
    ggplot(variable_section_df %>% filter(algo%in%c('S3_logistic','S3_probit','SKINNY')), 
           aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), shape=algo)) + 
    geom_point(size=3, color='black', stroke=0.5) + xlab(xlabel) + 
    ylab(TeX('Marginal posterior probabilities')) + scale_x_continuous(trans='log10')+
    scale_shape_manual(name=TeX('Sampler'), breaks = c('S3_logistic','S3_probit','SKINNY'),
                       labels=unname(TeX(c('$S^3$ Logistic', '$S^3$ Probit', 'Skinny Gibbs Logistic'))), 
                       values = c(4, 1, 3)) +
    scale_x_continuous(trans='log10') + theme_classic(base_size = 12) +
    theme(legend.position = 'bottom',legend.text=element_text(size=10))+
    guides(shape=guide_legend(nrow=2,byrow=TRUE))
  # z_variable_section
  
  #### Time per iteration plot: binary
  multi_datasets_df_time <-
    multi_datasets_df %>% filter(algo%in%c('S3_logistic','S3_probit','SOTA_logistic','SOTA_probit','SKINNY'))
  level_order <- 
    factor(multi_datasets_df_time$algo, 
           level = c('S3_logistic','SOTA_logistic','S3_probit','SOTA_probit','SKINNY'))
  label_title= TeX('')
  label_breaks <- c('S3_logistic','SOTA_logistic','S3_probit','SOTA_probit','SKINNY')
  label_names <- unname(TeX(c('$S^3$ Logistic', 'SOTA Logistic', '$S^3$ Probit','SOTA Probit','Skinny Gibbs Logistic')))
  
  # Both logistic and probit
  time_comparison <- 
    ggplot(multi_datasets_df_time, aes(x=level_order, y=time_mean*1000, color=algo, shape=algo)) + 
    geom_point(size=3) + 
    xlab(TeX('Sampler')) + ylab(TeX('Time per iteration (ms)')) +
    scale_color_manual(name=label_title, breaks=label_breaks, labels=label_names,
                       values = c('Black','Gray','Black','Gray','Black')) +
    scale_shape_manual(name=label_title, breaks=label_breaks, labels=label_names,
                       values = c(4,4,1,1,3)) +
    scale_x_discrete(name=TeX('Sampler'), breaks=label_breaks, labels=label_names) +
    geom_errorbar(aes(ymax=(time_mean+time_sd/no_chains)*1000, ymin=(time_mean-time_sd/no_chains)*1000),
                  position=position_dodge(.9), show.legend = FALSE) +
    # scale_y_continuous(limits =c(0,40)) +
    scale_y_continuous(trans='log10') +
    theme_classic(base_size = 12) +
    theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
    guides(color=guide_legend(nrow=3,byrow=TRUE)) + 
    theme(axis.text.x = element_blank())
  # time_comparison 
  shared_legend <- ggpubr::get_legend(time_comparison, position = NULL)
  
  runtimes_per_iteration <- multi_datasets_df %>% group_by(dataset,algo) %>% select(time_mean)
  kfold_cv_df <- 
    rmse_df %>% 
    left_join(runtimes_per_iteration, by='algo') %>%
    mutate(time_elasped=iteration*time_mean)
  
  kfold_cv_plot <-
    ggplot(kfold_cv_df %>% filter(iteration%%50==0), aes(x=time_elasped, y=rmse, color=algo, linetype=algo)) +
    geom_line(size=1) + 
    scale_color_manual(name=label_title, breaks=label_breaks, labels=label_names,
                       values = c('Black','Gray','Black','Gray','Black')) +
    scale_linetype_manual(name=label_title, breaks=label_breaks, labels=label_names,
                       values = c('solid','solid','dashed','dashed','dotted')) +
    # xlab(TeX('Time elasped (s)')) +
    scale_y_continuous(trans='log10', name = 'Cross-validation RMSE') +
    scale_x_continuous(limits = c(0,400),name = TeX('Time elasped (s)')) +
    theme_classic(base_size = 12) +
    theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  # kfold_cv_plot
}

if(dataset_response_type == 'continuous'){
  z_variable_section <- 
    ggplot(variable_section_df %>% filter(algo%in%c('S3_linear')), 
           aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), shape=algo)) + 
    geom_point(size=3, color='black', stroke=0.5) + xlab(xlabel) + 
    ylab(TeX('Marginal posterior probabilities')) + scale_x_continuous(trans='log10')+
    scale_shape_manual(name=TeX('Sampler'), breaks = c('S3_linear'),
                       labels=unname(TeX(c('$S^3$ Linear'))), values = c(4)) +
    scale_x_continuous(trans='log10') + theme_classic(base_size = 12) +
    theme(legend.position = 'bottom',legend.text=element_text(size=10))
  # z_variable_section
  
  #### Time per iteration plot: binary
  multi_datasets_df_time <-
    multi_datasets_df %>% filter(algo%in%c('S3_linear','SOTA_linear'))
  level_order <- 
    factor(multi_datasets_df_time$algo, level = c('S3_linear','SOTA_linear'))
  label_name <- TeX('')
  label_breaks <- c('S3_linear','SOTA_linear')
  label_names <- unname(TeX(c('$S^3$ Linear', 'SOTA Linear')))
  
  # Linear
  time_comparison <- 
    ggplot(multi_datasets_df_time, aes(x=level_order, y=time_mean*1000, color=algo)) + 
    geom_point(size=3, shape=4) + 
    xlab(TeX('Sampler')) + ylab(TeX('Time per iteration (ms)')) +
    scale_color_manual(name=label_name, breaks=label_breaks, labels=label_names,values = c('Black','Gray')) +
    scale_x_discrete(name=TeX(''), breaks=label_breaks, labels=label_names) +
    geom_errorbar(aes(ymax=(time_mean+time_sd/no_chains)*1000, ymin=(time_mean-time_sd/no_chains)*1000),
                  position=position_dodge(.9), show.legend = FALSE) +
    # scale_y_continuous(limits =c(0,40)) +
    theme_classic(base_size = 12) +
    theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
    guides(color=guide_legend(nrow=1,byrow=TRUE)) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, colour = 'white'))
  time_comparison
  shared_legend <- ggpubr::get_legend(time_comparison, position = NULL)
  
  runtimes_per_iteration <- multi_datasets_df %>% group_by(dataset,algo) %>% select(time_mean)
  kfold_cv_df <- 
    rmse_df %>% 
    left_join(runtimes_per_iteration, by='algo') %>%
    mutate(time_elasped=iteration*time_mean)
  
  kfold_cv_plot <-
    ggplot(kfold_cv_df, aes(x=time_elasped, y=rmse, color=algo, shape=algo)) +
    geom_line(size=1) + 
    scale_color_manual(name=label_name, breaks=label_breaks, labels=label_names,
                       values = c('Black','Gray')) +
    scale_shape_manual(name=label_name, breaks=label_breaks, labels=label_names,
                       values = c(4,1)) +
    # xlab(TeX('Time elasped (s)')) +
    scale_y_continuous(name = 'Cross-validation RMSE') +
    scale_x_continuous(limits = c(0,5000),name = TeX('Time elasped (s)')) +
    theme_classic(base_size = 12) +
    theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
  # kfold_cv_plot
}

# Combined plot
plot_combined <- ggarrange(z_variable_section, time_comparison, kfold_cv_plot, legend = "bottom", nrow=1)
# plot_combined
plot_name <- paste("inst/dataset_simulations/",dataset_name,"_plot.pdf",sep = '')
# ggsave(filename = plot_name, plot = plot_combined, width = 12, height = 4)

