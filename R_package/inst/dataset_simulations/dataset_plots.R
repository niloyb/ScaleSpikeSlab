# Malware plot

dataset <- 'malware'
filename <- paste("inst/dataset_simulations/",dataset,'_sims.RData',sep = '')
load(filename)
xlab_title <- TeX('Malware Covariates')


# Time comparison
level_order <- 
  factor(sss_sota_comparison_df$algo, 
         level = c('S3_logistic','SOTA_logistic','S3_probit','SOTA_probit'))
label= TeX('Sampler')
label_breaks <- c('S3_logistic','SOTA_logistic','S3_probit','SOTA_probit')
label_names <- unname(TeX(c('S^3 Logistic', 'SOTA Logistic',
                            'S^3 Probit', 'SOTA Probit')))

# sss_sota_comparison_df %>% group_by(algo) %>% summarise(time_mean = mean(time_mean))



# Both logistic and probit
time_comparison <- 
  ggplot(sss_sota_comparison_df, aes(x=level_order, y=time_mean*1000, color=algo, shape=algo)) + 
  geom_point(size=4) + 
  xlab(TeX('Sampler')) + ylab(TeX('Time per iteration (ms)')) +
  scale_color_manual(name=label, breaks=label_breaks, labels=label_names,
                     values = c('Black','Black','Gray','Gray')) +
  scale_shape_manual(name=label, breaks=label_breaks, labels=label_names,
                     values = c(4,1,4,1)) +
  scale_x_discrete(name=TeX('Sampler'), 
                   breaks=c('S3_logistic','SOTA_logistic','S3_probit','SOTA_probit'),
                   labels=unname(TeX(c('S^3 Logistic', 'SOTA Logistic',
                                       'S^3 Probit', 'SOTA Probit')))) +
  geom_errorbar(aes(ymax=(time_mean+time_sd)*1000,
                    ymin=(time_mean-time_sd)*1000),position=position_dodge(.9),
                show.legend = FALSE) +
  # scale_y_continuous(limits =c(0,40)) +
  theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
time_comparison 

shared_legend <- get_legend(time_comparison, position = NULL)

# Variable selection comparison
variable_section_df <- 
  sss_sota_comparison_df %>% select(z_avg, n, p) %>%
  summarise(prob_covariate_slab = rowMeans(do.call(cbind, z_avg)),n=mean(n),p=mean(p)) %>%
  mutate(prob_covariate_slab=prob_covariate_slab, covariate_index=1:n(), n=n, p=p) %>%
  arrange(desc(prob_covariate_slab)) %>%
  mutate(xaxis =1:n())

# variable_section <- 
#   ggplot(variable_section_df, 
#          aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo)) + 
#   geom_point(shape=4) + xlab(xlab_title) + ylab(TeX('Marginal posterior probabilities')) +
#   scale_color_manual(name=TeX('Sampler'), breaks = c("S3_logistic", "S3_probit"),
#                      labels=unname(TeX(c('S^3 Logistic','S^3 Probit'))), values = c('Black', 'Gray')) +
#   scale_x_continuous(trans='log10') + theme_classic(base_size = 9) +
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"),legend.text=element_text(size=10)) +
#   guides(color=guide_legend(nrow=1,byrow=TRUE))
# shared_legend <- get_legend(variable_section, position = NULL)


variable_section_logistic <- 
  ggplot(variable_section_df %>% filter(algo=="S3_logistic"), 
         aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo)) + 
  geom_point(shape=4) + xlab(xlab_title) + 
  ylab(TeX('Marginal posterior probabilities')) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("S3_logistic", "S3_probit"),
                     labels=unname(TeX(c('S^3 Logistic','S^3 Probit'))), values = c('Black', 'Gray')) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"),legend.text=element_text(size=10)) +
  guides(color='none')
variable_section_logistic
variable_section_probit <- 
  ggplot(variable_section_df %>% filter(algo=="S3_probit"), 
         aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo)) + 
  geom_point(shape=4) + xlab(xlab_title) + ylab(TeX('Marginal posterior probabilities')) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("S3_logistic", "S3_probit"),
                     labels=unname(TeX(c('S^3 Logistic','S^3 Probit'))), values = c('Black', 'Gray')) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"),legend.text=element_text(size=10)) +
  guides(color='none')
variable_section_probit

# variable_section_plot_combined <- 
#   ggarrange(variable_section_logistic, variable_section_probit, common.legend = TRUE, legend.grob =  shared_legend, legend = "bottom", nrow=1)
# variable_section_plot_combined



# Combined plot
plot_combined <- 
  ggarrange(variable_section_logistic, variable_section_probit,
            time_comparison, common.legend=TRUE, legend.grob=shared_legend, 
            legend = "bottom", nrow=1)
plot_combined
plot_name <- paste("inst/dataset_simulations/",dataset,'_plot_new.pdf',sep = '')
# ggsave(filename = plot_name, plot = plot_combined, width = 10, height = 4)




################################################################################
# Multiple datasets plot
# filename2 <- paste("inst/dataset_simulations/multiple_dataset_sims_new.RData",sep = '')
# load(filename2)

filename2 <- paste("inst/dataset_simulations/datasets_sims.RData",sep = '')
load(filename2)

sss_sota_multi_dataset_time_comparison_df <- 
  sss_sota_multi_dataset_time_comparison_df %>%
  mutate(algo=if_else((algo=='S3_linear')|(algo=='S3_probit'),'S3','SOTA'))
sss_sota_multi_dataset_time_comparison_df$dataset[sss_sota_multi_dataset_time_comparison_df$dataset=='Synthetic_Binary'] <- 'Synthetic Binary'
sss_sota_multi_dataset_time_comparison_df$dataset[sss_sota_multi_dataset_time_comparison_df$dataset=='Synthetic_Continuous'] <- 'Synthetic Continuous'

# Time comparison
time_comparison <- 
  ggplot(sss_sota_multi_dataset_time_comparison_df, aes(x=dataset, y=time_mean*1000, color=algo)) + 
  geom_point(size=3, shape=4) + xlab(TeX('Datasets')) + ylab(TeX('Time per iteration (ms)')) +
  scale_x_discrete(limits=sss_sota_multi_dataset_time_comparison_df$dataset) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("S3", "SOTA"), 
                     labels=unname(TeX(c('S^3','SOTA'))),
                     values = c('Black', 'Gray')) +
  geom_errorbar(aes(ymax=(time_mean+2*time_sd/sqrt(no_chains))*1000, 
                    ymin=(time_mean-2*time_sd/sqrt(no_chains))*1000),
                position=position_dodge(.0)) +
  scale_y_continuous(trans='log10') + theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
time_comparison

plot_name2 <- paste("inst/dataset_simulations/multi_dataset_time_plot_new.pdf",sep = '')
# ggsave(filename = plot_name2, plot = time_comparison, width = 5, height = 4)




