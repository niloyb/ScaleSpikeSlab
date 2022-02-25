# Malware plot

data <- 'malware'
filename <- paste("/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScalableSpikeSlab/inst/dataset_simulations/",data,'_sims.RData',sep = '')
load(filename)

# if(data=='synthetic_linreg'){xlab_title <- TeX('Synthetic Linear Regression Covariates')}
# if(data=='synthetic_logreg'){xlab_title <- TeX('Synthetic Logistic Regression Covariates')}
# if(data=='riboflavin'){xlab_title <- TeX('Riboflavin Covariates')}
# if(data=='maize'){xlab_title <- TeX('Maize Covariates')}
# if(data=='pcr'){xlab_title <- TeX('PCR Covariates')}
# if(data=='malware'){xlab_title <- TeX('Malware Covariates')}
# if(data=='lymph'){xlab_title <- TeX('Lymph Node Covariates')}

xlab_title <- TeX('Malware Covariates')

# Time comparison
time_comparison <- 
  ggplot(sss_sota_comparison_df, aes(x=as.factor(algo), y=time_mean*1000, color=algo)) + 
  geom_point(size=4, shape=4) + xlab(TeX('Sampler')) + ylab(TeX('Time per iteration (ms)')) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),labels=unname(TeX(c('S^3', 'SOTA'))),
                     values = c('Black', 'Gray')) +
  scale_x_discrete(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),labels=unname(TeX(c('S^3', 'SOTA')))) +
  geom_errorbar(aes(ymax=(time_mean+time_sd)*1000, ymin=(time_mean-time_sd)*1000),position=position_dodge(.9)) +
  scale_y_continuous(limits =c(0,40)) +
  theme_classic(base_size = 9) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
time_comparison

# Variable selection comparison
variable_section_df <- 
  sss_sota_comparison_df %>% select(z_avg, n, p) %>%
  summarise(prob_covariate_slab = rowMeans(do.call(cbind, z_avg)),n=mean(n),p=mean(p)) %>%
  mutate(prob_covariate_slab=prob_covariate_slab, covariate_index=1:n(), n=n, p=p) %>%
  arrange(desc(prob_covariate_slab)) %>%
  mutate(xaxis =1:n())

variable_section_comparison <- 
  ggplot(variable_section_df, aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo)) + 
  geom_point(shape=4) + xlab(xlab_title) + ylab(TeX('Marginal posterior probabilities')) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),
                     labels=unname(TeX(c('S^3', 'SOTA'))), values = c('Black', 'Gray')) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 9) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"),legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
shared_legend <- get_legend(variable_section_comparison, position = NULL)

variable_section_comparison <- 
  ggplot(variable_section_df %>% filter(algo == 'ScalableSpikeSlab'), aes(x=xaxis, y=sort(prob_covariate_slab, decreasing = TRUE), color=algo)) + 
  geom_point(shape=4) + xlab(xlab_title) + ylab(TeX('Marginal posterior probabilities')) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),
                     labels=unname(TeX(c('S^3', 'SOTA'))), values = c('Black', 'Gray')) +
  scale_x_continuous(trans='log10') + theme_classic(base_size = 9) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"),legend.text=element_text(size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
variable_section_comparison

# Combined plot
plot_combined <- 
  ggarrange(variable_section_comparison, 
            time_comparison, common.legend = TRUE, legend.grob =  shared_legend, legend = "bottom", nrow=1)
plot_combined
plot_name <- paste("/Users/niloybiswas/Dropbox/Apps/Overleaf/fast_spike_slab/images/",data,'_plot.pdf',sep = '')
# ggsave(filename = plot_name, plot = plot_combined, width = 4, height = 4.5)






# Multiple datasets plot
filename2 <- paste("/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScalableSpikeSlab/inst/dataset_simulations/multiple_dataset_sims.RData",sep = '')
load(filename2)

# Time comparison
time_comparison <- 
  ggplot(sss_sota_multi_dataset_time_comparison_df, aes(x=dataset, y=time_mean*1000, color=algo)) + 
  geom_point(size=2, shape=4) + xlab(TeX('Datasets')) + ylab(TeX('Time per iteration (ms)')) +
  scale_x_discrete(limits=sss_sota_multi_dataset_time_comparison_df$dataset) +
  scale_color_manual(name=TeX('Sampler'), breaks = c("ScalableSpikeSlab", "Sota"),labels=unname(TeX(c('S^3', 'SOTA'))),values = c('Black', 'Gray')) +
  # geom_errorbar(aes(ymax=(time_mean+time_sd)*1000, ymin=(time_mean-time_sd)*1000),position=position_dodge(.9)) +
  scale_y_continuous(trans='log10') + theme_classic(base_size = 10) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=10)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
time_comparison

plot_name2 <- paste("/Users/niloybiswas/Dropbox/Apps/Overleaf/fast_spike_slab/images/multi_dataset_time_plot.pdf",sep = '')
# ggsave(filename = plot_name2, plot = time_comparison, width = 4, height = 5.5)


