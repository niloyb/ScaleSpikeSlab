# Synthetic dataset simulations
library(ScaleSpikeSlab)

source('inst/comparisons/comparison_functions.R')
no_cores <- 2
# no_cores <- detectCores()-1
registerDoParallel(cores = no_cores)


library(dplyr)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(ggridges)
library(ggpubr)


######## Fix n vary p: time per iteration simulations ########
s0_seq <- c(seq(20,100,20),seq(200,400,200))
time_comparison_grid <- data.frame(n=10*s0_seq, p=100*s0_seq, error_std=2, s0=s0_seq)
# s0_seq <- seq(10,100,20)
# time_comparison_grid <- data.frame(n=4000, p=c(40000), error_std=2, s0=s0_seq)
time_comparison_list <- split(time_comparison_grid, 1:nrow(time_comparison_grid))
comparison_vary_p_df1 <-
  comparison_sims(time_comparison_list, chain_length=100,burnin=0,no_repeats=1,
                  algos=c('S3_logistic','S3_probit','SKINNY'),
                  store=FALSE, verbose=TRUE)
# Running Sota sampler for a smaller number of iterations to calculate average
# time per iteration as it is slower and the per iteration cost does not depend 
# on chain state
comparison_vary_p_df2 <- 
  comparison_sims(time_comparison_list,chain_length=5,burnin=0, no_repeats=1,
                  algos=c('SOTA_logistic', 'SOTA_probit'),store=FALSE, verbose=TRUE)
comparison_vary_p_df <- rbind(comparison_vary_p_df1, comparison_vary_p_df2)
# save(comparison_vary_p_df, file = 'inst/comparisons/comparison_vary_p_df_new.Rdata')
# load('/inst/comparisons/comparison_vary_p_df_new.Rdata')

# Run-time comparison plots
time_df <- 
  comparison_vary_p_df %>% 
  group_by(algo,n,p) %>% 
  summarise(time=mean(time))
# (time_df%>%filter(algo=='SKINNY'))$time/(time_df%>%filter(algo=='S3_logistic'))$time
# (time_df%>%filter(algo=='SOTA_logistic'))$time/(time_df%>%filter(algo=='S3_logistic'))$time
# (time_df%>%filter(algo=='SOTA_probit'))$time/(time_df%>%filter(algo=='S3_probit'))$time

legend_name <- TeX('Sampler')
legend_breaks <- c('S3_logistic', 'SOTA_logistic', 'SKINNY', 'S3_probit','SOTA_probit')
legend_labels <- unname(TeX(c('S^3 Logistic', 'SOTA Logistic', 'Skinny Gibbs Logistic', 
                              'S^3 Probit', 'SOTA Probit')))
linetype_values <- c('solid','dashed','dotted','solid','dashed')
color_values <- c('black','black', 'black', 'light gray', 'light gray')
time_comparison <- 
  ggplot(time_df, aes(x=p, y=time*1000, linetype=algo, col=algo)) + 
  geom_line(size=1) + xlab(TeX('Dimension $p$')) + ylab(TeX('Time per iteration (ms)')) +
  # scale_y_continuous(trans='log10') +
  labs(color=TeX('Sampler'), linetype=TeX('Sampler')) +
  scale_linetype_manual(name=legend_name,breaks=legend_breaks,labels=legend_labels,
                        values = linetype_values) +
  scale_color_manual(name=legend_name,breaks=legend_breaks,labels=legend_labels,
                     values=color_values) +
  theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
time_comparison
# ggsave(filename = "/inst/comparisons/time_comparison_new.pdf", plot = time_comparison, width = 7, height = 4)

(time_df%>%filter(algo=='SKINNY'))$time/(time_df%>%filter(algo=='S3_logistic'))$time
(time_df%>%filter(algo=='SOTA_logistic'))$time/(time_df%>%filter(algo=='S3_logistic'))$time
(time_df%>%filter(algo=='SOTA_probit'))$time/(time_df%>%filter(algo=='S3_probit'))$time


######## Fix n vary p: stat estimation comparison simulations ########
stat_comparison_grid <- data.frame(n=100, p=seq(200,1000,200), error_std=2, s0=5)
# s0_seq <- seq(2,10,2)
# stat_comparison_grid <- data.frame(n=10*s0_seq, p=100*s0_seq, error_std=2, s0=s0_seq)
stat_comparison_list <- split(stat_comparison_grid, 1:nrow(stat_comparison_grid))
stat_comparison_vary_p_df <- 
  comparison_sims(stat_comparison_list, chain_length=5e3, burnin=1e3, 
                  algos=c('S3_probit','S3_logistic','SKINNY'), 
                  no_repeats=10, signal = 'decay')
# save(stat_comparison_vary_p_df, file = 'inst/comparisons/stat_comparison_vary_p_df_new.Rdata')
# load('inst/comparisons/stat_comparison_vary_p_df_new.Rdata')

stat_performance_df <-
  stat_comparison_vary_p_df %>%
  select(algo, tpr, fdr, mse, n, p, iteration) %>%
  melt(id=c('algo','n','p','iteration')) %>%
  group_by(algo,n,p,variable) %>%
  summarise(mean=mean(value), sd =sd(value), no_repeats=n())
tpr_comparison <- 
  ggplot(stat_performance_df %>% filter(variable=='tpr'), aes(x=p, y=mean, linetype=algo)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean-sd/sqrt(no_repeats), ymax=mean+sd/sqrt(no_repeats),
                  fill=algo),alpha=0.2, colour = NA) +
  xlab(TeX('Dimension $p$')) + ylab(TeX('TPR')) +
  scale_linetype_manual(name=TeX('Sampler'),
                        breaks = c('S3_logistic','S3_probit','SKINNY'),
                        labels=unname(TeX(c('S^3 Logistic', 'S^3 Probit', 
                                            'Skinny Gibbs Logistic'))),
                        values = c('solid','dashed','dotted')) + 
  scale_fill_manual(values = rep('black',3),guide='none') +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
tpr_comparison

fdr_comparison <- 
  ggplot(stat_performance_df %>% filter(variable=='fdr'), aes(x=p, y=(mean), linetype=algo)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=pmax(mean-sd/sqrt(no_repeats),0), ymax=(mean+sd/sqrt(no_repeats)),
                  fill=algo),alpha=0.2, colour = NA) +
  xlab(TeX('Dimension $p$')) + ylab(TeX('FDR')) +
  scale_linetype_manual(name=TeX('Sampler'),
                        breaks = c('S3_logistic','S3_probit','SKINNY'),
                        labels=unname(TeX(c('S^3 Logistic', 'S^3 Probit', 
                                            'Skinny Gibbs Logistic'))),
                        values = c('solid','dashed','dotted')) + 
  scale_fill_manual(values = rep('black',3),guide='none') +
  scale_y_continuous(limits = c(0,0.04), labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 14) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
fdr_comparison

tpr_fdr_comparison <- 
  ggarrange(tpr_comparison, fdr_comparison, common.legend = TRUE, legend = "bottom")
tpr_fdr_comparison
# ggsave(filename = "/inst/comparisons/tpr_fdr_comparison_new.pdf", plot = tpr_fdr_comparison, width = 7, height = 4)

