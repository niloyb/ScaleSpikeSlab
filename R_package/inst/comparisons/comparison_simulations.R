rm(list=ls())
# Synthetic dataset simulations
library(ScaleSpikeSlab)

source('inst/comparisons/comparison_functions.R')
# no_cores <- 2
no_cores <- detectCores()-1
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
# load('inst/comparisons/comparison_vary_p_df_new.Rdata')

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
# ggsave(filename = "inst/comparisons/time_comparison_arxiv.pdf", plot = time_comparison, width = 7, height = 4)
# time_comparison <- 
#   time_comparison + theme_classic(base_size = 16) + theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"))
# ggsave(filename = "inst/comparisons/time_comparison_icml.pdf", plot = time_comparison, width = 7.5, height = 4)

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
# ggsave(filename = "inst/comparisons/tpr_fdr_comparison_new.pdf", plot = tpr_fdr_comparison, width = 7, height = 4)



######## Fix n, fix p: stat estimation comparison over time ########
n <- 200
p <- 1000
s0 <- 10
logreg_sync_data <-
  ScaleSpikeSlab:::synthetic_data(n, p, s0, type = 'logistic', signal='decay')

stat_comparison_over_time_sss_logistic <- 
  comparison_sims_over_time(logreg_sync_data$X,logreg_sync_data$y,logreg_sync_data$true_beta,
                            chain_length=15e3,burnin=1e3,
                            algos=c('S3_logistic','SOTA_logistic'),
                            no_chains=10)
stat_comparison_over_time_skinny_logistic <- 
  comparison_sims_over_time(logreg_sync_data$X,logreg_sync_data$y,logreg_sync_data$true_beta,
                            chain_length=15e3,burnin=1e3,
                            algos=c('SKINNY'),
                            no_chains=10)
stat_comparison_over_time_probit <- 
  comparison_sims_over_time(logreg_sync_data$X,logreg_sync_data$y,logreg_sync_data$true_beta,
                            chain_length=5e4,burnin=1e3, algos=c('S3_probit','SOTA_probit'),
                            no_chains=10)
stat_comparison_over_time_df <- 
  rbind(stat_comparison_over_time_sss_logistic,
        stat_comparison_over_time_skinny_logistic,
        stat_comparison_over_time_probit)

stat_performance_df <-
  stat_comparison_over_time_df %>%
  select(algo, tpr, fdr, iteration, time_per_iteration) %>%
  melt(id=c('algo','iteration','time_per_iteration')) %>%
  group_by(algo,variable,iteration) %>%
  summarise(mean=mean(value), sd =sd(value), no_repeats=n(), time_per_iteration=mean(time_per_iteration))

# save(stat_performance_df, file = 'inst/comparisons/stat_comparison_over_time.Rdata')
# load('inst/comparisons/stat_comparison_over_time.Rdata')


stat_performance_df <-
  stat_performance_df %>%
  mutate(elasped_clock_time=time_per_iteration*iteration)


label_name <- TeX('Sampler')
label_breaks <- c('S3_logistic','S3_probit','SKINNY','SOTA_logistic','SOTA_probit')
label_names <- unname(TeX(c('S^3 Logistic', 'S^3 Probit',
                            'Skinny Gibbs Logistic','SOTA Logistic', 'SOTA Probit')))

tpr_over_iterations <-
  ggplot(stat_performance_df %>% filter(variable=='tpr',algo%in%c('S3_logistic','S3_probit','SKINNY'),iteration%%500==0), 
         aes(x=iteration, y=mean, linetype=algo)) +
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean-sd/sqrt(no_repeats), ymax=mean+sd/sqrt(no_repeats),
                  fill=algo),alpha=0.2, colour = NA) +
  xlab(TeX('Chain length')) + ylab(TeX('TPR')) +
  scale_linetype_manual(name=label_name, breaks = label_breaks, labels=label_names,
                        values = c('solid','dashed','dotted','solid','dashed')) +
  scale_fill_manual(name=label_name, breaks = label_breaks, labels=label_names,
                    values = c('black','black','black','gray','gray')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(0,15e3)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,ncol=3,byrow=TRUE))
tpr_over_iterations

fdr_over_iterations <- 
  ggplot(stat_performance_df %>% filter(variable=='fdr',algo%in%c('S3_logistic','S3_probit','SKINNY'),iteration%%500==0), 
         aes(x=iteration, y=(mean), linetype=algo)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=pmax(mean-sd/sqrt(no_repeats),0), ymax=(mean+sd/sqrt(no_repeats)),
                  fill=algo),alpha=0.2, colour = NA) +
  xlab(TeX('Chain length')) + ylab(TeX('FDR')) +
  scale_linetype_manual(name=label_name, breaks = label_breaks, labels=label_names,
                        values = c('solid','dashed','dotted','solid','dashed')) +
  scale_fill_manual(name=label_name, breaks = label_breaks, labels=label_names,
                    values = c('black','black','black','gray','gray')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(limits = c(0,15e3)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,ncol=3,byrow=TRUE))
fdr_over_iterations

tpr_fdr_over_iterations <- 
  ggarrange(tpr_over_iterations, fdr_over_iterations, common.legend = TRUE, legend = "bottom")
tpr_fdr_over_iterations
# ggsave(filename = "inst/comparisons/tpr_fdr_over_iterations.pdf", plot = tpr_fdr_over_iterations, width = 7, height = 4)

tpr_over_time <-
  ggplot(stat_performance_df %>% filter(variable=='tpr',iteration%%500==0), aes(x=elasped_clock_time, y=mean, linetype=algo, color=algo)) +
  geom_line(size=1) +
  # geom_ribbon(aes(ymin=pmax(mean-sd/sqrt(no_repeats),0), ymax=pmin(mean+sd/sqrt(no_repeats),1)),alpha=0.2, colour = NA) +
  xlab(TeX('Time elasped (s)')) + ylab(TeX('TPR')) +
  scale_linetype_manual(name=label_name, breaks = label_breaks, labels=label_names,
                        values = c('solid','dashed','dotted','solid','dashed')) +
  scale_color_manual(name=label_name, breaks = label_breaks, labels=label_names,
                     values = c('black','black','black','darkgray','darkgray')) +
  # scale_color_manual(values = c('black','black','black','gray','gray'),guide='none') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(0,90)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
tpr_over_time

fdr_over_time <- 
  ggplot(stat_performance_df %>% filter(variable=='fdr',iteration%%500==0), aes(x=elasped_clock_time, y=mean, linetype=algo, color=algo)) +
  geom_line(size=1) +
  # geom_ribbon(aes(ymin=pmax(mean-sd/sqrt(no_repeats),0), ymax=pmin(mean+sd/sqrt(no_repeats),1)),alpha=0.2, colour = NA) +
  xlab(TeX('Time elasped (s)')) + ylab(TeX('FDR')) +
  scale_linetype_manual(name=label_name, breaks = label_breaks, labels=label_names,
                        values = c('solid','dashed','dotted','solid','dashed')) +
  scale_color_manual(name=label_name, breaks = label_breaks, labels=label_names,
                    values = c('black','black','black','darkgray','darkgray')) +
  # scale_color_manual(values = c('black','black','black','gray','gray'),guide='none') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(limits = c(0,90)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
fdr_over_time

tpr_fdr_over_time <-
  ggarrange(tpr_over_time, fdr_over_time, common.legend = TRUE, legend = "bottom",nrow=1)
tpr_fdr_over_time

tpr_fdr_combined <-
  ggarrange(tpr_fdr_over_iterations, tpr_fdr_over_time, common.legend = FALSE, legend = "bottom",ncol=2,nrow=1)
# tpr_fdr_combined
# ggsave(filename = "inst/comparisons/tpr_fdr_over_time_combined.pdf", plot = tpr_fdr_combined, width = 12, height = 5)

# tpr_fdr_over_time_combined <- 
#   ggarrange(tpr_over_iterations, fdr_over_iterations, tpr_over_time, fdr_over_time, 
#             common.legend = TRUE, legend = "bottom",ncol=4,nrow=1)
# tpr_fdr_over_time_combined
# ggsave(filename = "inst/comparisons/tpr_fdr_over_time_combined.pdf", plot = tpr_fdr_over_time_combined, width = 12, height = 5)



