# Synthetic dataset simulations
rm(list=ls())
library(ScaleSpikeSlab)
source('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/comparison_functions.R')

library(dplyr)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(ggridges)
library(ggpubr)


######## Fix n vary p: time per iteration simulations ########
s0_seq <- c(seq(10,100,10),seq(200,400,100))
# s0_seq <- seq(10,50,10)
time_comparison_grid <- data.frame(n=10*s0_seq, p=100*s0_seq, error_std=2, s0=s0_seq)
time_comparison_list <- split(time_comparison_grid, 1:nrow(time_comparison_grid))
comparison_vary_p_df1 <-
  comparison_sims(time_comparison_list, chain_length=1e3,burnin=0,no_repeats=1,
                  algos=c('ScalableSpikeSlab', 'SkinnyGibbs'),store=FALSE)

# Running Sota sampler for a smaller number of iterations to calculate average
# time per iteration as it is slower and the per iteration cost does not depend 
# on chain state
comparison_vary_p_df2 <- 
  comparison_sims(time_comparison_list,chain_length=1e2,burnin=0, no_repeats=1,
                  algos=c('Sota'),store=FALSE)
comparison_vary_p_df <- rbind(comparison_vary_p_df1, comparison_vary_p_df2)
# save(comparison_vary_p_df, file = '/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/comparison_vary_p_df.Rdata')
# load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/comparison_vary_p_df.Rdata')
# comparison_vary_p_df1
# Run-time comparison plots
time_df <- 
  comparison_vary_p_df %>% 
  group_by(algo,p) %>% 
  summarise(time=mean(time))
# (time_df%>%filter(algo=='SkinnyGibbs'))$time/(time_df%>%filter(algo=='ScalableSpikeSlab'))$time
# (time_df%>%filter(algo=='Sota'))$time/(time_df%>%filter(algo=='ScalableSpikeSlab'))$time

time_comparison <- ggplot(time_df, aes(x=p, y=time*1000, linetype=algo)) + 
  geom_line(size=1) + xlab(TeX('Dimension $p$')) + ylab(TeX('Time per iteration (ms)')) +
  # scale_y_continuous(limits = c(0,600)) +
  # scale_y_continuous(trans='log10') +
  scale_linetype_manual(name=TeX('Sampler'),
                        breaks = c("ScalableSpikeSlab", "Sota", "SkinnyGibbs"),
                        labels=unname(TeX(c('S^3', 'SOTA', 'Skinny Gibbs'))),
                        values = c('solid','dotted','dashed')) +
  theme_classic(base_size = 10) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
time_comparison
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/fast_spike_slab/images/time_comparison.pdf", plot = time_comparison, width = 4, height = 4.5)


######## Fix n vary p: stat estimation comparison simulations ########
stat_comparison_grid <- data.frame(n=100, p=seq(100,1000,100), error_std=2, s0=5)
stat_comparison_list <- split(stat_comparison_grid, 1:nrow(stat_comparison_grid))
stat_comparison_vary_p_df <- 
  comparison_sims(stat_comparison_list, chain_length=5e3, burnin=1e3, algos=c('ScalableSpikeSlab','SkinnyGibbs'), no_repeats=10, signal = 'decay')
# save(stat_comparison_vary_p_df, file = '/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/stat_comparison_vary_p_df.Rdata')
# load('/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/github/ScaleSpikeSlab/R_package/inst/comparisons/stat_comparison_vary_p_df.Rdata')

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
                        breaks = c("ScalableSpikeSlab", "SkinnyGibbs"),
                        labels=unname(TeX(c('$S^3$','Skinny Gibbs'))),
                        values = c('solid','dashed')) +
  scale_fill_manual(values = c('black','black', 'black'),guide='none') +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  # scale_x_continuous(labels = scales::scientific_format(accuracy = 1)) +
  theme_classic(base_size = 10) +
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
                        breaks = c("ScalableSpikeSlab", "SkinnyGibbs"),
                        labels=unname(TeX(c('$S^3$','Skinny Gibbs'))),
                        values = c('solid','dashed')) +
  scale_fill_manual(values = c('black','black', 'black'),guide='none') +
  scale_y_continuous(limits = c(0,0.04), labels = scales::percent_format(accuracy = 1)) +
  # scale_x_continuous(labels = scales::scientific_format(accuracy = 1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
fdr_comparison

tpr_fdr_comparison <- 
  ggarrange(tpr_comparison, fdr_comparison, common.legend = TRUE, legend = "bottom")
tpr_fdr_comparison
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/fast_spike_slab/images/tpr_fdr_comparison.pdf", plot = tpr_fdr_comparison, width = 4, height = 4.5)

