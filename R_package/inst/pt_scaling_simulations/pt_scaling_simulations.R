# Simulation study on synthetic datasets 
library(ScaleSpikeSlab)
library(doParallel)
no_cores <- 2
# no_cores <- detectCores()-1
registerDoParallel(cores = no_cores)
library(foreach)

library(dplyr)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(ggridges)
library(ggpubr)


###### Function to for run-time simulations ######
generate_sims <- function(n_p_error_s0_list,chain_length=1e4,burnin=5e3,no_repeats=1){
  foreach(n_p_error_s0 = n_p_error_s0_list, .combine = rbind)%:%
    foreach(i = c(1:no_repeats), .combine = rbind)%dopar%{
      n <- n_p_error_s0$n
      p <- n_p_error_s0$p
      error_std <- n_p_error_s0$error_std
      s0 <- n_p_error_s0$s0
      
      syn_data <- synthetic_data(n=n,p=p,s0=s0,error_std=error_std,type='linear',signal='constant')
      X <- syn_data$X
      y <- syn_data$y
      Xt <- t(X)
      signal_indices <- syn_data$true_beta!=0
      
      params <- spike_slab_params(n, p)
      
      sss_time_taken <-
        system.time(
          sss_chain <- 
            spike_slab_linear(chain_length=chain_length, X=X,Xt=Xt,y=y,
                            tau0=params$tau0, tau1=params$tau1, q=params$q, 
                            a0=params$a0,b0=params$b0, rinit=NULL, verbose=TRUE))
      
      delta_t <- c(NA,rowSums(sss_chain$z[c(1:(chain_length-1)),]!=sss_chain$z[c(2:chain_length),]))
      no_active <- rowSums(sss_chain$z[c(1:chain_length),])
      p_t <- pmin(delta_t, no_active, na.rm = TRUE)
      
      sss_tpr <- mean(colMeans(sss_chain$z[c(burnin:chain_length),signal_indices,drop=FALSE])>0.5)
      sss_fdr <- mean(colMeans(sss_chain$z[c(burnin:chain_length),!signal_indices,drop=FALSE])>0.5)
      sss_mse <- mean((colMeans(sss_chain$beta[c(burnin:chain_length),])-syn_data$true_beta)^2)
      
      print(c(n,p,s0,error_std))
      
      return(data.frame(algo='ScalableSpikeSlab', time=as.double(sss_time_taken[1])/chain_length, tpr=sss_tpr, fdr=sss_fdr, mse=sss_mse,
                        t = c(burnin:chain_length), delta_t_values = delta_t[burnin:chain_length], no_active_values = no_active[burnin:chain_length],
                        p_t_values = p_t[burnin:chain_length], n=n, p=p, s0=s0, iteration=i))
    }
}

## Fix n vary p ##
n_p_error_s0_grid <- data.frame(n=100, p=seq(100,3000,100), error_std=2, s0=10)
n_p_error_s0_list <- split(n_p_error_s0_grid, 1:nrow(n_p_error_s0_grid))
vary_p_df <- generate_sims(n_p_error_s0_list)
vary_p_summary <-
  vary_p_df %>% 
  dplyr::select(t, p_t_values, n, p, s0) %>%
  reshape2::melt(id=c('t','n', 'p')) %>%
  group_by(n,p,variable) %>% 
  summarise(mean=mean(value), sd=sd(value), length=n())
# save(vary_p_summary, file = 'inst/pt_scaling_simulations/vary_p_summary.Rdata')

## Fix p vary n ##
n_p_error_s0_grid <- data.frame(n=seq(100,1000,100), p=1000, error_std=2, s0=10)
n_p_error_s0_list <- split(n_p_error_s0_grid, 1:nrow(n_p_error_s0_grid))
vary_n_df <- generate_sims(n_p_error_s0_list)
vary_n_summary <-
  vary_n_df %>% 
  dplyr::select(t, p_t_values, n, p, s0) %>%
  reshape2::melt(id=c('t','n', 'p')) %>%
  group_by(n,p,variable) %>% 
  summarise(mean=mean(value), sd=sd(value), length=n())
# save(vary_n_summary, file = 'inst/pt_scaling_simulations/vary_n_summary.Rdata')

## Vary sparsity ##
s0_seq <- c(1:20)
# n_p_error_s0_grid <- data.frame(n=10*s0_seq, p=100*s0_seq, error_std=2, s0=s0_seq)
n_p_error_s0_grid <- data.frame(n=100, p=1000, error_std=2, s0=s0_seq)
n_p_error_s0_list <- split(n_p_error_s0_grid, 1:nrow(n_p_error_s0_grid))
vary_s_df <- generate_sims(n_p_error_s0_list)
vary_s_summary <-
  vary_s_df %>% 
  dplyr::select(t, p_t_values, n, p, s0) %>%
  dplyr::mutate(sparsity=s0) %>%
  reshape2::melt(id=c('t','n','p','s0')) %>%
  group_by(n,p,s0,variable) %>% 
  summarise(mean=mean(value), sd=sd(value), length=n())
# save(vary_s_summary, file = 'inst/pt_scaling_simulations/vary_s_summary.Rdata')


###### Plots ######
# load('inst/pt_scaling_simulations/vary_p_summary.Rdata')
vary_p_plot <- 
  ggplot(vary_p_summary, aes(x=p, y=mean, linetype=variable)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd,
                  fill=variable),alpha=0.2, colour = NA) +
  xlab(TeX('Dimension $p$')) + ylab(TeX('Average cost parameter $p_t$')) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("p_t_values", "s0"),
                        labels=unname(TeX(c('$p_t$','Sparsity $s$'))),
                        values = c('solid','dotted')) +
  scale_fill_manual(values = c('black','black', 'black'),guide='none') +
  scale_y_continuous(limits = c(0,12), breaks = seq(0,20,4)) +
  scale_x_continuous(limits = c(0,3e3), breaks = seq(0,5e3,1000)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=10)) +
  # theme(plot.margin=unit(c(1,1,2,0.5),"cm")) +
  # theme(legend.position = c(0.5, -0.1), legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
vary_p_plot
# ggsave(filename = "inst/pt_scaling_simulations/vary_p_plot.pdf", plot = vary_p_plot, width = 4, height = 2.5)


# load('inst/pt_scaling_simulations/vary_n_summary.Rdata')
vary_n_plot <- 
  ggplot(vary_n_summary, aes(x=n, y=mean, linetype=variable)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd,
                  fill=variable),alpha=0.2, colour = NA) +
  xlab(TeX('No. observations $n$')) + ylab(TeX('')) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("p_t_values", "s0"),
                        labels=unname(TeX(c('$p_t$','Sparsity $s$'))),
                        values = c('solid', 'dotted')) +
  scale_fill_manual(values = c('black','black', 'black'),guide='none') +
  scale_y_continuous(limits = c(0,12), breaks = seq(0,20,4)) +
  scale_x_continuous(breaks = seq(100,1000,300)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=10)) +
  # theme(plot.margin=unit(c(1,1,2,0.5),"cm")) +
  # theme(legend.position = c(0.5, -0.1), legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
vary_n_plot
# ggsave(filename = "inst/pt_scaling_simulations/vary_n_plot.pdf", plot = vary_n_plot, width = 4, height = 2.5)


# load('inst/pt_scaling_simulations/vary_s_summary.Rdata')
vary_s_plot <- 
  ggplot(vary_s_summary, aes(x=s0, y=mean, linetype=variable)) + 
  geom_line(size=0.5) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd,
                  fill=variable),alpha=0.2, colour = NA) +
  xlab(TeX('Sparsity $s$')) + ylab(TeX('')) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("p_t_values", "sparsity"),
                        labels=unname(TeX(c('$p_t$','Sparsity $s$'))),
                        values = c('solid', 'dotted')) +
  scale_fill_manual(values = c('black','black', 'black'),guide='none') +
  scale_y_continuous(limits = c(0,25), breaks = seq(0,30,5)) +
  scale_x_continuous(limits = c(0,20), breaks = seq(0,20,10)) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=10)) +
  # theme(plot.margin=unit(c(1,1,2,0.5),"cm")) +
  # theme(legend.position = c(0.5, -0.1), legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=1,byrow=TRUE))
vary_s_plot
# ggsave(filename = "inst/pt_scaling_simulations/vary_s_plot.pdf", plot = vary_s_plot, width = 4, height = 2.5)

scaling_plot <- 
  ggarrange(vary_p_plot, vary_n_plot, vary_s_plot, common.legend = TRUE, 
            legend = "bottom", nrow=1)
scaling_plot
# ggsave(filename = "inst/pt_scaling_simulations/scaling_plot_arxiv.pdf", plot = scaling_plot, width = 8, height = 4)


# scaling_plot <-
#   ggarrange(vary_p_plot + theme_classic(base_size = 16) + theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=16)),
#             vary_n_plot + theme_classic(base_size = 16) + theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=16)),
#             vary_s_plot + theme_classic(base_size = 16) + theme(legend.position = 'bottom', legend.key.width=unit(1,"cm"), legend.text = element_text(size=16)),
#             common.legend = TRUE,
#             legend = "bottom", nrow=1)
# ggsave(filename = "inst/pt_scaling_simulations/scaling_plot_icml.pdf", plot = scaling_plot, width = 8, height = 5)

# scaling_plot2 <- 
#   ggarrange(vary_p_plot + ylab(TeX('Average cost parameter $p_t$')), 
#             vary_n_plot + ylab(TeX('Average cost parameter $p_t$')), 
#             vary_s_plot + ylab(TeX('Average cost parameter $p_t$')), common.legend = TRUE, 
#             legend = "bottom", nrow=3)
# scaling_plot2
# ggsave(filename = "inst/pt_scaling_simulations/scaling_plot2.pdf", plot = scaling_plot2, width = 4, height = 6)



