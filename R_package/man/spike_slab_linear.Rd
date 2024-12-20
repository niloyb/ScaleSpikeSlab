% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{spike_slab_linear}
\alias{spike_slab_linear}
\title{spike_slab_linear}
\usage{
spike_slab_linear(
  chain_length,
  X,
  y,
  tau0,
  tau1,
  q,
  a0 = 1,
  b0 = 1,
  rinit = NULL,
  verbose = FALSE,
  burnin = 0,
  store = TRUE,
  Xt = NULL,
  XXt = NULL,
  tau0_inverse = NULL,
  tau1_inverse = NULL
)
}
\arguments{
\item{chain_length}{Markov chain length}

\item{X}{matrix of length n by p}

\item{y}{Response}

\item{tau0}{prior hyperparameter (non-negative real)}

\item{tau1}{prior hyperparameter (non-negative real)}

\item{q}{prior hyperparameter (strictly between 0 and 1)}

\item{a0}{prior hyperparameter (non-negative real)}

\item{b0}{prior hyperparameter (non-negative real)}

\item{rinit}{initial distribution of Markov chain (default samples from the prior)}

\item{verbose}{print iteration of the Markov chain (boolean)}

\item{burnin}{chain burnin (non-negative integer)}

\item{store}{store chain trajectory (boolean)}

\item{Xt}{Pre-calculated transpose of X}

\item{XXt}{Pre-calculated matrix X*transpose(X) (n by n matrix)}

\item{tau0_inverse}{Pre-calculated matrix inverse(I + tau0^2*XXt) (n by n matrix)}

\item{tau1_inverse}{Pre-calculated matrix inverse(I + tau1^2*XXt) (n by n matrix)}
}
\value{
Output from Markov chain targeting the posterior corresponding to
Bayesian linear regression with spike and slab priors
}
\description{
Generates Markov chain targeting the posterior corresponding to
Bayesian linear regression with spike and slab priors
}
\examples{
# Synthetic dataset
syn_data <- synthetic_data(n=100,p=200,s0=5,error_std=2,type='linear')
X <- syn_data$X
y <- syn_data$y

# Hyperparamters
params <- spike_slab_params(n=nrow(X),p=ncol(X))

# Run S^3
sss_chain <- spike_slab_linear(chain_length=4e3,burnin=1e3,X=X,y=y,
tau0=params$tau0,tau1=params$tau1,q=params$q,a0=params$a0,b0=params$b0,
verbose=FALSE,store=FALSE)

# Use posterior probabilities for variable selection
sss_chain$z_ergodic_avg[1:10]
}
