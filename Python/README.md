| title       | author | date |
|-----------|--------|--------|
| ScaleSpikeSlab    | Niloy Biswas | 24 April 2021 |


# ScaleSpikeSlab (S^3)

This folder contains algorithms in Python for *Scalable Spike-and-Slab* (S^3),
a scalable Gibbs sampling implementation for high-dimensional Bayesian 
regression with the continuous spike-and-slab prior.

It is based on the article "Scalable Spike-and-Slab" (https://arxiv.org/abs/2204.01668), 
by Niloy Biswas, Lester Mackey and Xiao-Li Meng.

## A tutorial with Riboflavin GWAS data

#### Import S^3 functions
```python
from mcmc_functions import *
from helper_functions import *
```

#### Import data 

```python
# Riboflavin linear regression dataset of Buhlmann et al. (2014)
import pandas as pd
riboflavin = pd.read_csv('riboflavin.csv')
y = riboflavin['y'].to_numpy()
X = riboflavin.drop(['y','Unnamed: 0'], axis=1).to_numpy()

```

#### Run MCMC with S^3. 

```python
# Select hyperparamters
hyperparams = spike_slab_params(n=X.shape[0],p=X.shape[1])

chain_length=5000
burnin=1000
chain_output = spike_slab_linear(chain_length=chain_length, X=X, y=y, 
                                 rinit=None, verbose=False, 
                                 burnin=burnin, store=False, **hyperparams)
# Get the results
z_averages = chain_output['z_ergodic_avg']
```
#### Plot Spike-and-Slab marginal posterior probabilities for variable selection
```python
import numpy as np
import matplotlib.pyplot as plt


sorted_indices = np.argsort(np.mean(z_averages, axis=0))[::-1]
z_means = np.mean(z_averages, axis=0)[sorted_indices]
z_sds = (np.var(z_averages, axis=0)[sorted_indices]/no_chains)**0.5

plt.errorbar(x=np.arange(len(z_means)), y=z_means, yerr=z_sds,fmt='o')
plt.xscale('log')
plt.xlabel('Riboflavin Covariates')
plt.ylabel('Marginal posterior probabilities')
plt.show()
```
