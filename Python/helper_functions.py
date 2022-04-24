# Helper functions to generate synthetic data and choose hyperparameters #

import numpy as np
from scipy.stats import logistic, binom


# Generate synthetic dataset
def synthetic_data(n, p, s0, error_std=1, type='linear', scale=True, signal='constant'):
    true_beta = np.zeros(shape=[p])
    s0 = min(p, s0)
    if s0>0:
        if signal=='constant':
            true_beta[range(s0)] = 2
        elif signal=='decay':
            true_beta[range(s0)] = 2**(-(np.arange(s0)+1-9)/4)
        else:
            return("Error: input parameter signal must be 'constant' or 'decay'")

    X = np.random.normal(size=[n,p])
    if scale:
        X = (X - X.mean(axis=0)) / X.std(axis=0)

    X_truebeta = X@true_beta

    if type=='linear':
        error_terms = np.random.normal(size=n)*error_std
        y = X_truebeta + error_terms
    elif type=='probit':
        true_aug_y = X_truebeta + np.random.normal(size=n)
        y = np.where(true_aug_y > 0, 1, 0)
    elif type=='logistic':
        true_aug_y = logistic.rvs(loc=X_truebeta)
        y = np.where(true_aug_y > 0, 1, 0)
    else:
        return("Error: input parameter 'type' must be 'linear' or 'logistic'")
    return({'X':X,'y':y,'true_beta':true_beta})

def spike_slab_params(n,p,type='linear'):
    K = max(10,np.log(5))
    q_seq = np.arange(1/p,(1-1/p),1/p)
    probs = abs(binom.cdf(k=K,n=p,p=q_seq)-0.9)
    q = min(q_seq[probs == np.min(probs)])
    tau0 = 1/(n**0.5)
    tau1 = 1
    # tau1 <- sqrt(max(1, p^(2.1)/(100*n))) # Alternative choice for tau1
    a0 = 1
    b0 = 1
    if type=='linear':
        return({'q':q,'tau0':tau0,'tau1':tau1,'a0':a0,'b0':b0})
    elif type=='probit':
        return({'q': q,'tau0':tau0,'tau1':tau1})
    elif type=='logistic':
        return({'q': q,'tau0':tau0,'tau1':tau1})
    else:
        return ("Error: input parameter 'type' must be 'linear' or 'logistic'")


