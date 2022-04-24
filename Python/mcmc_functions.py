# Functions for Scalable Spike-and-Slab

import numpy as np
import scipy.linalg
from scipy.stats import norm, gamma, truncnorm


# Helper functions for MCMC
def matrix_full(Xt,d):
    n = Xt.shape[1]
    d_Xt = Xt / (d**0.5)[:, np.newaxis]
    return(d_Xt.T @ d_Xt + np.eye(n))

def matrix_precompute(Xt, d, prev_d, prev_matrix, XXt, tau0, tau1):
    n = Xt.shape[1]
    p = Xt.shape[0]
    no_swaps = np.sum(d != prev_d)
    slab_indices = (d != 1/tau0**2)
    spike_indices = (d == 1 / tau0 ** 2)
    no_slabs = np.sum(slab_indices)
    no_spikes = p-no_slabs

    if no_swaps==0:
        return(prev_matrix)
    if no_slabs==0:
        return(np.eye(n) + XXt * (tau0 ** 2))

    if no_swaps<min(no_slabs,no_spikes):
        pos_indices = (1/d - 1/prev_d) > 0
        if np.sum(pos_indices)==0:
            pos_part = 0
        else:
            pos_d_update = (1 / d - 1 / prev_d)[pos_indices,np.newaxis]
            pos_Xt_update = Xt[pos_indices,:]
            pos_Xt_d_update = pos_Xt_update * (pos_d_update**0.5)
            pos_part = pos_Xt_d_update.T @ pos_Xt_d_update

        neg_indices = (1/d - 1/prev_d) < 0
        if np.sum(neg_indices)==0:
            neg_part = 0
        else:
            neg_d_update = -(1 / d - 1 / prev_d)[neg_indices,np.newaxis]
            neg_Xt_update = Xt[neg_indices, :]
            neg_Xt_d_update = neg_Xt_update * (neg_d_update ** 0.5)
            neg_part = neg_Xt_d_update.T @ neg_Xt_d_update
        return(pos_part - neg_part + prev_matrix)
    elif no_slabs<no_spikes:
        Xt_update = Xt[slab_indices,:]
        update_part = (Xt_update.T @ Xt_update)*(tau1**2-tau0**2)
        return(np.eye(n) + XXt*(tau0**2) + update_part)
    else:
        Xt_update = Xt[spike_indices,:]
        update_part = (Xt_update.T @ Xt_update)*(tau0**2-tau1**2)
        return(np.eye(n) + XXt*(tau1**2) + update_part)

def inverse_precompute(Xt, d, prev_d, prev_inverse, tau0_inverse, tau1_inverse, tau0, tau1):
    n = Xt.shape[1]
    p = Xt.shape[0]
    swap_indices = (d!=prev_d)
    no_swaps = np.sum(swap_indices)
    slab_indices = (d!=1/tau0**2)
    spike_indices = (d==1/tau0**2)
    no_slabs = np.sum(slab_indices)
    no_spikes = p - no_slabs

    if no_swaps==0:
        return(prev_inverse)
    if no_slabs==0:
        return(tau0_inverse)
    if no_spikes==0:
        return(tau1_inverse)

    if no_swaps<min(no_slabs,no_spikes):
        diag_diff = (1/d - 1/prev_d)[swap_indices]
        Xt_update = Xt[swap_indices,:]
        Xt_update_prev_inverse = Xt_update @ prev_inverse
        M_matrix = np.diag(1/diag_diff) + Xt_update_prev_inverse@Xt_update.T
        M_matrix_inverse = np.linalg.inv((M_matrix+M_matrix.T)/2)
        inverse_update = Xt_update_prev_inverse.T @ (M_matrix_inverse @ Xt_update_prev_inverse)
        return(prev_inverse-inverse_update)
    elif no_slabs<no_spikes:
        Xt_update = Xt[slab_indices,:]
        Xt_update_prev_inverse = Xt_update @ tau0_inverse
        M_matrix = np.eye(no_slabs)/(tau1**2-tau0**2) + (Xt_update_prev_inverse @ Xt_update.T)
        # Use (M_matrix+t(M_matrix))/2 to avoid numerical error, as M_matrix symmetric
        M_matrix_inverse = np.linalg.inv((M_matrix + M_matrix.T)/2)
        inverse_update = Xt_update_prev_inverse.T @ (M_matrix_inverse @ Xt_update_prev_inverse)
        return(tau0_inverse-inverse_update)
    else:
        Xt_update = Xt[spike_indices,:]
        Xt_update_prev_inverse = Xt_update @ tau1_inverse
        M_matrix = np.eye(no_spikes)/(tau0**2-tau1**2) + (Xt_update_prev_inverse @ Xt_update.T)
        # Use (M_matrix+t(M_matrix))/2 to avoid numerical error, as M_matrix symmetric
        M_matrix_inverse = np.linalg.inv((M_matrix + M_matrix.T) / 2)
        inverse_update = Xt_update_prev_inverse.T @ (M_matrix_inverse @ Xt_update_prev_inverse)
        return(tau1_inverse - inverse_update)



## Spike slab linear regression functions ##
# beta update
def update_beta(z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                XXt, tau0_inverse, tau1_inverse, tau0, tau1, u=None, xi=None):
    p = len(z)
    n = len(y)
    swap_indices = (z != prev_z)
    no_swaps = np.sum(swap_indices)
    slab_indices = (z != 0)
    no_slabs = np.sum(slab_indices)
    no_spikes = p-no_slabs

    d = z/(tau1**2) + (1-z)/(tau0**2)
    prev_d = prev_z/(tau1**2) + (1-prev_z)/(tau0**2)

    if u is None:
        u = np.random.normal(size=p)
    u = u/(d**0.5)
    if xi is None:
        xi = np.random.normal(size=n)
    v = X@u + xi
    M_matrix = matrix_precompute(Xt, d, prev_d, prev_matrix, XXt, tau0, tau1)

    if min(no_swaps,no_spikes,no_slabs)<n:
        M_matrix_inverse = inverse_precompute(Xt, d, prev_d, prev_inverse, tau0_inverse, tau1_inverse, tau0, tau1)
    else:
        M_matrix_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(M_matrix, lower=True), np.eye(n))
        # M_matrix_inverse = np.linalg.inv(M_matrix)

    v_star = M_matrix_inverse@(y/(sigma2**0.5)-v)
    beta = (sigma2**0.5)*(u + (Xt@v_star)/d)
    return({'beta':beta,'matrix':M_matrix,'matrix_inverse':M_matrix_inverse})

# z update
def update_z(beta, sigma2, tau0, tau1, q, u_crn=None):
    p = len(beta)
    if u_crn is None:
        u_crn = np.random.uniform(size=p)
    log_prob1 = np.log(q) + norm.logpdf(beta,scale=tau1*(sigma2**0.5))
    log_prob2 = np.log(1-q) + norm.logpdf(beta,scale=tau0*(sigma2**0.5))
    probs = 1/(1+np.exp(log_prob2-log_prob1))
    z = np.where(u_crn<probs,1,0)
    return(z)

# sigma2 update
def update_sigma2(beta, z, tau0, tau1, a0, b0, X, y, u_crn=None):
    n = len(y)
    p = len(beta)
    rss = np.sum((y-X@beta)**2)
    d = z/(tau1**2) + (1-z)/(tau0**2)
    beta2_d = np.sum((beta**2)*d)

    if u_crn is None:
        sigma2 = 1/np.random.gamma(shape=((a0+n+p)/2),scale=(2/(b0+rss+beta2_d)),size=1)
    else:
        log_sigma2 = -np.log(gamma.ppf(q=u_crn,a=((a0+n+p)/2),scale=(2/(b0+rss+beta2_d))))
        sigma2 = np.exp(log_sigma2)
    return(sigma2)

# spike_slab_linear_kernel
def spike_slab_linear_kernel(beta, z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                             XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, a0, b0, random_samples):
    beta_output = update_beta(z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                              XXt, tau0_inverse, tau1_inverse, tau0, tau1,
                              u=random_samples['beta_u'], xi=random_samples['beta_xi'])
    mat = beta_output['matrix']
    mat_inverse = beta_output['matrix_inverse']
    beta_new = beta_output['beta']
    z_new = update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples['z_u'])
    sigma2_new = update_sigma2(beta_new, z_new, tau0, tau1, a0, b0, X, y, u_crn=random_samples['sigma2_u'])
    return ({'beta': beta_new, 'z': z_new, 'sigma2':sigma2_new, 'prev_z':z,
             'prev_matrix': mat, 'prev_inverse':mat_inverse})

# spike_slab_linear
def spike_slab_linear(chain_length, X, y, tau0, tau1, q, a0, b0, burnin=0, rinit=None, verbose=False, store=True):
    n = X.shape[0]
    p = X.shape[1]
    Xt = X.T
    XXt = X@Xt
    tau0_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(np.eye(n) + (tau0 ** 2) * XXt, lower=True), np.eye(n))
    tau1_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(np.eye(n) + (tau1 ** 2) * XXt, lower=True), np.eye(n))

    # Initializing from the prior
    if rinit is None:
        z = np.random.binomial(1,q,p)
        sigma2 = 1/np.random.gamma(shape=(a0/2),scale=(2/b0),size=1)
        beta = np.random.normal(size=p)*(z*tau1+(1-z)*tau0)*(sigma2**0.5)
        prev_z = z
        # prev_d = prev_z/(tau1**2) + (1-prev_z)/(tau0 ** 2)
        # prev_matrix = matrix_full(Xt, prev_d)
        Xt_update = Xt[prev_z == 1, :]
        prev_matrix = np.eye(n) + XXt * (tau0 ** 2) + (Xt_update.T @ Xt_update) * (tau1 ** 2 - tau0 ** 2)
        prev_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(prev_matrix, lower=True), np.eye(n))

    if store:
        beta_samples = np.zeros(shape=[(chain_length-burnin),p])
        z_samples = np.zeros(shape=[(chain_length-burnin),p])
        sigma2_samples = np.zeros(shape=[(chain_length-burnin),1])
    else:
        beta_ergodic_sum = np.zeros(shape=[p])
        z_ergodic_sum = np.zeros(shape=[p])
        sigma2_ergodic_sum = np.zeros(shape=[1])

    for t in range(chain_length):
        random_samples = {'beta_u':np.random.normal(size=p), 'beta_xi':np.random.normal(size=n),
                          'z_u':np.random.uniform(size=p), 'sigma2_u':np.random.uniform(size=1)}
        new_state = spike_slab_linear_kernel(beta, z, sigma2, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                                             XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, a0, b0, random_samples)
        beta = new_state['beta']
        z = new_state['z']
        sigma2 = new_state['sigma2']
        prev_z = new_state['prev_z']
        prev_matrix = new_state['prev_matrix']
        prev_inverse = new_state['prev_inverse']

        if t>=burnin:
            if store:
                beta_samples[(t-burnin),] = beta
                z_samples[(t - burnin),] = z
                sigma2_samples[(t - burnin),] = sigma2
            else:
                beta_ergodic_sum += beta
                z_ergodic_sum += z
                sigma2_ergodic_sum += sigma2

        if verbose:
            print(t)

    if store:
        return({'beta':beta_samples, 'z':z_samples, 'sigma2':sigma2_samples})
    else:
        return({'beta_ergodic_avg': beta_ergodic_sum/(chain_length-burnin),
                'z_ergodic_avg': z_ergodic_sum/(chain_length-burnin),
                'sigma2_ergodic_avg': sigma2_ergodic_sum/(chain_length-burnin)})


## Spike slab probit regression functions ##
def update_e(beta, sigma2, y, X, u_crn=None):
    n = len(y)
    if u_crn is None:
        u_crn = np.random.uniform(size=n)
    means = X@beta
    sds = sigma2**0.5
    ubs = np.zeros(n)+np.inf
    ubs[y==0] = 0
    lbs = np.zeros(n)-np.inf
    lbs[y==1] = 0
    # e = truncnorm.ppf(u_crn, a=lbs, b=ubs, loc=means, scale=sds)
    e = truncnorm.ppf(u_crn, a=(lbs-means)/sds, b=(ubs-means)/sds)*sds + means
    return(e)

def spike_slab_probit_kernel(beta, z, e, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                             XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, random_samples):
    sigma2 = 1
    beta_output = update_beta(z, sigma2, X, Xt, e, prev_z, prev_matrix, prev_inverse, XXt, tau0_inverse, tau1_inverse, tau0, tau1,
                              u=random_samples['beta_u'], xi=random_samples['beta_xi'])
    mat = beta_output['matrix']
    mat_inverse = beta_output['matrix_inverse']
    beta_new = beta_output['beta']
    z_new = update_z(beta_new, sigma2, tau0, tau1, q, u_crn=random_samples['z_u'])
    e_new = update_e(beta_new, sigma2, y, X, u_crn=random_samples['e_u'])
    return({'beta': beta_new, 'z': z_new, 'e': e_new,
            'prev_z': z, 'prev_matrix': mat, 'prev_inverse': mat_inverse})

def spike_slab_probit(chain_length, X, y, tau0, tau1, q, burnin=0, rinit=None, verbose=False, store=True):
    n = X.shape[0]
    p = X.shape[1]
    Xt = X.T
    XXt = X@Xt
    tau0_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(np.eye(n) + (tau0 ** 2) * XXt, lower=True), np.eye(n))
    tau1_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(np.eye(n) + (tau1 ** 2) * XXt, lower=True), np.eye(n))

    # Initializing from the prior
    if rinit is None:
        z = np.random.binomial(1, q, p)
        beta = np.random.normal(size=p) * (z * tau1 + (1 - z) * tau0)
        e = X @ beta + np.random.normal(size=n)
        prev_z = z
        # prev_d = prev_z/(tau1**2) + (1-prev_z)/(tau0 ** 2)
        # prev_matrix = matrix_full(Xt, prev_d)
        Xt_update = Xt[prev_z == 1, :]
        prev_matrix = np.eye(n) + XXt * (tau0 ** 2) + (Xt_update.T @ Xt_update) * (tau1 ** 2 - tau0 ** 2)
        prev_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(prev_matrix, lower=True), np.eye(n))

    if store:
        beta_samples = np.zeros(shape=[(chain_length-burnin),p])
        z_samples = np.zeros(shape=[(chain_length-burnin),p])
        e_samples = np.zeros(shape=[(chain_length - burnin), n])
    else:
        beta_ergodic_sum = np.zeros(shape=[p])
        z_ergodic_sum = np.zeros(shape=[p])
        e_ergodic_sum = np.zeros(shape=[n])

    for t in range(chain_length):
        random_samples = {'beta_u':np.random.normal(size=p), 'beta_xi':np.random.normal(size=n),
                          'z_u':np.random.uniform(size=p), 'e_u':np.random.uniform(size=n)}

        new_state = spike_slab_probit_kernel(beta, z, e, X, Xt, y, prev_z, prev_matrix, prev_inverse,
                                             XXt, tau0_inverse, tau1_inverse, tau0, tau1, q, random_samples)
        beta = new_state['beta']
        z = new_state['z']
        e = new_state['e']
        prev_z = new_state['prev_z']
        prev_matrix = new_state['prev_matrix']
        prev_inverse = new_state['prev_inverse']

        if t>=burnin:
            if store:
                beta_samples[(t-burnin),] = beta
                z_samples[(t - burnin),] = z
                e_samples[(t - burnin),] = e
            else:
                beta_ergodic_sum += beta
                z_ergodic_sum += z
                e_ergodic_sum += e

        if verbose:
            print(t)

    if store:
        return({'beta':beta_samples, 'z':z_samples, 'e':e_samples})
    else:
        return({'beta_ergodic_avg': beta_ergodic_sum/(chain_length-burnin),
                'z_ergodic_avg': z_ergodic_sum/(chain_length-burnin),
                'e_ergodic_avg': e_ergodic_sum/(chain_length-burnin)})

## Spike slab logistic regression functions ##

# Output D%*%M%*%D for symmetric matrix M and a diagonal matrix D
def diagonal_inner_outer_prod(symetric_mat, diag_vec):
    return((symetric_mat*diag_vec[:, np.newaxis]).T*diag_vec[:, np.newaxis])

def matrix_full_logistic(Xt, d, sigma2tilde):
    n = Xt.shape[1]
    d_Xt = Xt/(d**0.5)[:, np.newaxis]
    return(diagonal_inner_outer_prod(d_Xt.T@d_Xt,1/(sigma2tilde**0.5))+np.eye(n))

def matrix_precompute_logistic(Xt, d, sigma2tilde,
                               prev_d, prev_matrix, prev_sigma2tilde,
                               XXt, tau0, tau1):
    p = Xt.shape[0]
    swap_indices = (d != prev_d)
    no_swaps = np.sum(swap_indices)
    slab_indices = (d != 1/tau0**2)
    spike_indices = (d == 1/tau0**2)
    no_slabs = np.sum(slab_indices)
    no_spikes = p - no_slabs

    if no_swaps == 0:
        np.fill_diagonal(prev_matrix, prev_matrix.diagonal() - 1)
        prev_matrix = diagonal_inner_outer_prod(prev_matrix, prev_sigma2tilde**0.5)
        new_matrix = diagonal_inner_outer_prod(prev_matrix, 1/sigma2tilde**0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal() + 1)
        return(new_matrix)
    if no_slabs == 0:
        new_matrix = diagonal_inner_outer_prod((tau0**2)*XXt, 1/sigma2tilde**0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal() + 1)
        return(new_matrix)
    if no_spikes == 0:
        new_matrix = diagonal_inner_outer_prod((tau1**2)*XXt, 1/sigma2tilde**0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal() + 1)
        return(new_matrix)

    if no_swaps<min(no_slabs,no_spikes):
        pos_indices = (1/d - 1/prev_d) > 0
        if np.sum(pos_indices)==0:
            pos_part = 0
        else:
            pos_d_update = (1 / d - 1 / prev_d)[pos_indices,np.newaxis]
            pos_Xt_update = Xt[pos_indices,:]
            pos_Xt_d_update = pos_Xt_update * (pos_d_update**0.5)
            pos_part = pos_Xt_d_update.T @ pos_Xt_d_update

        neg_indices = (1/d - 1/prev_d) < 0
        if np.sum(neg_indices)==0:
            neg_part = 0
        else:
            neg_d_update = -(1 / d - 1 / prev_d)[neg_indices,np.newaxis]
            neg_Xt_update = Xt[neg_indices, :]
            neg_Xt_d_update = neg_Xt_update * (neg_d_update ** 0.5)
            neg_part = neg_Xt_d_update.T @ neg_Xt_d_update
        np.fill_diagonal(prev_matrix, prev_matrix.diagonal() - 1)
        prev_matrix = diagonal_inner_outer_prod(prev_matrix, prev_sigma2tilde ** 0.5)
        prev_matrix = pos_part - neg_part + prev_matrix
        new_matrix = diagonal_inner_outer_prod(prev_matrix, 1 / sigma2tilde ** 0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal() + 1)
    elif no_slabs<no_spikes:
        Xt_update = Xt[slab_indices,:]
        update_part = (Xt_update.T @ Xt_update)*(tau1**2-tau0**2)
        new_matrix = diagonal_inner_outer_prod(XXt*(tau0**2)+update_part, 1/sigma2tilde**0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal()+1)
    else:
        Xt_update = Xt[spike_indices,:]
        update_part = (Xt_update.T @ Xt_update)*(tau0**2-tau1**2)
        new_matrix = diagonal_inner_outer_prod(XXt*(tau1**2)+update_part, 1/sigma2tilde**0.5)
        np.fill_diagonal(new_matrix, new_matrix.diagonal() + 1)
    return(new_matrix)


def update_beta_logistic(z, sigma2tilde, X, Xt, y, prev_z, prev_sigma2tilde, prev_matrix,
                         XXt, tau0, tau1, u=None, xi=None):
    p = len(z)
    n = len(y)
    swap_indices = (z != prev_z)
    no_swaps = np.sum(swap_indices)
    slab_indices = (z != 0)
    no_slabs = np.sum(slab_indices)
    no_spikes = p - no_slabs

    d = z / (tau1 ** 2) + (1 - z) / (tau0 ** 2)
    prev_d = prev_z / (tau1 ** 2) + (1 - prev_z) / (tau0 ** 2)

    if u is None:
        u = np.random.normal(size=p)
    u = u / (d ** 0.5)
    if xi is None:
        xi = np.random.normal(size=n)
    v = u@(Xt/(sigma2tilde**0.5)) + xi
    M_matrix = matrix_precompute_logistic(Xt, d, sigma2tilde, prev_d, prev_matrix, prev_sigma2tilde, XXt, tau0, tau1)
    M_matrix_inverse = scipy.linalg.cho_solve(scipy.linalg.cho_factor(M_matrix, lower=True), np.eye(n))
    # M_matrix_inverse = np.linalg.inv(M_matrix)

    v_star = M_matrix_inverse @ (y/(sigma2tilde**0.5) - v)
    beta = u + (Xt @ (v_star/(sigma2tilde**0.5))) / d
    return ({'beta': beta, 'matrix': M_matrix})


def update_sigma2tilde_logistic(beta, z, tau0, tau1, a0, b0, X, y, u_crn=None):
    n = len(y)
    rss = (y-X@beta)**2
    d = z / (tau1 ** 2) + (1 - z) / (tau0 ** 2)

    if u_crn is None:
        sigma2tilde = 1 / np.random.gamma(shape=((a0 + 1) / 2), scale=(2 / (b0 + rss)))
    else:
        log_sigma2tilde = -np.log(gamma.ppf(q=u_crn, a=((a0 + 1) / 2), scale=(2 / (b0 + rss))))
        sigma2tilde = np.exp(log_sigma2tilde)
    return(sigma2tilde)

def spike_slab_logistic_kernel(beta, z, e, sigma2tilde, X, Xt, y, prev_z, prev_sigma2tilde, prev_matrix,
                               XXt, tau0, tau1, q, random_samples):
    beta_output = update_beta_logistic(z, sigma2tilde, X, Xt, e, prev_z, prev_sigma2tilde, prev_matrix,
                                       XXt, tau0, tau1, u=random_samples['beta_u'], xi=random_samples['beta_xi'])
    mat = beta_output['matrix']
    beta_new = beta_output['beta']
    z_new = update_z(beta_new, 1, tau0, tau1, q, u_crn=random_samples['z_u'])
    e_new = update_e(beta_new, sigma2tilde, y, X, u_crn=random_samples['e_u'])
    nu = 7.3
    w2 = (np.pi**2)*(nu-2)/(3*nu)
    a0 = nu
    b0 = w2*nu
    sigma2tilde_new = update_sigma2tilde_logistic(beta_new, z_new, tau0, tau1, a0, b0, X, e_new, u_crn=random_samples['sigma2tilde_u'])
    return({'beta':beta_new,'z':z_new,'e':e_new,'sigma2tilde': sigma2tilde_new,
            'prev_z': z, 'prev_sigma2tilde': sigma2tilde, 'prev_matrix': mat})


def spike_slab_logistic(chain_length, X, y, tau0, tau1, q, burnin=0, rinit=None, verbose=False, store=True):
    n = X.shape[0]
    p = X.shape[1]
    Xt = X.T
    XXt = X@Xt

    # Initializing from the prior
    if rinit is None:
        z = np.random.binomial(1,q,p)
        nu = 7.3
        w2 = (np.pi ** 2) * (nu - 2) / (3 * nu)
        a0 = nu
        b0 = w2 * nu
        sigma2tilde = 1/np.random.gamma(shape=(a0/2),scale=(2/b0),size=n)
        beta = np.random.normal(size=p)*(z*tau1+(1-z)*tau0)
        e = X@beta + (sigma2tilde**0.5)*np.random.normal(size=n)
        prev_z = z
        prev_sigma2tilde = sigma2tilde
        Xt_update = Xt[prev_z == 1, :]
        update_part = (Xt_update.T @ Xt_update) * (tau1 ** 2 - tau0 ** 2)
        prev_matrix = diagonal_inner_outer_prod(XXt * (tau0 ** 2) + update_part, 1 / sigma2tilde ** 0.5)
        np.fill_diagonal(prev_matrix, prev_matrix.diagonal() + 1)

    if store:
        beta_samples = np.zeros(shape=[(chain_length-burnin),p])
        z_samples = np.zeros(shape=[(chain_length-burnin),p])
        e_samples = np.zeros(shape=[(chain_length - burnin), n])
        sigma2tilde_samples = np.zeros(shape=[(chain_length-burnin),n])
    else:
        beta_ergodic_sum = np.zeros(shape=[p])
        z_ergodic_sum = np.zeros(shape=[p])
        e_ergodic_sum = np.zeros(shape=[n])
        sigma2tilde_ergodic_sum = np.zeros(shape=[n])

    for t in range(chain_length):
        random_samples = {'beta_u':np.random.normal(size=p), 'beta_xi':np.random.normal(size=n),
                          'z_u':np.random.uniform(size=p), 'e_u':np.random.uniform(size=n),
                          'sigma2tilde_u':np.random.uniform(size=n)}
        new_state = spike_slab_logistic_kernel(beta, z, e, sigma2tilde, X, Xt, y, prev_z, prev_sigma2tilde, prev_matrix,
                                               XXt, tau0, tau1, q, random_samples)
        beta = new_state['beta']
        z = new_state['z']
        e = new_state['e']
        sigma2tilde = new_state['sigma2tilde']
        prev_z = new_state['prev_z']
        prev_matrix = new_state['prev_matrix']
        prev_sigma2tilde = new_state['prev_sigma2tilde']

        if t>=burnin:
            if store:
                beta_samples[(t-burnin),] = beta
                z_samples[(t - burnin),] = z
                e_samples[(t - burnin),] = e
                sigma2tilde_samples[(t - burnin),] = sigma2tilde
            else:
                beta_ergodic_sum += beta
                z_ergodic_sum += z
                e_ergodic_sum += e
                sigma2tilde_ergodic_sum += sigma2tilde

        if verbose:
            print(t)

    if store:
        return({'beta':beta_samples, 'z':z_samples, 'e':e_samples, 'sigma2tilde':sigma2tilde_samples})
    else:
        return({'beta_ergodic_avg': beta_ergodic_sum/(chain_length-burnin),
                'z_ergodic_avg': z_ergodic_sum/(chain_length-burnin),
                'e_ergodic_avg': e_ergodic_sum/(chain_length-burnin),
                'sigma2tilde_ergodic_avg': sigma2tilde_ergodic_sum/(chain_length-burnin)})