options(warn=2) # Turn warnings into errors

# Simulate data --------------------------------------

n = 300
w = c(0.3,0.4,0.2, 0.1)
y = sort(c(rnorm(w[1]*n, -10), rnorm(w[2]*n, 3), rnorm(w[3]*n, 10), rnorm(w[4]*n, 20)))
hist(y)
sigma2 = 1
a = 0
tau_02 = 1
v0 = 1
lambda0 = 1
sc = NULL
k_star = NULL
rho = NULL
set.seed(001)
m = 1 # auxiliary component
alpha = 1 # concentration parameter (DP)
# shape = 1; rate = 1 # hyperparameters in Ga(a, b) hyperior for alpha (s.t. E(alpha) = a/b)
# c = sample(1:10,n, replace=TRUE) # all observations into many clusters at the start
c = rep(1, n)
mu_j = tapply(y, c, mean)
sigma2_j = tapply(y, c, var)
MCMCiter = 1000
m0 = mean(y)
s2 = runif(1,0,2)
set.seed(001)
c_map = matrix(NA, nrow=MCMCiter, ncol=n)# for posteriori similarity matrix
mu_store = matrix(NA, nrow=MCMCiter, ncol=n) 

# Full conditionals --------------------------------------

simulate_mu_j = function(y, c, sigma2_j, m, s2) {
  
  sumY = tapply(y, c, sum)
  nj = table(c)
  denom =(nj/sigma2_j + 1/s2)
  
  new_mu = rnorm(length(nj),
                 (sumY/sigma2_j + m/s2)/denom,
                 1/sqrt(denom))
  
  return(new_mu)
}

simulate_sigma2_j = function(y, c, mu, s2) {
  
  sumYminusMuSq = tapply(y-mu, c, function(x) sum(x^2))
  nj = table(c)
  
  new_sig = 1/rgamma(length(nj),
                     nj/2,
                     sumYminusMuSq/2)
  
  return(new_sig)
}

update_m0 = function(mu_j, s2, length_rho) {
  sumMu_j = sum(mu_j)
  new_m = rnorm(1, sumMu_j/length_rho, sqrt(s2/length_rho))
  return(new_m)
}

update_s2 = function(mu_j, m0, length_rho) {
  sumMuminusmSq = sum((mu_j - m0)^2)
  new_s2 = 1/rgamma(1,
                    length_rho/2,
                    sumMuminusmSq/2)
  return(new_s2)
}

update_alpha = function(alpha, # current value of alpha 
                        shape, # hyperparam a in alpha ~ Ga(a,b)
                        rate,  # hyperparam b in alpha ~ Ga(a,b)
                        G,     # number of occupied clusters
                        n      # number of observations
) {    # Escobar & West
  shape2 <- shape + G - 1
  rate2 <- rate - log(rbeta(1, alpha + 1, n))
  weight <- shape2/(shape2 + n * rate2)
  return(weight * rgamma(1, shape = shape2 + 1, rate = rate2) + (1 - weight) * rgamma(1, shape = shape2, rate = rate2))
}
# Getting started with PPM --------- 

for(mcmc in 1:MCMCiter) {
  
  length_rho = length(unique(c)) # number of clusters
  mu_j = simulate_mu_j(y, c, sigma2_j, m, s2) # sample mean cluster
  mu = mu_j[c] # assign mean cluster to all obs
  sigma2_j = simulate_sigma2_j(y, c, mu, s2) # sample variance cluster
  m0 = update_m0(mu_j, s2, length_rho) # update hyperparameter
  s2 = update_s2(mu_j, m0, length_rho) # update hyperparameter
  # alpha = update_alpha(alpha, shape, rate, length_rho, n)
  
  for(i in 1:n) {
    labels_minus_i = unique(c[-i])
    n_ic = table(c[-i]) # number of obs in each cluster (without obs i)
    k = length(labels_minus_i) # number of distinct cluster (without obs i)
    h = k + m # number of new (empty clusters)
    c_i = c[i] # cluster label for i
    
    # Check whether i is a singleton (i.e., check whether the cluster has ONLY one observation)
    if(c_i %in% labels_minus_i) { # Not a singleton
      
      temp_mu_j = c(mu_j, rnorm(m, m0, s2)) # temporary mu_j
      temp_sigma2_j = c(sigma2_j, runif(m)) # temporary sigma2_j
      
      logp1 = log(n_ic) - log(n - 1 + alpha) + dnorm(y[i], temp_mu_j[1:k], temp_sigma2_j[1:k], log=TRUE) # log[p(c_i = c| -)] for 1 =< c <= k
      logp2 = log(alpha) - log(m) - log(n - 1 + alpha) + dnorm(y[i], temp_mu_j[-c(1:k)], temp_sigma2_j[-c(1:k)], log=TRUE) # log[p(c_i = c| -)] for k < c <= h
      logp2max = which.max(logp2) # auxiliary: picks the best one among new clusters
      log_binv = matrixStats::logSumExp(logp1, logp2[logp2max]) # (log) normalising constant
      log_pci = c(logp1, logp2[logp2max]) - log_binv # normalising the log probabilities
      new_c_i = which.max(log_pci - log(rexp(length(log_pci)))) # update cluster membership (Gumbel-max trick)
      c[i] = new_c_i
      
      #p1 = (n_ic/(n - 1 + alpha)) * dnorm(y[i], temp_mu_j[1:k], temp_sigma2_j[1:k]) # p(c_i = c| -) for 1 =< c <= k
      #p2 = ((alpha/m)/(n - 1 + alpha)) * dnorm(y[i], temp_mu_j[-c(1:k)], temp_sigma2_j[-c(1:k)]) # p(c_i = c| -) for k < c <= h
      #b_inv = sum(p1,p2) # b^(-1) normalising constant
      #pc_i = c(p1, p2)/b_inv # normalsing the probs so that sum to 1
      #new_c_i = extraDistr::rcat(1, prob = pc_i) # update cluster membership
      #c[i] = new_c_i
      
      if(new_c_i == (length(mu_j) + 1)) { # Was the obs allocated into the newest cluster?
        mu_j = temp_mu_j[c(1:k, k+logp2max)]
        sigma2_j = temp_sigma2_j[c(1:k, k+logp2max)]
      }
      check = 'Not a singleton' # don't need this. I put it when I was debugging
      
    } else { # it's a singleton
      temp_mu_j = mu_j[-c_i]
      temp_sigma2_j = sigma2_j[-c_i]
      logp1 = log(n_ic) - log(n - 1 + alpha) + dnorm(y[i], temp_mu_j, temp_sigma2_j, log=TRUE) # log[p(c_i = c| -)] for 1 =< c <= k
      logp2 = log(alpha) - log(m) - log(n - 1 + alpha) + dnorm(y[i], mu_j[c_i], sigma2_j[c_i], log=TRUE) # log[p(c_i = c| -)] for k < c <= h (probability of keeping obs i by itself in the current cluster)
      aux = rep(NA, length(mu_j))
      aux[c_i] = logp2 # allocating the probabilities correctly
      aux[-c_i] = logp1
      log_binv = matrixStats::logSumExp(aux) # (log) normalising constant
      log_pci = aux - log_binv # normalsing the probs so that they sum to 1
      new_c_i = which.max(log_pci - log(rexp(length(log_pci)))) # update cluster membership (Gumbel-max trick)
      c[i] = new_c_i
      
      #p1 = (n_ic/(n - 1 + alpha)) * dnorm(y[i], temp_mu_j, temp_sigma2_j) # p(c_i = c| -) for 1 =< c <= k
      #p2 = ((alpha/m)/(n - 1 + alpha)) * dnorm(y[i], mu_j[c_i], sigma2_j[c_i]) # p(c_i = c| -) for k < c <= h (probability of keeping obs i by itself in the current cluster)
      #aux = rep(NA, length(mu_j))
      #aux[c_i] = p2 # allocating the probabilities correctly
      #aux[-c_i] = p1
      #b_inv = sum(aux) # b^(-1) normalising constant
      #pc_i = aux/b_inv # normalsing the probs so that sum to 1
      #new_c_i = extraDistr::rcat(1, prob = pc_i) # update cluster membership
      #c[i] = new_c_i
      
      if(new_c_i != c_i) {
        mu_j = temp_mu_j
        sigma2_j = temp_sigma2_j
      }
      check = 'Singleton' # don't need this. I put it when I was debugging
    }
    
    # no-gap restriction ------
    cluster_labels = unique(c) # get labels
    
    while(length(cluster_labels) != max(cluster_labels)) {
      aux = seq_len(max(cluster_labels))
      label_gap = which(!(aux %in% c))
      
      if(length(label_gap) == 0) { break }
      obs_gap = which(c > label_gap) # identify observations that present gap in their label
      c[obs_gap] = c[obs_gap] - 1 # adjust the gap
      cluster_labels = unique(c) # get labels
    }
  }
  mu_store[mcmc,] = mu
  c_map[mcmc, ] = c
  print(mcmc)
}

plot(y)
points(y, col=c)

# Compute posterior similarity matrix -------- 

Zsimilarity  <- function(zs) {
  has.pkg    <- suppressMessages(requireNamespace("mcclust", quietly=TRUE))
  if(!has.pkg)                         stop("'mcclust' package not installed", call.=FALSE)
  if(!is.matrix(zs))                   stop("'zs' must be a matrix with rows corresponding to the number of observations and columns corresponding to the number of iterations", call.=FALSE)
  if(anyNA(zs))                        stop("Missing values are not allowed in 'zs'", call.=FALSE)
  zsim       <- mcclust::comp.psm(zs)
  mse.z      <- vapply(seq_len(nrow(zs)), function(i, x=mcclust::cltoSim(zs[i,]) - zsim) tryCatch(suppressWarnings(mean(x * x)), error=function(e) Inf), numeric(1L))
  Z.avg      <- zs[which.min(mse.z),]
  attr(Z.avg, "MSE")  <- min(mse.z)
  return(list(z.sim  = zsim, z.avg = Z.avg, MSE.z = mse.z))
}
simat = Zsimilarity(c_map)
image(simat$z.sim)
plot(y, col=simat$z.avg)

