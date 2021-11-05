# Mixture of normals
k = 2 # number of clusters
nMCMC = 10000
N = 200
q = c(0.3, 0.7)
n1 = floor(N*q[1])
n2 = floor(N*q[2])
mu = c(-10, 10)
tau = c(1,1)
x = c(rnorm(n1,mu[1],1/tau[1]), rnorm(n2,mu[2],1/tau[2]))
hist(x, freq=FALSE)
z = sample(1:k, N, replace=TRUE)
# z = c(rep(1,n1), rep(2,n2))
vj = 1
mj = 0
aj = 1
bj = 1
r = rep(1/k,k)

log_sum_exp = function(x) {
  # m = max(x);
  m = apply(x,1,max);
  # return(m + log(sum(exp(x - m))));
  return(m + log(rowSums(exp(x-m))));
}

# sample_z = function(m) max.col(m - log(rexp(length(m)))) # keefe's idea
sample_z = function(m) max.col(m -log(-log(runif(length(m))))) # keefe's idea
# store

sq = matrix(NA, ncol=k, nrow=nMCMC)
sz = matrix(NA, ncol=N, nrow=nMCMC)
smu = matrix(NA, ncol=k, nrow=nMCMC)
stau = matrix(NA, ncol=k, nrow=nMCMC)

for (i in 1:nMCMC){
  
  # Auxiliary quantities ------ 
  nj = table(z) # number of observation in each cluster
  sumX = aggregate(x, by = list(z), function(x) sum(x))[,2] # sufficient stats for full conditionals
  sumSqX = aggregate(x, by = list(z), function(x) sum(x^2))[,2] # sufficient stats for full conditionals
  
  # Sampling parameters -----
  Aj = nj/2 + aj
  Mj = (nj + 1/vj)^(-1) * (sumX + mj/vj)
  Bj = bj + mj^2/(2*vj) + sumSqX/2 - 0.5*(nj + 1/vj)*Mj^2
  tau = rgamma(k,Aj, Bj) # sample the variances of the Normals
  
  Vj = vj/((nj * vj + 1/vj)*tau)
  mu = rnorm(k,Mj, Vj)  # sample the means of the Normals
  
  alpha_post = matrix(r + nj, ncol=k)
  q = extraDistr::rdirichlet(1,alpha = alpha_post) # sample the weights of the mixtures
  
  pz = as.numeric(q)*mapply(function(m,s) dnorm(x,m,s), mu, sqrt(1/tau)) # p(Z_i = j| -) # naive calculation
  log_pz = log(as.numeric(q)) + mapply(function(m,s) dnorm(x,m,s, log=TRUE), mu, sqrt(1/tau)) # keefe's idea
  
  # omega = pz/apply(pz,1,sum)
  # omega = exp(log_pz - log_sum_exp(log_pz))
  z = sample_z(log_pz) # sampling z based on keefe's idea
  # z = extraDistr::rcat(N, prob = omega) # that works fine
  
  # Identifiability (q[1] > q[2] > ... > q[k])--- 
  # q_sorted = sort(q, index=TRUE)$ix # return the indices rather than the values
  # q = q[q_sorted]
  # mu = mu[q_sorted]
  # tau = tau[q_sorted]
  # should omega be changed as well?
  
  # Store -----------
  sq[i,] = as.numeric(q)
  sz[i,] = z
  smu[i,] = mu
  stau[i,] = tau
}
means = apply(smu,2,mean)
vars = apply(stau,2,mean)
qs = apply(sq,2,mean)

curve(dnorm(x,means[1],1/vars[1])*qs[1] + dnorm(x,means[2],1/vars[2])*qs[2],-10,10, add=TRUE, col=2)

plot(sq[,1], type='l');abline(h=Q[1], col=2)
plot(sq[,2], type='l');abline(h=Q[2], col=2)


