library(Amelia)
library(mice)
library(imi)
library(MASS)
library(matlib)

# Settings
K = 100            # sample size
k = 5              # length of random vetor
M = 500            # number of imputations
sigma_sq = 1       # c(1, 16, 64)
rho = 0.2          # c(0.2, 0.5, 0.8)
mu = 0 * c(1:k)
eps = 0.005        # c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05)
k0 = 1             # c(1, 3, 5)
reps = 100         # scenario repetitions
mis_prop = 0.1     # c(0.1, 0.7) Proportion of missing Values
M0 = 2
tau = rho * sigma_sq / (1-rho)





################################################################################
##### Covariance Matrices ######################################################
################################################################################

# Covariance Matrix Sigma_cs for Compound Symmetry
Sigma_cs = function(rho, sigma_sq, k, tau){
  return(sigma_sq * diag(k) + tau * matrix(data=1, nrow=k, ncol=k))
}


# Covariance Matrix Sigma_ar for AR(1)
Sigma_ar = function(sigma_sq, k, rho){
  C = matrix(data=NA, nrow=k, ncol=k)
  for (i in 1:k){
    for (j in 1:k){
      C[i,j] = rho ** abs(i-j)
    }
  }
  return(sigma_sq * C)
}





################################################################################
##### Data sets to construct the simulation data ###############################
################################################################################

# Vector with random numbers distributed as N(mu, Sigma_)
get_random_y = function(mu, sigma, k){
  return(mvrnorm(k, mu=mu, sigma))
}


# The Dataset for the Multivariate Normal case with sample size K
get_data_set = function(k, K, mu, sigma){
  data_set = matrix(data=NA, nrow=K, ncol=k)
  for (i in 1:K){
    data_set[i, 1:k] = get_random_y(mu, sigma, k)
  }
  return(data_set)
}





################################################################################
##### Simulation Data: Compound Symmetry #######################################
################################################################################

# Parameter Estimation for CS based on Hermans (2019a)
get_theta_hat_cs_model = function(Y){
  N = dim(Y)[1]
  n = dim(Y)[2]
  J = matrix(data=1, nrow=n, nco=n)

  mu_hat = mean(as.matrix(Y))
  
  Z = c(1:N)*0
  for (i in 1:N){
    Z[i] = t(as.numeric(Y[i,])-mu) %*% (as.numeric(Y[i,]) - mu)
  }
  Z = sum(Z)
  
  Z_J = c(1:N)*0
  for (i in 1:N){
    Z_J[i] = t(as.numeric(Y[i,])-mu) %*% J %*% (as.numeric(Y[i,]) - mu)
  }
  Z_J = sum(Z_J)
  
  sigma_sq_hat = 1/(N*n*(n-1)) * (n*Z - Z_J)
  
  d_hat = 1/(N*n*(n-1)) * (Z_J - Z)
  
  return(c(mu_hat, sigma_sq_hat, d_hat))
}


get_W_hat = function(sigma_sq, tau, k, K){
  W_hat = matrix(data=NA, nrow=3, ncol=3)
  W_hat[1,] = c((sigma_sq + k * tau) / (K*k), 0, 0)
  prefactor = (2*sigma_sq**2) / (k*K*(k-1))
  W_hat[2,] = c(0, prefactor * k, prefactor * (-1))
  W_hat[3,] = c(0, prefactor * (-1), prefactor * (sigma_sq**2 + 2*(k-1)*tau*sigma_sq + k*(k-1)*tau**2 ) / sigma_sq**2)
  
  return(W_hat)
}

get_B_hat = function(theta_hat, m){
  theta_bar = colMeans(theta_hat[1:m,])
  B_hat_m = matrix(data=0, nrow=3, ncol=3)
  for (j in 1:m){
    B_hat_m = B_hat_m + (theta_hat[j,] - theta_bar) %*% t(theta_hat[j,] - theta_bar)
  }
  B_hat = 1/(m-1) * B_hat_m  
  return(B_hat)
}

get_parameters_for_all_imputations = function(k, K, mu, rho, sigma_sq, M, mis_prop){
  Y = get_random_y(mu, Sigma_cs(rho, sigma_sq, k, tau), 100)

  Y_mis = ampute(Y, prop=mis_prop)$amp
  imputed_datasets = amelia(Y_mis, m=M)$imputations
  
  theta_hat = matrix(data=NA, nrow=M, ncol=3)
  i = 1
  for (dataset in imputed_datasets){
    theta_hat[i, 1:3] = get_theta_hat_cs_model(dataset)
    i = i + 1
  }

  return(theta_hat)
}




################################################################################
##### Simulation Data: Logistic Regression #####################################
################################################################################
get_logit_reg_data = function(){
  X_logit = get_data_set(2,K, mu, Sigma_cs)
  beta = c(0.2, -2, 0.5)
  Y_logit = 1 / (1+ exp(-(beta[1] + X_logit %*% beta[2:3])))
  
  glm(Y_logit ~ X_logit[,1] + X_logit[,2], family = binomial(link = "logit"))$coefficients
  
  
  X_logit_missing = ampute(X_logit, prop=0.1)$amp
  
  imputed_datasets = amelia(X_logit_missing, m=M)$imputations
  
  theta = matrix(data=NA, nrow=M, ncol=3)
  i = 1
  for (dataset in imputed_datasets){
    theta[i, 1:3] = glm(Y_logit ~ dataset[,1] + dataset[,2], family = binomial(link = "logit"))$coefficients
    i = i + 1
  }
}




################################################################################
##### Distances ################################################################
################################################################################

# Euclidean Distance
d_euc = function(theta_m1, theta_m, S=NA){
  return( sqrt(t(theta_m1 - theta_m) %*% (theta_m1 - theta_m)) )
}


# Maximums Distance
d_max = function(theta_m1, theta_m, S=NA){
  return(max(abs(theta_m1 - theta_m)))
}


# Mahalanobis Distance
d_mah = function(theta_m1, theta_m, S){
  return( sqrt(abs(t(theta_m1 - theta_m) %*% S**(-1) %*% (theta_m1 - theta_m))))
}





################################################################################
##### Iterative Multiple Imputation (imi) ######################################
################################################################################

imi_cs = function(theta_hat, W_hat, M0, eps, k0){
  # theta_hat: #imputations x #parameters
  # M0: initial number of imputed data sets
  # distance: function of a norm with parameters theta_m1, theta_m, S (optional)
  # eps: accuracy to be achieved
  # k0: number of successive steps that the stopping rule should be validated
  # S: #imputations x #parameters x #parameters
  
  theta_m = colMeans(theta_hat[1:M0,])
  d = c(1:1000)*0+eps
  for (m in (M0):1000){
    theta_m1 = (m * theta_m + theta_hat[(m+1),]) / (m+1)
    B_hat_m1 = get_B_hat(theta_hat, m+1)
    V_m1 = B_hat_m1 * (1+1/m) + W_hat
    d[(m+1)] = d_mah(theta_m1, theta_m, V_m1)
    
    if (m+2-k0 > 0 && sum(d[(m+2-k0) : (m+1)] < eps) == k0) {
      return(c((m+2-k0), colMeans(theta_hat[1:(m+2-k0),])))
    }
  }
}






theta_hat = get_parameters_for_all_imputations(k, K, mu, rho, sigma_sq, M, mis_prop)
W_hat = get_W_hat(sigma_sq, tau, k, K)
W_hat_list = list(W_hat)
for (i in 1:(499)){
  W_hat_list = append(W_hat_list, list(W_hat))
}
V_hat = combine.mi(theta_hat, W_hat_list)
output = imi_cs(theta_hat, W_hat, M0, eps, k0)




