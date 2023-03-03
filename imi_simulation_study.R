library(Amelia)
library(mice)
library(imi)
library(MASS)
library(openxlsx)
library(stringr)  


# Settings
K = 100            # sample size
k = 5              # length of random vector
M = 500            # number of imputations
reps = 100         # scenario repetitions
mis_prop = 0.1     # c(0.1, 0.7) Proportion of missing Values

# CS spezifische Settings
sigma_sq = 0.25   
rho = 0.1          # c(0.1, 0.9)
mu = 0 * c(1:k)
get_tau = function(rho, sigma_sq){rho * sigma_sq / (1-rho)}

# imi specific Settings
eps = 0.05        # c(0.005, 0.05)
k0 = 1            # c(1, 3, 5)
M0 = 2

# Model selection
CS = TRUE         # FALSE for Logreg



################################################################################
##### Covariance Matrices ######################################################
################################################################################

# Covariance Matrix Sigma_cs for Compound Symmetry
Sigma_cs = function(rho, sigma_sq, k, tau){
  return(sigma_sq * diag(k) + tau * matrix(data=1, nrow=k, ncol=k))
}





################################################################################
##### Data sets to construct the simulation data ###############################
################################################################################

# Vector with random numbers distributed as N(mu, Sigma_)
get_random_y = function(mu, sigma, K){
  return(mvrnorm(K, mu=mu, sigma))
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

# Creates covariance matrix for the within variability
get_W_hat = function(sigma_sq, tau, k, K){
  W_hat = matrix(data=NA, nrow=3, ncol=3)
  W_hat[1,] = c((sigma_sq + k * tau) / (K*k), 0, 0)
  prefactor = (2*sigma_sq**2) / (k*K*(k-1))
  W_hat[2,] = c(0, prefactor * k, prefactor * (-1))
  W_hat[3,] = c(0, prefactor * (-1), prefactor * (sigma_sq**2 + 2*(k-1)*tau*sigma_sq + k*(k-1)*tau**2 ) / sigma_sq**2)
  
  return(W_hat)
}

# Creates covariance matrix for the between variability for each iteration
get_B_hat = function(theta_hat, m){
  theta_bar = colMeans(theta_hat[1:m,])
  B_hat_m = matrix(data=0, nrow=3, ncol=3)
  for (j in 1:m){
    B_hat_m = B_hat_m + (theta_hat[j,] - theta_bar) %*% t(theta_hat[j,] - theta_bar)
  }
  B_hat = 1/(m-1) * B_hat_m  
  return(B_hat)
}

# Creates a matrix with the parameter estimates for each iteration
get_parameters_for_all_imputations = function(k, K, mu, rho, sigma_sq, M, mis_prop, tau){
  Y = get_random_y(mu, Sigma_cs(rho, sigma_sq, k, tau), K)

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

# Creates the M imputed data sets + the parameter estimates for each iteration
# + covariance matrix for the within variability for each iteration
get_logit_reg_data = function(rho, sigma_sq, k, tau, K){
  k = 2
  mu = 0 * c(1:k)
  X_logit = get_random_y(mu, Sigma_cs(rho, sigma_sq, k, tau), K)
  beta = c(0.2, -2, 0.5)
  Y_logit = 1 / (1+ exp(-(beta[1] + X_logit %*% beta[2:3])))
  
  
  X_logit_missing = ampute(X_logit, prop=0.1)$amp
  
  imputed_datasets = amelia(X_logit_missing, m=M)$imputations
  
  theta_hat_list = matrix(data=NA, nrow=M, ncol=3)
  W_hat_list = list()
  i = 1
  for (dataset in imputed_datasets){
    logreg_obj = glm(Y_logit ~ dataset[,1] + dataset[,2], family = binomial(link = "logit"))
    W_hat_list = append(W_hat_list, list(vcov(logreg_obj)))
    theta_hat_list[i, 1:3] = logreg_obj$coefficients
    i = i + 1
  }
  return(list('theta_hat_list' = theta_hat_list, 'W_hat_list' = W_hat_list))
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
  return( sqrt(t(theta_m1 - theta_m) %*% solve(S) %*% (theta_m1 - theta_m)))
}





################################################################################
##### Iterative Multiple Imputation (imi) ######################################
################################################################################

# Implementation of the imi algorithm (Nassiri et al., 2020)
imi = function(theta_hat, W_hat, M0, eps, k0){
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
    V_m1 = B_hat_m1 * (1+1/m) + W_hat[m][[1]]
    d[(m+1)] = d_mah(theta_m1, theta_m, V_m1)
    
    if (m+2-k0 > 0 && sum(d[(m+2-k0) : (m+1)] < eps) == k0) {
      return(c((m+2-k0), d[(m+2-k0)], colMeans(theta_hat[1:(m+2-k0),])))
    }
  }
}

# Modification to return a matrix with the distances for each iteration.
# Used for the convergence plot.
imi_for_distances = function(theta_hat, W_hat, M0, eps, k0){
  theta_m = colMeans(theta_hat[1:M0,])
  d = c(1:M)*0+eps
  for (m in (M0):(M-1)){
    theta_m1 = (m * theta_m + theta_hat[(m+1),]) / (m+1)
    B_hat_m1 = get_B_hat(theta_hat, m+1)
    V_m1 = B_hat_m1 * (1+1/m) + W_hat[m][[1]]
    d[(m+1)] = d_mah(theta_m1, theta_m, V_m1)
  }
  return(d)
}





################################################################################
##### Additional Help Functions ################################################
################################################################################

# Converts a float to a string rounded to 'dec' decimal positions
fts = function(input, dec){
  gsub(",", ".", gsub("\\.", "", sprintf(paste('%0.', dec, 'f', sep=''), input)))
}

# Creates the name of the file based on the settings of the simulation
get_path = function(mis_prop, sigma_sq, rho, eps, k0, M0, CS){
  if (CS){
    model = 'cs'
  } else {
    model = 'logreg'
  }
  return(paste(model, '-sim-missing_', 
               fts(mis_prop,1), 
               '-sigma_', fts(sigma_sq, 2), 
               '-rho_', fts(rho, 1),  
               '-eps_', fts(eps,3), 
               '-k0_', k0, 
               '-M0_', M0, 
               '.xlsx', sep=''))
}










################################################################################
##### Main #####################################################################
################################################################################

output_df = data.frame(matrix(data=NA, nrow=reps, ncol=5))
colnames(output_df) = c('#imputations', 'Final distance', 'Mu', 'Sigma', 'Tau')

if (CS) {
  for (r in 1:reps){
    tau = get_tau(rho, sigma_sq)
    theta_hat = get_parameters_for_all_imputations(k, K, mu, rho, sigma_sq, M, mis_prop, tau)
    W_hat = get_W_hat(sigma_sq, tau, k, K)
    W_hat_list = list(W_hat)
    for  (i in 1:(M-1)){
      W_hat_list = append(W_hat_list, list(W_hat))
    }
    output = imi(theta_hat, W_hat_list, M0, eps, k0)
    output_df[r,] = output
  }
} else {
  for (r in 1:reps){
    tau = get_tau(rho, sigma_sq)
    logreg_data = get_logit_reg_data(rho, sigma_sq, k, tau, K)
    theta_hat = logreg_data[1]$theta_hat_list
    W_hat_list = logreg_data[2]$W_hat_list
    output = imi(theta_hat, W_hat_list, M0, eps, k0)
    output_df[r,] = output
  }
}

path = get_path(mis_prop, sigma_sq, rho, eps, k0, M0, CS)
save_path = paste('/Users/gedeonvogt/Desktop/', path, sep='')
write.xlsx(output_df, save_path, asTable=TRUE)
