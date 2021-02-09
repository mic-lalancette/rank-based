library("mev")
library("parallel")
library("pbapply")
library("colorspace")
library("plotrix")

source("Functions.R")



#######################################
######## BIVARIATE SIMULATIONS ########
#######################################

# This file generates data sets according to the three models presented in Section 5.1 of the paper (with the additive Pareto noise described therein)
# and calcualtes M-estimators for each model and each simulation repetition, for a range of different parameters k
# The results and tuning parameters are then saved to a RData file to be called later to produce the plots in Section 5.1



# Some simulation parameters

n = 5000 # sample size
nrep = 4 # number of samples generated
k_values = seq(100, n/3, by=100) # set of values k to be used for estimation
K = length(k_values)
R1 = c(0, 0, 1, 1)
R2 = c(0, 0, 2, 2)
R3 = c(0.5, 0.5, 1.5, 1.5)
R4 = c(0, 0, 1, 3)
R5 = c(0, 0, 3, 1)
# each line of the following matrix corresponds to the coordinates of a rectangle (given as (x_min, y_min, x_max, y_max))
# over which the function c will be integrated
Rec = rbind(R1, R2, R3, R4, R5)
g_dim = nrow(Rec)
# Noise type needs to be "a" (additive) or "m" (mixture)
noise_type <- "a"
p_noise <- 0.1
tail_noise <- 4


######################## INVERTED HUSLER-REISS MODEL #####################


theta_values_IHR <- seq(0.525, 0.975, by = 0.025) # parameters from which to simlate
theta_values <- theta_values_IHR
L = length(theta_values)
Mat = cbind(c(0, 1), c(1, 0)) # needed for the Husler-Reiss distribution
rp = 0.6 # a reference point. Integrals are divided by the integrals at that parameter value

params = seq(0.51, 0.99, by = 0.01)
npar = length(params)
F_rectangle <- function(t) F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
  F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2])
true_Int = t(sapply(params, F_rectangle) / F_rectangle(rp))

dn = list(theta_values, paste0("k_", k_values), paste0("rep_", 1:nrep))
estimators = array(-10, dim = c(L, K, nrep), dimnames = dn)

for(i in 1:L){
  theta = theta_values[i]
  for(r in 1:nrep){
    
    Y = -1 / log( 1 - exp( -1 / rmev(n = n, d = 2, sigma = qnorm(theta)^2 * Mat, model = "hr")) )
    # Z = cbind(runif(n, 0, 0.5), runif(n, 0, 0.5))
    Z = (-log(cbind(runif(n), runif(n))))^(-1/tail_noise)
    if(noise_type == "m") {
      U = rbinom(n, 1, p_noise)
      X = cbind((1-U)*Y[,1] + U*Z[,1], (1-U)*Y[,2] + U*Z[,2])
    }
    else X = Y + Z
    r1 = rank(-X[,1])
    r2 = rank(-X[,2])
    R = r2[sort(r1, index.return = TRUE)$ix]
    
    for(j in 1:K){
      
      k = k_values[j]
      Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2]) ) / F_rectangle(rp)
      emp_Int = t(array(Int, dim = c(length(Int), npar)))
      estimators[i, j, r] = M_E(params, emp_Int, true_Int)$arg[1]
      
    }
    
  }
  print(paste(i, "/", L))
  
}

all_estimators_IHR_n5000_ind = estimators


######################## INVERTED ASYMMETRIC LOGISTIC MODEL #####################


alph <- 2
theta_values_IAlog <- alog_par_upper(alph, 1, 15)
theta_values <- theta_values_IAlog
L = nrow(theta_values)
rp = c(0.6, 0.6) # a reference point. Integrals are divided by the integrals at that parameter value

x = rep(1:100, 1:100)
y = 100
for(i in 99:1) y = c(y, i:100)
params = cbind(x, y)/100
npar = nrow(params)
F_rectangle <- function(t) F_IMS_ind(t, a = Rec[,3], b = Rec[,4]) + F_IMS_ind(t, a = Rec[,1], b = Rec[,2]) -
  F_IMS_ind(t, a = Rec[,1], b = Rec[,4]) - F_IMS_ind(t, a = Rec[,3], b = Rec[,2])
true_Int = t(apply(params, 1, F_rectangle) / F_rectangle(rp))

dn = list(paste0("theta_", 1:L), paste0("k_", k_values), paste0("rep_", 1:nrep), paste0("component_", 1:2))
estimators = array(-10, dim = c(L, K, nrep, 2))

for(i in 1:L){
  nu = theta_values[i,1]
  phi = theta_values[i,2]
  for(r in 1:nrep){
    
    Y = -1 / log( 1 - exp( -1 / rmev(n = n, d = 2, param = alph, asy = list(1-nu, 1-phi, c(nu, phi)), model = "alog") ) )
    # Z = cbind(runif(n, 0, 0.5), runif(n, 0, 0.5))
    Z = (-log(cbind(runif(n), runif(n))))^(-1/tail_noise)
    if(noise_type == "m") {
      U = rbinom(n, 1, p_noise)
      X = cbind((1-U)*Y[,1] + U*Z[,1], (1-U)*Y[,2] + U*Z[,2])
    }
    else X = Y + Z
    r1 = rank(-X[,1])
    r2 = rank(-X[,2])
    R = r2[sort(r1, index.return = TRUE)$ix]
    
    for(j in 1:K){
      
      k = k_values[j]
      Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2]) ) / F_rectangle(rp)
      emp_Int = t(array(Int, dim = c(length(Int), npar)))
      estimators[i, j, r,] = M_E(params, emp_Int, true_Int)$arg[1:2]
      
    }
    
  }
  print(paste(i, "/", L))
  
}

all_estimators_IAlog_n5000_ind = estimators


######################## PARETO RANDOM SCALE MODEL #####################


theta_values_P <- 1:30/10
theta_values <- theta_values_P # parameters from which to simlate
L = length(theta_values)
rp = 1

params = 2*(1:200)/200
npar = length(params)
F_rectangle <- function(t) F_P_ind(t, a = Rec[,3], b = Rec[,4]) + F_P_ind(t, a = Rec[,1], b = Rec[,2]) -
  F_P_ind(t, a = Rec[,1], b = Rec[,4]) - F_P_ind(t, a = Rec[,3], b = Rec[,2])
true_Int = t(apply(cbind(params), 1, F_rectangle) / F_rectangle(rp))

dn = list(theta_values, paste0("k_", k_values), paste0("rep_", 1:nrep))
estimators = array(-10, dim = c(L, K, nrep), dimnames = dn)

for(i in 1:L){
  theta = theta_values[i]
  for(r in 1:nrep){
    
    R = runif(n)^(-1/theta) #n iid Pareto(theta)
    w1 = runif(n)^(-1) #n iid Pareto(1)
    w2 = runif(n)^(-1) #n iid Pareto(1)
    Y = R * cbind(w1, w2)
    Z = cbind(runif(n), runif(n))^(-1/tail_noise)
    if(noise_type == "m") {
      U = rbinom(n, 1, p_noise)
      X = cbind((1-U)*Y[,1] + U*Z[,1], (1-U)*Y[,2] + U*Z[,2])
    }
    else X = Y + Z
    r1 = rank(-X[,1])
    r2 = rank(-X[,2])
    R = r2[sort(r1, index.return = TRUE)$ix]
    
    for(j in 1:K){
      
      k = k_values[j]
      Int = (int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
               int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2])) / F_rectangle(rp)
      emp_Int = t(array(Int, dim = c(length(Int), npar)))
      estimators[i, j, r] = M_E(params, emp_Int, true_Int)$arg[1]
      
    }
    
  }
  
  print(paste(i, "/", L))
}

all_estimators_P_n5000_ind = estimators


save(n, nrep, k_values, K, Rec, theta_values_IHR, alph, theta_values_IAlog, theta_values_P,
     all_estimators_IHR_n5000_ind, all_estimators_IAlog_n5000_ind, all_estimators_P_n5000_ind,
     file = "BivEstimators.RData")



