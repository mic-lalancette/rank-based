library("mev")
library("parallel")
library("pbapply")
library("colorspace")
library("plotrix")

source("Functions.R")



##########################################################
######## COMPARISON OF DIFFERENT WEIGHT FUNCTIONS ########
##########################################################

# This file generates data sets according to the three models presented in Section 5.1 of the paper (with the additive Pareto noise described therein)
# and calcualtes M-estimators for each model and each simulation repetition, for a range of different parameters k and weight functions g
# The results and tuning parameters are then saved to a RData file to be called later to produce the plots in Section S6.1.1



n = 5000 # sample size
nrep = 4 # number of samples generated
k_values = seq(100, n/3, by=100)
K = length(k_values)
R1 = c(0, 0, 1, 1)
R2 = c(0, 0, 2, 2)
R3 = c(0.5, 0.5, 1.5, 1.5)
R4 = c(0, 0, 1, 3)
R5 = c(0, 0, 3, 1)
Rec_list = list(rbind(R1, R2, R3, R4, R5), rbind(R1, R2), rbind(R1, R3), rbind(R1, R4, R5),
                rbind(R1, R2, R3), rbind(R1, R2, R4, R5), rbind(R1, R3, R4, R5))
W = length(Rec_list)
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
F_rectangle_list <- lapply(Rec_list, function(Rec) {
  function(t) F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
    F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2])
})
true_Int_list = lapply(F_rectangle_list, function(F_) t(sapply(params, F_) / F_(rp)))

estimators = array(-10, dim = c(L, K, W, nrep))

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
      
      for(w in 1:W){
        
        Rec = Rec_list[[w]]
        Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                  int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2])
        ) / F_rectangle_list[[w]](rp)
        emp_Int = t(array(Int, dim = c(length(Int), npar)))
        estimators[i, j, w, r] = M_E(params, emp_Int, true_Int_list[[w]])$arg[1]
        
      }
      
    }
    
  }
  print(paste(i, "/", L))
  
}

all_estimators_IHR_wf = estimators



######################## INVERTED ASYMMETRIC LOGISTIC MODEL #####################


alph <- 2
theta_values_IAlog <- alog_par_upper(alph, 1, 15)
theta_values <- theta_values_IAlog
L = nrow(theta_values)
rp = c(0.6, 0.6) # a reference point. Integrals are divided by the integrals at that parameter value

ncl <- 2 # Number of clusters used for parallel computation

cl <-makeCluster(ncl)
invisible(clusterEvalQ(cl, {
  
  library("mev")
  source("Functions.R")
  
}))
clusterExport(cl = cl, varlist = c("n", "k_values", "K", "Rec_list", "W", "noise_type", "p_noise",
                                   "tail_noise", "alph", "theta_values", "L", "rp"))

estimators <- pbreplicate(nrep, {
  
  x = rep(1:100, 1:100)
  y = 100
  for(i in 99:1) y = c(y, i:100)
  params = cbind(x, y)/100
  npar = nrow(params)
  F_rectangle_list <- lapply(Rec_list, function(Rec) {
    function(t) F_IMS_ind(t, a = Rec[,3], b = Rec[,4]) + F_IMS_ind(t, a = Rec[,1], b = Rec[,2]) -
      F_IMS_ind(t, a = Rec[,1], b = Rec[,4]) - F_IMS_ind(t, a = Rec[,3], b = Rec[,2])
  })
  true_Int_list = lapply(F_rectangle_list, function(F_) t(apply(params, 1, F_) / F_(rp)))
  
  e = array(-10, dim = c(L, K, W, 2))
  
  for(i in 1:L){
    nu = theta_values[i,1]
    phi = theta_values[i,2]
    
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
      
      for(w in 1:W){
        
        Rec = Rec_list[[w]]
        Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                  int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2]) ) / F_rectangle_list[[w]](rp)
        emp_Int = t(array(Int, dim = c(length(Int), npar)))
        e[i, j, w,] = M_E(params, emp_Int, true_Int_list[[w]])$arg[1:2]  
        
      }
      
    }
    
  }
  
  return(e)
  
}, cl=cl)
stopCluster(cl)

estimators <- aperm(estimators, c(1, 2, 3, 5, 4))

all_estimators_IAlog_wf = estimators



######################## PARETO RANDOM SCALE MODEL #####################


theta_values_P <- 1:30/10
theta_values <- theta_values_P # parameters from which to simlate
L = length(theta_values)
rp = 1

params = 2*(1:200)/200
npar = length(params)
F_rectangle_list <- lapply(Rec_list, function(Rec) {
  function(t) F_P_ind(t, a = Rec[,3], b = Rec[,4]) + F_P_ind(t, a = Rec[,1], b = Rec[,2]) -
    F_P_ind(t, a = Rec[,1], b = Rec[,4]) - F_P_ind(t, a = Rec[,3], b = Rec[,2])
})
true_Int_list = lapply(F_rectangle_list, function(F_) t(sapply(params, F_) / F_(rp)))

estimators = array(-10, dim = c(L, K, W, nrep))

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
      
      for(w in 1:W){
        
        Rec = Rec_list[[w]]
        Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                  int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2])
        ) / F_rectangle_list[[w]](rp)
        emp_Int = t(array(Int, dim = c(length(Int), npar)))
        estimators[i, j, w, r] = M_E(params, emp_Int, true_Int_list[[w]])$arg[1]
        
      }
      
    }
    
  }
  print(paste(i, "/", L))
  
}

all_estimators_P_wf = estimators



save(n, nrep, k_values, K, Rec_list, W, theta_values_IHR, alph, theta_values_IAlog, theta_values_P,
     all_estimators_IHR_wf, all_estimators_IAlog_wf, all_estimators_P_wf,
     file = "WfEstimators.RData")



