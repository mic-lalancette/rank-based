library("parallel")
library("pbapply")



#####################################
######## SPATIAL SIMULATIONS ########
#####################################

# This file generates spatial data sets according to the Brown-Resnick model presented in Section 5.2 of the paper (with the additive Pareto noise described therein)
# and calcualtes M-estimators for each simulation repetition, for a range of different parameters m
# The results and tuning parameters are then saved to a RData file to be called later to produce the plots in Section 5.2



# Simulation parameters
n = 5000 # sample size
nrep = 4 # number of data sets generated
ncl = 2 #number of clusters used for parallel computing
sv <- function(dist, alpha, scale) 0.5 * (dist/scale)^alpha
theta_f <- function(dist, a, s) pnorm(0.5 * (dist/s)^(a/2))
true_a = 1
true_s = 3
# Noise type needs to be "a" (additive) or "m" (mixture)
noise_type <- "a"
p_noise <- 0.1
tail_noise <- 4

# Estimation parameters
m_values = seq(25, 500, by=25)
M = length(m_values)
R1 = c(0, 0, 1, 1)
R2 = c(0, 0, 2, 2)
R3 = c(0.5, 0.5, 1.5, 1.5)
R4 = c(0, 0, 1, 3)
R5 = c(0, 0, 3, 1)
# each line of the following matrix corresponds to the coordinates of a rectangle (given as (x_min, y_min, x_max, y_max))
# over which the function c will be integrated
Rec = rbind(R1, R2, R3, R4, R5)
g_dim = nrow(Rec)
rp = 0.6

cl <- makeCluster(ncl)
invisible(clusterEvalQ(cl, {
  
  library("mev")
  library("pbapply")
  load("SpInfo.Rdata")
  source("Functions.R")
  
}))
clusterExport(cl = cl, varlist = c("n", "sv", "theta_f", "true_a", "true_s", "noise_type", "p_noise", "tail_noise",
                                   "m_values", "M", "Rec", "g_dim", "rp"))

estimators <- pbreplicate(nrep, {
  
  Y = -1 / log( 1 - exp(-1 / rmev(n = n, vario = sv, coord = gr, model = "br", alpha = true_a, scale = true_s)) )
  Z = (-log(matrix(runif(d*n), ncol = d)))^(-1/tail_noise)
  if(noise_type == "m") {
    U = rbinom(n, 1, p_noise)
    sim_data <- sapply(1:d, function(i) (1-U)*Y[,i] + U*Z[,i])
  }
  else sim_data <- Y + Z
  
  # Preparation
  params = seq(0.51, 0.99, by = 0.01)
  npar = length(params)
  F_rectangle <- function(t) F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
    F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2])
  Frp = F_rectangle(rp)
  F_rectangle <- function(t) ( F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
                                 F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2]) ) / Frp
  true_Int = t(sapply(params, F_rectangle))
  
  IntEst <- vapply(1:np, function(i) {
    
    X = sim_data[, pairs[i,]]
    r1 = rank(-X[,1])
    r2 = rank(-X[,2])
    R = r2[sort(r1, index.return = TRUE)$ix]
    k_values = sort(apply(cbind(r1, r2), 1, max))[m_values]
    
    A2 <- sapply(1:M, function(j) {
      
      k = k_values[j]
      Int = ( int_ind(R, k, a = Rec[,3], b = Rec[,4]) + int_ind(R, k, a = Rec[,1], b = Rec[,2]) -
                int_ind(R, k, a = Rec[,1], b = Rec[,4]) - int_ind(R, k, a = Rec[,3], b = Rec[,2]) ) / Frp
      emp_Int = t(array(Int, dim = c(length(Int), npar)))
      A1 = c(Int, M_E(params, emp_Int, true_Int)$arg[1])
      return(A1)
      
    })
    
    return(A2)
    
  }, FUN.VALUE = array(0, dim = c(g_dim+1, M)))
  IntEst <- aperm(IntEst, 3:1)
  
  estimators_sp_1 <- sapply(1:M, function(j) {
    
    Q1 <- function(as) sum( (IntEst[,j,g_dim+1] - theta_f(distance, a=as[1], s=as[2]))^2 )
    return(optim(par = c(1, 1), fn = Q1, method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(1.99, 10))$par)
    
  })
  
  estimators_sp_2 <- sapply(1:M, function(j) {
    
    Q2 <- function(as) {
      
      F_array = t(sapply(theta_f(distance, a=as[1], s=as[2]), F_rectangle))
      sigma = rowSums(IntEst[,j,1:g_dim]) / rowSums(F_array)
      return( sum((IntEst[,j,1:g_dim] - sigma*F_array)^2) )
      
    }
    return(optim(par = c(1, 1), fn = Q2, method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(1.99, 10))$par)
    
  })
  
  # M columns represent the m_values. Rows 1:np are pairwise estimatores,
  # rows np+(1:2) are alpha and beta "least squares" estimates
  # and rows np+(3+4) are alpha and beta "global" estimates
  e <- rbind(IntEst[,,g_dim+1], estimators_sp_1, estimators_sp_2)
  return(e)
  
}, cl=cl)
stopCluster(cl)

all_estimators_sp <- estimators

save(n, nrep, theta_f, true_a, true_s, p_noise, tail_noise, m_values, M, Rec, all_estimators_sp,
     file = "SpEstimators.Rdata")



