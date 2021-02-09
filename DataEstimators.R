load("SpInfo.Rdata")
load("RainData.Rdata")
source("Functions.R")



###############################
######## DATA ANALYSIS ########
###############################

# This file calcualtes M-estimators based on the rainfall data set as described in Section 6, for a range of different parameters m
# The results and tuning parameters are then saved to a RData file to be called later to produce the plots in Section 6



# Estimation parameters
m_values = seq(100, 2000, by=100)
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

theta_f <- function(dist, a, s) pnorm(0.5 * (dist/s)^(a/2))

# Preparation
params = seq(0.501, 0.999, by = 0.001)
npar = length(params)
F_rectangle <- function(t) F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
  F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2])
Frp = F_rectangle(rp)
F_rectangle <- function(t) ( F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,4]) + F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,2]) -
                               F_IMS_ind(c(t, t), a = Rec[,1], b = Rec[,4]) - F_IMS_ind(c(t, t), a = Rec[,3], b = Rec[,2]) ) / Frp
true_Int = t(sapply(params, F_rectangle))

IntEst <- vapply(1:np, function(i) {
  
  X = data[, pairs[i,]]
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

estimators_1 <- sapply(1:M, function(j) {
  
  Q1 <- function(as) sum( (IntEst[,j,g_dim+1] - theta_f(distance, a=as[1], s=as[2]))^2 )
  return(optim(par = c(1, 1), fn = Q1)$par)
  
})

estimators_2 <- sapply(1:M, function(j) {
  
  Q2 <- function(as) {
    
    F_array = t(sapply(theta_f(distance, a=as[1], s=as[2]), F_rectangle))
    sigma = rowSums(IntEst[,j,1:g_dim]) / rowSums(F_array)
    return( sum((IntEst[,j,1:g_dim] - sigma*F_array)^2) )
    
  }
  return(optim(par = c(1, 1), fn = Q2)$par)
  
})

# M columns represent the m_values. Rows 1:np are pairwise estimatores,
# rows np+(1:2) are alpha and beta "least squares" estimates
# and rows np+(3+4) are alpha and beta "global" estimates
estimators <- rbind(IntEst[,,g_dim+1], estimators_1, estimators_2)

# #Test plot (supposed to be same scatterplot as in Fig9 of the paper)
# plot(distance, estimators[1:np,4])

all_estimators_data <- estimators

save(theta_f, m_values, M, Rec, all_estimators_data, file = "DataEstimators.Rdata")



