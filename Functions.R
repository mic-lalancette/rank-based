#########################
# PACKAGE INSTALLATIONS #
#########################

# install.packages(c("mev", "parallel", "pbapply", "colorspace", "plotrix", "oz", "abind"))



######################################
# CALCULATION OF EMPIRICAL INTEGRALS #
######################################

# Calculation of the integral of m * the nonparametric estimator of c on [0, a] * [0, b], based on a (bivariate) sample X
# Assumes that input is already ranked and sorted, should be preceded by the following commands:
# r1 = n + 1 - rank(X[,1]) #Ranks of U's
# r2 = rank(X[,2]) #Ranks of V's
# R = r2[sort(r1, index.return = TRUE)$ix]

int_ind <- function(R, k, a, b){
  
  L = length(a)
  I = rep(-10, L)
  
  for(l in 1:L){
    
    if(a[l] == 0 || b[l] == 0) I[l] = 0
    
    else{
      
      n_terms = floor(a[l]*k)
      v1 = (n_terms - 1):1
      v2 = pmax(floor(b[l]*k) - R[1:(n_terms-1)], 0)
      I[l] = sum(v1*v2) / (k^2)
      
    }
    
  }
  
  return(I)
  
}



################
# M-ESTIMATION #
################

# Those two functions return the integral of c_theta(x, y) on [0, a] x [0, b],
# where c_theta is given by the inverted max-stable model and the Pareto random scale model, respectively

F_IMS_ind <- function(theta, a, b) (a^(theta[1]+1) * b^(theta[2]+1)) / ((theta[1]+1) * (theta[2]+1))

F_P_ind <- function(lambda, a, b){
  
  L = length(a)
  I = rep(-10, L)
  for(l in 1:L){
    
    AB = max(a[l], b[l])
    a[l] = min(a[l], b[l])
    b[l] = AB
    
    if(a[l] == 0 || b[l] == 0) I[l] = 0
    
    else{
      
      if(lambda > 0 && lambda < 1 && lambda != 0.5) I[l] = (
        
        ( (2-lambda) * (a[l]^2*b[l]/2 - a[l]^3/6) 
          - lambda * ( (1 - 2/lambda)*a[l]^3 / (3*(1/lambda + 1)*(2 - 1/lambda)) + a[l]^(1/lambda + 1)*b[l]^(2 - 1/lambda) / ((1/lambda + 1)*(2 - 1/lambda)) )
        ) / (2*(1 - lambda))
        
      )
      
      else if(lambda == 0.5) I[l] = 3*a[l]^2*b[l]/4 - 13*a[l]^3/36 - a[l]^3*log(b[l]/a[l])/6
      
      else if(lambda == 1) I[l] = a[l]^2 * ( a[l] + 9*b[l] + 6*b[l]*log(b[l]/a[l]) ) / 24
      
      else if(lambda > 1 && lambda < 2) I[l] = (
        
        ( lambda * ( (lambda - 2)*a[l]^(lambda + 2) / (2*lambda*(lambda + 2)) + a[l]^2*b[l]^lambda / (2*lambda) ) 
          - (2 - lambda) * ( -lambda*a[l]^(lambda + 2) / ((lambda + 1)*(lambda + 2)) + a[l]^(lambda + 1)*b[l] / (lambda + 1) )
        ) / (2*(lambda - 1))
        
      )
      
      else if(lambda == 2) I[l] = a[l]^2*b[l]^2/4
      
      else stop("Incorrect value of lambda")
      
    }
    
  }
  
  return(I)
  
}

# Finding the row where min_{sigma} ||emp_Int - sigma*true_Int|| is minimized

M_E <- function(params, emp_Int, true_Int){
  
  params = cbind(params)
  npar = nrow(params)
  sigma = rowSums(emp_Int) / rowSums(true_Int)
  sigma_array = array(sigma, dim = c(npar, ncol(emp_Int)))
  
  val = rowSums((sigma_array * true_Int - emp_Int)^2)
  wimi = which.min(val)
  
  return(list(arg = c(params[wimi,], sigma[wimi]), value = val[wimi]))
  
}



#################
#Alog SIMULATION#
#################

# Generating a grid of valid parameters for the asymmetric logistic model of Tawn,
# With the associated parameters theta_1, theta_2 of our representation
# alpha is the shape parameter (r in the original paper by Tawn)

alog_par_upper <- function(alpha, skip, gridsize) {
  
  theta = rep(0, 4)
  
  for(i in 2:gridsize){
    for(j in 1:(i-1)){
      
      a = i/(gridsize+1)
      b = j/(gridsize+1)
      
      test = c( 1 - a + a^alpha * (a^alpha + b^alpha)^(1/alpha - 1), 1 - b + b^alpha * (a^alpha + b^alpha)^(1/alpha - 1) )
      if(test[1] <= 0.95 && test[2] <= 0.95) theta = rbind(theta, c(a, b, test), c(b, a, test[2:1]))
      # Keep only parameters where theta_1 and theta_2 are under 0.95
      
    }
    
    a = i/(gridsize+1)
    b = a
    
    test = c( 1 - a + a^alpha * (a^alpha + b^alpha)^(1/alpha - 1), 1 - b + b^alpha * (a^alpha + b^alpha)^(1/alpha - 1) )
    if(test[1] <= 0.95 && test[2] <= 0.95) theta = rbind(theta, c(a, b, test), c(b, a, test[2:1]))
    
  }
  
  # Uncomment the next line for a plot of the generated parameters (in our parametrization)
  # plot(theta[3, seq(2, ncol(theta), by = skip)], theta[4, seq(2, ncol(theta), by = skip)], xlim = c(0, 1), ylim = c(0, 1))
  
  return(theta[seq(2, nrow(theta), by = skip),])
} #Columns 1 and 2 of this matrix are theta and phi (from Tawn's representation, used in rmev) and lines 3 and 4 are our theta (to be estimated)



