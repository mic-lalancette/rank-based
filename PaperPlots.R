library("mev")
library("colorspace")
library("plotrix")
library("oz")
library("pbapply")
library("abind")



###############################################
######## BIVARIATE SIMULATION PLOTS ###########
###############################################

# This file uses the estimators contained in various RData files to produce all the plots found in the paper



# source("BivEstimators.R") # Run once, to generate the data and produce the estimators
load("BivEstimators.RData")

########################
# SAMPLES SCATTERPLOTS #
########################

n_sample = 1000
Mat = cbind(c(0, 1), c(1, 0)) # needed for the Husler-Reiss distribution

for(theta in c(0.6, 0.75, 0.9)){
  
  setEPS()
  postscript(file = paste0("./Plots/sampleIHR", as.character(100*theta), ".eps"), width = 5, height = 5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(1 / rmev(n=n_sample, d=2, sigma = qnorm(theta)^2 * Mat, model = "hr"), xlab = "", ylab = "", main = "")
  dev.off()
  
}

theta_values = theta_values_IAlog
for(i in c(9, 128, 144)){
  
  nu = theta_values[i,1]
  phi = theta_values[i,2]
  setEPS()
  postscript(file = paste0("./Plots/sample_IAlog_", i, "_",
                    paste(as.character(1000*round(theta_values[i,3:4], 3)), collapse = "_"), ".eps"), width = 5, height = 5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(1 / rmev(n = n_sample, d = 2, param = alph, asy = list(1-nu, 1-phi, c(nu, phi)), model = "alog"), xlab="", ylab="", main="")
  dev.off()
  
}

for(lambda in c(0.4, 1, 1.6)) {
  
  R = runif(n_sample)^(-1/lambda) #n iid Pareto(alpha)
  w1 = runif(n_sample)^(-1) #n iid Pareto(1)
  w2 = runif(n_sample)^(-1) #n iid Pareto(1)
  X = R * cbind(w1, w2)
  U = apply(X, 2, rank) / (n_sample+1)
  setEPS()
  postscript(file = paste0("./Plots/sampleP", as.character(10*lambda), ".eps"), width = 5, height = 5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(-log(1-U), xlab = "", ylab = "", main = "")
  dev.off()
  
}

#####################
# CHOICE OF k PLOTS #
#####################

theta_values <- theta_values_IHR

for(i in c(4, 10, 16)){
  
  bias_M = rowMeans(all_estimators_IHR_n5000_ind[i,,]) - rep(theta_values[i], K)
  RMSE_M = sqrt(rowMeans( (all_estimators_IHR_n5000_ind[i,,] - array(theta_values[i], dim = c(K, nrep)))^2 ))
  setEPS()
  postscript(file = paste0("./Plots/ck_IHR_n5000_ind_", 1000*theta_values[i], ".eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(k_values, abs(bias_M), ylim = c(0, 0.12), type = "l", col = "blue", main = "", xlab = "k", ylab = "Bias and RMSE")
  lines(k_values, RMSE_M, lty = 2, col = "blue")
  dev.off()
  
}

theta_values <- theta_values_IAlog[,3:4]
L = nrow(theta_values)
theta_array = array(-10, dim = c(L, K, nrep, 2))
for(j in 1:K){
  for(r in 1:nrep) theta_array[, j, r,] = theta_values
}
f1 <- function(x) sqrt(sum((colMeans(x)^2)))
f2 <- function(x) sqrt(2*mean(x))
bias_M = apply(all_estimators_IAlog_n5000_ind - theta_array, c(1, 2), f1)
RMSE_M = apply((all_estimators_IAlog_n5000_ind - theta_array)^2, c(1, 2), f2)

for(i in c(9, 128, 144)){
  
  setEPS()
  postscript(file = paste0("./Plots/ck_IAlog_n5000_ind_", i, "_", paste(as.character(1000*round(theta_values[i,], 3)), collapse = "_"), ".eps"),
      width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(k_values, abs(bias_M[i,]), ylim = c(0, 0.2), type = "l", col = "blue", main = "", xlab = "k", ylab = "Bias and RMSE")
  lines(k_values, RMSE_M[i,], lty = 2, col = "blue")
  dev.off()
  
}

theta_values <- theta_values_P

for(i in c(4, 10, 16)){
  
  bias_M = rowMeans(all_estimators_P_n5000_ind[i,,]) - rep(theta_values[i], K)
  RMSE_M = sqrt(rowMeans( (all_estimators_P_n5000_ind[i,,] - array(theta_values[i], dim = c(K, nrep)))^2 ))
  setEPS()
  postscript(file = paste0("./Plots/ck_P_n5000_", 10*theta_values[i], ".eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(k_values, abs(bias_M), ylim = c(0, 0.3), type = "l", col = "blue", main = "", xlab = "k", ylab = "Bias and RMSE")
  lines(k_values, RMSE_M, lty = 2, col = "blue")
  dev.off()
  
}

############
# BOXPLOTS #
############

theta_values <- theta_values_IHR
L = length(theta_values)
j=8
setEPS()
postscript(file = paste0("./Plots/Boxplot_IHR_n5000_ind.eps"), width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
boxplot(t(all_estimators_IHR_n5000_ind[2*(1:9),j,]), col = "lightblue",
        at = theta_values[2*(1:9)], xlim = c(0.5, 1), boxwex = 0.015, xlab = expression(theta), ylab = "Estimators")
lines(theta_values, theta_values)
dev.off()

theta_values <- theta_values_P
L = length(theta_values)
ptp = 2*(1:10)
j=4
setEPS()
postscript(file = paste0("./Plots/Boxplot_P_n5000_ind.eps"), width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
boxplot(t(all_estimators_P_n5000_ind[ptp,j,]), col = "lightblue",
        at = theta_values[ptp], xlim = c(0, 2), boxwex = 0.06, xlab = expression(lambda), ylab = "Estimators")
lines(theta_values, theta_values)
dev.off()

###############
# COLOR PLOTS #
###############

theta_values <- theta_values_IAlog[,3:4]
estimators = all_estimators_IAlog_n5000_ind
L = dim(estimators)[1]
p = dim(estimators)[4]
theta_array = array(-10, dim = dim(estimators))
for(j in 1:K){
  for(r in 1:nrep) theta_array[, j, r,] = theta_values
}

f1 <- function(x) sqrt(sum((colMeans(x)^2)))
f2 <- function(x) sqrt(2*mean(x))

errors = array(-10, dim = c(L, K, p))
errors[,,1] = apply(estimators - theta_array, c(1, 2), f1)
errors[,,2] = apply((estimators - theta_array)^2, c(1, 2), f2)

mypal = colorRampPalette(colors = c("cornflowerblue", "palegreen", "orangered"), space = "Lab", interpolate = "linear")
ncolor = 100
j = 8
max_color = max(errors[,j,])

setEPS()
postscript(file = "./Plots/Color_IALog_n5000_ind.eps", width = 12, height = 5.5)
par(mfrow = c(1, 2), oma = c(0, 0, 0, 10), xpd = TRUE, pty = "s")
plot( theta_values, pch = 16, xlab = expression(theta * scriptstyle(1)), ylab = expression(theta * scriptstyle(2)), main = "",
      col = mypal(100)[pmax(1, ceiling( ncolor*errors[,j,1]/max_color ))] )
plot( theta_values, pch = 16, xlab = expression(theta * scriptstyle(1)), ylab = expression(theta * scriptstyle(2)), main = "",
      col = mypal(100)[pmax(1, ceiling( ncolor*errors[,j,2]/max_color ))] )
color.legend(0.99, 0.72, 1.02, 0.95,
             legend = round(seq(0, max_color, length.out = 5), digits = 2), rect.col = mypal(100),
             cex = 1, gradient = "y")
dev.off()



#########################################################
######### COMPARISON OF WEIGHT FUNCTIONS PLOT ###########
#########################################################


# source("WfEstimators.R") # Run once, to generate the data and produce the estimators
load("WfEstimators.RData")


# IHR model
theta_values <- theta_values_IHR
L = length(theta_values)

j = 8
RMSE = apply((all_estimators_IHR_wf - array(theta_values, dim = c(L, K, W, nrep)))^2, c(1, 2, 3), function(x) sqrt(mean(x)))
setEPS()
postscript(file = paste0("./Plots/wf_IHR.eps"), width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot( theta_values, RMSE[,j,1], type = "l",
      ylim = c(0, max(RMSE[,j,])),
      xlab = expression(theta), ylab = "RMSE",
      main = "" )
for(w in 2:7){
  lines(theta_values, RMSE[,j,w], col = w)
}
dev.off()


# IAlog model
theta_values <- theta_values_IAlog
L = nrow(theta_values)
p = 2 # dimension of the parameter
f2 <- function(x) sqrt(2*mean(x))

j = 8
RMSE = apply((all_estimators_IAlog_wf - aperm(array(theta_values[,3:4], dim = c(L, p, K, W, nrep)), c(1, 3, 4, 5, 2)))^2,
             c(1, 2, 3), f2)

crit = apply(theta_values[,3:4], 1, sum)
par_order = sort(crit, index.return = TRUE)$ix
setEPS()
postscript(file = paste0("./Plots/wf_IAlog.eps"), width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot( crit[par_order],
      RMSE[par_order, j, 1], type = "l",
      ylim = c(0, max(RMSE[,j,c(1, 4, 6, 7)])),
      xlab = expression(paste(theta[1] , "+" , theta[2])), ylab = "RMSE",
      main = "" )
for(w in c(4, 6, 7)){
  lines(crit[par_order], RMSE[par_order, j, w], col = w)
}
dev.off()


# P model
theta_values <- theta_values_P
L = length(theta_values)

j = 4
RMSE = apply((all_estimators_P_wf - array(theta_values, dim = c(L, K, W, nrep)))^2, c(1, 2, 3), function(x) sqrt(mean(x)))
setEPS()
postscript(file = paste0("./Plots/wf_P.eps"), width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot( theta_values[1:20], RMSE[1:20,j,1], type = "l",
      ylim = c(0, max(RMSE[1:20,j,c(1, 3:7)])),
      xlab = expression(lambda), ylab = "RMSE",
      main = "" )
for(w in 3:7){
  lines(theta_values[1:20], RMSE[1:20,j,w], col = w)
}
legend(x = 1.2, y = max(RMSE[1:20,j,c(1, 3:7)]),
       legend = c(expression(g^(1)), expression(g^(2)), expression(g^(3)), expression(g^(4)),
                  expression(g^(5)), expression(g^(6)), expression(g^(7))),
       col=1:7, lty=1, lwd=1, ncol=2)
dev.off()



#############################################
######## SPATIAL SIMULATION PLOTS ###########
#############################################


# source("SpEstimators.R") # Run once, to generate the data and produce the estimators
load("SpEstimators.Rdata")
load("SpInfo.Rdata")

for(i in 1:2){
  
  e = all_estimators_sp[np+2*(i-1)+(1:2),,]
  sup = pbapply(e, c(2, 3), function(as) max(abs(sapply((1:60)/20, theta_f, a=as[1], s=as[2]) -
                                                   sapply((1:60)/20, theta_f, a=true_a, s=true_s))))
  assign(paste0("e", i), abind(e, sup, along=1))
  
}

J = 3:M

for(i in 1:2){
  
  e = get(paste0("e", i))
  
  # Estimation error of alpha, beta and theta (the function)
  setEPS()
  postscript(file = paste0("./Plots/ck_sp_", i, "_d", d, "alpha.eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  bias = apply(e[1,,], 1, function(x) abs(mean(x) - true_a))
  RMSE = apply(e[1,,], 1, function(x) sqrt(mean((x - true_a)^2)))
  plot(m_values[J], bias[J], type="l", col="blue", ylim = c(min(bias[J]), max(RMSE[J])),
       xlab="m", ylab="Bias and RMSE")
  lines(m_values[J], RMSE[J], col="blue", lty=2)
  dev.off()
  
  setEPS()
  postscript(file = paste0("./Plots/ck_sp_", i, "_d", d, "beta.eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  bias = apply(e[2,,], 1, function(x) abs(mean(x) - true_s))
  RMSE = apply(e[2,,], 1, function(x) sqrt(mean((x - true_s)^2)))
  plot(m_values[J], bias[J], type="l", col="blue", ylim = c(min(bias[J]), max(RMSE[J])),
       xlab="m", ylab="Bias and RMSE")
  lines(m_values, RMSE, col="blue", lty=2)
  dev.off()
  
  setEPS()
  postscript(file = paste0("./Plots/ck_sp_", i, "_d", d, "theta.eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(m_values[J], rowMeans(e[3,J,]), type="l", col="blue", xlab="m", ylab="Mean supremum error")
  dev.off()
  
  # Spatial boxplots
  ptp = c(252, 567, 544, 753, 521)
  dtp = distance[ptp]
  j = 6
  setEPS()
  postscript(file = paste0("./Plots/Boxplot_sp_", i, "_2d.eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  boxplot(t(all_estimators_sp[ptp, j, ]),
          col = "lightgreen", at = dtp-0.05, boxwex = 0.1,
          xlab = expression(Delta), ylab = expression(paste("Estimators of ", theta(Delta))),
          xaxt='n')
  boxplot(t(apply(e[,j,], 2, function(x) sapply(dtp, theta_f, a=x[1], s=x[2]))),
          col = "lightblue", at = dtp+0.05, boxwex = 0.1,
          xaxt='n', add=T)
  for(dist in dtp){
    lines(dist + c(-0.1, 0.1), rep(theta_f(dist, true_a, true_s), 2),
          lw = 4, col = "darkgoldenrod1")
  }
  axis(side = 1, at = dtp, labels = round(dtp, 1))
  dev.off()
  
  # Curves
  stp = 1:3
  setEPS()
  postscript(file = paste0("./Plots/Curves_sp_", i, "_d", d, ".eps"), width = 5, height = 5.5)
  par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
  plot(NULL, xlim = c(0, 3), ylim = c(0.5, 0.75), xlab=expression(Delta),
       ylab = expression(paste(theta, "(", Delta, "; ",  hat(alpha), ",", hat(beta), ")")))
  for(s in stp){
    lines((1:400)/100, sapply((1:400)/100, theta_f, a=e[1,j,s], s=e[2,j,s]), lw=0.5)
  }
  lines((1:400)/100, sapply((1:400)/100, theta_f, a=true_a, s=true_s), col="blue", lw=4)
  dev.off()
  
}



########################################
######## DATA ANALYSIS PLOTS ###########
########################################

If the data set RainData.Rdata is available, run DataEstimators.R and then uncomment and run
the following code to produce the "real data" plots from Section 6 of the paper


# source("DataEstimators.R") # Run once, to produce the estimators
load("DataEstimators.Rdata")
load("SpInfo.Rdata")


# Scatter plot with curve
j = 4
setEPS()
postscript(file = paste0("./Plots/Estimators_m", m_values[j], ".eps"), width = 5, height = 5.5)
par(cex = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(distance, all_estimators_data[1:np,j], col = "blue",
     xlab = "Distance (units of latitude)", ylab = expression(paste("Estimated ", theta)))
lines(sort(distance), sapply(sort(distance), theta_f,
                                 a = all_estimators_data[np+1,j], s = all_estimators_data[np+2,j]))
# lines(sort(distance), sapply(sort(distance), theta_f,
#                              a = all_estimators_data[np+3,j], s = all_estimators_data[np+4,j]), lty = 2)
dev.off()

# Bunch of curves
J = 2*(1:5)
setEPS()
postscript(file = "./Plots/Curves1.eps", width = 5, height = 5.5)
par(cex = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(NULL, xlim = sort(distance)[c(1, np)], ylim = c(0.5, 0.75),
     xlab = expression(paste("Distance in units of latitude ", (Delta))),
     ylab = expression(paste(theta, "(", Delta, "; ",  hat(alpha), ",", hat(beta), ")")))
for(j in J){
  lines(sort(distance), sapply(sort(distance), theta_f,
                               a = all_estimators_data[np+1,j], s = all_estimators_data[np+2,j]), col = j/2)
}
legend(x=0, y=0.75, legend=paste("m = ", m_values[J]), lty=1, col=1:5)
dev.off()

# Map
setEPS()
postscript(file = "./Plots/Map.eps", width = 5, height = 5.5)
par(cex = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
plot(gr, pch=16, col="red", xlim = c(140.5, 147.5), ylim = c(-40, -33),
     xlab = "Longitude (degrees)", ylab = "Latitude (degrees)")
oz(add = TRUE)
dev.off()

# Histogram of distances
setEPS()
postscript(file = "./Plots/Distances.eps", width = 5, height = 5.5)
par(cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pty="s", mar = c(5,5,4,2) +.1)
hist(distance, xlab = "Distance (units of latitude)", main = "")
dev.off()


