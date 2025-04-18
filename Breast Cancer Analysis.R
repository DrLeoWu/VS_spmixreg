
################# load libraries
library(mixtools)
library(Hmisc)
library(graphics)
#library(mixreg)  
library(parallel)
#install.packages("devtools")
#devtools::install_github("thiyangt/fformpp")
library(fformpp)
library(VariableScreening)
library(MixSemiRob)
library(ncvreg)

#####  define relevant functions

###################################################################################################
#####  define the penalty functions (consider Lasso SCAD and HARD)
##### The linear quatratic approximation method only needs the derivative of peinstall.packages("mixreg")nalty functions 
pLasso_deri <- function(x,gama,n){
  return(gama*sqrt(n)*sign(x))
}

pHard_deri <- function(x,gama,n){
  return(-2*(sqrt(n)*abs(x)-gama)*sqrt(n)*sign(x)*1*(sqrt(n)*abs(x)<gama))
}

pSCAD_deri <- function(x,gama,n,a){
  return(gama*sqrt(n)*1*(sqrt(n)*abs(x)<gama)+sqrt(n)*(a*gama-sqrt(n)*abs(x))*1*(sqrt(n)*abs(x)<a*gama)/(a-1)*1*(sqrt(n)*abs(x)>gama))
  
########################################################################################################################################
############################################  define EM algorithms for the oracle setup
############################################
#######################################################################################################################################
  
  semipnormmixEM1_orac <- function(z,x1,x2,y,N,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter,epsilon,bw){
    #u <- seq(min(z),max(z),l=N)
    p1 <- ncol(x1)
    p2 <- ncol(x2)
    l <- 0
    beta1_EM <- initl_beta1  ## n by p
    beta2_EM <- initl_beta2  ## n by p
    sgm21_EM <- initl_sgm1^2    ## n by 1
    sgm22_EM <- initl_sgm2^2   ## n by 1
    pi1_EM <- initl_pi1      ## n by 1
    pi2_EM <- initl_pi2      ## n by 1
    #pi_EM <- cbind(initl_pi1,initl_pi2) ## n by 2
    theta_old <- rep(0,(p1+p2+4)*N)
    theta_diff <- epsilon 
    
    while(l < maxiter & theta_diff >= epsilon){
      l <- l + 1
      pi_EM <- cbind(pi1_EM,pi2_EM)
      sgm2_EM <- cbind(sgm21_EM,sgm22_EM)
      r <- pi_EM*dnorm(cbind(y,y),mean=cbind(diag(beta1_EM%*%t(x1)),diag(beta2_EM%*%t(x2))),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(diag(beta1_EM%*%t(x1)),diag(beta2_EM%*%t(x2))),sd=sqrt(sgm2_EM)))
      
      for(j in 1:N){
        pi1_EM[j] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
        pi2_EM[j] <- 1-pi1_EM[j]
        #pi2[j] <- r[,2]%*%EpaKernel(x-u[j],bw)/sum(EpaKernel(x-u[j],bw))
        beta1_EM[j,] <- c(solve(t(x1)%*%diag(r[,1]*EpaKernel(z-z[j],bw))%*%x1)%*%t(x1)%*%diag(r[,1]*EpaKernel(z-z[j],bw))%*%y)
        beta2_EM[j,] <- c(solve(t(x2)%*%diag(r[,2]*EpaKernel(z-z[j],bw))%*%x2)%*%t(x2)%*%diag(r[,2]*EpaKernel(z-z[j],bw))%*%y)
        #m1[j] <- (r[,1]*EpaKernel(x-u[j],bw))%*%y/(r[,1]%*%EpaKernel(x-u[j],bw))
        #m2[j] <- (r[,2]*EpaKernel(x-u[j],bw))%*%y/(r[,2]%*%EpaKernel(x-u[j],bw))
        sgm21_EM[j] <- (r[,1]*EpaKernel(z-z[j],bw))%*%((y-x1%*%beta1_EM[j,])^2)/(r[,1]%*%EpaKernel(z-z[j],bw))
        sgm22_EM[j] <- (r[,2]*EpaKernel(z-z[j],bw))%*%((y-x2%*%beta2_EM[j,])^2)/(r[,2]%*%EpaKernel(z-z[j],bw))
      }
      
      theta_new <- c(pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
      theta_diff <- max(abs(theta_new - theta_old))
      theta_old <- theta_new
      
    }
    #if (l >= maxiter) {print("The algorithm does not converge")}
    return(list(pi1=pi1_EM,pi2=pi2_EM,beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
  }
  
  semipnormmixEM2_orac <- function(z,x1,x2,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter,epsilon){
    p1 <- ncol(x1)
    p2 <- ncol(x2)
    N <- nrow(x1)
    l <- 0
    beta1_EM <- initl_beta1  ## p by 1
    beta2_EM <- initl_beta2  ## p by 1
    sgm21_EM <- initl_sgm21    
    sgm22_EM <- initl_sgm22   
    pi_EM <- cbind(pi1z,pi2z)
    theta_old <- rep(0,p1+p2+2)
    theta_diff <- epsilon 
    
    while(l < maxiter & theta_diff >= epsilon){
      l <- l + 1
      pi_EM <- cbind(pi1z,pi2z)
      sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
      r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x1%*%beta1_EM,x2%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x1%*%beta1_EM,x2%*%beta2_EM),sd=sqrt(sgm2_EM)))
      beta1_EM <- c(solve(t(x1)%*%diag(r[,1])%*%x1)%*%t(x1)%*%diag(r[,1])%*%y)
      beta2_EM <- c(solve(t(x2)%*%diag(r[,2])%*%x2)%*%t(x2)%*%diag(r[,2])%*%y)
      sgm21_EM <- r[,1]%*%((y-x1%*%beta1_EM)^2)/sum(r[,1])
      sgm22_EM <- r[,2]%*%((y-x2%*%beta2_EM)^2)/sum(r[,2])
      
      theta_new <- c(beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
      theta_diff <- max(abs(theta_new - theta_old))
      theta_old <- theta_new
      
    }
    if (l >= maxiter) {print("The algorithm does not converge")}
    return(list(beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
  }
  
  semipnormmixEM3_orac <- function(z,x1,x2,y,N,beta1,beta2,sgm21,sgm22,initl_pi1,initl_pi2,maxiter,epsilon,bw){
    #u <- seq(min(z),max(z),l=N)
    p1 <- ncol(x1)
    p2 <- ncol(x2)
    l <- 0
    pi1_EM <- initl_pi1      ## n by 1
    pi2_EM <- initl_pi2      ## n by 1
    sgm2_EM <- cbind(rep(sgm21,N),rep(sgm22,N))
    #pi_EM <- cbind(initl_pi1,initl_pi2) ## n by 2
    theta_old <- rep(0,N*2)
    theta_diff <- epsilon 
    
    while(l < maxiter & theta_diff >= epsilon){
      l <- l + 1
      pi_EM <- cbind(pi1_EM,pi2_EM)
      r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x1%*%beta1,x2%*%beta2),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x1%*%beta1,x2%*%beta2),sd=sqrt(sgm2_EM)))
      
      for(j in 1:N){
        pi1_EM[j] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
        pi2_EM[j] <- 1-pi1_EM[j]
      }
      
      theta_new <- c(pi1_EM,pi2_EM)
      theta_diff <- max(abs(theta_new - theta_old))
      theta_old <- theta_new
      
    }
    #if (l >= maxiter) {print("The algorithm does not converge")}
    return(list(pi1=pi1_EM,pi2=pi2_EM,iter=l))
  }
}

######################################################################### accelerating EM algorithms 

semipnormmixEM_L1 <- function(z,x,y,initl_pi1,initl_pi2,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon,bw){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22  
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*p+2+2*N)
  theta_diff <- epsilon
  pi_EM <- cbind(pi1_EM,pi2_EM)
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    r1 <- 
      beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])<=1){break}
    Sig1 <- diag(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])<=1){break}
    Sig2 <- diag(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
    if (is.singular(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)|is.singular(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)){singular <- TRUE; break}
    else {singular <- FALSE}
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    for(j in 1:N){
      pi_EM[j,1] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
      pi_EM[j,2] <- 1-pi_EM[j,1]
    }
    
    
    theta_new <- c(pi_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  if (singular){
    pi_EM <- cbind(rep(NA,length(pi1_EM)),rep(NA,length(pi2_EM)))
    beta1_EM <- beta2_EM <-rep(NA,p);sgm21_EM <- sgm22_EM <- NA
    l <- maxiter
  }
  return(list(pi1= pi_EM[,1],pi2= pi_EM[,2],beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM_HARD <- function(z,x,y,initl_pi1,initl_pi2,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon,bw){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22  
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*p+2+2*N)
  theta_diff <- epsilon
  pi_EM <- cbind(pi1_EM,pi2_EM)
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])<=1){break}
    Sig1 <- diag(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])<=1){break}
    Sig2 <- diag(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
    if (is.singular(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)|is.singular(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)){singular <- TRUE;break}
    else {singular <- FALSE}
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    for(j in 1:N){
      pi_EM[j,1] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
      pi_EM[j,2] <- 1-pi_EM[j,1]
    }
    
    
    theta_new <- c(pi_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  if (singular){
    pi_EM <- cbind(rep(NA,length(pi1_EM)),rep(NA,length(pi2_EM)))
    beta1_EM <- beta2_EM <-rep(NA,p);sgm21_EM <- sgm22_EM <- NA
    l <- maxiter
  }
  return(list(pi1= pi_EM[,1],pi2= pi_EM[,2],beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM_SCAD <- function(z,x,y,initl_pi1,initl_pi2,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon,bw){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22  
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*p+2+2*N)
  theta_diff <- epsilon
  pi_EM <- cbind(pi1_EM,pi2_EM)
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])<=1){break}
    Sig1 <- diag(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])
    if (length(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])<=1){break}
    Sig2 <- diag(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])
    if (is.singular(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)|is.singular(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)){singular <- TRUE;break}
    else {singular <- FALSE}
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    for(j in 1:N){
      pi_EM[j,1] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
      pi_EM[j,2] <- 1-pi_EM[j,1]
    }
    
    
    theta_new <- c(pi_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  if (singular){
    pi_EM <- cbind(rep(NA,length(pi1_EM)),rep(NA,length(pi2_EM)))
    beta1_EM <- beta2_EM <-rep(NA,p);sgm21_EM <- sgm22_EM <- NA
    l <- maxiter
  }
  return(list(pi1= pi_EM[,1],pi2= pi_EM[,2],beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}
############################################################################################################################################################
###############   Need to load the breast cancer data (BCdata)

X <- BCdata[,-c(1:4)]
y <- BCdata$size
z <- BCdata$age
N <- nrow(X)
##### perform feature screening first 
rank <- screenIID(X = data.reduce[,-c(1:3)], Y = BCdata$size, method="DC-SIS")
X_temp <- X[,rank$rank<=ceiling(N/log(N))]
Xd_temp <- cbind(rep(1,nrow(X_temp)),X_temp)

##### determine the initial values 
lmr <- mixreg(X_temp,y)
initl_beta1 <- lmr$beta[1,] 
initl_beta2 <- lmr$beta[2,] 
initl_sgm21 <- lmr$sigma[1]^2
initl_sgm22 <- lmr$sigma[2]^2
initl_pi1 <- rep(lmr$pi[1],N)
initl_pi2 <- 1 - initl_pi1

##### cross validation on squared prediction errors to determine optimal tuning parameters
gama1_grid <- seq(0.2,1.4,0.2)
gama2_grid <- seq(0.2,1.4,0.2)
h_grid <- 10 
max_tpar <- expand.grid(gama1_grid,gama2_grid,h_grid)  ## the matrix for all combinations of tuning parameters
k_fold <- 10
p <- ncol(Xd_temp)

set.seed(7767)
ind_cv <- sample(x = c(rep(1:k_fold, each = floor(nrow(X_temp) / k_fold)),rep(1,nrow(X_temp)-floor(nrow(X_temp) / k_fold)*k_fold)), size = nrow(X_temp))

cv_spmixreg_fold <- function(i_test){
  library(fformpp)
  cv_error_SCAD <- rep(0,nrow(max_tpar))
  for(i_para in 1:nrow(max_tpar)){
    
    y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
    z_test <- z[ind_cv==i_test] ; z_train <- z[ind_cv!=i_test]
    xd_test <- Xd_temp[ind_cv==i_test,] ; xd_train <- Xd_temp[ind_cv!=i_test,]
    pi1z <- rep(lmr$pi[1],length(z_train)) ; pi2z <- 1 - pi1z
    
    out_acSCAD <- semipnormmixEM_SCAD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
    l_ini <- 0
    while(out_acSCAD$iter==1000 & l_ini < 20) {
      l_ini <- l_ini + 1
      initl_beta1_new <- initl_beta1 + rnorm(p,sd=0.2)
      initl_beta2_new <- initl_beta2 + rnorm(p,sd=0.2)
      out_acSCAD <- semipnormmixEM_SCAD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
      print(l_ini)
    }
    
    r1_SCAD <- out_acSCAD$pi1*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1))/(out_acSCAD$pi1*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1))+out_acSCAD$pi2*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta2,sd=sqrt(out_acSCAD$sigma2_2)))
    r2_SCAD <- 1 - r1_SCAD
    
    pi1_SCADtest <- sapply(z_test,function(q) r1_SCAD%*%EpaKernel(z_train-q,max_tpar[i_para,3])/sum(EpaKernel(z_train-q,max_tpar[i_para,3])))
    pi2_SCADtest <- 1 - pi1_SCADtest 
    
    m1_SCADtest <- xd_test%*%out_acSCAD$beta1
    m2_SCADtest <- xd_test%*%out_acSCAD$beta2
    
    r1_SCADtest <- pi1_SCADtest*dnorm(y_test,mean=m1_SCADtest,sd=sqrt(out_acSCAD$sigma2_1))/(pi1_SCADtest*dnorm(y_test,mean=m1_SCADtest,sd=sqrt(out_acSCAD$sigma2_1))+pi2_SCADtest*dnorm(y_test,mean=m2_SCADtest,sd=sqrt(out_acSCAD$sigma2_2)))
    r2_SCADtest <- 1 - r1_SCADtest
    
    cv_error_SCAD[i_para] <- sum((y_test - (r1_SCADtest*m1_SCADtest+r2_SCADtest*m2_SCADtest))^2)
    
  }
  return(list(error = cv_error_SCAD))
}

the_cluster_cvfold <- makeCluster(10)
clusterSetRNGStream(the_cluster_cvfold,49)
clusterExport(the_cluster_cvfold,c("pLasso_deri","pHard_deri","pSCAD_deri","semipnormmixEM_L1","semipnormmixEM_HARD","semipnormmixEM_SCAD","EpaKernel",
                                   "p","N","z","k_fold","max_tpar","initl_beta1","initl_beta2","initl_sgm21","initl_sgm22","ind_cv","y","lmr","X_temp","Xd_temp"))
#tpar_cv <- parLapply(the_cluster_cv,1:30,cv_spmixreg)
cv_error <- clusterCall(cl = the_cluster_cvfold, cv_spmixreg_fold,1:10)
stopCluster(the_cluster_cvfold)

errors <- NULL 
for(i in 1:10){
  errors <- rbind(errors,cv_error[[i]]$error)
}

tpar_opt <- max_tpar[apply(errors,2,mean)==min(apply(errors,2,mean)),]
print(tpar_opt)  ### the optimal gama1 is 1.2 and optimal gama2 is 0.4

out_acSCAD <- semipnormmixEM_SCAD(z,x=Xd_temp,y,initl_pi1,initl_pi2,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=1.2,gama2=0.4,maxiter=1000,epsilon=0.001,bw=10)
set.seed(666)
l_ini <- 0
while(out_acSCAD$iter==1000 & l_ini < 20) {
  l_ini <- l_ini + 1
  initl_beta1_new <- initl_beta1 + rnorm(p,sd=0.2)
  initl_beta2_new <- initl_beta2 + rnorm(p,sd=0.2)
  out_acSCAD <- semipnormmixEM_SCAD(z=z,x=Xd_temp,y,initl_pi1,initl_pi2,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=1.2,gama2=0.4,maxiter=1000,epsilon=0.001,bw=10)
  print(l_ini)
}

##### print the selected genes 
probe.select1 <- colnames(X_temp)[out_acSCAD$beta1[-1]!=0]
probe.select2 <- colnames(X_temp)[out_acSCAD$beta2[-1]!=0]

print(probe_gene[probe.select1,2])
print(probe_gene[probe.select2,2])

length(probe.select1)
length(probe.select2)

##### plot the estimated mixing proportion function
plot(z[order(z)],out_acSCAD$pi1[order(z)],type = "l",ylim=c(0,0.6),main="Estimated Proportion Function",xlab = "Age",ylab="T2-type Group Proportion")

##### plot the estimated densities
post_pi <- (out_acSCAD$pi1*dnorm(y,mean=Xd_temp%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1)))/(out_acSCAD$pi1*dnorm(y,mean=Xd_temp%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1))+(1-out_acSCAD$pi1)*dnorm(y,mean=Xd_temp%*%out_acSCAD$beta2,sd=sqrt(out_acSCAD$sigma2_2)))
y1 <- y[post_pi>0.5]
y2 <- y[post_pi<=0.5]
dy1 <- density(y1,bw=0.6)
dy2 <- density(y2,bw=0.4)

plot(dy2, lwd = 2, col = "red", main = "Conditional Densities of Two Components",xlim=c(0,5),ylim=c(0,0.6),xlab="Tumor Size")
lines(dy1, lwd = 2,col="blue")
legend(x = 3.5,y=0.6, legend=c("Component 1", "Component 2"),  
       fill = c("blue","red"),bty = "n"
)

################################################## cross validation for penalized linear regression without mixture (SCAD)
library(ncvreg)
lambda.grid <- seq(0.1,0.005,by=-0.005)
k_fold <- 5
set.seed(7767)
ind_cv <- sample(x = c(rep(1:k_fold, each = floor(nrow(X_temp) / k_fold)),rep(1,nrow(X_temp)-floor(nrow(X_temp) / k_fold)*k_fold)), size = nrow(X_temp))
error <- rep(0,length(lambda.grid))


for(i_test in 1:k_fold){
  
  y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
  X_test <- X_temp[ind_cv==i_test,] ; X_train <- X_temp[ind_cv!=i_test,]
  Xd_test <- cbind(rep(1,nrow(X_test)),X_test)
  ## fit the model for the training set
  fit_train <- ncvreg(X_train,y_train,penalty="SCAD",lambda=lambda.grid)
  
  for(i_grid in 1:length(lambda.grid)){
    error[i_grid] <- error[i_grid] + sum((y_test - c(Xd_test%*%fit_train$beta[,i_grid]))^2)
  }
  
}

print(error)
opt.lambda <- lambda.grid[error==min(error)]
print(opt.lambda)
# the optimal lambda is 0.05
fit <- ncvreg(X_temp,y,penalty="SCAD",lambda=c(opt.lambda,0))
length(fit$beta[,1][fit$beta[,1]!=0])
probe_1 <- fit$beta[,1][fit$beta[,1]!=0]
probe_1 <- probe_1[-1]
probe_1 <- names(probe_1)
###### print the selected genes from the penalized linear regression 
print(probe_gene[probe_1,2])


###########################################################################################################################################################
#########################   Compare the performance of semiparametric mixtures and linear regression on selected genes

############# First compare their in-sample performance based on squared prediction errors 

X_LM <- X_temp[,probe_1]
Xd_LM <- cbind(rep(1,nrow(X_LM)),X_LM)

X1_SPMIX <- X_temp[,probe.select1]
Xd1_SPMIX <- cbind(rep(1,nrow(X1_SPMIX)),X1_SPMIX)
X2_SPMIX <- X_temp[,probe.select2]
Xd2_SPMIX <- cbind(rep(1,nrow(X2_SPMIX)),X2_SPMIX)

initl_beta1 <- t(matrix(rep(out_acSCAD$beta1[out_acSCAD$beta1!=0],N),ncol=N))
initl_beta2 <- t(matrix(rep(out_acSCAD$beta2[out_acSCAD$beta2!=0],N),ncol=N))
initl_sgm1 <- rep(sqrt(out_acSCAD$sigma2_1),N)
initl_sgm2 <- rep(sqrt(out_acSCAD$sigma2_2),N)
initl_pi1 <- out_acSCAD$pi1
initl_pi2 <- out_acSCAD$pi2

testout <- semipnormmixEM1_orac(z,x1=Xd1_SPMIX,x2=Xd2_SPMIX,y,N,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter=1000,epsilon=0.001,bw=20)

initl_beta1 <- apply(testout$beta1,2,mean)
initl_beta2 <- apply(testout$beta2,2,mean)
initl_sgm21 <- mean(testout$sigma2_1)
initl_sgm22 <- mean(testout$sigma2_2)
pi1z <- testout$pi1
pi2z <- testout$pi2

testout2 <- semipnormmixEM2_orac(z=z,x1=Xd1_SPMIX,x2=Xd2_SPMIX,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter=1000,epsilon=0.001)
testout3 <- semipnormmixEM3_orac(z=z,x1=Xd1_SPMIX,x2=Xd2_SPMIX,y,N,beta1=testout2$beta1,beta2=testout2$beta2,sgm21=testout2$sigma2_1,sgm22=testout2$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=8)

post_pi <- (testout3$pi1*dnorm(y,mean=Xd1_SPMIX%*%testout2$beta1,sd=sqrt(testout2$sigma2_1)))/(testout3$pi1*dnorm(y,mean=Xd1_SPMIX%*%testout2$beta1,sd=sqrt(testout2$sigma2_1))+(1-testout3$pi1)*dnorm(y,mean=Xd2_SPMIX%*%testout2$beta2,sd=sqrt(testout2$sigma2_2)))

eLM <- (y - Xd_LM%*%solve(t(Xd_LM)%*%Xd_LM)%*%t(Xd_LM)%*%y)^2
eSPMIX <- (y - post_pi*Xd1_SPMIX%*%testout2$beta1 - (1- post_pi)*Xd2_SPMIX%*%testout2$beta2)^2


#############  compare their out-of-sample performance based on squared prediction errors
set.seed(677)
N_train <- floor(0.9*nrow(X_temp))
train_ind <- sample.int(n = nrow(X_temp), size = N_train)


initl_beta1 <- t(matrix(rep(out_acSCAD$beta1[out_acSCAD$beta1!=0],N_train),ncol=N_train))
initl_beta2 <- t(matrix(rep(out_acSCAD$beta2[out_acSCAD$beta2!=0],N_train),ncol=N_train))
initl_sgm1 <- rep(sqrt(out_acSCAD$sigma2_1),N_train)
initl_sgm2 <- rep(sqrt(out_acSCAD$sigma2_2),N_train)
initl_pi1 <- out_acSCAD$pi1[train_ind]
initl_pi2 <- out_acSCAD$pi2[train_ind]
z_train <- z[train_ind];z_test <- z[-train_ind]
y_train <- y[train_ind];y_test <- y[-train_ind]
Xd1_SPMIX_train <- Xd1_SPMIX[train_ind,];Xd1_SPMIX_test <- Xd1_SPMIX[-train_ind,]
Xd2_SPMIX_train <- Xd2_SPMIX[train_ind,];Xd2_SPMIX_test <- Xd2_SPMIX[-train_ind,]

testout_train <- semipnormmixEM1_orac(z=z_train,x1=Xd1_SPMIX_train,x2=Xd2_SPMIX_train,y=y_train,N=N_train,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter=1000,epsilon=0.001,bw=20)

initl_beta1 <- apply(testout_train$beta1,2,mean)
initl_beta2 <- apply(testout_train$beta2,2,mean)
initl_sgm21 <- mean(testout_train$sigma2_1)
initl_sgm22 <- mean(testout_train$sigma2_2)
pi1z <- testout_train$pi1
pi2z <- testout_train$pi2

testout2_train <- semipnormmixEM2_orac(z=z_train,x1=Xd1_SPMIX_train,x2=Xd2_SPMIX_train,y=y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter=1000,epsilon=0.001)
testout3_train <- semipnormmixEM3_orac(z=z_train,x1=Xd1_SPMIX_train,x2=Xd2_SPMIX_train,y=y_train,N=N_train,beta1=testout2_train$beta1,beta2=testout2_train$beta2,sgm21=testout2_train$sigma2_1,sgm22=testout2_train$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=4)

r1_train <- (testout3_train$pi1*dnorm(y_train,mean=Xd1_SPMIX_train%*%testout2_train$beta1,sd=sqrt(testout2_train$sigma2_1)))/(testout3_train$pi1*dnorm(y_train,mean=Xd1_SPMIX_train%*%testout2_train$beta1,sd=sqrt(testout2_train$sigma2_1))+(1-testout3_train$pi1)*dnorm(y_train,mean=Xd2_SPMIX_train%*%testout2_train$beta2,sd=sqrt(testout2_train$sigma2_2)))
r2_train <- 1 - r1_train

pi1_test <- sapply(z_test,function(q) r1_train%*%EpaKernel(z_train-q,10)/sum(EpaKernel(z_train-q,10)))
pi2_test <- 1 - pi1_test 

m1_test <- Xd1_SPMIX_test%*%testout2_train$beta1
m2_test <- Xd2_SPMIX_test%*%testout2_train$beta2

r1_test <- pi1_test*dnorm(y_test,mean=m1_test,sd=sqrt(testout2_train$sigma2_1))/(pi1_test*dnorm(y_test,mean=m1_test,sd=sqrt(testout2_train$sigma2_1))+pi2_test*dnorm(y_test,mean=m2_test,sd=sqrt(testout2_train$sigma2_2)))
r2_test <- 1 - r1_test

eSPMIX_test <- (y_test - (r1_test*m1_test+r2_test*m2_test))^2
eLM_test <- (y_test - Xd_LM_test%*%solve(t(Xd_LM_train)%*%Xd_LM_train)%*%t(Xd_LM_train)%*%y_train)^2

####### Give the boxplots 
par(mfrow=c(1,2))

method <- c(rep("Linear",N),rep("Mixture",N))
err <- c(eLM,eSPMIX)
boxplot(err~method,ylab = "Squared Prediction Error",main="Comparison of In-Sample Performance",xlab="Method")

N_test <- N -  N_train
method <- c(rep("Linear",N_test),rep("Mixture",N_test))
err_test <- c(eLM_test,eSPMIX_test)
boxplot(err_test~method,ylab = "Squared Prediction Error",main="Comparison of Out-of-Sample Performance",xlab="Method")

par(mfrow=c(1,1))



























