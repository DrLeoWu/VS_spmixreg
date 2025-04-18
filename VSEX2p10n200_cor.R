##########################################################################################################
##########################################################################################################
############## Simulation for variable selection (p=10, n=200, EX2 correlated covariates)
##########################################################################################################
##########################################################################################################

library(mixtools)
library(Hmisc)
library(graphics)
#library(mixreg)  
library(parallel)

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
}

####################################################################################################
##### define a function to implement the modified EM agorithm for 2.3 (for the number of component c =2)
##### N is the number of grid points, set N=n as the sample size is not large 
semipnormmixEM1 <- function(z,x,y,N,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter,epsilon,bw){
  #u <- seq(min(z),max(z),l=N)
  p <- ncol(x)
  l <- 0
  beta1_EM <- initl_beta1  ## n by p
  beta2_EM <- initl_beta2  ## n by p
  sgm21_EM <- initl_sgm1^2    ## n by 1
  sgm22_EM <- initl_sgm2^2   ## n by 1
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  #pi_EM <- cbind(initl_pi1,initl_pi2) ## n by 2
  theta_old <- rep(0,2*(p+2)*N)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1_EM,pi2_EM)
    sgm2_EM <- cbind(sgm21_EM,sgm22_EM)
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(diag(beta1_EM%*%t(x)),diag(beta2_EM%*%t(x))),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(diag(beta1_EM%*%t(x)),diag(beta2_EM%*%t(x))),sd=sqrt(sgm2_EM)))
    
    for(j in 1:N){
      pi1_EM[j] <- r[,1]%*%EpaKernel(z-z[j],bw)/sum(EpaKernel(z-z[j],bw))
      pi2_EM[j] <- 1-pi1_EM[j]
      #pi2[j] <- r[,2]%*%EpaKernel(x-u[j],bw)/sum(EpaKernel(x-u[j],bw))
      beta1_EM[j,] <- c(solve(t(x)%*%diag(r[,1]*EpaKernel(z-z[j],bw))%*%x)%*%t(x)%*%diag(r[,1]*EpaKernel(z-z[j],bw))%*%y)
      beta2_EM[j,] <- c(solve(t(x)%*%diag(r[,2]*EpaKernel(z-z[j],bw))%*%x)%*%t(x)%*%diag(r[,2]*EpaKernel(z-z[j],bw))%*%y)
      #m1[j] <- (r[,1]*EpaKernel(x-u[j],bw))%*%y/(r[,1]%*%EpaKernel(x-u[j],bw))
      #m2[j] <- (r[,2]*EpaKernel(x-u[j],bw))%*%y/(r[,2]%*%EpaKernel(x-u[j],bw))
      sgm21_EM[j] <- (r[,1]*EpaKernel(z-z[j],bw))%*%((y-x%*%beta1_EM[j,])^2)/(r[,1]%*%EpaKernel(z-z[j],bw))
      sgm22_EM[j] <- (r[,2]*EpaKernel(z-z[j],bw))%*%((y-x%*%beta2_EM[j,])^2)/(r[,2]%*%EpaKernel(z-z[j],bw))
    }
    
    theta_new <- c(pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=pi1_EM,pi2=pi2_EM,beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

####################################################################################################
####################################################################################################
##### define a function to implement the modified EM agorithm for 2.4 (for the number of component c =2)
##### N is the number of grid points, set N=n as the sample size is not large 

semipnormmixEM2 <- function(z,x,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22   
  pi_EM <- cbind(pi1z,pi2z)
  theta_old <- rep(0,2*p+2)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1z,pi2z)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM <- c(solve(t(x)%*%diag(r[,1])%*%x)%*%t(x)%*%diag(r[,1])%*%y)
    beta2_EM <- c(solve(t(x)%*%diag(r[,2])%*%x)%*%t(x)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    theta_new <- c(beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM2_L1 <- function(z,x,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22   
  pi_EM <- cbind(pi1z,pi2z)
  theta_old <- rep(0,2*p+2)
  theta_diff <- epsilon
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1z,pi2z)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    theta_new <- c(beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM2_HARD <- function(z,x,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22   
  pi_EM <- cbind(pi1z,pi2z)
  theta_old <- rep(0,2*p+2)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1z,pi2z)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    theta_new <- c(beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM2_SCAD <- function(z,x,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1,gama2,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta1_EM <- initl_beta1  ## p by 1
  beta2_EM <- initl_beta2  ## p by 1
  sgm21_EM <- initl_sgm21    
  sgm22_EM <- initl_sgm22   
  pi_EM <- cbind(pi1z,pi2z)
  theta_old <- rep(0,2*p+2)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1z,pi2z)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])
    if (length(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])
    beta1_EM[beta1_EM!=0] <- c(solve(t(x1_0out)%*%diag(r[,1])%*%x1_0out-Sig1)%*%t(x1_0out)%*%diag(r[,1])%*%y)
    beta2_EM[beta2_EM!=0] <- c(solve(t(x2_0out)%*%diag(r[,2])%*%x2_0out-Sig2)%*%t(x2_0out)%*%diag(r[,2])%*%y)
    sgm21_EM <- r[,1]%*%((y-x%*%beta1_EM)^2)/sum(r[,1])
    sgm22_EM <- r[,2]%*%((y-x%*%beta2_EM)^2)/sum(r[,2])
    
    theta_new <- c(beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    #print(beta1_EM)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

####################################################################################################
##### define a function to implement the modified EM agorithm for 2.5 (for the number of component c =2)
##### N is the number of grid points, set N=n as the sample size is not large 

semipnormmixEM3 <- function(z,x,y,N,beta1,beta2,sgm21,sgm22,initl_pi1,initl_pi2,maxiter,epsilon,bw){
  #u <- seq(min(z),max(z),l=N)
  p <- ncol(x)
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
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1,x%*%beta2),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1,x%*%beta2),sd=sqrt(sgm2_EM)))
    
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

########################################################################################################################################
############################################  define EM algorithms for generalized cross validation
############################################
#######################################################################################################################################

semipnormmixEM2_L1p <- function(z,x,y,pi_mle,initl_beta_pmle,beta_mle,initl_sgm2_pmle,sgm2_mle,gama,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta_EM <- initl_beta_pmle  ## p by 1
  sgm2_EM <- initl_sgm2_pmle
  theta_old <- rep(0,p+1)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    r <- (1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))/((1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))+pi_mle*dnorm(y,mean=x%*%beta_mle,sd=sqrt(sgm2_mle)))
    beta_EM[abs(beta_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    x_0out <- x[,beta_EM!=0] 
    Sig <- diag(pLasso_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    #Sig <- diag(pHard_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    #Sig <- diag(pSCAD_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    beta_EM[beta_EM!=0] <- c(solve(t(x_0out)%*%diag(r)%*%x_0out-Sig)%*%t(x_0out)%*%diag(r)%*%y)
    sgm2_EM <- c(r%*%((y-x%*%beta_EM)^2))/sum(r)
    
    theta_new <- c(beta_EM,sgm2_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pmle=beta_EM,sigma2_pmle=sgm2_EM,iter=l))
}

semipnormmixEM2_HARDp <- function(z,x,y,pi_mle,initl_beta_pmle,beta_mle,initl_sgm2_pmle,sgm2_mle,gama,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta_EM <- initl_beta_pmle  ## p by 1
  sgm2_EM <- initl_sgm2_pmle
  theta_old <- rep(0,p+1)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    r <- (1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))/((1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))+pi_mle*dnorm(y,mean=x%*%beta_mle,sd=sqrt(sgm2_mle)))
    beta_EM[abs(beta_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    x_0out <- x[,beta_EM!=0] 
    #Sig <- diag(pLasso_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    Sig <- diag(pHard_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    #Sig <- diag(pSCAD_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    beta_EM[beta_EM!=0] <- c(solve(t(x_0out)%*%diag(r)%*%x_0out-Sig)%*%t(x_0out)%*%diag(r)%*%y)
    sgm2_EM <- c(r%*%((y-x%*%beta_EM)^2))/sum(r)
    
    theta_new <- c(beta_EM,sgm2_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pmle=beta_EM,sigma2_pmle=sgm2_EM,iter=l))
}

semipnormmixEM2_SCADp <- function(z,x,y,pi_mle,initl_beta_pmle,beta_mle,initl_sgm2_pmle,sgm2_mle,gama,maxiter,epsilon){
  p <- ncol(x)
  N <- nrow(x)
  l <- 0
  beta_EM <- initl_beta_pmle  ## p by 1
  sgm2_EM <- initl_sgm2_pmle
  theta_old <- rep(0,p+1)
  theta_diff <- epsilon 
  
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    r <- (1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))/((1-pi_mle)*dnorm(y,mean=x%*%beta_EM,sd=sqrt(sgm2_EM))+pi_mle*dnorm(y,mean=x%*%beta_mle,sd=sqrt(sgm2_mle)))
    beta_EM[abs(beta_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    x_0out <- x[,beta_EM!=0] 
    #Sig <- diag(pLasso_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    #Sig <- diag(pHard_deri(beta_EM[beta_EM!=0],gama,N)/beta_EM[beta_EM!=0])
    Sig <- diag(pSCAD_deri(beta_EM[beta_EM!=0],gama,N,a=3.7)/beta_EM[beta_EM!=0])
    beta_EM[beta_EM!=0] <- c(solve(t(x_0out)%*%diag(r)%*%x_0out-Sig)%*%t(x_0out)%*%diag(r)%*%y)
    sgm2_EM <- c(r%*%((y-x%*%beta_EM)^2))/sum(r)
    
    theta_new <- c(beta_EM,sgm2_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    #print(theta_diff)
    theta_old <- theta_new
    
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pmle=beta_EM,sigma2_pmle=sgm2_EM,iter=l))
}

######################################################################### Fully iterative EM algorithms  

semipnormmixEM2_fiSCAD <- function(z,x,y,N,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,initl_pi1,initl_pi2,gama1,gama2,maxiter,maxitf=100,epsilon,bw){
  
  p <- ncol(x)
  l <- 0
  beta1_EM <- initl_beta1  ## dim = p
  beta2_EM <- initl_beta2  ## dim = p
  sgm21_EM <- initl_sgm21    ## dim = 1
  sgm22_EM <- initl_sgm22   ## dim = 1
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*(N+1+p))
  theta_diff <- epsilon 
  
  while(l < maxitf & theta_diff >= epsilon){
    l <- l + 1
    step1_EM <- semipnormmixEM2_SCAD(z,x,y,pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM,gama1,gama2,maxiter,epsilon)
    step2_EM <-  semipnormmixEM3(z,x,y,N,beta1=step1_EM$beta1,beta2=step1_EM$beta2,sgm21=step1_EM$sigma2_1,sgm22=step1_EM$sigma2_2,initl_pi1=pi1_EM,initl_pi2=pi2_EM,maxiter,epsilon,bw)
    
    beta1_EM <- step1_EM$beta1  
    beta2_EM <- step1_EM$beta2  
    sgm21_EM <- step1_EM$sigma2_1    
    sgm22_EM <- step1_EM$sigma2_2  
    pi1_EM <- step2_EM$pi1      
    pi2_EM <- step2_EM$pi1      
    
    theta_new <- c(pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    print(l)
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=pi1_EM,pi2=pi2_EM,beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM2_fiL1 <- function(z,x,y,N,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,initl_pi1,initl_pi2,gama1,gama2,maxiter,maxitf=100,epsilon,bw){
  
  p <- ncol(x)
  l <- 0
  beta1_EM <- initl_beta1  ## dim = p
  beta2_EM <- initl_beta2  ## dim = p
  sgm21_EM <- initl_sgm21    ## dim = 1
  sgm22_EM <- initl_sgm22   ## dim = 1
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*(N+1+p))
  theta_diff <- epsilon 
  
  while(l < maxitf & theta_diff >= epsilon){
    l <- l + 1
    step1_EM <- semipnormmixEM2_L1(z,x,y,pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM,gama1,gama2,maxiter,epsilon)
    step2_EM <-  semipnormmixEM3(z,x,y,N,beta1=step1_EM$beta1,beta2=step1_EM$beta2,sgm21=step1_EM$sigma2_1,sgm22=step1_EM$sigma2_2,initl_pi1=pi1_EM,initl_pi2=pi2_EM,maxiter,epsilon,bw)
    
    beta1_EM <- step1_EM$beta1  
    beta2_EM <- step1_EM$beta2  
    sgm21_EM <- step1_EM$sigma2_1    
    sgm22_EM <- step1_EM$sigma2_2  
    pi1_EM <- step2_EM$pi1      
    pi2_EM <- step2_EM$pi1      
    
    theta_new <- c(pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    print(l)
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=pi1_EM,pi2=pi2_EM,beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

semipnormmixEM2_fiHARD <- function(z,x,y,N,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,initl_pi1,initl_pi2,gama1,gama2,maxiter,maxitf=100,epsilon,bw){
  
  p <- ncol(x)
  l <- 0
  beta1_EM <- initl_beta1  ## dim = p
  beta2_EM <- initl_beta2  ## dim = p
  sgm21_EM <- initl_sgm21    ## dim = 1
  sgm22_EM <- initl_sgm22   ## dim = 1
  pi1_EM <- initl_pi1      ## n by 1
  pi2_EM <- initl_pi2      ## n by 1
  theta_old <- rep(0,2*(N+1+p))
  theta_diff <- epsilon 
  
  while(l < maxitf & theta_diff >= epsilon){
    l <- l + 1
    step1_EM <- semipnormmixEM2_HARD(z,x,y,pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM,gama1,gama2,maxiter,epsilon)
    step2_EM <-  semipnormmixEM3(z,x,y,N,beta1=step1_EM$beta1,beta2=step1_EM$beta2,sgm21=step1_EM$sigma2_1,sgm22=step1_EM$sigma2_2,initl_pi1=pi1_EM,initl_pi2=pi2_EM,maxiter,epsilon,bw)
    
    beta1_EM <- step1_EM$beta1  
    beta2_EM <- step1_EM$beta2  
    sgm21_EM <- step1_EM$sigma2_1    
    sgm22_EM <- step1_EM$sigma2_2  
    pi1_EM <- step2_EM$pi1      
    pi2_EM <- step2_EM$pi1      
    
    theta_new <- c(pi1_EM,pi2_EM,beta1_EM,beta2_EM,sgm21_EM,sgm22_EM)
    theta_diff <- max(abs(theta_new - theta_old))
    theta_old <- theta_new
    print(l)
  }
  #if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(pi1=pi1_EM,pi2=pi2_EM,beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
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
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1_EM,pi2_EM)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pLasso_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pLasso_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
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
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1_EM,pi2_EM)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pHard_deri(beta1_EM[beta1_EM!=0],gama1,N)/beta1_EM[beta1_EM!=0])
    if (length(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pHard_deri(beta2_EM[beta2_EM!=0],gama2,N)/beta2_EM[beta2_EM!=0])
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
  while(l < maxiter & theta_diff >= epsilon){
    l <- l + 1
    pi_EM <- cbind(pi1_EM,pi2_EM)
    sgm2_EM <- cbind(rep(sgm21_EM,N),rep(sgm22_EM,N))
    r <- pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM))/rowSums(pi_EM*dnorm(cbind(y,y),mean=cbind(x%*%beta1_EM,x%*%beta2_EM),sd=sqrt(sgm2_EM)))
    beta1_EM[abs(beta1_EM)<0.001] <- 0  ## if the current beta is close to 0, then set it to be 0
    beta2_EM[abs(beta2_EM)<0.001] <- 0  
    x1_0out <- x[,beta1_EM!=0] 
    x2_0out <- x[,beta2_EM!=0] 
    if (length(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])==1){break}
    Sig1 <- diag(pSCAD_deri(beta1_EM[beta1_EM!=0],gama1,N,a=3.7)/beta1_EM[beta1_EM!=0])
    if (length(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])==1){break}
    Sig2 <- diag(pSCAD_deri(beta2_EM[beta2_EM!=0],gama2,N,a=3.7)/beta2_EM[beta2_EM!=0])
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
  return(list(pi1= pi_EM[,1],pi2= pi_EM[,2],beta1=beta1_EM,beta2=beta2_EM,sigma2_1=sgm21_EM,sigma2_2=sgm22_EM,iter=l))
}

######################################  Parallel Computing for simulation
##### define a function to generate samples from the nonparametric normal mixture regression
rnpnormmix <- function(p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm)
  }
  return(y)
}

#### define a Epanechnikov kernel function with bandwith h
EpaKernel <- function(u,h){
  return(3*(1-(u/h)^2)/4*1*(abs(u)<=h)*(1/h))
}

#### define a function to create autocorrelation matrix
autocorr.mat <- function(p, rho = 0.5) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

pi_1 <- function(x){
  return(0.2+0.6*sin(pi*x)) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x,beta_1){
  return(c(x%*%beta_1))
}

m_2 <- function(x,beta_2){
  return(c(x%*%beta_2))
}

beta_1 <- c(2,0,0,2,1,-2,3,0,0,0,0)
beta_2 <- c(-1,0,0,2,0,-2,3,0,0,1,0)
p <- length(beta_1)

sigma2_1 <- 0.25
sigma2_2 <- 0.64


set.seed(257)
N <- 200
n_sim <- 500  ### the number of simulated data sets
z <- runif(N)
x <- rmvnorm(N,rep(0,p-1),autocorr.mat(p-1,0.5))
xd <- cbind(rep(1,N),x)
xd1_orac <- xd[,beta_1!=0]
xd2_orac <- xd[,beta_2!=0]
pro <- cbind(pi_1(z),pi_2(z))
mu <- cbind(m_1(xd,beta_1),m_2(xd,beta_2))
sgm <- c(sqrt(sigma2_1),sqrt(sigma2_2))

##################################################################################################################
######################################   10 fold cross validation to select the optimal bandwidth based on likelihood
h_grid <- c(0.21,0.23,0.25,0.27,0.29)
k_fold <- 10
#set.seed(248)
cv_spmixreg <- function(j){
  library(mixtools)
  y <- rnpnormmix(pro,mu,sgm)
  
  ind_cv <- sample(x = rep(1:k_fold, each = N / k_fold), size = N)
  cv_lh <- rep(0,length(h_grid))
  
  for(i_h in 1:length(h_grid)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      z_test <- z[ind_cv==i_test] ; z_train <- z[ind_cv!=i_test]
      xd_test <- xd[ind_cv==i_test,] ; xd_train <- xd[ind_cv!=i_test,]
      x_train <- xd_train[,-1]
      
      lmr <- regmixEM(y_train,x_train,lambda=apply(pro,2,mean),beta=cbind(beta_1,beta_2),sigma=sgm)
      N <- length(y_train)
      initl_sgm1 <- rep(sqrt(lmr$sigma[1]),N)
      initl_sgm2 <- rep(sqrt(lmr$sigma[2]),N)
      initl_beta1 <- t(matrix(rep(lmr$beta[,1],N),ncol=N))
      initl_beta2 <- t(matrix(rep(lmr$beta[,2],N),ncol=N))
      initl_pi1 <- rep(lmr$lambda[1],N)
      initl_pi2 <- 1 - initl_pi1
      
      out1 <- semipnormmixEM1(z=z_train,x=xd_train,y=y_train,N,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter=1000,epsilon=0.001,bw=h_grid[i_h])
      
      initl_beta1 <- apply(out1$beta1,2,mean)
      initl_beta2 <- apply(out1$beta2,2,mean)
      initl_sgm21 <- mean(out1$sigma2_1)
      initl_sgm22 <- mean(out1$sigma2_2)
      pi1z <- out1$pi1
      pi2z <- out1$pi2
      
      out2 <- semipnormmixEM2(z=z_train,x=xd_train,y=y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter=1000,epsilon=0.001)
      
      out3 <- semipnormmixEM3(z=z_train,x=xd_train,y_train,N,beta1=out2$beta1,beta2=out2$beta2,sgm21=out2$sigma2_1,sgm22=out2$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=h_grid[i_h])
      
      
      r1 <- out3$pi1*dnorm(y_train,mean=xd_train%*%out2$beta1,sd=sqrt(out2$sigma2_1))/(out3$pi1*dnorm(y_train,mean=xd_train%*%out2$beta1,sd=sqrt(out2$sigma2_1))+out3$pi2*dnorm(y_train,mean=xd_train%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
      r2 <- 1 - r1
      
      pi1_test <- sapply(z_test,function(q) r1%*%EpaKernel(z_train-q,h_grid[i_h])/sum(EpaKernel(z_train-q,h_grid[i_h])))
      pi2_test <- 1 - pi1_test
      
      cv_lh_test <- sum(log(pi1_test*dnorm(y_test,mean=xd_test%*%out2$beta1,sd=sqrt(out2$sigma2_1))+pi2_test*dnorm(y_test,mean=xd_test%*%out2$beta2,sd=sqrt(out2$sigma2_2))))
      
      cv_lh[i_h] <- cv_lh[i_h] + cv_lh_test
      
    }
    
  }
  
  return(list(h_opt = h_grid[cv_lh==max(cv_lh)]))
}

#cv_out <- cv_spmixreg(1)
#h_cv <- mclapply(1:10,cv_spmixreg,mc.cores=16)
the_cluster_cv <- makeCluster(10)
clusterSetRNGStream(the_cluster_cv,35)
clusterExport(the_cluster_cv,c("pLasso_deri","pHard_deri","pSCAD_deri","semipnormmixEM1","semipnormmixEM2","semipnormmixEM3","semipnormmixEM2_L1",
                               "semipnormmixEM2_HARD","semipnormmixEM2_SCAD","semipnormmixEM1_orac","semipnormmixEM2_orac","semipnormmixEM3_orac",
                               "semipnormmixEM2_L1p","semipnormmixEM2_HARDp","semipnormmixEM2_SCADp","semipnormmixEM2_fiL1","semipnormmixEM2_fiHARD","semipnormmixEM2_fiSCAD","rnpnormmix","EpaKernel","autocorr.mat","pi_1",
                               "pi_2","m_1","m_2","beta_1","beta_2","p","sigma2_1","sigma2_2","N","z","x","xd","xd1_orac","xd2_orac","pro","mu",
                               "sgm","h_opt","h_grid","k_fold"))
#h_cv <- parLapply(the_cluster,1:n_sim,sim_spmixreg_VS)
h_cv <- clusterCall(cl = the_cluster_cv, cv_spmixreg,1:10)
stopCluster(the_cluster_cv)

h_opt <- mean(unlist(h_cv))  # the optimal bw is 0.224
print(h_opt) 
##################################################################################################################
gama1_grid =seq(0.2,3,0.2)
gama2_grid =seq(0.2,3,0.2)
sim_spmixreg_VS <- function(i){
  library(mixtools)
  y <- rnpnormmix(pro,mu,sgm)
  
  ## set initial values for the first EM
  lmr <- regmixEM(y,x,lambda=apply(pro,2,mean),beta=cbind(beta_1,beta_2),sigma=sgm)
  
  initl_sgm1 <- rep(sqrt(lmr$sigma[1]),N)
  initl_sgm2 <- rep(sqrt(lmr$sigma[2]),N)
  initl_beta1 <- t(matrix(rep(lmr$beta[,1],N),ncol=N))
  initl_beta2 <- t(matrix(rep(lmr$beta[,2],N),ncol=N))
  initl_pi1 <- rep(lmr$lambda[1],N)
  initl_pi2 <- 1 - initl_pi1
  
  out1 <- semipnormmixEM1(z=z,x=xd,y,N,initl_beta1,initl_beta2,initl_sgm1,initl_sgm2,initl_pi1,initl_pi2,maxiter=1000,epsilon=0.001,bw=h_opt)
  
  ## set initial values for the second EM
  ##initl_beta1 <- apply(out1$beta1,2,mean)
  ##initl_beta2 <- apply(out1$beta2,2,mean)
  ##initl_beta1 <- lmr$beta[,1]
  ##initl_beta2 <- lmr$beta[,2]
  initl_beta1 <- beta_1 + rnorm(p,sd=0.1)
  initl_beta2 <- beta_2 + rnorm(p,sd=0.1)
  initl_sgm21 <- mean(out1$sigma2_1)
  initl_sgm22 <- mean(out1$sigma2_2)
  pi1z <- out1$pi1
  pi2z <- out1$pi2
  out2 <- semipnormmixEM2(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,maxiter=1000,epsilon=0.001)
  out3 <- semipnormmixEM3(z=z,x=xd,y,N=length(y),beta1=out2$beta1,beta2=out2$beta2,sgm21=out2$sigma2_1,sgm22=out2$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=h_opt)
  
  GCV1_L1p <- GCV1_HARDp <- GCV1_SCADp <- rep(0,length(gama1_grid))
  GCV2_L1p <- GCV2_HARDp <- GCV2_SCADp <- rep(0,length(gama2_grid))
  ## compute the GCV criteria for the two components 
  w1 <- out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1))/(out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
  w2 <- out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2))/(out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
  ## the length of gama1_grid is the same as the length of gama2_grid
  for(i_gama in 1:length(gama1_grid)){
    out_L1p_1 <- semipnormmixEM2_L1p(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=out2$beta1,beta_mle=out2$beta2,initl_sgm2_pmle=out2$sigma2_1,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    l_ini <- 0 
    while(out_L1p_1$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_1 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_1
      out_L1p_1 <- semipnormmixEM2_L1p(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta2,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    
    l_ini <- 0
    out_HARDp_1 <- semipnormmixEM2_HARDp(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=out2$beta1,beta_mle=out2$beta2,initl_sgm2_pmle=out2$sigma2_1,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    while(out_HARDp_1$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_1 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_1
      out_HARDp_1 <- semipnormmixEM2_HARDp(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta2,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    
    l_ini <- 0
    out_SCADp_1 <- semipnormmixEM2_SCADp(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=out2$beta1,beta_mle=out2$beta2,initl_sgm2_pmle=out2$sigma2_1,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    while(out_SCADp_1$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_1 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_1
      out_SCADp_1 <- semipnormmixEM2_SCADp(z,xd,y,pi_mle=out3$pi2,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta2,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_2,gama1_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    
    l_ini <- 0
    out_L1p_2 <- semipnormmixEM2_L1p(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=out2$beta2,beta_mle=out2$beta1,initl_sgm2_pmle=out2$sigma2_2,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    while(out_L1p_2$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_2 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_2
      out_L1p_2 <- semipnormmixEM2_L1p(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta1,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    
    l_ini <- 0
    out_HARDp_2 <- semipnormmixEM2_HARDp(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=out2$beta2,beta_mle=out2$beta1,initl_sgm2_pmle=out2$sigma2_2,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    while(out_HARDp_2$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_2 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_2
      out_HARDp_2 <- semipnormmixEM2_HARDp(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta1,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    
    l_ini <- 0
    out_SCADp_2 <- semipnormmixEM2_SCADp(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=out2$beta2,beta_mle=out2$beta1,initl_sgm2_pmle=out2$sigma2_2,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    while(out_SCADp_2$iter==1000 & l_ini<10) {
      l_ini <- l_ini + 1
      print(l_ini)
      initl_beta_pmle_new <- beta_2 + rnorm(p,sd=0.1); initl_sgm2_pmle_new <- sigma2_2
      out_SCADp_2 <- semipnormmixEM2_SCADp(z,xd,y,pi_mle=out3$pi1,initl_beta_pmle=initl_beta_pmle_new,beta_mle=out2$beta1,initl_sgm2_pmle=initl_sgm2_pmle_new,sgm2_mle=out2$sigma2_1,gama2_grid[i_gama],maxiter=1000,epsilon=0.001)
    }
    ## find the second derivatives of log-likelihood
    xd1_L1p_0out <- xd[,out_L1p_1$beta_pmle!=0]; xd2_L1p_0out <- xd[,out_L1p_2$beta_pmle!=0]
    xd1_HARDp_0out <- xd[,out_HARDp_1$beta_pmle!=0]; xd2_HARDp_0out <- xd[,out_HARDp_2$beta_pmle!=0]
    xd1_SCADp_0out <- xd[,out_SCADp_1$beta_pmle!=0]; xd2_SCADp_0out <- xd[,out_SCADp_2$beta_pmle!=0]
    
    p_L1p_1 <- ncol(xd1_L1p_0out); p_L1p_2 <- ncol(xd2_L1p_0out)
    p_HARDp_1 <- ncol(xd1_HARDp_0out); p_HARDp_2 <- ncol(xd2_HARDp_0out)
    p_SCADp_1 <- ncol(xd1_SCADp_0out); p_SCADp_2 <- ncol(xd2_SCADp_0out)
    
    ll_L1p_1 <- matrix(rep(0,p_L1p_1^2),nrow=p_L1p_1);ll_L1p_2 <- matrix(rep(0,p_L1p_2^2),nrow=p_L1p_2)
    ll_HARDp_1 <- matrix(rep(0,p_HARDp_1^2),nrow=p_HARDp_1);ll_HARDp_2 <- matrix(rep(0,p_HARDp_2^2),nrow=p_HARDp_2)
    ll_SCADp_1 <- matrix(rep(0,p_SCADp_1^2),nrow=p_SCADp_1);ll_SCADp_2 <- matrix(rep(0,p_SCADp_2^2),nrow=p_SCADp_2)
    
    c1_L1p_1 <- -out3$pi1*dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle))^2*(y-c(xd%*%out_L1p_1$beta_pmle))^2/(out_L1p_1$sigma2_pmle*(out3$pi1*dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))^2)
    c2_L1p_1 <- (dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle))*(y-c(xd%*%out_L1p_1$beta_pmle))^2/out_L1p_1$sigma2_pmle-dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle)))/(out3$pi1*dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
    c1_L1p_2 <- -out3$pi2*dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle))^2*(y-c(xd%*%out_L1p_2$beta_pmle))^2/(out_L1p_2$sigma2_pmle*(out3$pi2*dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))^2)
    c2_L1p_2 <- (dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle))*(y-c(xd%*%out_L1p_2$beta_pmle))^2/out_L1p_2$sigma2_pmle-dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle)))/(out3$pi2*dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))
    
    c1_HARDp_1 <- -out3$pi1*dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle))^2*(y-c(xd%*%out_HARDp_1$beta_pmle))^2/(out_HARDp_1$sigma2_pmle*(out3$pi1*dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))^2)
    c2_HARDp_1 <- (dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle))*(y-c(xd%*%out_HARDp_1$beta_pmle))^2/out_HARDp_1$sigma2_pmle-dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle)))/(out3$pi1*dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
    c1_HARDp_2 <- -out3$pi2*dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle))^2*(y-c(xd%*%out_HARDp_2$beta_pmle))^2/(out_HARDp_2$sigma2_pmle*(out3$pi2*dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))^2)
    c2_HARDp_2 <- (dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle))*(y-c(xd%*%out_HARDp_2$beta_pmle))^2/out_HARDp_2$sigma2_pmle-dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle)))/(out3$pi2*dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))
    
    c1_SCADp_1 <- -out3$pi1*dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle))^2*(y-c(xd%*%out_SCADp_1$beta_pmle))^2/(out_SCADp_1$sigma2_pmle*(out3$pi1*dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))^2)
    c2_SCADp_1 <- (dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle))*(y-c(xd%*%out_SCADp_1$beta_pmle))^2/out_SCADp_1$sigma2_pmle-dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle)))/(out3$pi1*dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle))+out3$pi2*dnorm(y,mean=xd%*%out2$beta2,sd=sqrt(out2$sigma2_2)))
    c1_SCADp_2 <- -out3$pi2*dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle))^2*(y-c(xd%*%out_SCADp_2$beta_pmle))^2/(out_SCADp_2$sigma2_pmle*(out3$pi2*dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))^2)
    c2_SCADp_2 <- (dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle))*(y-c(xd%*%out_SCADp_2$beta_pmle))^2/out_SCADp_2$sigma2_pmle-dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle)))/(out3$pi2*dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle))+out3$pi1*dnorm(y,mean=xd%*%out2$beta1,sd=sqrt(out2$sigma2_1)))
    
    for(i_n in 1:N){
      
      ll_L1p_1 <- ll_L1p_1 + out3$pi1[i_n]/out_L1p_1$sigma2_pmle*xd1_L1p_0out[i_n,]%*%t((c1_L1p_1[i_n]+c2_L1p_1[i_n])*xd1_L1p_0out[i_n,])
      ll_L1p_2 <- ll_L1p_2 + out3$pi2[i_n]/out_L1p_2$sigma2_pmle*xd2_L1p_0out[i_n,]%*%t((c1_L1p_2[i_n]+c2_L1p_2[i_n])*xd2_L1p_0out[i_n,])
      ll_HARDp_1 <- ll_HARDp_1 + out3$pi1[i_n]/out_HARDp_1$sigma2_pmle*xd1_HARDp_0out[i_n,]%*%t((c1_HARDp_1[i_n]+c2_HARDp_1[i_n])*xd1_HARDp_0out[i_n,])
      ll_HARDp_2 <- ll_HARDp_2 + out3$pi2[i_n]/out_HARDp_2$sigma2_pmle*xd2_HARDp_0out[i_n,]%*%t((c1_HARDp_2[i_n]+c2_HARDp_2[i_n])*xd2_HARDp_0out[i_n,])
      ll_SCADp_1 <- ll_SCADp_1 + out3$pi1[i_n]/out_SCADp_1$sigma2_pmle*xd1_SCADp_0out[i_n,]%*%t((c1_SCADp_1[i_n]+c2_SCADp_1[i_n])*xd1_SCADp_0out[i_n,])
      ll_SCADp_2 <- ll_SCADp_2 + out3$pi2[i_n]/out_SCADp_2$sigma2_pmle*xd2_SCADp_0out[i_n,]%*%t((c1_SCADp_2[i_n]+c2_SCADp_2[i_n])*xd2_SCADp_0out[i_n,])
      
    }
    
    D1_L1p <- sum(w1*(dnorm(y,mean=y,sd=sqrt(out_L1p_1$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_L1p_1$beta_pmle,sd=sqrt(out_L1p_1$sigma2_pmle),log=TRUE)))
    D2_L1p <- sum(w2*(dnorm(y,mean=y,sd=sqrt(out_L1p_2$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_L1p_2$beta_pmle,sd=sqrt(out_L1p_2$sigma2_pmle),log=TRUE)))
    D1_HARDp <- sum(w1*(dnorm(y,mean=y,sd=sqrt(out_HARDp_1$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_HARDp_1$beta_pmle,sd=sqrt(out_HARDp_1$sigma2_pmle),log=TRUE)))
    D2_HARDp <- sum(w2*(dnorm(y,mean=y,sd=sqrt(out_HARDp_2$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_HARDp_2$beta_pmle,sd=sqrt(out_HARDp_2$sigma2_pmle),log=TRUE)))
    D1_SCADp <- sum(w1*(dnorm(y,mean=y,sd=sqrt(out_SCADp_1$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_SCADp_1$beta_pmle,sd=sqrt(out_SCADp_1$sigma2_pmle),log=TRUE)))
    D2_SCADp <- sum(w2*(dnorm(y,mean=y,sd=sqrt(out_SCADp_2$sigma2_pmle),log=TRUE)-dnorm(y,mean=xd%*%out_SCADp_2$beta_pmle,sd=sqrt(out_SCADp_2$sigma2_pmle),log=TRUE)))
    
    e1_L1p <- sum(diag(solve(ll_L1p_1-diag(pLasso_deri(out_L1p_1$beta_pmle[out_L1p_1$beta_pmle!=0],gama1_grid[i_gama],N)/out_L1p_1$beta_pmle[out_L1p_1$beta_pmle!=0]))%*%ll_L1p_1))
    e2_L1p <- sum(diag(solve(ll_L1p_2-diag(pLasso_deri(out_L1p_2$beta_pmle[out_L1p_2$beta_pmle!=0],gama2_grid[i_gama],N)/out_L1p_2$beta_pmle[out_L1p_2$beta_pmle!=0]))%*%ll_L1p_2))
    e1_HARDp <- sum(diag(solve(ll_HARDp_1-diag(pHard_deri(out_HARDp_1$beta_pmle[out_HARDp_1$beta_pmle!=0],gama1_grid[i_gama],N)/out_HARDp_1$beta_pmle[out_HARDp_1$beta_pmle!=0]))%*%ll_HARDp_1))
    e2_HARDp <- sum(diag(solve(ll_HARDp_2-diag(pHard_deri(out_HARDp_2$beta_pmle[out_HARDp_2$beta_pmle!=0],gama2_grid[i_gama],N)/out_HARDp_2$beta_pmle[out_HARDp_2$beta_pmle!=0]))%*%ll_HARDp_2))
    e1_SCADp <- sum(diag(solve(ll_SCADp_1-diag(pSCAD_deri(out_SCADp_1$beta_pmle[out_SCADp_1$beta_pmle!=0],gama1_grid[i_gama],N,a=3.7)/out_SCADp_1$beta_pmle[out_SCADp_1$beta_pmle!=0]))%*%ll_SCADp_1))
    e2_SCADp <- sum(diag(solve(ll_SCADp_2-diag(pSCAD_deri(out_SCADp_2$beta_pmle[out_SCADp_2$beta_pmle!=0],gama2_grid[i_gama],N,a=3.7)/out_SCADp_2$beta_pmle[out_SCADp_2$beta_pmle!=0]))%*%ll_SCADp_2))
    
    if (out_L1p_1$iter[1]==1000) {GCV1_L1p[i_gama] <- 1000} ## if the algorithm does not converge for 1000 iterations, then discard the results
    else {GCV1_L1p[i_gama] <- D1_L1p/(N*(1-e1_L1p/N)^2)}
    if (out_L1p_2$iter[1]==1000) {GCV2_L1p[i_gama] <- 1000}
    else {GCV2_L1p[i_gama] <- D2_L1p/(N*(1-e2_L1p/N)^2)}
    if (out_HARDp_1$iter[1]==1000) {GCV1_HARDp[i_gama] <- 1000}
    else {GCV1_HARDp[i_gama] <- D1_HARDp/(N*(1-e1_HARDp/N)^2)}
    if (out_HARDp_2$iter[1]==1000) {GCV2_HARDp[i_gama] <- 1000}
    else {GCV2_HARDp[i_gama] <- D2_HARDp/(N*(1-e2_HARDp/N)^2)}
    if (out_SCADp_1$iter[1]==1000) {GCV1_SCADp[i_gama] <- 1000}
    else {GCV1_SCADp[i_gama] <- D1_SCADp/(N*(1-e1_SCADp/N)^2)}
    if (out_SCADp_2$iter[1]==1000) {GCV2_SCADp[i_gama] <- 1000}
    else {GCV2_SCADp[i_gama] <- D2_SCADp/(N*(1-e2_SCADp/N)^2)}
    
  }
  optgama1_L1 <- gama1_grid[GCV1_L1p==min(GCV1_L1p)]
  optgama2_L1 <- gama2_grid[GCV2_L1p==min(GCV2_L1p)]
  optgama1_HARD <- gama1_grid[GCV1_HARDp==min(GCV1_HARDp)]
  optgama2_HARD <- gama2_grid[GCV2_HARDp==min(GCV2_HARDp)]
  optgama1_SCAD <- gama1_grid[GCV1_SCADp==min(GCV1_SCADp)]
  optgama2_SCAD <- gama2_grid[GCV2_SCADp==min(GCV2_SCADp)]
  
  out_acL1 <- semipnormmixEM_L1(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  l_ini <- 0
  while(out_acL1$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out_acL1 <- semipnormmixEM_L1(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  }
  
  out_acHARD <- semipnormmixEM_HARD(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  l_ini <- 0
  while(out_acHARD$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out_acHARD <- semipnormmixEM_HARD(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  }  
  
  out_acSCAD <- semipnormmixEM_SCAD(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  l_ini <- 0
  while(out_acSCAD$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out2_acSCAD <- semipnormmixEM_SCAD(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001,bw=h_opt)
  }
  
  initl_fbeta1 <- initl_beta1
  initl_fbeta2 <- initl_beta2
  out2_L1 <- semipnormmixEM2_L1(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001)
  l_ini <- 0
  while(out2_L1$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out2_L1 <- semipnormmixEM2_L1(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.001)
    if(out2_L1$iter < 1000) {initl_fbeta1 <- initl_beta1_new; initl_fbeta2 <- initl_beta2_new}
  }
  out2_fL1 <- semipnormmixEM2_fiL1(z=z,x=xd,y,N,initl_fbeta1,initl_fbeta2,initl_sgm21,initl_sgm22,pi1z,pi2z,gama1=optgama1_L1,gama2=optgama2_L1,maxiter=1000,epsilon=0.0001,bw=h_opt)
  
  initl_fbeta1 <- initl_beta1
  initl_fbeta2 <- initl_beta2
  out2_HARD <- semipnormmixEM2_HARD(z=z,x=xd,y,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=optgama1_HARD,gama2=optgama2_HARD,maxiter=1000,epsilon=0.001)
  l_ini <- 0
  while(out2_HARD$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out2_HARD <- semipnormmixEM2_HARD(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_HARD,gama2=optgama2_HARD,maxiter=1000,epsilon=0.001)
    if(out2_HARD$iter < 1000) {initl_fbeta1 <- initl_beta1_new; initl_fbeta2 <- initl_beta2_new}
  }
  out2_fHARD <- semipnormmixEM2_fiHARD(z=z,x=xd,y,N,initl_fbeta1,initl_fbeta2,initl_sgm21,initl_sgm22,pi1z,pi2z,gama1=optgama1_HARD,gama2=optgama2_HARD,maxiter=1000,epsilon=0.0001,bw=h_opt)
  
  initl_fbeta1 <- initl_beta1
  initl_fbeta2 <- initl_beta2
  out2_SCAD <- semipnormmixEM2_SCAD(z=z,x=xd,y,pi1z,pi2z,initl_fbeta1,initl_fbeta2,initl_sgm21,initl_sgm22,gama1=optgama1_SCAD,gama2=optgama2_SCAD,maxiter=1000,epsilon=0.001)
  l_ini <- 0
  while(out2_SCAD$iter==1000 & l_ini < 10) {
    l_ini <- l_ini + 1
    initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
    initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
    out2_SCAD <- semipnormmixEM2_SCAD(z=z,x=xd,y,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=optgama1_SCAD,gama2=optgama2_SCAD,maxiter=1000,epsilon=0.001)
    if(out2_SCAD$iter < 1000) {initl_fbeta1 <- initl_beta1_new; initl_fbeta2 <- initl_beta2_new}
  }
  out2_fSCAD <- semipnormmixEM2_fiSCAD(z=z,x=xd,y,N,initl_fbeta1,initl_fbeta2,initl_sgm21,initl_sgm22,pi1z,pi2z,gama1=optgama1_SCAD,gama2=optgama2_SCAD,maxiter=1000,epsilon=0.0001,bw=h_opt)
  
  
  #######################################################################################################################################################
  out3_L1 <- semipnormmixEM3(z,xd,y,N,beta1=out2_L1$beta1,beta2=out2_L1$beta2,sgm21=out2_L1$sigma2_1,sgm22=out2_L1$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=h_opt)
  out3_HARD <- semipnormmixEM3(z,xd,y,N,beta1=out2_HARD$beta1,beta2=out2_HARD$beta2,sgm21=out2_HARD$sigma2_1,sgm22=out2_HARD$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=h_opt)
  out3_SCAD <- semipnormmixEM3(z,xd,y,N,beta1=out2_SCAD$beta1,beta2=out2_SCAD$beta2,sgm21=out2_SCAD$sigma2_1,sgm22=out2_SCAD$sigma2_2,initl_pi1=pi1z,initl_pi2=pi2z,maxiter=1000,epsilon=0.001,bw=h_opt)
  
  ######################################################################################################################################################
  ## set initial values for oracle set up
  initl_sgm1_o <- rep(sqrt(sigma2_1),N)
  initl_sgm2_o <- rep(sqrt(sigma2_2),N)
  initl_beta1_o <- t(matrix(rep(beta_1[beta_1!=0],N),ncol=N))
  initl_beta2_o <- t(matrix(rep(beta_2[beta_2!=0],N),ncol=N))
  initl_pi1_o <- rep(mean(pi_1(z)),N)
  initl_pi2_o <- 1 - initl_pi1
  
  out1_orac <- semipnormmixEM1_orac(z=z,x1=xd1_orac,x2=xd2_orac,y,N,initl_beta1_o,initl_beta2_o,initl_sgm1_o,initl_sgm2_o,initl_pi1_o,initl_pi2_o,maxiter=1000,epsilon=0.0001,bw=0.14)
  ## set initial values for the second EM for oracle set up
  initl_beta1_o <- apply(out1_orac$beta1,2,mean)
  initl_beta2_o <- apply(out1_orac$beta2,2,mean)
  initl_sgm21_o <- mean(out1_orac$sigma2_1)
  initl_sgm22_o <- mean(out1_orac$sigma2_2)
  pi1z_o <- out1_orac$pi1
  pi2z_o <- out1_orac$pi2
  out2_orac <- semipnormmixEM2_orac(z=z,x1=xd1_orac,x2=xd2_orac,y,pi1z_o,pi2z_o,initl_beta1_o,initl_beta2_o,initl_sgm21_o,initl_sgm22_o,maxiter=1000,epsilon=0.0001)
  out3_orac <- semipnormmixEM3_orac(z=z,x1=xd1_orac,x2=xd2_orac,y,N=length(y),beta1=out2_orac$beta1,beta2=out2_orac$beta2,sgm21=out2_orac$sigma2_1,sgm22=out2_orac$sigma2_2,initl_pi1=pi1z_o,initl_pi2=pi2z_o,maxiter=1000,epsilon=0.0001,bw=0.15)
  
  ##### calculate the errors 
  SE_beta1_orac <- sum((beta_1[beta_1!=0]-out2_orac$beta1)^2) ;SE_beta2_orac <- sum((beta_2[beta_2!=0]-out2_orac$beta2)^2) 
  SE_beta1_mle <- sum((beta_1-out2$beta1)^2) ; SE_beta2_mle <- sum((beta_2-out2$beta2)^2)
  SE_beta1_L1 <- sum((beta_1-out2_L1$beta1)^2) ; SE_beta2_L1 <- sum((beta_2-out2_L1$beta2)^2)
  SE_beta1_HARD <- sum((beta_1-out2_HARD$beta1)^2) ; SE_beta2_HARD <- sum((beta_2-out2_HARD$beta2)^2)
  SE_beta1_SCAD <- sum((beta_1-out2_SCAD$beta1)^2) ; SE_beta2_SCAD <- sum((beta_2-out2_SCAD$beta2)^2)
  SE_beta1_fL1 <- sum((beta_1-out2_fL1$beta1)^2) ; SE_beta2_fL1 <- sum((beta_2-out2_fL1$beta2)^2)
  SE_beta1_fHARD <- sum((beta_1-out2_fHARD$beta1)^2) ; SE_beta2_fHARD <- sum((beta_2-out2_fHARD$beta2)^2)
  SE_beta1_fSCAD <- sum((beta_1-out2_fSCAD$beta1)^2) ; SE_beta2_fSCAD <- sum((beta_2-out2_fSCAD$beta2)^2)
  SE_beta1_acL1 <- sum((beta_1-out_acL1$beta1)^2) ; SE_beta2_acL1 <- sum((beta_2-out_acL1$beta2)^2)
  SE_beta1_acHARD <- sum((beta_1-out_acHARD$beta1)^2) ; SE_beta2_acHARD <- sum((beta_2-out_acHARD$beta2)^2)
  SE_beta1_acSCAD <- sum((beta_1-out_acSCAD$beta1)^2) ; SE_beta2_acSCAD <- sum((beta_2-out_acSCAD$beta2)^2)
  
  SE_sigma2_1_orac <- (sigma2_1-out2_orac$sigma2_1)^2 ; SE_sigma2_2_orac <- (sigma2_2-out2_orac$sigma2_2)^2 
  SE_sigma2_1_mle <- (sigma2_1-out2$sigma2_1)^2 ; SE_sigma2_2_mle <- (sigma2_2-out2$sigma2_2)^2
  SE_sigma2_1_L1 <- (sigma2_1-out2_L1$sigma2_1)^2 ; SE_sigma2_2_L1 <- (sigma2_2-out2_L1$sigma2_2)^2 
  SE_sigma2_1_HARD <- (sigma2_1-out2_HARD$sigma2_1)^2 ; SE_sigma2_2_HARD <- (sigma2_2-out2_HARD$sigma2_2)^2 
  SE_sigma2_1_SCAD <- (sigma2_1-out2_SCAD$sigma2_1)^2 ; SE_sigma2_2_SCAD <- (sigma2_2-out2_SCAD$sigma2_2)^2
  SE_sigma2_1_fL1 <- (sigma2_1-out2_fL1$sigma2_1)^2 ; SE_sigma2_2_fL1 <- (sigma2_2-out2_fL1$sigma2_2)^2 
  SE_sigma2_1_fHARD <- (sigma2_1-out2_fHARD$sigma2_1)^2 ; SE_sigma2_2_fHARD <- (sigma2_2-out2_fHARD$sigma2_2)^2 
  SE_sigma2_1_fSCAD <- (sigma2_1-out2_fSCAD$sigma2_1)^2 ; SE_sigma2_2_fSCAD <- (sigma2_2-out2_fSCAD$sigma2_2)^2
  SE_sigma2_1_acL1 <- (sigma2_1-out_acL1$sigma2_1)^2 ; SE_sigma2_2_acL1 <- (sigma2_2-out_acL1$sigma2_2)^2 
  SE_sigma2_1_acHARD <- (sigma2_1-out_acHARD$sigma2_1)^2 ; SE_sigma2_2_acHARD <- (sigma2_2-out_acHARD$sigma2_2)^2 
  SE_sigma2_1_acSCAD <- (sigma2_1-out_acSCAD$sigma2_1)^2 ; SE_sigma2_2_acSCAD <- (sigma2_2-out_acSCAD$sigma2_2)^2
  
  RASE_orac <- sqrt(sum((out3_orac$pi1-pi_1(z))^2)/N)
  RASE_mle <- sqrt(sum((out3$pi1-pi_1(z))^2)/N)
  RASE_L1 <- sqrt(sum((out3_L1$pi1-pi_1(z))^2)/N)
  RASE_HARD <- sqrt(sum((out3_HARD$pi1-pi_1(z))^2)/N)
  RASE_SCAD <- sqrt(sum((out3_SCAD$pi1-pi_1(z))^2)/N)
  RASE_fL1 <- sqrt(sum((out2_fL1$pi1-pi_1(z))^2)/N)
  RASE_fHARD <- sqrt(sum((out2_fHARD$pi1-pi_1(z))^2)/N)
  RASE_fSCAD <- sqrt(sum((out2_fSCAD$pi1-pi_1(z))^2)/N)
  RASE_acL1 <- sqrt(sum((out_acL1$pi1-pi_1(z))^2)/N)
  RASE_acHARD <- sqrt(sum((out_acHARD$pi1-pi_1(z))^2)/N)
  RASE_acSCAD <- sqrt(sum((out_acSCAD$pi1-pi_1(z))^2)/N)
  
  SE <- cbind(c(SE_beta1_orac,SE_beta1_mle,SE_beta1_L1,SE_beta1_HARD,SE_beta1_SCAD,SE_beta1_fL1,SE_beta1_fHARD,SE_beta1_fSCAD,SE_beta1_acL1,SE_beta1_acHARD,SE_beta1_acSCAD),c(SE_beta2_orac,SE_beta2_mle,SE_beta2_L1,SE_beta2_HARD,SE_beta2_SCAD,SE_beta2_fL1,SE_beta2_fHARD,SE_beta2_fSCAD,SE_beta2_acL1,SE_beta2_acHARD,SE_beta2_acSCAD),
              c(SE_sigma2_1_orac,SE_sigma2_1_mle,SE_sigma2_1_L1,SE_sigma2_1_HARD,SE_sigma2_1_SCAD,SE_sigma2_1_fL1,SE_sigma2_1_fHARD,SE_sigma2_1_fSCAD,SE_sigma2_1_acL1,SE_sigma2_1_acHARD,SE_sigma2_1_acSCAD),c(SE_sigma2_2_orac,SE_sigma2_2_mle,SE_sigma2_2_L1,SE_sigma2_2_HARD,SE_sigma2_2_SCAD,
                                                                                                                                                         SE_sigma2_2_fL1,SE_sigma2_2_fHARD,SE_sigma2_2_fSCAD,SE_sigma2_2_acL1,SE_sigma2_2_acHARD,SE_sigma2_2_acSCAD))
  RASE <- c(RASE_orac,RASE_mle,RASE_L1,RASE_HARD,RASE_SCAD,RASE_fL1,RASE_fHARD,RASE_fSCAD,RASE_acL1,RASE_acHARD,RASE_acSCAD)
  
  return(list(beta1_mle=out2$beta1,beta2_mle=out2$beta2,beta1_L1=out2_L1$beta1,beta2_L1=out2_L1$beta2,beta1_HARD=out2_HARD$beta1,beta2_HARD=out2_HARD$beta2,
              beta1_SCAD=out2_SCAD$beta1,beta2_SCAD=out2_SCAD$beta2,beta1_orac=out2_orac$beta1,beta2_orac=out2_orac$beta2,beta1_fL1=out2_fL1$beta1,beta2_fL1=out2_fL1$beta2,beta1_fHARD=out2_fHARD$beta1,beta2_fHARD=out2_fHARD$beta2,
              beta1_fSCAD=out2_fSCAD$beta1,beta2_fSCAD=out2_fSCAD$beta2,beta1_acL1=out_acL1$beta1,beta2_acL1=out_acL1$beta2,beta1_acHARD=out_acHARD$beta1,beta2_acHARD=out_acHARD$beta2,beta1_acSCAD=out_acSCAD$beta1,beta2_acSCAD=out_acSCAD$beta2,sigma2_1_mle=out2$sigma2_1,sigma2_2_mle=out2$sigma2_2,sigma2_1_L1=out2_L1$sigma2_1,sigma2_2_L1=out2_L1$sigma2_2,sigma2_1_HARD=out2_HARD$sigma2_1,sigma2_2_HARD=out2_HARD$sigma2_2,
              sigma2_1_SCAD=out2_SCAD$sigma2_1,sigma2_2_SCAD=out2_SCAD$sigma2_2,sigma2_1_orac=out2_orac$sigma2_1,sigma2_2_orac=out2_orac$sigma2_2,sigma2_1_fL1=out2_fL1$sigma2_1,sigma2_2_fL1=out2_fL1$sigma2_2,sigma2_1_fHARD=out2_fHARD$sigma2_1,sigma2_2_fHARD=out2_fHARD$sigma2_2,
              sigma2_1_fSCAD=out2_fSCAD$sigma2_1,sigma2_2_fSCAD=out2_fSCAD$sigma2_2,sigma2_1_acL1=out_acL1$sigma2_1,sigma2_2_acL1=out_acL1$sigma2_2,sigma2_1_acHARD=out_acHARD$sigma2_1,sigma2_2_acHARD=out_acHARD$sigma2_2,sigma2_1_acSCAD=out_acSCAD$sigma2_1,sigma2_2_acSCAD=out_acSCAD$sigma2_2,pi1_mle=out3$pi1,pi1_L1=out3_L1$pi1,pi1_HARD=out3_HARD$pi1,pi1_SCAD=out3_SCAD$pi1,pi1_orac=out3_orac$pi1,pi1_fL1=out2_fL1$pi1,pi1_fHARD=out2_fHARD$pi1,pi1_fSCAD=out2_fSCAD$pi1,pi1_acL1=out_acL1$pi1,
              pi1_acHARD=out_acHARD$pi1,pi1_acSCAD=out_acSCAD$pi1,SE=SE,RASE=RASE))
  
}
#sim1 <- sim_spmixreg_VS(1)

the_cluster <- makeCluster(15)
clusterSetRNGStream(the_cluster,27)
clusterExport(the_cluster,c("pLasso_deri","pHard_deri","pSCAD_deri","semipnormmixEM1","semipnormmixEM2","semipnormmixEM3","semipnormmixEM2_L1",
                            "semipnormmixEM2_HARD","semipnormmixEM2_SCAD","semipnormmixEM1_orac","semipnormmixEM2_orac","semipnormmixEM3_orac",
                            "semipnormmixEM2_L1p","semipnormmixEM2_HARDp","semipnormmixEM2_SCADp","semipnormmixEM2_fiL1","semipnormmixEM2_fiHARD","semipnormmixEM2_fiSCAD","rnpnormmix","EpaKernel","autocorr.mat","pi_1",
                            "pi_2","m_1","m_2","beta_1","beta_2","p","sigma2_1","sigma2_2","N","z","x","xd","xd1_orac","xd2_orac","pro","mu",
                            "sgm","h_opt","gama1_grid","gama2_grid","semipnormmixEM_L1","semipnormmixEM_HARD","semipnormmixEM_SCAD"))
test_sim <- parLapply(the_cluster,1:n_sim,sim_spmixreg_VS)
#test_sim <- clusterCall(cl = the_cluster, sim_spmixreg_VS,1:3)
stopCluster(the_cluster)

#set.seed(579)
##test_sim <- mclapply(1:n_sim,sim_spmixreg_VS,mc.cores=8)
beta1s_L1 <- beta2s_L1 <- beta1s_HARD <- beta2s_HARD <- beta1s_SCAD <- beta2s_SCAD <- matrix(rep(0,n_sim*ncol(x)),ncol=n_sim)
beta1s_fL1 <- beta2s_fL1 <- beta1s_fHARD <- beta2s_fHARD <- beta1s_fSCAD <- beta2s_fSCAD <- matrix(rep(0,n_sim*ncol(x)),ncol=n_sim)
beta1s_acL1 <- beta2s_acL1 <- beta1s_acHARD <- beta2s_acHARD <- beta1s_acSCAD <- beta2s_acSCAD <- matrix(rep(0,n_sim*ncol(x)),ncol=n_sim)
SEs <- matrix(rep(0,11*4),nrow=11) ; RASEs <- matrix(rep(0,n_sim*11),ncol=n_sim)

for(i_sim in 1:n_sim){
  
  beta1s_L1[,i_sim] <- test_sim[[i_sim]]$beta1_L1[-1]
  beta2s_L1[,i_sim] <- test_sim[[i_sim]]$beta2_L1[-1]
  beta1s_HARD[,i_sim] <- test_sim[[i_sim]]$beta1_HARD[-1]
  beta2s_HARD[,i_sim] <- test_sim[[i_sim]]$beta2_HARD[-1]
  beta1s_SCAD[,i_sim] <- test_sim[[i_sim]]$beta1_SCAD[-1]
  beta2s_SCAD[,i_sim] <- test_sim[[i_sim]]$beta2_SCAD[-1]
  beta1s_fL1[,i_sim] <- test_sim[[i_sim]]$beta1_fL1[-1]
  beta2s_fL1[,i_sim] <- test_sim[[i_sim]]$beta2_fL1[-1]
  beta1s_fHARD[,i_sim] <- test_sim[[i_sim]]$beta1_fHARD[-1]
  beta2s_fHARD[,i_sim] <- test_sim[[i_sim]]$beta2_fHARD[-1]
  beta1s_fSCAD[,i_sim] <- test_sim[[i_sim]]$beta1_fSCAD[-1]
  beta2s_fSCAD[,i_sim] <- test_sim[[i_sim]]$beta2_fSCAD[-1]
  beta1s_acL1[,i_sim] <- test_sim[[i_sim]]$beta1_acL1[-1]
  beta2s_acL1[,i_sim] <- test_sim[[i_sim]]$beta2_acL1[-1]
  beta1s_acHARD[,i_sim] <- test_sim[[i_sim]]$beta1_acHARD[-1]
  beta2s_acHARD[,i_sim] <- test_sim[[i_sim]]$beta2_acHARD[-1]
  beta1s_acSCAD[,i_sim] <- test_sim[[i_sim]]$beta1_acSCAD[-1]
  beta2s_acSCAD[,i_sim] <- test_sim[[i_sim]]$beta2_acSCAD[-1]
  
  SEs <- SEs + test_sim[[i_sim]]$SE
  RASEs[,i_sim] <- test_sim[[i_sim]]$RASE
}

### find sensitivity and specificity 

sen_L1_1 <- sum(beta1s_L1[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)   
spe_L1_1 <- sum(beta1s_L1[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_L1_2 <- sum(beta2s_L1[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_L1_2 <- sum(beta2s_L1[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim)  

sen_HARD_1 <- sum(beta1s_HARD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_HARD_1 <- sum(beta1s_HARD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim) 
sen_HARD_2 <- sum(beta2s_HARD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_HARD_2 <- sum(beta2s_HARD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

sen_SCAD_1 <- sum(beta1s_SCAD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_SCAD_1 <- sum(beta1s_SCAD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_SCAD_2 <- sum(beta2s_SCAD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_SCAD_2 <- sum(beta2s_SCAD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

sen_fL1_1 <- sum(beta1s_fL1[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)   
spe_fL1_1 <- sum(beta1s_fL1[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_fL1_2 <- sum(beta2s_fL1[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_fL1_2 <- sum(beta2s_fL1[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim)  

sen_fHARD_1 <- sum(beta1s_fHARD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_fHARD_1 <- sum(beta1s_fHARD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim) 
sen_fHARD_2 <- sum(beta2s_fHARD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_fHARD_2 <- sum(beta2s_fHARD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

sen_fSCAD_1 <- sum(beta1s_acSCAD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_fSCAD_1 <- sum(beta1s_acSCAD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_fSCAD_2 <- sum(beta2s_acSCAD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_fSCAD_2 <- sum(beta2s_acSCAD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

sen_acL1_1 <- sum(beta1s_acL1[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)   
spe_acL1_1 <- sum(beta1s_acL1[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_acL1_2 <- sum(beta2s_acL1[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_acL1_2 <- sum(beta2s_acL1[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim)  

sen_acHARD_1 <- sum(beta1s_acHARD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_acHARD_1 <- sum(beta1s_acHARD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim) 
sen_acHARD_2 <- sum(beta2s_acHARD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_acHARD_2 <- sum(beta2s_acHARD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

sen_acSCAD_1 <- sum(beta1s_acSCAD[beta_1[-1]!=0,]!=0)/((length(beta_1[beta_1!=0])-1)*n_sim)
spe_acSCAD_1 <- sum(beta1s_acSCAD[beta_1[-1]==0,]==0)/(length(beta_1[beta_1==0])*n_sim)  
sen_acSCAD_2 <- sum(beta2s_acSCAD[beta_2[-1]!=0,]!=0)/((length(beta_2[beta_2!=0])-1)*n_sim)
spe_acSCAD_2 <- sum(beta2s_acSCAD[beta_2[-1]==0,]==0)/(length(beta_2[beta_2==0])*n_sim) 

result_sen_spe <- cbind(c(sen_L1_1,sen_HARD_1,sen_SCAD_1),c(sen_L1_2,sen_HARD_2,sen_SCAD_2),c(spe_L1_1,spe_HARD_1,spe_SCAD_1),c(spe_L1_2,spe_HARD_2,spe_SCAD_2))
resultf_sen_spe <- cbind(c(sen_fL1_1,sen_fHARD_1,sen_fSCAD_1),c(sen_fL1_2,sen_fHARD_2,sen_fSCAD_2),c(spe_fL1_1,spe_fHARD_1,spe_fSCAD_1),c(spe_fL1_2,spe_fHARD_2,spe_fSCAD_2))
resultac_sen_spe <- cbind(c(sen_acL1_1,sen_acHARD_1,sen_acSCAD_1),c(sen_acL1_2,sen_acHARD_2,sen_acSCAD_2),c(spe_acL1_1,spe_acHARD_1,spe_acSCAD_1),c(spe_acL1_2,spe_acHARD_2,spe_acSCAD_2))
print(result_sen_spe)
print(resultf_sen_spe) 
print(resultac_sen_spe) 
### find MSEs and RASEs
MSEs <- sqrt(SEs/n_sim)
print(MSEs)
apply(RASEs,1,mean)
apply(RASEs,1,sd)




















