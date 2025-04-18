##########################################################################################################
##########################################################################################################
############## Simulation for variable selection (p=50, n=200, EX2 independent covariates)
##########################################################################################################
##########################################################################################################

library(mixtools)
library(Hmisc)
library(graphics)
#library(mixreg)  
library(parallel)
#install.packages("devtools")
#devtools::install_github("thiyangt/fformpp")
library(fformpp)

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

beta_1 <- c(2,0,3,2,0,-2,3,0,4,0,-3,rep(0,90))
beta_2 <- c(-1,2,-2,3,0,-3,2,0,2,3,0,rep(0,90))
p <- length(beta_1)

sigma2_1 <- 0.25
sigma2_2 <- 0.49

set.seed(257)
N <- 200
n_sim <- 500  ### the number of simulated data sets
z <- runif(N)
x <- rmvnorm(N,rep(0,p-1),autocorr.mat(p-1,0))
xd <- cbind(rep(1,N),x)
xd1_orac <- xd[,beta_1!=0]
xd2_orac <- xd[,beta_2!=0]
pro <- cbind(pi_1(z),pi_2(z))
mu <- cbind(m_1(xd,beta_1),m_2(xd,beta_2))
sgm <- c(sqrt(sigma2_1),sqrt(sigma2_2))

gama1_grid <- seq(0.2,1.2,0.2)
gama2_grid <- seq(0.2,1.2,0.2)
h_grid <- c(0.14,0.16,0.18,0.20,0.22)
max_tpar <- expand.grid(gama1_grid,gama2_grid,h_grid)  ## the matrix for all combinations of tuning parameters

#### conduct k-fold cross validation for a simulated dataset
set.seed(978)
k_fold <- 10

initl_beta1 <- beta_1 + rnorm(p,sd=0.3)
initl_beta2 <- beta_2 + rnorm(p,sd=0.3)
initl_sgm21 <- sigma2_1
initl_sgm22 <- sigma2_2

cv_spmixreg <- function(j){
  library(mixtools)
  library(fformpp)
  y <- rnpnormmix(pro,mu,sgm)
  
  ind_cv <- sample(x = rep(1:k_fold, each = N / k_fold), size = N) 
  ### set the criteria for cross validation
  cv_error_L1 <- matrix(rep(0,nrow(max_tpar)*k_fold),ncol=k_fold)   
  cv_error_HARD <- matrix(rep(0,nrow(max_tpar)*k_fold),ncol=k_fold) 
  cv_error_SCAD <- matrix(rep(0,nrow(max_tpar)*k_fold),ncol=k_fold) 
  
  #cv_lh_L1 <- rep(0,nrow(max_tpar))    
  #cv_lh_HARD <- rep(0,nrow(max_tpar))
  #cv_lh_SCAD <- rep(0,nrow(max_tpar))
  
  for(i_para in 1:nrow(max_tpar)){
    
    for(i_test in 1:k_fold){
      
      y_test <- y[ind_cv==i_test] ; y_train <- y[ind_cv!=i_test]
      z_test <- z[ind_cv==i_test] ; z_train <- z[ind_cv!=i_test]
      xd_test <- xd[ind_cv==i_test,] ; xd_train <- xd[ind_cv!=i_test,]
      pi1z <- pi_1(z_train)  ; pi2z <- 1 - pi1z
      
      out_acL1 <- semipnormmixEM_L1(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
      l_ini <- 0
      while(out_acL1$iter==1000 & l_ini < 10) {
        l_ini <- l_ini + 1
        initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
        initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
        out_acL1 <- semipnormmixEM_L1(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
        
      }
      
      out_acHARD <- semipnormmixEM_HARD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
      l_ini <- 0
      while(out_acHARD$iter==1000 & l_ini < 10) {
        l_ini <- l_ini + 1
        initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
        initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
        out_acHARD <- semipnormmixEM_HARD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
        
      }  
      
      out_acSCAD <- semipnormmixEM_SCAD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1,initl_beta2,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
      l_ini <- 0
      while(out_acSCAD$iter==1000 & l_ini < 10) {
        l_ini <- l_ini + 1
        initl_beta1_new <- beta_1 + rnorm(p,sd=0.5)
        initl_beta2_new <- beta_2 + rnorm(p,sd=0.5)
        out_acSCAD <- semipnormmixEM_SCAD(z=z_train,x=xd_train,y_train,pi1z,pi2z,initl_beta1_new,initl_beta2_new,initl_sgm21,initl_sgm22,gama1=max_tpar[i_para,1],gama2=max_tpar[i_para,2],maxiter=1000,epsilon=0.001,bw=max_tpar[i_para,3])
        
      }
      
      r1_L1 <- out_acL1$pi1*dnorm(y_train,mean=xd_train%*%out_acL1$beta1,sd=sqrt(out_acL1$sigma2_1))/(out_acL1$pi1*dnorm(y_train,mean=xd_train%*%out_acL1$beta1,sd=sqrt(out_acL1$sigma2_1))+out_acL1$pi2*dnorm(y_train,mean=xd_train%*%out_acL1$beta2,sd=sqrt(out_acL1$sigma2_2)))
      r2_L1 <- 1 - r1_L1
      r1_HARD <- out_acHARD$pi1*dnorm(y_train,mean=xd_train%*%out_acHARD$beta1,sd=sqrt(out_acHARD$sigma2_1))/(out_acHARD$pi1*dnorm(y_train,mean=xd_train%*%out_acHARD$beta1,sd=sqrt(out_acHARD$sigma2_1))+out_acHARD$pi2*dnorm(y_train,mean=xd_train%*%out_acHARD$beta2,sd=sqrt(out_acHARD$sigma2_2)))
      r2_HARD <- 1 - r1_HARD
      r1_SCAD <- out_acSCAD$pi1*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1))/(out_acSCAD$pi1*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta1,sd=sqrt(out_acSCAD$sigma2_1))+out_acSCAD$pi2*dnorm(y_train,mean=xd_train%*%out_acSCAD$beta2,sd=sqrt(out_acSCAD$sigma2_2)))
      r2_SCAD <- 1 - r1_SCAD
      
      pi1_L1test <- sapply(z_test,function(q) r1_L1%*%EpaKernel(z_train-q,max_tpar[i_para,3])/sum(EpaKernel(z_train-q,max_tpar[i_para,3])))
      pi2_L1test <- 1 - pi1_L1test
      pi1_HARDtest <- sapply(z_test,function(q) r1_HARD%*%EpaKernel(z_train-q,max_tpar[i_para,3])/sum(EpaKernel(z_train-q,max_tpar[i_para,3])))
      pi2_HARDtest <- 1 - pi1_HARDtest 
      pi1_SCADtest <- sapply(z_test,function(q) r1_SCAD%*%EpaKernel(z_train-q,max_tpar[i_para,3])/sum(EpaKernel(z_train-q,max_tpar[i_para,3])))
      pi2_SCADtest <- 1 - pi1_SCADtest 
      
      m1_L1test <- xd_test%*%out_acL1$beta1
      m2_L1test <- xd_test%*%out_acL1$beta2
      m1_HARDtest <- xd_test%*%out_acHARD$beta1
      m2_HARDtest <- xd_test%*%out_acHARD$beta2
      m1_SCADtest <- xd_test%*%out_acSCAD$beta1
      m2_SCADtest <- xd_test%*%out_acSCAD$beta2
      
      r1_L1test <- pi1_L1test*dnorm(y_test,mean=m1_L1test,sd=sqrt(out_acL1$sigma2_1))/(pi1_L1test*dnorm(y_test,mean=m1_L1test,sd=sqrt(out_acL1$sigma2_1))+pi2_L1test*dnorm(y_test,mean=m2_L1test,sd=sqrt(out_acL1$sigma2_2)))
      r2_L1test <- 1 - r1_L1test 
      r1_HARDtest <- pi1_HARDtest*dnorm(y_test,mean=m1_HARDtest,sd=sqrt(out_acHARD$sigma2_1))/(pi1_HARDtest*dnorm(y_test,mean=m1_HARDtest,sd=sqrt(out_acHARD$sigma2_1))+pi2_HARDtest*dnorm(y_test,mean=m2_HARDtest,sd=sqrt(out_acHARD$sigma2_2)))
      r2_HARDtest <- 1 - r1_HARDtest 
      r1_SCADtest <- pi1_SCADtest*dnorm(y_test,mean=m1_SCADtest,sd=sqrt(out_acSCAD$sigma2_1))/(pi1_SCADtest*dnorm(y_test,mean=m1_SCADtest,sd=sqrt(out_acSCAD$sigma2_1))+pi2_SCADtest*dnorm(y_test,mean=m2_SCADtest,sd=sqrt(out_acSCAD$sigma2_2)))
      r2_SCADtest <- 1 - r1_SCADtest
      
      #cv_lh_L1test <- sum(log(pi1_L1test*dnorm(y_test,mean=m1_L1test,sd=sqrt(out_acL1$sigma2_1))+pi2_L1test*dnorm(y_test,mean=m2_L1test,sd=sqrt(out_acL1$sigma2_2))))
      #cv_lh_HARDtest <- sum(log(pi1_HARDtest*dnorm(y_test,mean=m1_HARDtest,sd=sqrt(out_acHARD$sigma2_1))+pi2_HARDtest*dnorm(y_test,mean=m2_HARDtest,sd=sqrt(out_acHARD$sigma2_2))))
      #cv_lh_SCADtest <- sum(log(pi1_SCADtest*dnorm(y_test,mean=m1_SCADtest,sd=sqrt(out_acSCAD$sigma2_1))+pi2_SCADtest*dnorm(y_test,mean=m2_SCADtest,sd=sqrt(out_acSCAD$sigma2_2))))
      
      #cv_lh_L1[i_para] <- cv_lh_L1[i_para] + cv_lh_L1test
      #cv_lh_HARD[i_para] <- cv_lh_HARD[i_para] + cv_lh_HARDtest
      #cv_lh_SCAD[i_para] <- cv_lh_SCAD[i_para] + cv_lh_SCADtest
      
      cv_error_L1[i_para,i_test] <- sum((y_test - (r1_L1test*m1_L1test+r2_L1test*m2_L1test))^2)
      cv_error_HARD[i_para,i_test] <- sum((y_test - (r1_HARDtest*m1_HARDtest+r2_HARDtest*m2_HARDtest))^2)
      cv_error_SCAD[i_para,i_test] <- sum((y_test - (r1_SCADtest*m1_SCADtest+r2_SCADtest*m2_SCADtest))^2)
      
    }
    print(i_para)
  }
  
  mincv_L1 <- min(rowMeans(cv_error_L1,na.rm=TRUE),na.rm = TRUE)
  indexcv_L1 <- rowMeans(cv_error_L1,na.rm=TRUE)== mincv_L1
  indexcv_L1[is.na(indexcv_L1)] <- FALSE
  
  mincv_HARD <- min(rowMeans(cv_error_HARD,na.rm=TRUE),na.rm = TRUE)
  indexcv_HARD <- rowMeans(cv_error_HARD,na.rm=TRUE)== mincv_HARD
  indexcv_HARD[is.na(indexcv_HARD)] <- FALSE
  
  mincv_SCAD <- min(rowMeans(cv_error_SCAD,na.rm=TRUE),na.rm = TRUE)
  indexcv_SCAD <- rowMeans(cv_error_SCAD,na.rm=TRUE)==mincv_SCAD
  indexcv_SCAD[is.na(indexcv_SCAD)] <- FALSE
  
  return(list(tparL1_opt = max_tpar[indexcv_L1,],tparHARD_opt = max_tpar[indexcv_HARD,],tparSCAD_opt = max_tpar[indexcv_SCAD,]))
}

the_cluster_cv <- makeCluster(16)
clusterSetRNGStream(the_cluster_cv,35)
clusterExport(the_cluster_cv,c("pLasso_deri","pHard_deri","pSCAD_deri","semipnormmixEM_L1","semipnormmixEM_HARD","semipnormmixEM_SCAD","rnpnormmix","EpaKernel","autocorr.mat","pi_1",
                               "pi_2","m_1","m_2","beta_1","beta_2","p","sigma2_1","sigma2_2","N","z","x","xd","pro","mu","sgm","k_fold","max_tpar","initl_beta1","initl_beta2",
                               "initl_sgm21","initl_sgm22"))
tpar_cv <- parLapply(the_cluster_cv,1:30,cv_spmixreg)
#h_cv <- clusterCall(cl = the_cluster_cv, cv_spmixreg,1:10)
stopCluster(the_cluster_cv)

tparL1_opt <- tparHARD_opt <- tparSCAD_opt <- NULL
for(i in 1:30){
  tparL1_opt <- rbind(tparL1_opt,tpar_cv[[i]]$tparL1_opt)
  tparHARD_opt <- rbind(tparHARD_opt,tpar_cv[[i]]$tparHARD_opt)
  tparSCAD_opt <- rbind(tparSCAD_opt,tpar_cv[[i]]$tparSCAD_opt)
}

write.csv(tparL1_opt,"tparL1_opt_p100indep.csv")
write.csv(tparHARD_opt,"tparHARD_opt_p100indep.csv")
write.csv(tparSCAD_opt,"tparSCAD_opt_p100indep.csv")













