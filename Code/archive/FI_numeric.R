A <- function(rho,beta1,beta2,X=1,mu1=0,mu2=0 ){
  beta <- expand.grid(beta1,beta2)
  beta1 <- beta[,1]
  beta2 <- beta[,2]
  ((3*rho^2+1)/(1-rho^2)^3)*((beta1*X+mu2)^2+(beta2*X+mu2)^2)-
    ((2*rho)*(rho^2+3)/(1-rho^2)^3) * ((beta1*X+mu1)*(beta2*X+mu2)) + 
    1*((rho^2+1)/(1-rho^2))^2
}

Isg <- function(rho,beta1,beta2,X=1,mu1=0,mu2=0 ){
  beta <- expand.grid(beta1,beta2)
  beta1 <- beta[,1]
  beta2 <- beta[,2]
  a <- exp(X*beta1+mu1+X*beta2+mu2-rho) + exp(-X*beta1-mu1-X*beta2-mu2-rho)
  b <- exp(-X*beta1-mu1+X*beta2+mu2+rho) + exp(X*beta1+mu1-X*beta2-mu2+rho)
  
  (1-(a-b)^2/(a+b)^2)
}

Isg_FI <- function(rho,beta1,beta2){
  a <- exp(beta1+beta2-rho)
  b <- exp(beta1-beta2+rho)
  c <- exp(-beta1+beta2+rho)
  d <- exp(-beta1-beta2-rho)
  
  s <- a+b+c+d
  
  res <- matrix(0,3,3)
  
  res[1,2] <- -(a+c-b-d)/s+(a+b-c-d)*(a+d-b-c)/(s^2)
  res[1,3] <- -(a-c+b-d)/s+(a-b+c-d)*(a+d-b-c)/(s^2)
  res[2,3] <- (a+d-b-c)/(s)-(a+b-c-d)*(a-b+c-d)/(s^2)
  res <- res + t(res)
  res[1,1] <- 1-(a+d-b-c)^2/(s^2)
  res[2,2] <- 1-(a+b-c-d)^2/(s^2)
  res[3,3] <- 1-(a-b+c-d)^2/(s^2)
  return(res)
}

beta1 <- seq(-3,3,0.01)
rhos <- c(-0.8,-0.3,-0.1,0.8,0.3,0.1)
par(mfcol = c(3,2))
for(rho in rhos){
  ww <- A(rho,beta1,beta1,1,0,0)#+A(rho,beta1,beta1,0,0,0)
  contour(beta1,beta1,matrix(1/ww,length(beta1),length(beta1),byrow = F),xlab = "beta1",ylab = "beta2",main = paste0("rho=",rho))
  #abline(0,1)
  #abline(-2,-1)
}

beta1 <- seq(-10,10,0.1)
rhos <- c(-1,-0.5,-0,1,0.5,0)
par(mfcol = c(3,2))
for(rho in rhos){
  #cat(rho)
  ww <- Isg(rho,beta1,beta1,1,0,0)+Isg(rho,beta1,beta1,0,0,0)
  image(matrix(1/ww,length(beta1),length(beta1)))
  #abline(0,sign(rho))
  #contour(beta1,beta1,matrix(1/ww,length(beta1),length(beta1)),xlab = "beta1",ylab = "beta2",main = paste0("rho=",rho))
}
    





