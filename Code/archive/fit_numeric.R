loglikf <- function(rho,y,beta,X){
  ll <- 0
  Omega <- matrix(c(1,rho,rho,1),2,2)
  Sigma <- solve(Omega)
  for(i in 1:length(X)){
    mu <- Sigma %*% (beta * X[i])
    ll <- ll + (-0.5*log(det(Sigma))-0.5*t(y[i,]-mu) %*% Omega %*% (y[i,]-mu))
  }
  return(-ll)
}

loglikf_ub <- function(par,y,X){
  rho <- par[1]
  beta <- par[2:3]
  ll <- 0
  Omega <- matrix(c(1,rho,rho,1),2,2)
  Sigma <- solve(Omega)
  for(i in 1:length(X)){
    mu <- Sigma %*% (beta * X[i])
    ll <- ll + (-0.5*log(det(Sigma))-0.5*t(y[i,]-mu) %*% Omega %*% (y[i,]-mu))
  }
  return(-ll)
}

simu <- function(rho,beta1,X){
  res <- matrix(NA,length(X),2)
  Omega <- matrix(c(1,rho,rho,1),2,2)
  Sigma <- solve(Omega)
  for(i in 1:length(X)){
    mu <- Sigma%*%(beta1 * X[i])
    res[i,] <- MASS::mvrnorm(1,mu,Sigma)
  }
  return(res)
}

rho <- 0.5
beta_1 <- c(1,1)
beta_2 <- c(1,-1)
beta_3 <- c(1,0)
beta_4 <- c(0,0)


nrep <- 5000
n <- 10
res1 <- rep(NA,nrep)
res4 <- res3 <- res2 <- res1
res4_ub <- res3_ub <- res2_ub <- res1_ub <- res1
X <- 1*(runif(n)>-0.5)
for(i in 1:nrep){
  cat(i,"\n")
  data1 <- simu(rho,beta_1,X)
  data2 <- simu(rho,beta_2,X)
  data3 <- simu(rho,beta_3,X)
  data4 <- simu(rho,beta_4,X)
  fit1 <- optim(0.1,fn = loglikf,
                lower = -0.999,upper = 0.999,method = "L-BFGS-B",
                y = data1,beta = beta_1,X = X)
  fit2 <- optim(0.1,fn = loglikf,
                lower = -0.999,upper = 0.999,method = "L-BFGS-B",
                y = data2,beta = beta_2,X = X)
  fit3 <- optim(0.1,fn = loglikf,
                lower = -0.999,upper = 0.999,method = "L-BFGS-B",
                y = data3,beta = beta_3,X = X)
  
  fit4 <- optim(0.1,fn = loglikf,
                 lower = -0.999,upper = 0.999,method = "L-BFGS-B",
                 y = data4,beta = beta_4,X = X)
  
  fit1_ub <- optim(c(0.1,0,0),fn = loglikf_ub,
                   lower = c(-0.999,-Inf,-Inf),upper = c(0.999,Inf,Inf),
                   method = "L-BFGS-B",
                   y = data1,X = X)
  fit2_ub <- optim(c(0.1,0,0),fn = loglikf_ub,
                   lower = c(-0.999,-Inf,-Inf),upper = c(0.999,Inf,Inf),
                   method = "L-BFGS-B",
                   y = data2,X = X)
  fit3_ub <- optim(c(0.1,0,0),fn = loglikf_ub,
                   lower = c(-0.999,-Inf,-Inf),upper = c(0.999,Inf,Inf),
                   method = "L-BFGS-B",
                   y = data3,X = X)
  
  fit4_ub <- optim(c(0.1,0,0),fn = loglikf_ub,
                   lower = c(-0.999,-Inf,-Inf),upper = c(0.999,Inf,Inf),
                   method = "L-BFGS-B",
                   y = data4,X = X)
  
  res1[i] <- fit1$par
  res2[i] <- fit2$par
  res3[i] <- fit3$par
  res4[i] <- fit4$par
  
  res1_ub[i] <- fit1_ub$par[1]
  res2_ub[i] <- fit2_ub$par[1]
  res3_ub[i] <- fit3_ub$par[1]
  res4_ub[i] <- fit4_ub$par[1]
}

c1 <- rgb(0,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,0,0, max = 255, alpha = 100, names = "lt.pink")

par(mfrow = c(2,2))
hist(res1_ub,xlim = c(-1,1),main = "beta=(1,1)",
     xlab = "rho",freq = F,col = c2,ylim = c(0,4.75))
res1_plot <- hist(res1,plot = F)
plot(res1_plot,freq = F,add = T,xlim = c(-1,1),
     main = "beta=(1,1)",xlab = "rho",col = c1,)

hist(res2_ub,xlim = c(-1,1),main = "beta=(1,-1)",
     xlab = "rho",freq = F,col = c2,ylim = c(0,4.75))
res2_plot <- hist(res2,plot = F)
plot(res2_plot,freq = F,add = T,xlim = c(-1,1),
     main = "beta=(1,-1)",xlab = "rho",col = c1,)

hist(res3_ub,xlim = c(-1,1),main = "beta=(1,0)",
     xlab = "rho",freq = F,col = c2,ylim = c(0,4.75))
res3_plot <- hist(res3,plot = F)
plot(res3_plot,freq = F,add = T,xlim = c(-1,1),
     main = "beta=(1,0)",xlab = "rho",col = c1,)

hist(res4_ub,xlim = c(-1,1),main = "beta=(0,0)",
     xlab = "rho",freq = F,col = c2,ylim = c(0,4.75))
res4_plot <- hist(res4,plot = F)
plot(res4_plot,freq = F,add = T,xlim = c(-1,1),
     main = "beta=(0,0)",xlab = "rho",col = c1,)


var(res1)
var(res2)
var(res3)
var(res4)


var(res1_ub)
var(res2_ub)
var(res3_ub)
var(res4_ub)
  
