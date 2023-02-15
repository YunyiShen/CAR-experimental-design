library(rMGIG)
library(MASS)
library(CholWishart)
source("./Code/Graph_generator.R")
source("./Code/utils.R")

colnorm <- function(x){
  (x-mean(x))/sd(x)
}
colcent <- function(x){
  x-mean(x)
}
logit <- function(x) log(x/(1-x))
divided_by_last <- function(x){
  last <- x[,ncol(x)] 
  first <- x[,-ncol(x)] + 1
  apply(X=as.matrix(first),MARGIN = 2,FUN = function(xx,y){xx/y},y=last)
}

load("./data/mgp154.rda")

formulea <- paste0("Alistipes~",paste0(c(colnames(mgp154)[17:21]),collapse = "+"))
X <- model.matrix(as.formula(formulea), data = mgp154)[,-c(1)] |>
  apply(2,colnorm)
colnames(X) <- NULL

Y <- mgp154[,2:16] |> divided_by_last() |> log() |> apply(2,colcent) #|> logit()
init <- lm.fit(x=X,y=Y)$coefficients

n <- nrow(Y)
p <- ncol(X)
k <- ncol(Y)
lambda <- k+1

#B0 <- matrix(0,p,k) #init
B0 <- init
set.seed(42)
Yprime <- mvrnorm(n, mu = rep(0,k), cov(Y-X%*%init))
par(mfrow = c(1,2))
hist(Y)
hist(Yprime)

psi <- .01 * diag(k)
phi <- .01 * diag(k)
prior_cert <- 10^(-1:1)

wishart0 <- wishart <- list()
mgig0 <- mgig <- list()

q <- 2
qprime <- 7
df <- 300
df0 <- 200
for(i in 1:3){
  wishart_param <- get_posterior_para_Omega_Wishart(Y, X, B0%*%cov(Y), phi, 2*lambda, prior_cert[i]*diag(p))
  wishart_sample <- rWishart(5e4, wishart_param$lambda, solve(wishart_param$phi))[q,qprime,]
  wishart[[i]] <- wishart_sample
  
  wishart_param <- get_posterior_para_Omega_Wishart(Yprime, 0*X, B0%*%cov(Y), phi, 2*lambda, prior_cert[i]*diag(p))
  wishart_sample <- rWishart(5e4, wishart_param$lambda, solve(wishart_param$phi))[q,qprime,]
  wishart0[[i]] <- wishart_sample
  
  mgig_param <- get_posterior_para_Omega_MGIG(Y, X, B0, phi, psi,lambda,  prior_cert[i]*diag(p))
  mgig_sample <- rMGIG::rMGIG(n=5e4, mgig_param$nu,
                              mgig_param$phi,
                              mgig_param$psi, 
               df = df , maxit = 1e5)|>
    sapply(`[`,q,qprime)
  mgig[[i]] <- mgig_sample
  
  mgig_param <- get_posterior_para_Omega_MGIG(Yprime, 0*X, B0, phi, psi,lambda,  prior_cert[i]*diag(p))
  mgig_sample <- rMGIG::rMGIG(n=5e4, mgig_param$nu,
                              mgig_param$phi,
                              mgig_param$psi, 
                              df = df0, maxit = 1e5)|>
    sapply(`[`,q,qprime)
  mgig0[[i]] <- mgig_sample
}

png("gut-posterior.png", width = 10, height = 4, units = "in", res = 500)
par(mar = c(1.5,2.5,1.5,1.5), mgp = c(1.1, 0.5, 0), oma =c(0,1,0,0) )
par(mfcol = c(2,3))
for(i in 1:3){
  plot(density(mgig[[i]]), main = prior_cert[i], xlim = c(-1,.5), xlab = "", ylab = ifelse(i==1,"MGIG",""))
  lines(density(mgig0[[i]]), col = "red")
  if(i==1){
    legend("topleft", legend = c("experiment","none"), col = c("black","red"),lty = c(1,1))
  }
  plot(density(wishart[[i]]),main = "", xlim = c(-1,.5), xlab = "",ylab = ifelse(i==1,"Wishart",""))
  lines(density(wishart0[[i]]), col = "red")
}
dev.off()

