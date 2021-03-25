library(rMGIG)
library(MASS)
source("./Code/Graph_generator.R")
source("./Code/utils.R")


k <- 5
p <- 2*k

n <- 150



G <- g_model1(k)
Omega <- G$Omega
Sigma <- G$Sigma

set.seed(42)
X <- matrix(rnorm(n*p),n,p)
B <- matrix(rnorm(p*k),p,k) > 0

Y <- simu_data(X,B,Sigma)

lambda <- k
psi <- diag(k)
phi <- diag(k)
#phi <- Omega

Lambda <- 0.1 * diag(p)

invLambda <- solve(Lambda)
B0 <- B %*% Sigma
B0tL <- t(B0) %*% invLambda
YtX <- t(Y) %*% X

XtX <- t(X) %*% X

YtY <- t(Y) %*% Y

psi_hat <- psi 
phi_hat <- phi + YtY +B0tL%*%B0 - (B0tL+YtX) %*% solve(XtX+invLambda,t(B0tL+YtX))

#www <- rMGIG::rMGIG(1000,6 + n/2, psi = psi_hat,phi=phi_hat, df = 10000,list = T)
#Reduce("+",www)/1000

MAP_conf_prior <- solve(phi_hat) * (lambda+n-k-1-p)
CARlasso:::stein_loss(Omega, MAP_conf_prior)

getMAP_Omega_Wishart(Y, X, B0, phi, lambda, Lambda)
B_Wishart <- getMAP_B_Wishart(Y, X, B0, phi, lambda, Lambda, MAP_conf_prior)


XtY <- t(YtX)
sP <- (XtX) + invLambda
invLB <- solve(Lambda, B)

phi_hat <-phi + t(Y) %*% Y-t(XtY) %*%solve(sP, XtY) 
phi_hat <-( phi_hat + t(phi_hat))/2
psi_hat <- psi+ t(B)%*%invLB-t(invLB)%*%solve(sP,invLB)
psi_hat <- (psi_hat + t(psi_hat))/2



MAP_good_prior <- rMGIG::mMGIG(5 + n/2, phi_hat, 
                               psi_hat)

getMAP_Omega_MGIG(Y, X, B0, phi, lambda, Lambda)


CARlasso:::stein_loss(Omega, MAP_good_prior)

B_MGIG <- getMAP_B_MGIG(Y, X, B0, phi, lambda, Lambda, MAP_good_prior)

Design_D_opt(B, Sigma, phi, Lambda, lambda, 5)
