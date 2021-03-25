#if(!require(rMGIG)){
#  devtools::install_github("YunyiShen/rMGIG")
#}
library(rMGIG)
library(MASS)
source("./Code/Graph_generator.R")
source("./Code/utils.R")

n <- 1000
p <- 10
k <- 10

X <- matrix( rnorm(n * p) , n , p )
B <- matrix( rnorm(p * k) , p , k )

G <- g_model1(k)
Omega <- G$Omega
Sigma <- G$Sigma

Y <- matrix(NA, n,k)
XB <- X %*% B
XtX <- t(X) %*% X



for(i in 1:n){
  Y[i,] <- MASS::mvrnorm(1,Sigma %*% XB[i,] , Sigma)
}

psi <- 0.01 * diag(k)
phi <- 0.01 * diag(k)
Lambda <- 1 * diag(p)
sP <- (XtX) + solve(Lambda)
invLB <- solve(Lambda, B)
XtY <- t(X) %*% Y 

# + t(Y) %*% Y + t(B)%*%solve(Lambda, B)  + t(B)%*%t(X)%*%(X)%*%B
nn <- 100
posterior_less <- rMGIG::rMGIG(nn, 5+n/2, phi + t(Y) %*% Y-t(XtY) %*%solve(sP, XtY), 
                        psi+ t(B)%*%invLB-t(invLB)%*%solve(sP,invLB) , 
                        100*n, T,1e5)

post_mean <- Reduce("+",posterior_less)/nn

image(post_mean)

post_B <- lapply(posterior_less, sample_B_post, Y, X, B, Lambda)
post_mean_B <- Reduce("+",post_B)/nn
image(B)
image(post_mean_B)


CARlasso:::stein_loss_cpp(post_mean,Omega)

resp <- paste0("Y",1:k)
colnames(Y) <- resp
pred <- paste0("X",1:p)
colnames(X) <- pred
df <- as.data.frame(cbind(X,Y))
fomu <- paste(paste0(resp,collapse = "+"),paste0(pred,collapse = "+"),sep = "~")

www <- CARlasso::CARlasso(as.formula(fomu),data = df, adaptive = T)

car_Omega <- lapply(1:200, rMGIG:::get_mat, t(www$Omega), k)
post_mean_car <- Reduce("+",car_Omega)/200
