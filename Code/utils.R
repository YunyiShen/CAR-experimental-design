## This get the duplicated matrix D_k

getD <- function(n){
  temp <-upper.tri( diag(n),T)
  diag(n^2)[,c(temp)]
}


# laplacian approximation of the prior
getLA_prior <- function(B,Sigma,phi,Lambda, lambda){
  k <- nrow(Sigma)
  D_k <-c(upper.tri( diag(k),T))
  iLambda <- solve(Lambda)
  a <- (lambda-k-p-1)/2
  A <- a * kronecker(Sigma,Sigma) + kronecker(Sigma,Sigma %*% (t(B)%*%iLambda%*%B))
  A <- A[D_k,D_k]  
  
  D <- kronecker(Sigma, iLambda)
  
  C <- -kronecker(Sigma, Sigma%*% t(B)%*%iLambda)
  C <- C[D_k,]
  
  return(list(A=A,D=D,C=C))
}

# Fisher information 
getFI <- function(X,B,Sigma){
  k <- nrow(Omega)
  n <- nrow(X)
  D_k <-c(upper.tri( diag(k),T))
  Sigma <- solve(Sigma)
  XB <- X %*% B
  XtX <- t(X) %*% X 
  A <- (n/2) * kronecker(Sigma, Sigma) + kronecker(Sigma, Sigma %*% t(XB)%*%XB%*%Sigma)
  A <- A[D_k,D_k]
  D <- kronecker(Sigma, XtX)
  C <- - kronecker(Sigma, Sigma %*% t(B)%*%XtX)
  C <- C[D_k,]
  
  return(list(A=A,D=D,C=C))
}


D_optim_Omega <- function(X, B, Sigma, phi, Lambda, lambda,X_old){
  p <- nrow(B)
  X <- matrix(X, ncol = p)
  X1 <- rbind(X_old,X)
  prior <- getLA_prior(B,Sigma,phi,Lambda, lambda)
  FI <- getFI(X1,B,Sigma)
  A <- prior$A + FI$A
  D <- prior$D + FI$D
  C <- prior$C + FI$C
  
  Schur <- A - C %*% solve(D , t(C))
  return(-determinant(Schur)$modulus)
}


Design_D_opt <- function(B, Sigma, phi, Lambda, lambda, n, X_old){
  p <- nrow(B)
  temp <- optim(runif(n*p,-1,1),D_optim_Omega, B=B,Sigma=Sigma,Lambda=Lambda,lambda=lambda, X_old = X_old,method = "L-BFGS-B",lower = -1, upper = 1)
  matrix(temp$par,ncol = p)
}


simu_data <- function(X,B,Sigma){
  n <- nrow(X)
  k <- nrow(Sigma)
  Y <- matrix(NA, n,k)
  XB <- X %*% B
  XtX <- t(X) %*% X
  for(i in 1:n){
    Y[i,] <- MASS::mvrnorm(1,Sigma %*% XB[i,] , Sigma)
  }
  return(Y)
}

# posterior mode of Wishart prior
getMAP_Omega_Wishart <- function(Y, X, B0, phi, lambda, Lambda ){
  invLambda <- solve(Lambda)
  B0tL <- t(B0) %*% invLambda
  YtX <- t(Y) %*% X
  XtX <- t(X) %*% X
  YtY <- t(Y) %*% Y
  phi_hat <- phi + YtY +B0tL%*%B0 - (B0tL+YtX) %*% solve(XtX+invLambda,t(B0tL+YtX))
  
  return(solve(phi_hat) * (lambda+nrow(Y)-ncol(Y)-1-ncol(X)))
}


getMAP_B_Wishart <- function(Y, X, B0, Lambda, Omega){
  invLambda <- solve(Lambda)
  B0tL <- t(B0) %*% invLambda
  YtX <- t(Y) %*% X
  XtX <- t(X) %*% X
  YtY <- t(Y) %*% Y
  (solve((invLambda+XtX),t(B0tL+YtX)) %*% Omega)
}

# posterior Mode of MGIG prior
getMAP_Omega_MGIG <- function(Y, X, B0, phi, psi,lambda, Lambda ){
  invLambda <- solve(Lambda)
  XtY <- t(X) %*% Y
  XtX <- t(X) %*% X
  sP <- (XtX) + invLambda
  invLB <- solve(Lambda, B)
  
  phi_hat <-phi + t(Y) %*% Y-t(XtY) %*%solve(sP, XtY) 
  phi_hat <-( phi_hat + t(phi_hat))/2
  psi_hat <- psi+ t(B)%*%invLB-t(invLB)%*%solve(sP,invLB)
  psi_hat <- (psi_hat + t(psi_hat))/2
  rMGIG::mMGIG(lambda + nrow(Y)/2, phi_hat, psi_hat)
}

getMAP_B_MGIG <- function(Y, X, B0, Lambda, Omega){
  invLambda <- solve(Lambda)
  B0tL <- t(B0) %*% invLambda
  YtX <- t(Y) %*% X
  XtX <- t(X) %*% X
  YtY <- t(Y) %*% Y
  (solve((invLambda+XtX),t(solve(Omega, B0tL)+YtX)) %*% Omega)
}

linear_Omega <- function(Omega_list){
  temp <- lapply(Omega_list, function(w){c(w[upper.tri(w,T)])})
  Reduce(rbind,temp)
}

linear_B <- function(B_list){
  temp <- lapply(B_list, function(w){c(w)})
  Reduce(rbind,temp)
}

get_posterior_para_Omega_MGIG <- function(Y, X, B0, phi, psi,lambda, Lambda ){
  invLambda <- solve(Lambda)
  XtY <- t(X) %*% Y
  XtX <- t(X) %*% X
  sP <- (XtX) + invLambda
  invLB <- solve(Lambda, B)
  
  phi_hat <-phi + t(Y) %*% Y-t(XtY) %*%solve(sP, XtY) 
  phi_hat <-( phi_hat + t(phi_hat))/2
  psi_hat <- psi+ t(B)%*%invLB-t(invLB)%*%solve(sP,invLB)
  psi_hat <- (psi_hat + t(psi_hat))/2
  list(nu = lambda + nrow(Y)/2, phi = phi_hat, psi = psi_hat)
}


get_posterior_para_Omega_Wishart <- function(Y, X, B0, phi, lambda, Lambda ){
  invLambda <- solve(Lambda)
  B0tL <- t(B0) %*% invLambda
  YtX <- t(Y) %*% X
  XtX <- t(X) %*% X
  YtY <- t(Y) %*% Y
  phi_hat <- phi + YtY +B0tL%*%B0 - (B0tL+YtX) %*% solve(XtX+invLambda,t(B0tL+YtX))
  
  list(phi = phi_hat, lambda = lambda+nrow(Y))
}

KLdiv_wishart <- function(lambda, phi1, phi2){
  ph2invphi <- solve(phi2, phi1)
  (-lambda/2 ) * determinant(ph2invphi)$modulus + (lambda/2) * (sum(diag(ph2invphi))-nrow(phi1))
}
