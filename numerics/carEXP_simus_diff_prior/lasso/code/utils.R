rlaplace <- function(n, lambdasqr){
  tmp <- rexp(n, lambdasqr)
  signs <- sign(rnorm(n))
  tmp * signs
}

ispd <- function(x){
  eigen_val <- eigen(x, TRUE, TRUE)
  all(eigen_val$values>0)
}

rglasso <- function(n, lambda, k=30, factoring = 1000){
  res <- list()
  n_got <- 0
  for(i in 1:(factoring * n)){
    tmp <- matrix(0, k,k)
    diag(tmp) <- rexp(k, lambda)
    tmp[upper.tri(tmp)] <- rlaplace(k*(k-1)/2, lambda^2)
    tmp <- tmp + t(tmp)
    if(ispd(tmp)){
      n_got <- n_got + 1
      res[[n_got]] <- tmp
    }
    if(n_got >= n){
      res <- lapply(res, function(x){
        x[upper.tri(x, TRUE)]
      })
      res <- Reduce(rbind, res)
      return(res)
    }
  }
  res <- lapply(res, function(x){
    x[upper.tri(x, TRUE)]
  })
  res <- Reduce(rbind, res)
  return(res)
}





symtry <- function(x){
    return((x+t(x))/2)
}




log_sum_exp <- function(x){
    xstar <- max(x)
    xstar + log(sum(exp(x-xstar)))
}


## This get the duplicated matrix D_k

getD <- function(n){
  temp <-upper.tri( diag(n),T)
  diag(n^2)[,c(temp)]
}


simu_data <- function(X,B,Sigma){
  n <- nrow(X)
  k <- nrow(Sigma)
  Y <- matrix(NA, n,k)
  XB <- X %*% B
  XtX <- t(X) %*% X
  for(i in 1:n){
    Y[i,] <- MASS::mvrnorm(1,rep(0,k) , Sigma) + Sigma %*% XB[i,]
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
  list(nu = lambda + nrow(Y)/2, phi = symtry(phi_hat), psi = symtry(psi_hat))
}


get_posterior_para_Omega_Wishart <- function(Y, X, B0, phi, lambda, Lambda ){
  invLambda <- solve(Lambda)
  B0tL <- t(B0) %*% invLambda
  YtX <- t(Y) %*% X
  XtX <- t(X) %*% X
  YtY <- t(Y) %*% Y
  phi_hat <- phi + YtY +B0tL%*%B0 - (B0tL+YtX) %*% solve(XtX+invLambda,t(B0tL+YtX))
  
  list(phi = symtry(phi_hat), lambda = lambda+nrow(Y))
}

KLdiv_wishart <- function(lambda, phi1, phi2){
  ph2invphi <- solve(phi2, phi1)
  (-lambda/2 ) * determinant(ph2invphi)$modulus + (lambda/2) * (sum(diag(ph2invphi))-nrow(phi1))
}

KLdiv_wishart2 <- function(samples, lambda1, phi1, lambda2, phi2){
  logfoverg <- sapply(samples, function(X, lambda1,phi1, lambda2, phi2){
    CholWishart::dWishart(X, lambda1, phi1) - CholWishart::dWishart(X, lambda2, phi2)
  }, lambda1,phi1, lambda2, phi2)
  mean(logfoverg)
}

KLdiv_MGIG2 <- function(samples, nu1, phi1, psi1,nu2, phi2, psi2){
    logfoverg <- sapply(samples, function(X,nu1, phi1, psi1, nu2, phi2, psi2 ){
        fMGIG(X,nu1, phi1, psi1) - fMGIG(X,nu2, phi2, psi2)
    },nu1, phi1, psi1, nu2, phi2, psi2)

    term1 <- log_sum_exp(-logfoverg) - log(length(logfoverg))
    term2 <- mean(logfoverg)
    
    term1+term2
}

stein_loss <- function(X_hat,X){
  p <- nrow(X)
  inv_X <- solve(X)
  XhatinvX <- X_hat %*% inv_X
  sum(diag(XhatinvX)) - determinant(XhatinvX)$modulus[1] - p
}




