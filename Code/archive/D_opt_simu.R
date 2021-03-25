source("./Code/utils.R")
k <- 10
p <- k+1 
#p <- 10
n <- 20


#G <- g_model2(k)
#Omega <- G$Omega
#Sigma <- G$Sigma
#set.seed(42)
Omega <- rWishart(1, k+5, 0.1*diag(k))[,,1]
Sigma <- solve(Omega)
image(Omega)


B <- matrix( rnorm(p * k) , p , k )
#B <- (B >=0) * 1.0 
#B <- rbind(rnorm(k),diag(k))
Lambda <- .01 * diag(p)

Xs <- lapply(1:1000, function(i,n,p){cbind(1,matrix(sign( runif(n*(p-1),-1,1)),n,p-1))},n,p)

FIs <- lapply(Xs, getLA2, B, Omega, Lambda)

dets <- sapply(FIs, function(w){determinant(w)$modulus})

hist(dets)

D_opt <- optim(runif(n * (p-1),-1,1),fn = Dopt_Omega2, X_old = NULL,B = B+matrix(rnorm(k*p,0,0.5),p,k), Omega = solve( rWishart(1, 20,  Sigma/20)[,,1]), Lambda = Lambda, lower = -1, upper = 1, method = "L-BFGS-B")
D_opt_des <- cbind(1,matrix(D_opt$par,nrow = n))
    
LA_D_opt <- getLA2(D_opt_des, B, Omega, Lambda)
determinant(LA_D_opt)$modulus

abline(v = determinant(LA_D_opt)$modulus, col = "red")
