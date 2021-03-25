library(rMGIG)
library(MASS)
source("./Code/Graph_generator.R")
source("./Code/utils.R")

res_al <- "./Res/AL/Model1_AL_k5p10_kb.csv"
res_rd <- "./Res/AL/Model1_RD_k5p10_kb.csv"

#L2_al <- "./Res/AL/Model1_AL_L2_k3p6_kn.csv"
#L2_rd <- "./Res/AL/Model1_RD_L2_k3p6_kb.csv"




k <- 5
p <- 2*k 
#p <- 2*k
step_size <- 5 
n_steps <- 10
n_reps <- 20
n_init <- 10



#Omega <- G$Omega
#Sigma <- G$Sigma
set.seed(123)
Omega <- rWishart(1, k+5, 0.1*diag(k))[,,1]
Sigma <- solve(Omega)
image(Omega)

B <- matrix( rnorm(p * k,0,3) , p , k )
#B <- (B >=0) * 1.0 
#B <- rbind(0,diag(k))
#B <- rbind(rnorm(k),diag(k))

psi <- .01 * diag(k)
phi <- .01 * diag(k)
Lambda <- 1 * diag(p)
invLB <- solve(Lambda, B)

nn <- 500
a <- 1000 # troublesome tuning parameter for rMGIG
b <- 1000
c <- 1000
d <- 1000
maxit <- 10000


res_random_design <- data.frame(matrix(NA, nrow = n_reps, ncol = n_steps+1))
res_active_learn <- res_random_design

L2_random_design <- data.frame(matrix(NA, nrow = n_reps, ncol = n_steps+1))
L2_active_learn <- res_random_design



for(i_rep in 1:n_reps){
  ## initial samples
  n <- n_init
  X_init <- matrix( runif(n * p - n, -3,3 ) , n , p-1 )
  X_init <- cbind(1,X_init)
  Y_init <- simu_data(X_init,B,Sigma)
  
  
  XtX <- t(X_init) %*% X_init
  XtY <- t(X_init) %*% Y_init 
  sP <- (XtX) + solve(Lambda)

  cat("Round:", i_rep,"Initial estimation\n")
  posterior_Omega <- rMGIG::rMGIG(nn, 5+n/2, phi + t(Y_init) %*% Y_init, 
                                psi+ t(B)%*%XtX %*% B , 
                                a*k+b + c * n, T,maxit)

  post_mean_Omega <- Reduce("+",posterior_Omega)/nn
  
  res_active_learn[i_rep,1] <- res_random_design[i_rep,1] <- CARlasso:::stein_loss_cpp(post_mean_Omega,Omega)
  #L2_active_learn[i_rep,1] <- L2_random_design[i_rep,1] <- sum((post_mean_Omega[upper.tri(Omega)]-Omega[upper.tri(Omega)])^2)
  
  
  
  
  
  post_mean_Omega_design <- post_mean_Omega
  
  Y_design <- Y_rand <- Y_init
  X_design <- X_rand <- X_init
  for(i_step in 1:n_steps + 1){
    
    # get data
    n <- n + step_size
    cat("  step :" , i_step, " Designing..\n")
    X_temp_design <- Design_X_KB(B, post_mean_Omega_design, X_design, step_size)#,optimality = Dopt_full2)
    #X_temp_design <- Design_X_KB(B, Omega, X_design, step_size)#,optimality = Dopt_full2)
    
    #X_temp_design <- Design_X2(B, Omega, X_design,Lambda, step_size)
    
    
    X_temp_rand <- matrix( 3 * sign( runif(step_size * p - step_size,-3,3)), nrow = step_size)
    X_temp_rand <- cbind(1,X_temp_rand)
    
    Y_temp_design <- simu_data(X_temp_design, B, Sigma)
    Y_temp_rand <- simu_data(X_temp_rand, B, Sigma)
    
    X_rand <- rbind(X_rand, X_temp_rand)
    Y_rand <- rbind(Y_rand, Y_temp_rand)
    
    X_design <- rbind(X_design, X_temp_design)
    Y_design <- rbind(Y_design, Y_temp_design)
    
    # sample posterior rand
    XtX_rand <- t(X_rand) %*% X_rand
    XtY_rand <- t(X_rand) %*% Y_rand 
    
    
    cat("  samling random\n")
    phi_hat_rand <-phi + t(Y_rand) %*% Y_rand
    phi_hat_rand <-( phi_hat_rand + t(phi_hat_rand))/2
    psi_hat_rand <- psi+ t(B) %*% XtX_rand %*% B
    psi_hat_rand <- (psi_hat_rand + t(psi_hat_rand))/2
    
    posterior_Omega_rand <- rMGIG::rMGIG(nn, 5 + n/2, phi_hat_rand, 
                                    psi_hat_rand , 
                                    a*k+b + c * n, T,maxit)
    
    post_mean_Omega_rand <- Reduce("+",posterior_Omega_rand)/nn
    
    
    res_random_design[i_rep,i_step] <- CARlasso:::stein_loss_cpp(post_mean_Omega_rand,Omega)
    #L2_random_design[i_rep,i_step] <- sum((post_mean_Omega_rand[upper.tri(Omega)]-Omega[upper.tri(Omega)])^2)
    
    
    write.csv(res_random_design,res_rd,row.names = F)
    #write.csv(L2_random_design,L2_rd,row.names = F)
    
    
    ## sample the designed
    XtX_design <- t(X_design) %*% X_design
    XtY_design <- t(X_design) %*% Y_design 
    
    cat("  sampling active \n")
    phi_hat_design <-phi + t(Y_design) %*% Y_design
    phi_hat_design <-( phi_hat_design + t(phi_hat_design))/2
    psi_hat_design <- psi+ t(B)%*%XtX_design %*% B
    psi_hat_design <- (psi_hat_design + t(psi_hat_design))/2
    
    posterior_Omega_design <- rMGIG::rMGIG(nn, 5 + n/2, phi_hat_design, 
                                         psi_hat_design , 
                                         a*k+b + c * n + d * n, T,maxit)
    
    post_mean_Omega_design <- Reduce("+",posterior_Omega_design)/nn
    
    
    res_active_learn[i_rep,i_step] <- CARlasso:::stein_loss_cpp(post_mean_Omega_design,Omega)
    #L2_active_learn[i_rep,i_step] <- sum((post_mean_Omega_design[upper.tri(Omega)]-Omega[upper.tri(Omega)])^2)
    
    
    write.csv(res_active_learn,res_al,row.names = F)
    #write.csv(res_active_learn,L2_al,row.names = F)
    
    #par(mfrow = c(2,1))
    plot_out <- data.frame(step = 0:n_steps, 
                           random = t(as.matrix(res_random_design[i_rep,])), 
                           active = t(as.matrix(res_active_learn[i_rep,] )), row.names = NULL)
    matplot(plot_out[,1], plot_out[,2:3], type = "l",xlab = "steps", ylab = "Stein's loss")
    legend("topright", legend = c("random","active"),lty = 1:2,col = 1:2)
    
    #plot_out <- data.frame(step = 0:n_steps, 
    #                       random = t(as.matrix(L2_random_design[i_rep,])), 
    #                       active = t(as.matrix(L2_active_learn[i_rep,] )), row.names = NULL)
    #matplot(plot_out[,1], plot_out[,2:3], type = "l",xlab = "steps", ylab = "L2 loss")
    #legend("topright", legend = c("random","active"),lty = 1:2,col = 1:2)
    
  }

}

