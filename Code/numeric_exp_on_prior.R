library(rMGIG)
library(MASS)
source("./Code/Graph_generator.R")
source("./Code/utils.R")


k <- 50
p <- k

prior_cert <- c(1e-3,1e3)

step_size <- 50
n_steps <- 40
n_reps <- 100
n_init <- 200


set.seed(12345)
B <- diag(k)
B0 <- B  #+ matrix(rnorm(k*p,0,0.2),p,k)
G <- g_model6(k)
file_names_base <- "./Res/init_200_stepsize_50_steps_40_lambda_kp1_stein/Model6_B0_lambda_kp1_"

Omega <- G$Omega
Sigma <- G$Sigma



psi <- .01 * diag(k)
phi <- .01 * diag(k)

Lambda1 <- prior_cert[1] * diag(p) # less uncertainty
Lambda2 <- prior_cert[2] * diag(p) # more uncertainty 

lambda <- k+1

# object with results
MGIG_random_cert <- data.frame(matrix(NA, nrow = n_reps, ncol = n_steps+1))
MGIG_random_uncert <- MGIG_random_cert

Wishart_random_cert <- data.frame(matrix(NA, nrow = n_reps, ncol = n_steps+1))
Wishart_random_uncert <- Wishart_random_cert

### place to save results
prior_class <- c("MGIG_rand_","Wishart_rand_")
simu_designs <- as.matrix( expand.grid(prior_class,prior_cert))
design_name <- apply( as.matrix(simu_designs),1,paste0,collapse = "")
file_names_last <- paste0("_k",k,".csv" )
all_file_name <- paste0(file_names_base,design_name,file_names_last)

### start the simulation 

for(i_rep in 1:n_reps){
  ## initial samples
  n <- n_init
  X_init <- matrix( sign( runif(n * p, -1,1 )) , n , p )
  Y_init <- simu_data(X_init,B,Sigma)
  Y_no_exp <- simu_data(0 * X_init,B,Sigma)
  
  cat("Round:", i_rep,"Initial estimation\n")
  
  MAP_Omega_wrc <- getMAP_Omega_Wishart(Y_init, X_init, B0 %*% Sigma, phi, 2*lambda, Lambda1)
  MAP_Omega_wrc_noexp <- getMAP_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, 2*lambda, Lambda1)
  MAP_Omega_wru <- getMAP_Omega_Wishart(Y_init, X_init, B0 %*% Sigma, phi, 2*lambda, Lambda2)
  MAP_Omega_wru_noexp <- getMAP_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, 2*lambda, Lambda2)
  
  
  MAP_Omega_mrc <- getMAP_Omega_MGIG(Y_init, X_init, B0, phi, psi,lambda, Lambda1)
  MAP_Omega_mrc_noexp <- getMAP_Omega_MGIG(Y_no_exp, 0*X_init, B0, phi, psi,lambda, Lambda1)
  MAP_Omega_mru <- getMAP_Omega_MGIG(Y_init, X_init, B0, phi, psi,lambda, Lambda2)
  MAP_Omega_mru_noexp <- getMAP_Omega_MGIG(Y_no_exp, 0*X_init, B0, phi, psi,lambda, Lambda2)
  

  
  
  #MAP_Omega_nex <- getMAP_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, lambda, Lambda2)
  #MAP_Omega_nex <- solve(cov(Y_no_exp))
  #Lambda_hat_design <- Lambda_hat
  Y_r <- Y_mac  <- Y_mau <- Y_init
  X_r <- X_mac  <- X_mau <- X_init
  
  
  MGIG_random_cert[i_rep,1] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_mrc) ) - 
    log(CARlasso:::stein_loss(Omega,MAP_Omega_mrc_noexp))
  MGIG_random_uncert[i_rep,1] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_mru) ) - 
    log(CARlasso:::stein_loss(Omega,MAP_Omega_mru_noexp))
  
  Wishart_random_cert[i_rep,1]  <- log( CARlasso:::stein_loss(Omega,MAP_Omega_wrc) ) - 
    log( CARlasso:::stein_loss(Omega,MAP_Omega_wrc_noexp) )
  Wishart_random_uncert[i_rep,1] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_wru) ) - 
    log( CARlasso:::stein_loss(Omega,MAP_Omega_wru_noexp) )
  
  
  
  for(i_step in 1:n_steps + 1){
    
    # get data
    n <- n + step_size
    cat("  step :" , i_step, " Random Design..\n")
    
    ## Random, all random share the same design
    
    X_rand_temp <- matrix(runif(step_size * p,-1,1),step_size,p)
    Y_rand_temp <- simu_data(X_rand_temp, B, Sigma)
    Y_nexp_temp <- simu_data(0 * X_rand_temp, B, Sigma)
    
    X_r <- rbind(X_r, X_rand_temp)
    Y_r <- rbind(Y_r, Y_rand_temp)
    Y_no_exp <- rbind(Y_no_exp, Y_nexp_temp)
    
    # Result for MGIG
    MAP_Omega_wrc <- getMAP_Omega_Wishart(Y_r, X_r, B0 %*% Sigma, phi, 2*lambda, Lambda1)
    MAP_Omega_wrc_noexp <- getMAP_Omega_Wishart(Y_no_exp, 0*X_r, B0 %*% Sigma, phi, 2*lambda, Lambda1)
    MAP_Omega_wru <- getMAP_Omega_Wishart(Y_r, X_r, B0 %*% Sigma, phi, 2*lambda, Lambda2)
    MAP_Omega_wru_noexp <- getMAP_Omega_Wishart(Y_no_exp, 0*X_r, B0 %*% Sigma, phi, 2*lambda, Lambda2)
    
    # result for Wishart
    MAP_Omega_mrc <- getMAP_Omega_MGIG(Y_r, X_r, B0, phi, psi,lambda, Lambda1)
    MAP_Omega_mrc_noexp <- getMAP_Omega_MGIG(Y_no_exp, 0*X_r, B0, phi, psi,lambda, Lambda1)
    MAP_Omega_mru <- getMAP_Omega_MGIG(Y_r, X_r, B0, phi, psi,lambda, Lambda2)
    MAP_Omega_mru_noexp <- getMAP_Omega_MGIG(Y_no_exp, 0*X_r, B0, phi, psi,lambda, Lambda2)
    
    
    # save
    MGIG_random_cert[i_rep,i_step] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_mrc) ) - 
      log(CARlasso:::stein_loss(Omega,MAP_Omega_mrc_noexp))
    MGIG_random_uncert[i_rep,i_step] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_mru) ) - 
      log(CARlasso:::stein_loss(Omega,MAP_Omega_mru_noexp))
    
    Wishart_random_cert[i_rep,i_step]  <- log( CARlasso:::stein_loss(Omega,MAP_Omega_wrc) ) - 
      log( CARlasso:::stein_loss(Omega,MAP_Omega_wrc_noexp) )
    Wishart_random_uncert[i_rep,i_step] <- log( CARlasso:::stein_loss(Omega,MAP_Omega_wru) ) - 
      log( CARlasso:::stein_loss(Omega,MAP_Omega_wru_noexp) )
    
    
    ## active learning
    # cat("  step :" , i_step, " D-optimal Design..\n")
    # X_al_mau <- Design_D_opt(MAP_B_mau, solve(MAP_Omega_mau),
    #                          phi, Lambda2, lambda, step_size, X_mau)
    # X_al_mac <-  Design_D_opt(MAP_B_mac, solve(MAP_Omega_mac),
    #                           phi, Lambda1, lambda, step_size, X_mac)
    # 
    # Y_al_mau <- simu_data(X_al_mau, B, Sigma)
    # Y_al_mac <- simu_data(X_al_mac, B, Sigma)
    # 
    # X_mau <- rbind(X_mau, X_al_mau)
    # Y_mau <- rbind(Y_mau, Y_al_mau)
    # 
    # X_mac <- rbind(X_mac, X_al_mac)
    # Y_mac <- rbind(Y_mac, Y_al_mac)
    # 
    # 
    # MAP_Omega_mac <- getMAP_Omega_MGIG(Y_mac, X_mac, B0, phi, psi,lambda, Lambda1)
    # MAP_Omega_mau <- getMAP_Omega_MGIG(Y_mau, X_mac, B0, phi, psi,lambda, Lambda2)
    # 
    # MAP_B_mac <- getMAP_B_MGIG(Y_mac, X_mac, B0, Lambda1, MAP_Omega_mac)
    # MAP_B_mau <- getMAP_B_MGIG(Y_mac, X_mac, B0, Lambda2,MAP_Omega_mau)
    # 
    # 
    # MGIG_active_cert[i_rep,i_step] <- CARlasso:::stein_loss(Omega,MAP_Omega_mac)
    # MGIG_active_uncert[i_rep,i_step] <- CARlasso:::stein_loss(Omega,MAP_Omega_mau)
    
    
    
    #write.csv(MGIG_active_cert,all_file_name[1],row.names = F)
    write.csv(MGIG_random_cert,all_file_name[1],row.names = F)
    write.csv(Wishart_random_cert,all_file_name[2],row.names = F)
    #write.csv(MGIG_active_uncert,all_file_name[4],row.names = F)
    write.csv(MGIG_random_uncert,all_file_name[3],row.names = F)
    write.csv(Wishart_random_uncert,all_file_name[4],row.names = F)
    #write.csv(no_exp,all_file_name[7],row.names = F)
    
    
    plot_out <- data.frame(step = 0:n_steps * step_size + n_init, 
                           mrc <- t(as.matrix(MGIG_random_cert[i_rep, ])), 
                           #mac <- t(as.matrix(MGIG_active_cert[i_rep, ])),
                           
                           mru <- t(as.matrix(MGIG_random_uncert[i_rep, ])), 
                           #mau <- t(as.matrix(MGIG_active_uncert[i_rep, ])), 
                           
                           wrc <- t(as.matrix(Wishart_random_cert[i_rep, ])),
                           wru <- t(as.matrix(Wishart_random_uncert[i_rep, ])), 
                           
                           #nex <- t(as.matrix(no_exp[i_rep, ])), 
                           row.names = NULL)
    matplot(plot_out[,1], plot_out[,2:5], type = "l",xlab = "sample size", ylab = "log Stein's loss ratio to no experiment",lty = 1:7,col = 1:7)
    abline(h = 0, lty = 2)
    #legend("topright", legend = c("MGIG-cert","MGIG-uncert","Wishart-cert","Wishart-uncert"),lty = 1:4,col = 1:4)
    
  }
  
}






