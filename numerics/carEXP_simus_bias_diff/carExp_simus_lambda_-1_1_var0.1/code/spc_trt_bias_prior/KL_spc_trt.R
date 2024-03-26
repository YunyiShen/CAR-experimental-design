g_model1 <- function(k, rho=.7, tol = 1e-10){
  temp <- matrix(rep(1:k,k),ncol = k)
  Sigma <- rho ^ (abs(temp-t(temp)))
  Omega <- solve(Sigma)
  Omega <- Omega * (abs(Omega)>tol)
  return(list(Sigma = Sigma, Omega = Omega))
}


g_model2 <- function(k, rho = .5){
  temp <- matrix(rep(1:k,k),ncol = k)
  Omega <- rho ^ (abs(temp-t(temp))) * (abs(temp-t(temp)) <= 2)
  return(list(Omega = Omega, Sigma = solve(Omega)))
}

g_model3 <- function(k, rho = .5){
  row_ind <- matrix(rep(1:k,k),ncol = k)
  col_ind <- t(row_ind)
  
  Sigma <- diag(1,k,k)
  
  Sigma[(1<=row_ind ) &
          (1<=col_ind ) &
          (row_ind != col_ind) & 
          (row_ind <= k/2) & 
          (col_ind <= k/2)] <- rho
  Sigma[((1+k/2)<=row_ind ) & 
          ((1+k/2)<=col_ind ) &
          (row_ind != col_ind)] <- rho # & 
  #(row_ind <= 10) & 
  #(col_ind <= 10)] 
  Omega <- solve(Sigma)
  
  return(list(Sigma = Sigma, Omega = Omega))
}

g_model4 <- function(k, rho=.1){
  
  Omega <- diag(1,k,k)
  Omega[2:k,1] <- rho
  Omega[1,2:k] <- rho
  Sigma <- solve(Omega)
  return(list(Sigma = Sigma, Omega = Omega))
}

g_model5 <- function(k, rhos = c(2,1,.9)){
  row_ind <- matrix(rep(1:k,k),ncol = k)
  col_ind <- t(row_ind)
  Omega <- diag(rhos[1],k,k)
  Omega[abs(row_ind-col_ind)==1] <- rhos[2]
  Omega[1,k] <- Omega[k,1] <- rhos[3]
  return(list(Sigma = solve(Omega), Omega = Omega))
}

g_model6 <- function(k, rhos = c(2,1)){
  Omega <- matrix(rhos[2], k, k)
  diag(Omega) <- rhos[1]
  
  return(list(Sigma = solve(Omega),Omega = Omega))
}


args <- commandArgs(trailingOnly=TRUE)
library(rMGIG)
library(MASS)
library(CholWishart)
#source("./Graph_generator.R")
source("./utils.R")


k <- 50
p <- k

prior_cert <- c(1e-1,1e1)

step_size <- 50
n_steps <- 40
n_reps <- 1
n_init <- 200

args[3] <- as.character(floor(as.numeric(args[2])/100)+1)
B <- diag(k)
set.seed(42)
B0 <- B  + matrix(rnorm(k*p,0,.1),p,k)
G <- do.call(paste0("g_model",args[3]), list(k=k))
base_dir <- paste0("./res_",args[1],"-",args[2])
file_names_base <- paste0( base_dir,"/Model")
file_names_base <- paste0(file_names_base,args[3],"_B0_lambda_kp1_KLdiv_bias_prior_")

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


### some useful prior samples

sample_mgig <- rMGIG::rMGIG(n=5e3, lambda, phi, psi, df = 2*k, maxit = 1e5)
sample_wishart <- rWishart(5e3, 2*lambda, phi)
sample_wishart <- lapply(1:5e3, function(i){sample_wishart[,,i]})

### start the simulation 
set.seed(as.numeric(args[1])+floor(as.numeric(args[2])))
for(i_rep in 1:n_reps){
  ## initial samples
  n <- n_init
  X_init <- lapply(1:(floor(n_init/step_size)),function(i,k){return(diag(k))},k)
  X_init <- 3*Reduce(rbind, X_init)
  Y_init <- simu_data(X_init,B,Sigma)
  Y_init <- simu_data(X_init,B,Sigma)
  Y_no_exp <- simu_data(0 * X_init,B,Sigma)
  
  cat("Round:", i_rep,"Initial estimation\n")
  
  posterior_para_Omega_wrc <- get_posterior_para_Omega_Wishart(Y_init, X_init, B0 %*% Sigma, phi, 2*lambda, Lambda1)
  posterior_para_Omega_wrc_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, 2*lambda, Lambda1)
  posterior_para_Omega_wru <- get_posterior_para_Omega_Wishart(Y_init, X_init, B0 %*% Sigma, phi, 2*lambda, Lambda2)
  posterior_para_Omega_wru_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, 2*lambda, Lambda2)
  
  
  posterior_para_Omega_mrc <- get_posterior_para_Omega_MGIG(Y_init, X_init, B0, phi, psi,lambda, Lambda1)
  posterior_para_Omega_mrc_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X_init, B0, phi, psi,lambda, Lambda1)
  posterior_para_Omega_mru <- get_posterior_para_Omega_MGIG(Y_init, X_init, B0, phi, psi,lambda, Lambda2)
  posterior_para_Omega_mru_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X_init, B0, phi, psi,lambda, Lambda2)
  

  
  
  #posterior_para_Omega_nex <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X_init, B0 %*% Sigma, phi, lambda, Lambda2)
  #posterior_para_Omega_nex <- solve(cov(Y_no_exp))
  #Lambda_hat_design <- Lambda_hat
  Y_r <- Y_mac  <- Y_mau <- Y_init
  X_r <- X_mac  <- X_mau <- X_init
  
  
  MGIG_random_cert[i_rep,1] <- KLdiv_MGIG2(sample_mgig,lambda, phi, psi,
                                            posterior_para_Omega_mrc$nu,
                                            posterior_para_Omega_mrc$phi,
                                            posterior_para_Omega_mrc$psi
                                            ) /
                                KLdiv_MGIG2(sample_mgig,lambda, phi,psi,
                                            posterior_para_Omega_mrc_noexp$nu,
                                            posterior_para_Omega_mrc_noexp$phi,
                                            posterior_para_Omega_mrc_noexp$psi
                                            )
  MGIG_random_uncert[i_rep,1] <- KLdiv_MGIG2(sample_mgig,lambda, phi, psi,
                                            posterior_para_Omega_mru$nu,
                                            posterior_para_Omega_mru$phi,
                                            posterior_para_Omega_mru$psi
                                            ) /
                                KLdiv_MGIG2(sample_mgig,lambda, phi,psi,
                                            posterior_para_Omega_mru_noexp$nu,
                                            posterior_para_Omega_mru_noexp$phi,
                                            posterior_para_Omega_mru_noexp$psi
                                            )
  
  Wishart_random_cert[i_rep,1]  <- KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wrc$lambda,
                                                posterior_para_Omega_wrc$phi
                                                ) / 
                                    KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wrc_noexp$lambda,
                                                posterior_para_Omega_wrc_noexp$phi
                                                )
  Wishart_random_uncert[i_rep,1] <- KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wru$lambda,
                                                posterior_para_Omega_wru$phi
                                                ) / 
                                    KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wru_noexp$lambda,
                                                posterior_para_Omega_wru_noexp$phi
                                                )
  
  
  
  for(i_step in 1:n_steps + 1){
    
    # get data
    n <- n + step_size
    cat("  step :" , i_step, " specific treatment..\n")
    
    ## Random, all random share the same design
    
    X_rand_temp <- X_rand_temp <- 3*diag(k)
    Y_rand_temp <- simu_data(X_rand_temp, B, Sigma)
    Y_nexp_temp <- simu_data(0 * X_rand_temp, B, Sigma)
    
    X_r <- rbind(X_r, X_rand_temp)
    Y_r <- rbind(Y_r, Y_rand_temp)
    Y_no_exp <- rbind(Y_no_exp, Y_nexp_temp)
    
    # Result for MGIG
    posterior_para_Omega_wrc <- get_posterior_para_Omega_Wishart(Y_r, X_r, B0 %*% Sigma, phi, 2*lambda, Lambda1)
    posterior_para_Omega_wrc_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X_r, B0 %*% Sigma, phi, 2*lambda, Lambda1)
    posterior_para_Omega_wru <- get_posterior_para_Omega_Wishart(Y_r, X_r, B0 %*% Sigma, phi, 2*lambda, Lambda2)
    posterior_para_Omega_wru_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X_r, B0 %*% Sigma, phi, 2*lambda, Lambda2)
    
    # result for Wishart
    posterior_para_Omega_mrc <- get_posterior_para_Omega_MGIG(Y_r, X_r, B0, phi, psi,lambda, Lambda1)
    posterior_para_Omega_mrc_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X_r, B0, phi, psi,lambda, Lambda1)
    posterior_para_Omega_mru <- get_posterior_para_Omega_MGIG(Y_r, X_r, B0, phi, psi,lambda, Lambda2)
    posterior_para_Omega_mru_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X_r, B0, phi, psi,lambda, Lambda2)
    
    
    # save
    MGIG_random_cert[i_rep,i_step] <- KLdiv_MGIG2(sample_mgig,lambda, phi, psi,
                                            posterior_para_Omega_mrc$nu,
                                            posterior_para_Omega_mrc$phi,
                                            posterior_para_Omega_mrc$psi
                                            ) /
                                KLdiv_MGIG2(sample_mgig,lambda, phi,psi,
                                            posterior_para_Omega_mrc_noexp$nu,
                                            posterior_para_Omega_mrc_noexp$phi,
                                            posterior_para_Omega_mrc_noexp$psi
                                            )

    MGIG_random_uncert[i_rep,i_step] <- KLdiv_MGIG2(sample_mgig,lambda, phi, psi,
                                            posterior_para_Omega_mru$nu,
                                            posterior_para_Omega_mru$phi,
                                            posterior_para_Omega_mru$psi
                                            ) /
                                KLdiv_MGIG2(sample_mgig,lambda, phi,psi,
                                            posterior_para_Omega_mru_noexp$nu,
                                            posterior_para_Omega_mru_noexp$phi,
                                            posterior_para_Omega_mru_noexp$psi
                                            )
    
    Wishart_random_cert[i_rep,i_step] <- KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wrc$lambda,
                                                posterior_para_Omega_wrc$phi
                                                ) / 
                                    KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wrc_noexp$lambda,
                                                posterior_para_Omega_wrc_noexp$phi
                                                )

    Wishart_random_uncert[i_rep,i_step] <- KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wru$lambda,
                                                posterior_para_Omega_wru$phi
                                                ) / 
                                    KLdiv_wishart2(sample_wishart, 2*lambda, phi,
                                                posterior_para_Omega_wru_noexp$lambda,
                                                posterior_para_Omega_wru_noexp$phi
                                                )
    
        
    
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
    #matplot((plot_out[,1]), log(plot_out[,2:5]), type = "l",xlab = "sample size", ylab = "log KL to prior ratio",lty = 1:7,col = 1:7)
    #abline(h = 0, lty = 2)
    #legend("topright", legend = c("MGIG-cert","MGIG-uncert","Wishart-cert","Wishart-uncert"),lty = 1:4,col = 1:4)
    
  }
  
}
