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
library(MASS)
library(CholWishart)
library(FNN)
#source("./Code/Graph_generator.R")
source("./utils.R")


k <- 10
p <- k

prior_cert_level <- c(2,10)
scaling <- 1e10

step_size <- 50
n_steps <- 40
n_reps <- 1
n_init <- 200

args[3] <- as.character(floor(as.numeric(args[2])/100)+1)

B <- diag(k)
set.seed(42)
B0 <- B  
G <- do.call(paste0("g_model",args[3]), list(k=k))
base_dir <- paste0("./res_",args[1],"-",args[2],"-", args[4])
file_names_base <- paste0( base_dir,"/Model")
file_names_base <- paste0(file_names_base,args[3],"_KLdiv_")

Omega <- G$Omega
Sigma <- G$Sigma

# object with results
cglasso_random_cert <- data.frame(matrix(NA, nrow = n_reps, ncol = n_steps+1))
mlasso_random_uncert <- mlasso_random_cert <- cglasso_random_uncert <- cglasso_random_cert

mlasso_random_uncert_stein <- mlasso_random_cert_stein <- cglasso_random_uncert_stein <- cglasso_random_cert_stein <- cglasso_random_cert


### place to save results

if(as.numeric(args[4])==1){
  prior_class <- c("cglasso_spc_trt_", "mlasso_spc_trt_")
}else{
  prior_class <- c("cglasso_rand_", "mlasso_rand_")
}
simu_designs <- as.matrix( expand.grid(prior_class,prior_cert_level))
design_name <- sub(" ","", apply( as.matrix(simu_designs),1,paste0,collapse = ""))
file_names_last <- paste0("_k",k,".csv" )
all_file_name <- paste0(file_names_base,design_name,file_names_last)
all_file_name_stein <- sub("KLdiv", "stein", all_file_name)


### some useful prior samples
prior_uncert <- rglasso(1000, prior_cert_level[1], k = k, factoring = 1e7)
prior_cert <- rglasso(1000, prior_cert_level[2], k = k, factoring = 1e7)


### start the simulation 
set.seed(as.numeric(args[1])+floor(as.numeric(args[2])))
for(i_rep in 1:n_reps){
  ## initial samples
  n <- n_init
  
  if(as.numeric(args[4])==1){ # 1 for spc_trt
    X_init <- lapply(1:(floor(n_init/k)),function(i,k){return(diag(k))},k)
    X_init <- 3*Reduce(rbind, X_init)
    Y_init <- simu_data(X_init,B,Sigma)
    Y_init <- simu_data(X_init,B,Sigma)
    Y_no_exp <- simu_data(0 * X_init,B,Sigma)
    
  }else{
    X_init <- matrix( sign( runif(n * p, -1,1 )) , n , p )
    Y_init <- simu_data(X_init,B,Sigma)
    Y_no_exp <- simu_data(0 * X_init,B,Sigma)
  }
  
  cat("Round:", i_rep,"Initial estimation\n")
  
  #### cg uncert
  cglasso_uncert <- CARlasso:::CAR_LASSO_Cpp(Y_init, X_init, 1000, 1000, 1, 
                                         scaling * prior_cert_level[1]^2, scaling, 
                                         scaling * prior_cert_level[1], scaling,
                                         TRUE)
  cglasso_uncert_point <- CARlasso:::get_graph(cglasso_uncert, k)
  
  
  #### cg uncert w/o exp
  cglasso_uncert_noexp <- CARlasso:::CAR_LASSO_Cpp(Y_no_exp, 0 * X_init, 1000, 1000, 1, 
                                           scaling * prior_cert_level[1]^2, scaling, 
                                           scaling * prior_cert_level[1], scaling, 
                                           TRUE)
  cglasso_uncert_noexp_point <- CARlasso:::get_graph(cglasso_uncert_noexp, k)
  
  
  #### cg cert 
  cglasso_cert <- CARlasso:::CAR_LASSO_Cpp(Y_init, X_init, 1000, 1000, 1, 
                                           scaling * prior_cert_level[2]^2, scaling, 
                                           scaling * prior_cert_level[2], scaling, 
                                           TRUE)
  cglasso_cert_point <- CARlasso:::get_graph(cglasso_cert, k)
  
  
  #### cg cert w/o exp
  cglasso_cert_noexp <- CARlasso:::CAR_LASSO_Cpp(Y_no_exp, 0 * X_init, 1000, 1000, 1, 
                                                 scaling * prior_cert_level[2]^2, scaling, 
                                                 scaling * prior_cert_level[2], scaling, 
                                                 TRUE)
  cglasso_cert_noexp_point <- CARlasso:::get_graph(cglasso_cert_noexp, k)
  
  
  ############ srgs ###################
  #### mlasso uncert 
  mlasso_uncert <- CARlasso:::SRG_LASSO_Cpp(Y_init, X_init, 1000, 1000, 1, 
                                             scaling * prior_cert_level[1]^2, scaling, 
                                             scaling * prior_cert_level[1], scaling,
                                             TRUE)
  mlasso_uncert_point <- CARlasso:::get_graph(mlasso_uncert, k)
  
  
  
  #### mlasso uncert w/o exp
  mlasso_uncert_noexp <- CARlasso:::SRG_LASSO_Cpp(Y_no_exp, 0 * X_init, 1000, 1000, 1, 
                                                   scaling * prior_cert_level[1]^2, scaling, 
                                                   scaling * prior_cert_level[1], scaling, 
                                                   TRUE)
  mlasso_uncert_noexp_point <- CARlasso:::get_graph(mlasso_uncert_noexp, k)
  
  
  #### mlasso, cert w/ exp
  mlasso_cert <- CARlasso:::SRG_LASSO_Cpp(Y_init, X_init, 1000, 1000, 1, 
                                           scaling * prior_cert_level[2]^2, scaling, 
                                           scaling * prior_cert_level[2], scaling, 
                                           TRUE)
  mlasso_cert_point <- CARlasso:::get_graph(mlasso_cert, k)
  
  
  #### mlasso, certain w/o exp
  mlasso_cert_noexp <- CARlasso:::SRG_LASSO_Cpp(Y_no_exp, 0 * X_init, 1000, 1000, 1, 
                                                 scaling * prior_cert_level[2]^2, scaling, 
                                                 scaling * prior_cert_level[2], scaling, 
                                                 TRUE)
  mlasso_cert_noexp_point <- CARlasso:::get_graph(mlasso_cert_noexp, k)
  
  
  Y_r <- Y_cglasso_c  <- Y_cglasso_u <- Y_init
  X_r <- X_cglasso_c  <- X_cglasso_u <- X_init
  
  
  cglasso_random_cert[i_rep,1] <- mean( FNN::KL.divergence(cglasso_cert$Omega, prior_cert) / 
                                  FNN::KL.divergence(cglasso_cert_noexp$Omega, prior_cert))
  cglasso_random_uncert[i_rep,1] <- mean( FNN::KL.divergence(cglasso_uncert$Omega, prior_uncert) / 
                                            FNN::KL.divergence(cglasso_uncert_noexp$Omega, prior_uncert))
  
  cglasso_random_cert_stein[i_rep, 1] <- stein_loss(cglasso_cert_point, Omega)/stein_loss(cglasso_cert_noexp_point, Omega)
  cglasso_random_uncert_stein[i_rep, 1] <- stein_loss(cglasso_uncert_point, Omega)/stein_loss(cglasso_uncert_noexp_point, Omega)
  
  
  mlasso_random_cert[i_rep,1]  <- mean( FNN::KL.divergence(mlasso_cert$Omega, prior_cert) / 
                                          FNN::KL.divergence(mlasso_cert_noexp$Omega, prior_cert))
  
  
  mlasso_random_uncert[i_rep,1] <- mean( FNN::KL.divergence(mlasso_uncert$Omega, prior_uncert) / 
                                           FNN::KL.divergence(mlasso_uncert_noexp$Omega, prior_uncert))
  
  mlasso_random_cert_stein[i_rep, 1] <- stein_loss(mlasso_cert_point, Omega)/stein_loss(mlasso_cert_noexp_point, Omega)
  mlasso_random_uncert_stein[i_rep, 1] <- stein_loss(mlasso_uncert_point, Omega)/stein_loss(mlasso_uncert_noexp_point, Omega)
  
  
  for(i_step in 1:n_steps + 1){
    
    # get data
    n <- n + step_size
   
    
    ## Random, all random share the same design
    
    if(as.numeric(args[4])==1){
      cat("  step :" , i_step, " Specific Design..\n")
      X_rand_temp <- lapply(1:(floor(step_size/k)),function(i,k){return(diag(k))},k)
      X_rand_temp <- Reduce(rbind, X_rand_temp)
      Y_rand_temp <- simu_data(X_rand_temp, B, Sigma)
      Y_nexp_temp <- simu_data(0 * X_rand_temp, B, Sigma)
      
      X_r <- rbind(X_r, X_rand_temp)
      Y_r <- rbind(Y_r, Y_rand_temp)
      Y_no_exp <- rbind(Y_no_exp, Y_nexp_temp)
    }else{
      cat("  step :" , i_step, " Random Design..\n")
      X_rand_temp <- matrix(runif(step_size * p,-1,1),step_size,p)
      Y_rand_temp <- simu_data(X_rand_temp, B, Sigma)
      Y_nexp_temp <- simu_data(0 * X_rand_temp, B, Sigma)
      
      X_r <- rbind(X_r, X_rand_temp)
      Y_r <- rbind(Y_r, Y_rand_temp)
      Y_no_exp <- rbind(Y_no_exp, Y_nexp_temp)
    }
    
    
    
    #### cg uncert
    cglasso_uncert <- CARlasso:::CAR_LASSO_Cpp(Y_r, X_r, 1000, 1000, 1, 
                                               scaling * prior_cert_level[1]^2, scaling, 
                                               scaling * prior_cert_level[1], scaling,
                                               TRUE)
    cglasso_uncert_point <- CARlasso:::get_graph(cglasso_uncert, k)
    
    
    #### cg uncert w/o exp
    cglasso_uncert_noexp <- CARlasso:::CAR_LASSO_Cpp(Y_no_exp, 0 * X_r, 1000, 1000, 1, 
                                                     scaling * prior_cert_level[1]^2, scaling, 
                                                     scaling * prior_cert_level[1], scaling, 
                                                     TRUE)
    cglasso_uncert_noexp_point <- CARlasso:::get_graph(cglasso_uncert_noexp, k)
    
    
    #### cg cert 
    cglasso_cert <- CARlasso:::CAR_LASSO_Cpp(Y_r, X_r, 1000, 1000, 1, 
                                             scaling * prior_cert_level[2]^2, scaling, 
                                             scaling * prior_cert_level[2], scaling, 
                                             TRUE)
    cglasso_cert_point <- CARlasso:::get_graph(cglasso_cert, k)
    
    
    #### cg cert w/o exp
    cglasso_cert_noexp <- CARlasso:::CAR_LASSO_Cpp(Y_no_exp, 0 * X_r, 1000, 1000, 1, 
                                                   scaling * prior_cert_level[2]^2, scaling, 
                                                   scaling * prior_cert_level[2], scaling, 
                                                   TRUE)
    cglasso_cert_noexp_point <- CARlasso:::get_graph(cglasso_cert_noexp, k)
    
    
    ############ srgs ###################
    #### mlasso uncert 
    mlasso_uncert <- CARlasso:::SRG_LASSO_Cpp(Y_r, X_r, 1000, 1000, 1, 
                                              scaling * prior_cert_level[1]^2, scaling, 
                                              scaling * prior_cert_level[1], scaling,
                                              TRUE)
    mlasso_uncert_point <- CARlasso:::get_graph(mlasso_uncert, k)
    
    
    
    #### mlasso uncert w/o exp
    mlasso_uncert_noexp <- CARlasso:::SRG_LASSO_Cpp(Y_no_exp, 0 * X_r, 1000, 1000, 1, 
                                                    scaling * prior_cert_level[1]^2, scaling, 
                                                    scaling * prior_cert_level[1], scaling, 
                                                    TRUE)
    mlasso_uncert_noexp_point <- CARlasso:::get_graph(mlasso_uncert_noexp, k)
    
    
    #### mlasso, cert w/ exp
    mlasso_cert <- CARlasso:::SRG_LASSO_Cpp(Y_r, X_r, 1000, 1000, 1, 
                                            scaling * prior_cert_level[2]^2, scaling, 
                                            scaling * prior_cert_level[2], scaling, 
                                            TRUE)
    mlasso_cert_point <- CARlasso:::get_graph(mlasso_cert, k)
    
    
    #### mlasso, certain w/o exp
    mlasso_cert_noexp <- CARlasso:::SRG_LASSO_Cpp(Y_no_exp, 0 * X_r, 1000, 1000, 1, 
                                                  scaling * prior_cert_level[2]^2, scaling, 
                                                  scaling * prior_cert_level[2], scaling, 
                                                  TRUE)
    mlasso_cert_noexp_point <- CARlasso:::get_graph(mlasso_cert_noexp, k)  
    
    # save
    cglasso_random_cert[i_rep,i_step] <- mean( FNN::KL.divergence(cglasso_cert$Omega, prior_cert) / 
                                            FNN::KL.divergence(cglasso_cert_noexp$Omega, prior_cert))
    cglasso_random_uncert[i_rep,i_step] <- mean( FNN::KL.divergence(cglasso_uncert$Omega, prior_uncert) / 
                                              FNN::KL.divergence(cglasso_uncert_noexp$Omega, prior_uncert))
    
    cglasso_random_cert_stein[i_rep,i_step] <- stein_loss(cglasso_cert_point, Omega)/stein_loss(cglasso_cert_noexp_point, Omega)
    cglasso_random_uncert_stein[i_rep,i_step] <- stein_loss(cglasso_uncert_point, Omega)/stein_loss(cglasso_uncert_noexp_point, Omega)
    
    
    mlasso_random_cert[i_rep,i_step]  <- mean( FNN::KL.divergence(mlasso_cert$Omega, prior_cert) / 
                                            FNN::KL.divergence(mlasso_cert_noexp$Omega, prior_cert))
    
    
    mlasso_random_uncert[i_rep,i_step] <- mean( FNN::KL.divergence(mlasso_uncert$Omega, prior_uncert) / 
                                             FNN::KL.divergence(mlasso_uncert_noexp$Omega, prior_uncert))
    
    mlasso_random_cert_stein[i_rep,i_step] <- stein_loss(mlasso_cert_point, Omega)/stein_loss(mlasso_cert_noexp_point, Omega)
    mlasso_random_uncert_stein[i_rep,i_step] <- stein_loss(mlasso_uncert_point, Omega)/stein_loss(mlasso_uncert_noexp_point, Omega)
    
    
    ### write KL
    write.csv(cglasso_random_uncert,all_file_name[1],row.names = F)
    write.csv(mlasso_random_uncert,all_file_name[2],row.names = F)
    write.csv(cglasso_random_cert,all_file_name[3],row.names = F)
    write.csv(mlasso_random_cert,all_file_name[4],row.names = F)
    ### write stein
    write.csv(cglasso_random_uncert_stein,all_file_name_stein[1],row.names = F)
    write.csv(mlasso_random_uncert_stein,all_file_name_stein[2],row.names = F)
    write.csv(cglasso_random_cert_stein,all_file_name_stein[3],row.names = F)
    write.csv(mlasso_random_cert_stein,all_file_name_stein[4],row.names = F)
    
    
    plot_out <- data.frame(step = 0:n_steps * step_size + n_init, 
                           cgc <- t(as.matrix(cglasso_random_cert[i_rep, ])), 
                           #mac <- t(as.matrix(cglasso_active_cert[i_rep, ])),
                           
                           cgu <- t(as.matrix(cglasso_random_uncert[i_rep, ])), 
                           #mau <- t(as.matrix(cglasso_active_uncert[i_rep, ])), 
                           
                           mc <- t(as.matrix(mlasso_random_cert[i_rep, ])),
                           mu <- t(as.matrix(mlasso_random_uncert[i_rep, ])), 
                           
                           #nex <- t(as.matrix(no_exp[i_rep, ])), 
                           row.names = NULL)
    #matplot((plot_out[,1]), log(plot_out[,2:5]), type = "l",xlab = "sample size", ylab = "log KL to prior ratio",lty = 1:7,col = 1:7)
    #abline(h = 0, lty = 2)
    #legend("topright", legend = c("cg-cert","cg-uncert","m-cert","m-uncert"),lty = 1:4,col = 1:4)
    
  }
  
}
