library(rMGIG)
library(MASS)
library(CholWishart)
source("./Code/Graph_generator.R")
source("./Code/utils.R")


n <- 200
p <- 1
k <- 3

B0 <- matrix(0,p,k)
B0[p,k] <- 2
B <- B0

G <- g_model1(k)
Omega <- G$Omega
Sigma <- G$Sigma



psi <- .01 * diag(k)
phi <- .01 * diag(k)
prior_cert <- c(1e-5,1e5)

Lambda1 <- prior_cert[1] * diag(p) # less uncertainty
Lambda2 <- prior_cert[2] * diag(p) # more uncertainty 

lambda <- k+1


set.seed(12345)
X <- matrix( sign( runif(n * p, -1,1 )) , n , p )
Y <- simu_data(X,B,Sigma)
Y_no_exp <- simu_data(0 * X,B,Sigma)

posterior_para_Omega_wrc <- get_posterior_para_Omega_Wishart(Y, X, B0 %*% Sigma, phi, 2*lambda, Lambda1)
posterior_para_Omega_wrc_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X, B0 %*% Sigma, phi, 2*lambda, Lambda1)
posterior_para_Omega_wru <- get_posterior_para_Omega_Wishart(Y, X, B0 %*% Sigma, phi, 2*lambda, Lambda2)
posterior_para_Omega_wru_noexp <- get_posterior_para_Omega_Wishart(Y_no_exp, 0*X, B0 %*% Sigma, phi, 2*lambda, Lambda2)

sample_wishart_wrc <- rWishart(5e3, posterior_para_Omega_wrc$lambda, solve(posterior_para_Omega_wrc$phi))[1,2,]
sample_wishart_wrc_noexp <- rWishart(5e3, posterior_para_Omega_wrc_noexp$lambda, solve(posterior_para_Omega_wrc_noexp$phi))[1,2,]
sample_wishart_wru <- rWishart(5e3, posterior_para_Omega_wru$lambda, solve(posterior_para_Omega_wru$phi))[1,2,]
sample_wishart_wru_noexp <- rWishart(5e3, posterior_para_Omega_wru_noexp$lambda, solve(posterior_para_Omega_wru_noexp$phi))[1,2,]




posterior_para_Omega_mrc <- get_posterior_para_Omega_MGIG(Y, X, B0, phi, psi,lambda, Lambda1)
posterior_para_Omega_mrc_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X, B0, phi, psi,lambda, Lambda1)
posterior_para_Omega_mru <- get_posterior_para_Omega_MGIG(Y, X, B0, phi, psi,lambda, Lambda2)
posterior_para_Omega_mru_noexp <- get_posterior_para_Omega_MGIG(Y_no_exp, 0*X, B0, phi, psi,lambda, Lambda2)

sample_mgig_mrc <- rMGIG::rMGIG(n=5e3, posterior_para_Omega_mrc$nu,
                                posterior_para_Omega_mrc$phi,
                                posterior_para_Omega_mrc$psi, 
                                df = 5000, maxit = 1e5)|>
  sapply(`[`,2)
sample_mgig_mrc_noexp <- rMGIG::rMGIG(n=5e3, posterior_para_Omega_mrc_noexp$nu,
                                      posterior_para_Omega_mrc_noexp$phi,
                                      posterior_para_Omega_mrc_noexp$psi, 
                                df = 1000, maxit = 1e5) |>
  sapply(`[`,2)
sample_mgig_mru <- rMGIG::rMGIG(n=5e3, posterior_para_Omega_mru$nu,
                                posterior_para_Omega_mru$phi,
                                posterior_para_Omega_mru$psi, 
                                df = 1000, maxit = 1e5)|>
  sapply(`[`,2)
sample_mgig_mru_noexp <- rMGIG::rMGIG(n=5e3, posterior_para_Omega_mru_noexp$nu,
                                posterior_para_Omega_mru_noexp$phi,
                                posterior_para_Omega_mru_noexp$psi, 
                                df = 1000, maxit = 1e5)|>
  sapply(`[`,2)


png("toy-posterior.png", width = 8, height = 4, units = "in", res = 500)
par(mar = c(1.5,1.5,1.5,1.5), mgp = c(1.8, 0.5, 0))
par(mfrow = c(2,2))
plot(density(sample_mgig_mrc),xlim = c(-2.5,-0.5),xlab = "",main = "MGIG-certain")
abline(v = Omega[1,2])
lines(density(sample_mgig_mrc_noexp), col = "red", lty = 2)
legend("topleft", legend = c("experiment","none"), 
       lty = c(1,2), 
       col = c("black","red"))


plot(density(sample_mgig_mru),xlim = c(-2.5,-0.5),xlab = "",main = "MGIG-uncertain")
abline(v = Omega[1,2])
lines(density(sample_mgig_mru_noexp), col = "red", lty = 2)


plot(density(sample_wishart_wrc),xlim = c(-2.5,-0.5),xlab = "",main = "Wishart-certain")
abline(v = Omega[1,2])
lines(density(sample_wishart_wrc_noexp), col = "red", lty = 2)


plot(density(sample_wishart_wru),xlim = c(-2.5,-0.5),xlab = "",main = "Wishart-uncertain")
abline(v = Omega[1,2])
lines(density(sample_wishart_wru_noexp), col = "red", lty = 2)

dev.off()

