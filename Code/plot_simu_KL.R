library(ggplot2)
library(reshape2)

res_dir <- "./Res/init_200_stepsize_50_steps_40_lambda_kp1_kl"
res_files <- list.files(res_dir, pattern = "csv$", full.names = T)
#res_files <- res_files[-grep("Model5",res_files)]
res_files <- res_files[-grep("spc",res_files)]
res_df <- lapply(res_files,read.csv)

res_mean <- sapply(res_df, colMeans, na.rm = T)
colnames(res_mean) <- c("MGIG-cert-1", "MGIG-uncert-1",
                        "Wishart-cert-1","Wishart-uncert-1",
                        "MGIG-cert-2", "MGIG-uncert-2",
                        "Wishart-cert-2","Wishart-uncert-2",
                        "MGIG-cert-3", "MGIG-uncert-3",
                        "Wishart-cert-3","Wishart-uncert-3",
                        
                        "MGIG-cert-4", "MGIG-uncert-4",
                        "Wishart-cert-4","Wishart-uncert-4",
                        
                        "MGIG-cert-5", "MGIG-uncert-5",
                        "Wishart-cert-5","Wishart-uncert-5",
                        
                        "MGIG-cert-6", "MGIG-uncert-6",
                        "Wishart-cert-6","Wishart-uncert-6"
                        )

res_lw <- sapply(res_df, function(w){
  apply(as.matrix(w),2, quantile, probs = 0.025, na.rm = T)
})
colnames(res_lw) <- colnames(res_mean)

res_hi <- sapply(res_df, function(w){
  apply(as.matrix(w),2, quantile, probs = 0.925, na.rm = T)
})

colnames(res_hi) <- colnames(res_mean)

sample_size <- 200 + 0:39 * 50

res_mean <- data.frame(sample_size=sample_size, res_mean)
res_lw <- data.frame(sample_size=sample_size, res_lw)
res_hi <- data.frame(sample_size=sample_size, res_hi)

res_mean_long <- melt(res_mean, id.vars="sample_size", value.name = "mean", variable.name = "prior")
res_mean_long$model <- paste0("Model",substr(res_mean_long$prior, nchar(as.character(res_mean_long$prior))-1, 
                                             nchar(as.character(res_mean_long$prior))))
res_mean_long$prior <- substr(res_mean_long$prior, 1, 
                              nchar(as.character(res_mean_long$prior))-2)
res_lw_long <- melt(res_lw, id.vars="sample_size", value.name = "lw", variable.name = "prior")
res_lw_long$prior <- res_mean_long$prior
res_lw_long$model <- res_mean_long$model
res_hi_long <- melt(res_hi, id.vars="sample_size", value.name = "hi", variable.name = "prior")
res_hi_long$prior <- res_mean_long$prior
res_hi_long$model <- res_mean_long$model

plot_data <- merge(res_mean_long, res_hi_long)
plot_data <- merge(plot_data, res_lw_long)

ggplot(data=plot_data, aes(x=sample_size, y = log(mean), color = prior, shape = prior)) + 
  geom_hline(yintercept =  0, lty = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = log(lw), ymax = log(hi)), alpha = 0.3, width = 5) + 
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) + 
  xlab("Sample size") + 
  ylab("Log difference in KL divergence between prior and posterior (random vs null experiment)") + 
  facet_wrap(~model, nrow = 3, scales = "free_y")

ggsave("./Res/init_200_stepsize_50_steps_40_lambda_kp1_kl/log_KL_diff2_rand_exp.pdf",
       width = 10, height = 8, scale = 0.8)
