library(ggplot2)
library(reshape2)

res_dir <- "./Res/init_200_stepsize_20_steps_50"
res_files <- list.files(res_dir, pattern = "csv$", full.names = T)
res_df <- lapply(res_files,read.csv)

res_mean <- sapply(res_df, colMeans)
colnames(res_mean) <- c("MGIG-cert", "MGIG-uncert","Wishart-cert","Wishart-uncert")

res_lw <- sapply(res_df, function(w){
  apply(as.matrix(w),2, quantile, probs = 0.025)
})
colnames(res_lw) <- c("MGIG-cert", "MGIG-uncert","Wishart-cert","Wishart-uncert")

res_hi <- sapply(res_df, function(w){
  apply(as.matrix(w),2, quantile, probs = 0.925)
})

colnames(res_hi) <- c("MGIG-cert", "MGIG-uncert","Wishart-cert","Wishart-uncert")

sample_size <- 200 + 0:50 * 20

res_mean <- data.frame(sample_size=sample_size, res_mean)
res_lw <- data.frame(sample_size=sample_size, res_lw)
res_hi <- data.frame(sample_size=sample_size, res_hi)

res_mean_long <- melt(res_mean, id.vars="sample_size", value.name = "mean", variable.name = "prior")
res_lw_long <- melt(res_lw, id.vars="sample_size", value.name = "lw", variable.name = "prior")
res_hi_long <- melt(res_hi, id.vars="sample_size", value.name = "hi", variable.name = "prior")


plot_data <- merge(res_mean_long, res_hi_long)
plot_data <- merge(plot_data, res_lw_long)

ggplot(data=plot_data, aes(x=sample_size, y = mean, color = prior, shape = prior)) + 
  geom_hline(yintercept =  0, lty = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = lw, ymax = hi), alpha = 0.3, width = 5) + 
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) + 
  xlab("Sample size") + 
  ylab("log Stein's loss difference to no experiment")

ggsave("./Res/init_200_stepsize_20_steps_50/log_stein_diff.pdf",
       width = 10, height = 6, scale = 0.8)
