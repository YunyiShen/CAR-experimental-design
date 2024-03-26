library(ggplot2)
library(reshape2)
make_base_name_list <- function(bname){
  base_name_mc <- sub( "_2_", "_10_", sub("cglasso",
                                           "mlasso", bname) )
  base_name_mu <- sub("_10_","_2_", base_name_mc)
  base_name_cgc <- sub("mlasso","cglasso",base_name_mc)
  base_name_cgu <- sub("mlasso","cglasso",base_name_mu)
  basename_list <- c(base_name_mu, base_name_mc, base_name_cgu, base_name_cgc)
  basename_pat <- basename_list
  return(list(basename_list = basename_list, basename_pat = basename_pat))
}

##### set up what we want ####
cert_level <- list(cert = "_10_", uncert = "_2_")
priors <- c("cglasso", "mlasso")

for(mod in 1:6){
for(stein in c(TRUE, FALSE)){
##### get the file names with the results we want ####

exp_des <-c( "rand_exp", "spc_trt")

if(stein){
  exp_des <- paste0(exp_des, "_stein")
  
}
exp_des <- paste0("clean_", exp_des)

mod_pat <- paste0("Model",mod)
mod_name <- c("AR(1)", "AR(2)", "Block", "Star", "Circle", "Full")[mod]

fig_file <- paste0("./", mod_name,ifelse(stein, "_stein", "_KL"),".pdf")
sample_size <- 200 + 0:39 * 50 # sample sizes

#### read in data ####
# regexp for cglassocert
pattern_cglasso_cert <- paste0("Model",mod, ".+cglasso*.+", cert_level$cert)
pattern_cglasso_uncert <- paste0("Model",mod, ".+cglasso*.+", cert_level$uncert)
pattern_mlasso_cert <- paste0("Model",mod, ".+mlasso*.+", cert_level$cert)
pattern_mlasso_uncert <- paste0("Model",mod, ".+mlasso*.+", cert_level$uncert)
all_priors <- c(pattern_cglasso_cert, 
                   pattern_cglasso_uncert,
                   pattern_mlasso_cert, 
                   pattern_mlasso_uncert)
prior_names <- c("cglasso-cert", "cglasso-uncert", "mlasso-cert", "mlasso-uncert")

all_designs <- c(exp_des)
design_names <- c("random design", "specific design")

# a list of list, 4 designs then four priors, then reduce
all_res <- lapply(1:length(all_designs), 
                  function(i, all_designs, 
                           all_priors,sample_size,
                           prior_names,design_names){
                      tmp <- lapply(1:length(all_priors), 
                             function(j, i, all_designs, 
                                      all_priors, sample_size, 
                                      prior_names,design_names){
                               tmp <- list.files(all_designs[i], all_priors[j], full.names = T) |>
                                 read.csv()
                               data.frame(
                                 lw = apply(as.matrix(tmp),2, quantile, probs = 0.025, na.rm = T),
                                 hi = apply(as.matrix(tmp),2, quantile, probs = 0.975, na.rm = T),
                                 mean = colMeans(as.matrix(tmp), na.rm = T),
                                 sample_size = sample_size,
                                 prior = prior_names[j],
                                 design = design_names[i]
                               )
                             }, i, all_designs, all_priors,sample_size, prior_names,design_names)
                      Reduce(rbind, tmp)
}, all_designs, all_priors,sample_size, prior_names,design_names) |>
  Reduce(f = rbind)

ylabb <- paste( ifelse(stein, "log Stein's loss gain\n to null design", 
                "log KL divergence to prior gain\n to null design"),"in", mod_name, "model" )


ggplot(all_res, aes(x = sample_size, y = log(mean), 
                    color = prior, shape = prior, lty = prior)) + 
  geom_hline(yintercept =  0, lty = 2) + 
  geom_line(size = .75) + 
  geom_errorbar(aes(ymin = log(lw), ymax = log(hi)), alpha = 0.5, width = 5, size = .5) + 
  scale_color_manual(values = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")) + 
  xlab("Sample size") + 
  ylab(ylabb) + 
  facet_wrap(~design, nrow = 2)+
  theme_bw() + 
  theme(legend.position="top") + 
  theme(strip.background =element_rect(fill="white"))

ggsave(fig_file,
       width = 6.2, height = 3.5, scale = 0.9)  

}  
}
