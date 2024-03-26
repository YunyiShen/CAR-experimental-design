library(ggplot2)
library(reshape2)
make_base_name_list <- function(bname){
  base_name_wc <- sub( "1e\\+03", " 1e-03", sub("MGIG",
                                           "Wishart", bname) )
  base_name_wu <- sub("1e-03","1e+03", base_name_wc)
  base_name_mc <- sub("Wishart","MGIG",base_name_wc)
  base_name_mu <- sub("Wishart","MGIG",base_name_wu)
  basename_list <- c(base_name_wu, base_name_wc, base_name_mu, base_name_mc)
  basename_pat <- sub("1e\\+03", "1e\\\\+03", basename_list)
  return(list(basename_list = basename_list, basename_pat = basename_pat))
}

##### set up what we want ####
cert_level <- list(cert = "1e-03", uncert = "1e\\+03")
priors <- c("MGIG", "Wishart")

for(mod in 1:6){
for(stein in c(TRUE, FALSE)){
##### get the file names with the results we want ####

exp_des <-c( "rand_exp", "spc_trt")
biased_prior <- paste0(exp_des, "_bias_prior")
if(stein){
  exp_des <- paste0(exp_des, "_stein")
  biased_prior <- paste0(biased_prior, "_stein")
}
exp_des <- paste0("clean_", exp_des)
biased_prior <- paste0("clean_", biased_prior)
mod_pat <- paste0("Model",mod)
mod_name <- c("AR(1)", "AR(2)", "Block", "Star", "Circle", "Full")[mod]

fig_file <- paste0("./", mod_name,ifelse(stein, "_stein", "_KL"),".pdf")
sample_size <- 200 + 0:39 * 50 # sample sizes

#### read in data ####
# regexp for mgigcert
pattern_mgig_cert <- paste0("Model",mod, ".+MGIG*.+", cert_level$cert)
pattern_mgig_uncert <- paste0("Model",mod, ".+MGIG*.+", cert_level$uncert)
pattern_wishart_cert <- paste0("Model",mod, ".+Wishart*.+", cert_level$cert)
pattern_wishart_uncert <- paste0("Model",mod, ".+Wishart*.+", cert_level$uncert)
all_priors <- c(pattern_mgig_cert, 
                   pattern_mgig_uncert,
                   pattern_wishart_cert, 
                   pattern_wishart_uncert)
prior_names <- c("MGIG-cert", "MGIG-uncert", "Wishart-cert", "Wishart-uncert")

all_designs <- c(exp_des, biased_prior)
design_names <- c("random design", "specific design", 
                  "random design-biased prior", "specific design-biased prior")

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
