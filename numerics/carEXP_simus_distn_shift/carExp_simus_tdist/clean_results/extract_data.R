split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
make_base_name_list <- function(bname){
  base_name_wc <- sub( "1e\\+03", "1e-03", sub("MGIG",
                      "Wishart", bname) )
  base_name_wu <- sub("1e-03","1e+03", base_name_wc)
  
  base_name_mc <- sub("Wishart","MGIG",base_name_wc)
  
  base_name_mu <- sub("Wishart","MGIG",base_name_wu)
  
  basename_list <- c(base_name_wu, base_name_wc, base_name_mu, base_name_mc)
  basename_pat <- sub("1e\\+03", "1e\\\\+03", basename_list)
  
  return(list(basename_list = basename_list, basename_pat = basename_pat))
}



dirs <- list.dirs(".", recursive = FALSE)
clean_already <- gregexpr("clean", dirs) |> sapply(FUN = function(w){w[1]!=-1})
dirs <- dirs[!clean_already]

for(current_dir in dirs){
  setwd(current_dir)
  base_name <- list.files(".","csv", recursive = T)
  if(length(base_name) <= 10){
    setwd("..")
    next
  }
  base_name <- base_name[1] |>
    split_path()
  base_name <- sub("Model[1,2,3,4,5,6]_","Model1_",base_name[1])
  
  system(paste0("mkdir ../clean_", sub( "./","",current_dir)))
  target_folder <- paste0("../clean_", sub( "./","",current_dir),"/")
  
  basename_list <- make_base_name_list(base_name)
  for(mod in 1:6){
    for(i in 1:4){
      base_name <- basename_list$basename_list[i]
      pat_mat <- basename_list$basename_pat[i]
      out_name <- sub("Model1",paste0("Model",mod), base_name)
      #out_name <- sub("rand","spc_trt_bias_prior", out_name) # only needed for spc_trt_bias_prior setting
      pat_mat <- sub("Model1", paste0("Model",mod), pat_mat)
      #pat_mat <- sub("kp1","kp1_KLdiv_bias_prior", pat_mat)
      all_files <- list.files(path = ".",pattern = pat_mat, recursive = T, full.names = T)
      all_df <- lapply(all_files, read.csv, row.names = 1)
      all_df <- Reduce(rbind,all_df)
      write.csv(all_df, paste0(target_folder,out_name), row.names = F)
    }
  }
  setwd("..")
}
