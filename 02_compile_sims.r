library(tidyverse)

## collect arguments passed in from SLURM job script for number of sims and subdirectory name
args <- commandArgs(trailingOnly=TRUE)
num.sims <- as.numeric(args[1])
subdir_name <- as.character(args[2])

# set working directory
user_home_dir <- " " # enter path for top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")
setwd(dir_path)

# get all results and put in one dataframe
  list_of_all_results <- vector(mode = "list", length = num.sims)
  # read in all results from results directory
  for(i in 1:num.sims){
    print(i)
    list_of_all_results[[i]] <- read.csv(file = paste("results_", i, ".csv", sep = "")
    )[,-1] # remove first column
  }
  
  all_results <- bind_rows(list_of_all_results)
  setwd(user_home_dir)
  write.csv(all_results, file = paste("results_", subdir_name, ".csv", sep = ""))
  




