start_time <- Sys.time()
#setwd("D:/fbatslee/OneDrive - UGent/Modelling Thesis")
#defining scenarios to analyse; get it from input
##random, bimodal, bottom-up, CF, conspec-attraction, prev-exp-BU
args <- commandArgs(TRUE)
scenarios <- as.character(args[1])
file_path <- as.character(args[2])

# scenarios <- 'random'
# file_path <- 'D:/fbatslee/OneDrive - UGent/WP 2 Modelling Thesis/HPC/data/Outputs/random/20-11'

#setwd("D:/fbatslee/OneDrive - UGent/Modelling Thesis")
#load Spatstat
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

#open files

params_scen <- tibble(pf=factor(),scenario=factor(),active=factor(), CA_w_env=factor(),
                      CA_day=factor(), response_func=factor(),choosers_on_followers=numeric(),
                      bi_CA=factor(), weight_day=numeric())

#Separate function in a file was created to have a function that puts all parameters in a dataframe/tibble
source("./Analysis/Simple model/H_Data Params_s.R")
#looping over the different scenarios and putting the parameters values into one dataframe
#function to analyse spatial pattern, returns list with 3 tibbles: densities, ripley'sK values, mark corr values

for (scen in scenarios) {
  print(scen)
  b <- data_params(scen, file_path)
  params_scen <- rbind(params_scen, b)
  }


#making output file

#omit columns which contain NA values
params_scen = params_scen[,colSums(is.na(params_scen)) == 0]

#output table
write.table(params_scen, sprintf("./data/Outputs/stats/Params_%s_%s.txt", scenarios, Sys.Date()), sep="\t")

end_time <- Sys.time()
print(end_time-start_time)