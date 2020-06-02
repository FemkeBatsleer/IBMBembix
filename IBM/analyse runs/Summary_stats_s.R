start_time <- Sys.time()
#defining scenarios to analyse; get it from input
##random, bimodal, bottom-up, CF, conspec-attraction, prev-exp-BU
args <- commandArgs(TRUE)
scenarios <- as.character(args[1])
file_path <- as.character(args[2])
# scenarios <- 'bimodal'
# file_path <- 'D:/fbatslee/OneDrive - UGent/Modelling Thesis/Simple Bembix model'

#setwd("D:/fbatslee/OneDrive - UGent/Modelling Thesis")
#load Spatstat
library(spatial)
library(spatstat)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(igraph)
library(gridExtra)
library(readr)

#making raster of focal field
x_pol = c(20.8,101.6,94, 82.2,79, 79,76,62,57,49.8,46.2,46.2,41.2,38.6,35,32.6,27,4.6)
y_pol = c(3,32.6,50.4,52.2,52.2,47,47.6,49,49,45.6,45.6,46,47,47.6,52.2,52.2,50.8,30)
raster=owin(poly=list(x=x_pol, y=y_pol))


##SPATIAL PATTERN ANALYSIS##
#Initialising the dataframes to put in the values of spatial pattern analysis
#Mean and sd's densities of the runs
#densities <- tibble(scenario=factor(), mean=numeric(0), sd=numeric(0))
#Value of Ripley's K at different r (0, 2, 5, 10, 20, 30 m)
RKs <- tibble(scenario=factor(),file_path = factor(),  '2'=numeric(0), '5'=numeric(0),
              '10'=numeric(0),'15'=numeric(0), '20'=numeric(0), '30'=numeric(0), '40'=numeric(0))

#mean, sd and value of mark correlation function
MCs <- tibble(scenario=factor(), '0'=numeric(0),'2'=numeric(0), '4' = numeric(0),
              '6'=numeric(0), '8'=numeric(0),'10'=numeric(0), '12'=numeric(0), '14' = numeric(0))

source("./Analysis/Simple model/H_Analyse Point Pattern_s.R")

for (scen in scenarios) {
  print(scen)
  p <- analyse_ppp(scen, file_path)
  RKs <- bind_rows(RKs, p[1])
  MCs <- bind_rows(MCs, p[2])
}

##NETWORK ANALYSIS##
#initialise tibbles for network statistics
network_metrics <- tibble(scenario=character(), internal_loops=numeric(0), all_loops=numeric(0),
                          dens_undirected=numeric(0), dens_directed=numeric(0), reciproc=numeric(0),
                          transitivity_und=numeric(0),transitivity_dir=numeric(0))

set.seed(5)
number_of_clusters = 11
source("./Analysis/Simple model/H_Analyse Network_s.R")

for (scen in scenarios){
    print(scen)
    netw <- analyse_network(scen, file_path)
    network_metrics <- bind_rows(network_metrics, netw)
    
}

##OUTPUT SUMMARY STATISTICS##
RKs
summary_stats <- cbind(RKs[1:ncol(RKs)], MCs[2:ncol(MCs)], network_metrics[2:ncol(network_metrics)])
write.table(summary_stats, sprintf("./data/Outputs/stats/Summary_stats_%s_%s.txt", scenarios, Sys.Date()), sep="\t")

end_time = Sys.time()

print(end_time-start_time)
