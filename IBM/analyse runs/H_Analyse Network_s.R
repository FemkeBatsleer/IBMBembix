analyse_network <- function(scenario_run, file_path){
  #initialise tibbles for looping
  network_metrics <- tibble(scenario=character(), internal_loops=numeric(0), all_loops=numeric(0),
                            dens_undirected=numeric(0), dens_directed=numeric(0), reciproc=numeric(0),
                            transitivity_und=numeric(0),transitivity_dir=numeric(0))
  
  Output_files <- list.files(path=file_path, pattern=sprintf('Output %s *', scenario_run))
  Distances_files <- list.files(path=file_path, pattern=sprintf('Distances %s *', scenario_run))
  
  for (file_title_output in Output_files){
    #print(file_title_output)
    #cluster analysis
    title_full = paste(file_path, file_title_output, sep="/")
    data_cl = read.table(title_full, header=T)
    X<-data_cl$x
    Y<- (274*0.2-data_cl$y)
    day<-data_cl$day
    points <- ppp(X, Y, window=raster, marks=day)
    unitname(points) <- c('metres', 'metres')
    #plot(points, use.marks=F)
    
    kc = kmeans(cbind(X,Y), number_of_clusters, nstart=20)
    #add assigned cluster to the data
    data_cl$cluster = factor(kc$cluster)
    cl_ppp = ppp(X,Y, c(23231, 23333), c(197844, 197905), marks=data_cl$cluster,window=raster)
    #plot(cl_ppp, cex=0.5, chars=rep(1, number_of_clusters),cols=colfunc2(number_of_clusters), legend=F)
    
    #network analysis, read distance file with samen index as the current outputfile
    #format waspID, nest1, nest2, distance
    title_dist_full = paste(file_path, Distances_files[which(file_title_output == Output_files)], sep="/")
    data_net <- read.table(title_dist_full, header=T)
    
    #couple these
    cl_cl <- tibble(from=integer(0), to=integer(0))
    for (i in c(1:nrow(data_net))){
      #get nest from and to
      from_nest = as.character(data_net[i,]$nest1)
      to_nest = as.character(data_net[i,]$nest2)
      #seek corresponding cluster
      #print(from_nest %in% as.character(data_cl$nestID))
      from_cluster = as.integer(data_cl[as.character(data_cl$nestID)==as.character(from_nest),]$cluster)
      to_cluster = data_cl[as.character(data_cl$nestID)==as.character(to_nest),]$cluster
      #add them to the data_frame of cluster to cluster
      cl_cl = rbind(cl_cl, tibble(from=as.integer(from_cluster), to=as.integer(to_cluster)))
    }
    #give the links a weight (to combine multiple links later)
    cl_cl <- cbind(cl_cl, weight = rep(1,nrow(cl_cl)))
    
    #make networks
    net_d <- graph_from_data_frame(d=cl_cl, directed=T)
    net_ud <- graph_from_data_frame(d=cl_cl, directed=F)
    
    #amount loops, interal + external
    all_loops = gsize(net_ud)
    all_loops_d = gsize(net_d)
    
    #simplify by removing internal loops
    net_ud <- simplify(net_ud, remove.loops=T, remove.multiple = F)
    net_d <- simplify(net_d, remove.loops=T, remove.multiple = F)
    ##########################
    #amount of internal loops:
    ###########################
    int_l <- 1 - gsize(net_ud)/all_loops
    
    #simplify by removing multiple loops (sum them in weight)
    net_ud_sum <- simplify(net_ud, remove.loops=T, remove.multiple= T)
    net_d_sum <- simplify(net_d, remove.loops=T, remove.multiple= T)
    
    ###################
    #density/connectance
    ####################
    #self-loops not considered
    dens = edge_density(net_ud_sum, loops=F)
    #directed considered
    dens_d = edge_density(net_d_sum, loops=F)
    
    ##############
    #Reciprocity
    #############
    rp <- reciprocity(net_d_sum)
    
    ##############
    #Transitivity
    ################
    trans_ud = transitivity(net_ud_sum, type="global")
    trans_d = transitivity(net_d_sum, type="global")
    
    
    network_metrics <- add_row(network_metrics,  scenario=as.factor(scenario_run),
                               internal_loops=int_l, all_loops=all_loops,
                               dens_undirected=dens, dens_directed=dens_d,
                               reciproc=rp,
                               transitivity_und=trans_ud, transitivity_dir=trans_d)
  }
  return(network_metrics)
}