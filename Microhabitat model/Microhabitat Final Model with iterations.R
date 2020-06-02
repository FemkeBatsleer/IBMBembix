#############
#Loading packages
library(lattice)
library(ggplot2)
library(rgdal)
library(sp)
library(fields)
library(gstat)
library(ggmap)
library(reshape)
library(raster)
require(dismo)
library(INLA)
library(maps)
library(maptools)
library(mapdata)
library(rgeos)
library(GGally)
library(MASS)
library(ROCR)
library(plotly)
library(inlabru)
library(dplyr)
library(tidyr)
library(readxl)
library(rgdal)
library(caret)
library(reshape2)
library(imager)

#######################
####Loading data#######
#######################
#dataset
setwd("D:/Femke/OneDrive - UGent/WP 2 Field study/INLA")
Bembix_raw <- read_excel("../MicrohabitatmodelComplete.xlsx")

#dataset of regular raster-points for predictions
Bembix_reg_raw <- read.table("../RegPoints.txt", header=T)
################################
##### Format data correctly ####
###############################

colnames(Bembix_raw)[4] <- "NestID"
#only get the values on 2m
Bembix_2m <- Bembix_raw[,c('Presence', 'X', 'Y', "NestID", "RealNest", "NDVI2mean", "W2mean", "Sl2mean")]
site <- factor(rep("a",nrow(Bembix_2m))) 
Bembix <-  cbind(Bembix_2m, site)
### REAL NESTS
#select real nests, select same amount of random control points, concatenate the df, replace original dataset
Bembix_real_p <- Bembix[Bembix$RealNest==1,]
am_realnests <- nrow(Bembix_real_p)
Bembix_real_c <- sample_n(Bembix[which(Bembix$RealNest==0 & Bembix$Presence==0),], am_realnests)
Bembix_real <- rbind(Bembix_real_p, Bembix_real_c)
Bembix <- Bembix_real
#standardise data
standard_score <- function(x) {return ((x - mean(x))/sd(x))}
Bembix$NDVIstd <- standard_score(Bembix$NDVI2mean)
Bembix$Wstd    <- standard_score(Bembix$W2mean)
Bembix$Slstd   <- standard_score(Bembix$Sl2mean)

#formatting dataset regular points for predictions
Bembix_reg <- data.frame(Presence = rep(NA, nrow(Bembix_reg_raw)),
                          X = Bembix_reg_raw$x,
                          Y = Bembix_reg_raw$y,
                          NDVI2mean = Bembix_reg_raw$NDVI2mean,
                          W2mean = Bembix_reg_raw$W2mean,
                          Sl2mean = Bembix_reg_raw$Sl2mean)
###Standardisation (with the origianl dataset as reference!)
Bembix_reg$NDVIstd <- (Bembix_reg$NDVI2mean - mean(Bembix$NDVI2mean))/sd(Bembix$NDVI2mean)
Bembix_reg$Wstd    <- (Bembix_reg$W2mean-mean(Bembix$W2mean))/sd(Bembix$W2mean)


############################################
#### Making datasets ready to save data ####
############################################
# effect sizes -> per model
effect_sizes <- data.frame("rowname"=factor(), "effect"=factor(), "size"=double(), "model"=factor())
# variogramdata -> separate file
variogramsV0_df <- data.frame("np"=integer(), "dist"=double(), "gamma"=double(),
                              "dir.hor"=integer(), "dir.ver"=integer(), "id"=character(),
                              "model"=factor())
variogramsV1_df <- data.frame("np"=integer(), "dist"=double(), "gamma"=double(),
                              "dir.hor"=integer(), "dir.ver"=integer(), "id"=character(),
                              "model"=factor())
# spatial field? -> models and meshes are saved in separate files
#save(m1, file = "my_model1.rda") if(file.exists("my_model1.rda")) {load("my_model1.rda")
# posterior distributions -> separate file
posteriors_df <- data.frame("model"=factor(), "id"=factor(), "x"=double(), "y"=double())
# predictions for 30% test-data -> add to file with for each prediction-set X, Y, prediction. Each combination of X&Y --> unique point
Bembix_30_pred <- data.frame("Presence"=integer(), "X"=double(), "Y"=double(),
                            "mean"=double(), "selow"=double(), "seup"=double(), "model"=factor())
# AUC -> per model
Bembix_auc <- data.frame(model=character(), AUC=double())
# Predictions suitability map -> add to existing dataset Bembix_reg
Predictions_map_df <- data.frame("X"=double(), "Y"=double(), "mu"=double(), "selow"=double(),
                              "seup"=double(), "model"=double())

#############################
#### iterations of model ####
#############################
#it_numbers <- sample(1:999, 10, replace=F)
it_numbers <- sample(1:999, 10, replace=F)
#it_number <- it_numbers
for (it_number in it_numbers){
  set.seed(it_number) #to make sampling random, but remember it_number, so reproducible
  print(it_number)
  #First subset 70% training data and 30% test data
  Bembix70 <- Bembix[sample(nrow(Bembix), floor(0.7*nrow(Bembix))),]
  Bembix70_compl <- dplyr::setdiff(Bembix, Bembix70)
  
  ###################
  ### make mesh, run model ####
  ##################
  Loc <- cbind(Bembix70$X, Bembix70$Y)
  range <- 10
  MaxEdge <- range/5
  mesh <- inla.mesh.2d(Loc, max.edge = c(1,5)*MaxEdge, cutoff   = MaxEdge/5)
  A <- inla.spde.make.A(mesh, loc = Loc)
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(20, 0.05),prior.sigma = c(1, 0.05))
  w.index <- inla.spde.make.index(name    = 'w', n.spde  = spde$n.spde)
  #make a stack
  N <- nrow(Bembix70)
  X <- data.frame(NDVIstd= Bembix70$NDVIstd,Wstd= Bembix70$Wstd,Slstd= Bembix70$Slstd)
  StackFit<- inla.stack(tag  = "Fit", data = list(y = Bembix70$Presence),  
                        A    = list(1, 1, A),
                        effects = list(Intercept = rep(1, N),
                                       X = as.data.frame(X), #Covariates
                                       w = w.index))         #Spatial field
  
  f <- y ~ -1 + Intercept + NDVIstd + Wstd + f(w, model = spde)
  print("Model running for General model")
  
  begin_time <- Sys.time()
  M <- inla(f, family = "binomial", data = inla.stack.data(StackFit),
            control.compute = list(dic = TRUE, waic = TRUE, config=TRUE),
            control.predictor = list(A = inla.stack.A(StackFit)))
  end_time <- Sys.time()
  print(paste0("Duration was ", end_time-begin_time, " min"))
  
  #save the effect sizes
  post_beta1 <- as.data.frame(M$summary.fixed[, c("mean", "0.025quant", "0.975quant")]) %>% 
    tibble::rownames_to_column() %>%
    gather(key="effect", value="size", mean:'0.975quant')
  
  post_beta1$model <- rep(it_number, nrow(post_beta1))
  effect_sizes <- rbind(effect_sizes, post_beta1)
  
  #save model in rda
  save(M, file = paste0("./output_iterations/M", it_number, "-", Sys.Date(), ".rda"))
  #save mesh
  save(mesh, file= paste0("./output_iterations/mesh", it_number, "-", Sys.Date(), ".rda"))
  
  #############################
  ### making variograms ####
  ############################
  #model without spatial autocorrelation taken into account
  f0 <- y ~ -1 + Intercept + NDVIstd + Wstd
  
  print("Model running without spatial autocorrelation")
  begin_time <- Sys.time()
  M0 <- inla(f0, family = "binomial", data = inla.stack.data(StackFit),
             control.compute = list(dic = TRUE, waic = TRUE, config=TRUE),
             control.predictor = list(A = inla.stack.A(StackFit)))
  end_time <- Sys.time()
  print(paste0("Duration was ", end_time-begin_time, " min"))
  
  #Calculate pearson's residuals
  #M0-model
  Pi0 <- M0$summary.fitted.values[1:nrow(Bembix70), "mean"]
  ExpY0 <- Pi0
  VarY0 <- (1 - Pi0) * Pi0
  E0 <- (Bembix70$Presence - ExpY0) / sqrt(VarY0)
  #M-model
  Pi <- M$summary.fitted.values[1:nrow(Bembix70), "mean"]
  ExpY <- Pi
  VarY <- (1 - Pi) * Pi
  E1 <- (Bembix70$Presence - ExpY) / sqrt(VarY)
  
  data_variogram <- data.frame(E0 = E0, E1 = E1, Presence = Bembix70$Presence,
                               X = Bembix70$X, Y = Bembix70$Y)
  coordinates(data_variogram) <- c("X", "Y")
  #define variograms
  V0 <- variogram(E0 ~ X + Y, data = data_variogram, cressie = TRUE, cutoff = 50)
  V1 <- variogram(E1 ~ X + Y, data = data_variogram, cressie = TRUE, cutoff = 50)
  
  V0$model <- rep(it_number, nrow(V0))
  V1$model <- rep(it_number, nrow(V1))
  
  #add to dataset
  variogramsV0_df <- rbind(variogramsV0_df, V0)
  variogramsV1_df <- rbind(variogramsV1_df, V1)

  ########################################
  #### posterior distributions###############
  #########################################
  
  post_int <- as.data.frame(M$marginals.fixed$Intercept)
  post_ndvi <- as.data.frame(M$marginals.fixed$NDVIstd)
  post_insol <- as.data.frame(M$marginals.fixed$Wstd)
  posteriors <- bind_rows(list(Intercept=post_int, NDVI=post_ndvi, Insolation=post_insol), .id="id")
  posteriors$model <- rep(it_number, nrow(posteriors))
  
  posteriors_df <- rbind(posteriors_df, posteriors)
  
  ######################################################
  ### Calculate predictions for the 30% test data ####
  ####################################################
  #this is the Bembix70_compl dataset
  Xmm_td <- model.matrix(~ NDVIstd + Wstd,data = Bembix70_compl)
  Xp_td <- data.frame(NDVIstd  = Xmm_td[,2], Wstd = Xmm_td[,3])
  #make a stack for these predictions we want to make
  StackCov_td <- inla.stack(tag = "Covariates",data = list(y = NA),
                            A = list(1, 1), effects = list(Intercept = rep(1, nrow(Xp_td)), Xp_td = Xp_td))
  # This stack does not have a spatial field
  # because we, only want to have the predictions for the covariates and the intercept.
  #combine original and stack from predicted values, coarse and fine
  All.stacks_td <- inla.stack(StackFit, StackCov_td)
  #rerun inla with new stack
  
  print("running model for predictions test data")
  begin_time <- Sys.time()
  Pred_td <- inla(f,family = "binomial", data = inla.stack.data(All.stacks_td),
                  control.predictor=list(link=1, A = inla.stack.A(All.stacks_td)))
  end_time <- Sys.time()
  print(paste0("Duration was ", end_time-begin_time, " min"))
  
  #extract the indexes and correct rows, rename them
  index.Cov_td <- inla.stack.index(All.stacks_td, tag = "Covariates")$data
  Pred.pred_td <- Pred_td$summary.fitted.values[index.Cov_td, c(1,3,5)]
  Pred.pred_td <- plyr::rename(Pred.pred_td, c("mean" = "mu",
                                               "0.025quant" = "selow",
                                               "0.975quant" = "seup"))
  Bembix_pred_td <- cbind("Presence"=Bembix70_compl$Presence, "X"=Bembix70_compl$X,
                          "Y"=Bembix70_compl$Y, Pred.pred_td)
  Bembix_pred_td$model <- rep(it_number, nrow(Bembix_pred_td))

  #add predictions to dataset with all predictions of other models
  Bembix_30_pred <- rbind(Bembix_30_pred, Bembix_pred_td)
  
  ###################################################################
  ### Calculate AUC etc for these predictions of the test data###
  ################################################################
  auc <- performance(prediction(Bembix_pred_td$mu, Bembix_pred_td$Presence),
                     measure = "auc")@y.values[[1]]
  
  Bembix_auc <- rbind(Bembix_auc, data.frame(model=it_number, AUC=auc))
  
  #confusion matrix can be made with prediction data afterwards
  
  ##########################################
  ##### Predictions for suitability map ######
  ############################################
  ##Now lets predict some values!
  #make a design matrix
  Xmm_reg <- model.matrix(~ NDVIstd + Wstd, data = Bembix_reg)
  #head(Xmm_reg)
  Xp_reg <- data.frame(NDVIstd = Xmm_reg[,2], Wstd = Xmm_reg[,3])
  #make a stack for these predictions we want to make
  StackCov_reg <- inla.stack(tag = "Covariates", data = list(y = NA),
                             A = list(1, 1), effects = list(Intercept = rep(1, nrow(Xp_reg)), Xp = Xp_reg))
  # This stack does not have a spatial field, because we
  # only want to have the predictions for the covariates and the intercept
  #combine original and stack from predicted values, coarse and fine
  All.stacks_reg <- inla.stack(StackFit, StackCov_reg)
  #rerun inla with new stack
  #run model
  
  print("Running model for predictions suitability map")
  begin_time <- Sys.time()
  Pred_reg <- inla(f,family = "binomial", data = inla.stack.data(All.stacks_reg),
                      control.predictor=list(link=1, A= inla.stack.A(All.stacks_reg)))
  end_time <- Sys.time()
  print(paste0("Duration was ", end_time-begin_time, " min"))
  #extract the data
  index.Cov_reg <- inla.stack.index(All.stacks_reg,tag = "Covariates")$data
  Pred.pred_reg <- Pred_reg$summary.fitted.values[index.Cov_reg, c(1,3,5)]#predictions for the regular grid
  Pred.pred_reg <- plyr::rename(Pred.pred_reg, c("mean" = "mu",
                                                 "0.025quant" = "selow",
                                                 "0.975quant" = "seup"))
  Pred.pred_reg$model <- rep(it_number, nrow(Pred.pred_reg))
  
  #add predictions to Bembix_reg dataset
  Bembix_reg_pred <- cbind("X"=Bembix_reg$X, "Y"=Bembix_reg$Y, Pred.pred_reg)
  
  Predictions_map_df <- rbind(Predictions_map_df, Bembix_reg_pred)
  ###################################
}

###################################
#######Saving data#################
###################################

# effect sizes -> per model
write.table(effect_sizes, "./output_iterations/effect_sizes.txt", sep="\t")
# variogramdata -> separate file
write.table(variogramsV0_df, "./output_iterations/variogramsV0.txt", sep="\t")
write.table(variogramsV1_df, "./output_iterations/variogramsV1.txt", sep="\t")
# spatial field -> models and meshes are saved during iterations
# posterior distributions -> separate file
write.table(posteriors_df, "./output_iterations/posterior_distributions.txt", sep="\t")
# predictions for 30% test-data -> add to file with for each prediction-set X, Y, prediction. Each combination of X&Y --> unique point
write.table(Bembix_30_pred, "./output_iterations/predictions_testdata.txt", sep="\t")
#AUC
write.table(Bembix_auc, "./output_iterations/auc.txt", sep="\t")
#suitability map predictions
write.table(Predictions_map_df, "./output_iterations/predictions_reg.txt", sep="\t")
