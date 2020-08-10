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
library(maps)
library(maptools)
library(mapdata)
library(rgeos)
library(GGally)
library(MASS)
library(ROCR)
library(dplyr)
library(tidyr)
library(readxl)
library(rgdal)
library(caret)
library(reshape2)
library(imager)
library(INLA)

setwd("D:/fbatslee/OneDrive - UGent/WP 2 Field study/INLA")

scale <- '1'
######################
#### AUC #############
######################
#Loading data
auc1 <- read.table(paste0("./output_iterations/scale_", scale, "/auc.txt"), header=T)
mean(auc1$AUC)
sd(auc1$AUC)

############################
#### Plotting effect sizes##
############################
effect_sizes <- read.table(paste0("./output_iterations/scale_", scale, "/effect_sizes.txt"), header=T)
effect_sizes %>% group_by(rowname, effect) %>% summarise(avg=mean(size), sd=sd(size))
View(effect_sizes %>% group_by(rowname) %>% spread(effect, size))

posteriors <- read.table(paste0("./output_iterations/scale_", scale, "/posterior_distributions.txt"), header=T) %>%
  group_by(model) %>%  mutate(id = factor(id, levels=c("Intercept", "NDVI", "Insolation")))
posteriors$model <- as.factor(posteriors$model)
p_posteriors <- ggplot(data=posteriors, aes(x=x, y=y, colour=model)) + #color=model
  geom_line()+
  geom_vline(xintercept=0, color="blue")+
  theme(text = element_text(size = 15)) +
  facet_wrap( ~ id, nrow=3)+
  xlab("parameter value") + ylab("density") +
  #xlim(-20,10)+
  theme(legend.position="none")+
  theme_bw()
p_posteriors

####################################
#### Plotting variogram ############
#### for spatial autocorrelation ###
####################################

variogramV0 <- read.table(paste0("./output_iterations/scale_", scale, "/variogramsV0.txt"), header=T)
variogramV1 <- read.table(paste0("./output_iterations/scale_", scale, "/variogramsV1.txt"), header=T)

variogramV0$model <- as.factor(variogramV0$model)
variogramV1$model <- as.factor(variogramV1$model)

p_variogram <- ggplot() + 
  geom_point(data = variogramV0, aes(x = dist, y = gamma, group=model), color="blue") + 
  geom_line(data = variogramV0, aes(x = dist, y = gamma, group=model), color="blue")+
  geom_point(data = variogramV1,aes(x = dist, y = gamma, group=model)) + 
  geom_line(data = variogramV1,aes(x = dist, y = gamma, group=model))+
  scale_x_continuous(limits=c(0,50))+
  xlab("Distance (m)") +
  ylab("Variance") + 
  theme(text = element_text(size = 15))  + 
  theme(legend.position="none") +
  theme_classic()
#ylim(0,1.5)
p_variogram

######################################
#### Confusion matrix ###############
#### for sensitivty and specificity #
#####################################

pred_testdata <- read.table(paste0("./output_iterations/scale_", scale, "/predictions_testdata.txt"), header=T)
thresholds <- pred_testdata %>%  group_by(model) %>% count(Presence) %>%
  filter(Presence==1) %>% mutate(threshold=n/330)
sensitivity <- c()
specificity <- c()
accuracy    <- c()
#how to make confusion matrix for each prediction set?
#threshold_CM <- sum(Bembix_pred_td$Presence==1)/nrow(Bembix_pred_td)
for (i in c(1:nrow(thresholds))){
  model_i <- thresholds[[i,1]]
  threshold_i <- thresholds[[i,4]]
  pred_testdata_filter <- filter(pred_testdata, model==model_i)
  pred_testdata_filter$mupred <- ifelse(pred_testdata_filter$mu>threshold_i, 1, 0)
  CFmatrix <- confusionMatrix(ref=as.factor(pred_testdata_filter$Presence),
                              data=as.factor(pred_testdata_filter$mupred), positive="1")
  sensitivity <- c(sensitivity, CFmatrix$byClass[[1]])
  specificity <- c(specificity, CFmatrix$byClass[[2]])
  accuracy   <- c(accuracy, CFmatrix$byClass[[11]])
  
  print(CFmatrix)
}

c("sensitivity" ,mean(sensitivity), sd(sensitivity))#proportion presences correctly classified
c("specificity", mean(specificity), sd(specificity))#proportion absences correctly classified
c("accurcary", mean(accuracy), sd(accuracy)) #proportions overall correctly classified


#########################################
#### Maps for microhabitat suitaibility #
#### and spatial field ##################
#########################################

#loading dataset for suitability map
pred_reg_raw <- read.table(paste0("./output_iterations/scale_", scale, "/predictions_reg.txt"), header=T)
pred_reg_all <- pred_reg_raw %>% group_by(X,Y) %>% summarise(avg_mu = mean(mu), sd_mu=sd(mu))
Nests <- read.csv("../data/MicrohabitatModel_FullData_04062020.csv") %>%
  filter(RealNest==1)
#load predictions that are within study area
reg_area <- read.csv("./output_iterations/RegPoints_clipStudyArea.csv")
pred_reg <- reg_area %>% select(id, X, Y) %>%
  inner_join(pred_reg_all, by=c("X", "Y"))

#plot average predictions
plot_mu <- ggplot(pred_reg, aes(x = X, y = Y, z = avg_mu)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_raster(aes(fill = avg_mu)) + #, interpolate=TRUE
  theme_classic()+
  coord_fixed(ratio = 1)+
  labs(fill="Suitability")+
  xlab("x (m)")+
  ylab("y (m)")+
  xlim(c(23230, 23335))+
  ylim(c(197850, 197905))+
  scale_fill_distiller(palette="RdBu")+
  geom_point(data = Nests, aes(x = X, y = Y),
             shape = 16, size = 0.2, color = "black", inherit.aes=FALSE)
  #scale_fill_gradient(low="black", high="gray60")+
  #theme(legend.title = element_text("Suitability"))
plot_mu


# par(mfrow=c(1,2))
# plot_mu05
# plot_mu1

#plot of sd of predictions
plot_sd <- ggplot(pred_reg, aes(x = X, y = Y, z = sd_mu)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_raster(aes(fill = sd_mu)) + #, interpolate=TRUE
  theme_void()+
  coord_fixed(ratio = 1)#+
  #scale_fill_gradient(low="black", high="gray60")+
  #theme(legend.position="none")
plot_sd

#plot with nests
plot_sd_nests <- plot_sd +
  geom_point(data = Nests, aes(x = X, y = Y),
             shape = 16, size = 0.2, color = "yellow", inherit.aes=FALSE)
plot_sd_nests

plot_mu_nests <- plot_mu +
  geom_point(data = Nests, aes(x = X, y = Y),
             shape = 16, size = 0.2, color = "yellow", inherit.aes=FALSE)
plot_mu_nests

#load dataset for spatial field
load(paste0("output_iterations/scale_", scale, "/M67-2020-06-07.rda"))
load(paste0("output_iterations/scale_", scale, "/mesh67-2020-06-07.rda"))

#contour
DSN <- "Study_Area_thesis.shp"
boundary <- readOGR(dsn = DSN, layer = "Study_Area_thesis")
x1 <- c(23230, 23230, 23333, 23333) 
y1 <- c(197850, 197904, 197904, 197850) 
# Make a polygon of it.
AreaPoly <- Polygon(cbind(x1, y1), hole = FALSE)
AreaSP   <- SpatialPolygons(list(Polygons(list(AreaPoly), ID = '1')))
AreaSP@proj4string  <- boundary@proj4string
#define what outside the raster
Outraster <- gDifference(AreaSP, boundary)

w.pm <- M$summary.random$w$mean
w.proj <- inla.mesh.projector(mesh)
w.pmf <- inla.mesh.project(w.proj, w.pm)


xygrid <- expand.grid(w.proj$x, w.proj$y)
Data3D <- data.frame(x = xygrid[,1],
                     y = xygrid[,2],
                     z = melt(w.pmf)[,3])
names(Data3D) <- c("x", "y", "z")
Data3D <- Data3D %>% filter((x>=23231 & x<=23332) & (y>=197851 & y<=197902))

p_sp2 <- ggplot(Data3D, 
            aes(x, y, z = z),
            col = rgb(1, 0.5, 0.5, 0.7)) +
  stat_contour(geom="polygon", aes(fill = ..level..)) +
  geom_raster(aes(fill = z)) +
  labs(fill="w.pm")+
  # xlim(23231,23332)+
  # ylim(197851,197902)+
  geom_polygon(data=fortify(Outraster), aes(x=long, y=lat), fill="white", color="white", inherit.aes=FALSE)+
  #geom_polygon(data=fortify(Outraster), aes(x=lat, y=long, z=NULL), color="black", fill=NA)+
  #stat_contour(geom="polygon", aes(fill = ..level..)) +
  #stat_contour(aes(colour = ..level..))+
  xlab("x (m)")+
  ylab("y (m)")+
  xlim(c(23230, 23335))+
  ylim(c(197850, 197905))+
  scale_fill_distiller(palette="PuOr")+
  theme_classic()+
  coord_fixed(ratio = 1)+
  geom_point(data = Nests, aes(x = X, y = Y),
             shape = 16, size = 0.2, color = "black", inherit.aes=FALSE)
p_sp2

################################
####Saving grayscale for IBM####
################################
plot_g <- ggplot(pred_reg, aes(x = X, y = Y, z = avg_mu)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_raster(aes(fill = avg_mu)) + #, interpolate=TRUE
  theme_void()+
  coord_fixed(ratio = 1)+
  theme(legend.position = "none")
  #labs(fill="Suitability")+
  #xlab("x (m)")+
  #ylab("y (m)")+
  #scale_fill_distiller(palette="RdBu")+
  #geom_point(data = Nests, aes(x = X, y = Y),
  #           shape = 16, size = 0.2, color = "black", inherit.aes=FALSE)
#scale_fill_gradient(low="black", high="gray60")+
#theme(legend.title = element_text("Suitability"))
plot_g
ggsave('output_iterations/suitabilityINLA_10062020v2.png', dpi=50, width=10.5, height=5.4, units="in")
#this saving will result in image with 0.21m/px
z = load.image("output_iterations/suitabilityINLA_10062020v2.png")
y = grayscale(z)
save.image(im=y, "output_iterations/suitabilityINLAv2.png")
