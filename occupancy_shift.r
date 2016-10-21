## =================================
## get the occupancy model results
## Date 17.10.2016
## Author Reto Schmucki
## =================================

## by having more event, I am increasing my "certainty" about the outcome and the probability of at least 1 positive event happening that the grid size I present my results.
## similarly, if I have 1 site with very high prob, versus 3 sites with lower, my chance will not increase, but similarly.
## my 3 sites might also have very different sampling quality and thereby one of them might be very little information compare to 1 single site thoroughly sampled elsewhere.

R

rm(list=ls())

library(raster); library(sp); library(rgdal)

## library(rgdal); library(rgeos); library(SDMTools); library(maptools)

if(Sys.info()[4]=="WLD-3VHP992"){
	dir.create("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
	setwd("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
} else {
	dir.create("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
	setwd("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
}

## GET DATA
voltinism <- c("UNI","MULTI")

for (voltinism in voltinism){

l_f <- list.files(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67"))
species_list <- l_f[-grep(".",l_f,fixed = TRUE)]

for (sp in species_list[12:15]){

species_name <- sp
print(sp)

if(Sys.info()[4]=="WLD-3VHP992"){
	## on PC
	coord_zone <- read.csv(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIcoord_zone.csv"))
	occ_summary <- read.csv(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIoutdaytrimmed.csv"))
	occ_year <- read.csv(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIyear.csv"))
} else {

## on MAC /Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/Audusseau_et_al_revised0303 (002)Reto_edit.docx
	coord_zone <- read.csv(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIcoord_zone.csv"))
	occ_summary <- read.csv(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIoutdaytrimmed.csv"))
	occ_year <- read.csv(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_",voltinism,"_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIyear.csv"))
}


########################################
## Estimating distribution of occupancy
########################################

## set the OCC prob per site/year

occ_value <- occ_summary[substr(occ_summary$X,1,2)=="z[",]

occ_coord <- data.frame()

for (y in occ_year$Yearnumeric){
	y1 <-cbind(occ_value[grep(paste0(y,"]$"),occ_value$X, value=FALSE),],coord_zone,YEAR=occ_year$Year[y])
	occ_coord <- rbind(occ_coord,y1)
}

## build spatial points to overlay and work under the EPSG:3025 projection system

laea.proj <- "+init=epsg:3035"

coordinates(occ_coord) <- ~ longitude+latitude
proj4string(occ_coord) = CRS(laea.proj)

## build raster across Europe with the appropriate resolution for the grid cell, 50x50Km

gc.resolution <- 50000

full_europe <- raster(extent(2763687,5733766,1605279,5075336),crs=crs(laea.proj))
res(full_europe) <- gc.resolution
full_europe[] <- 1:ncell(full_europe)

## locate sampling site within the grid (i.e. allocate grid cell to each site)
point_cell <- extract(full_europe,occ_coord)
occ_coord.df <- as.data.frame(occ_coord)
occ_coord.df$grid_id <- point_cell

## generate X realization of occupancy based on the probability obtained from the mean probability of occupancy estimated from the hierarchical occupancy model.

nbr_iteration <- 1000

trial_iteration <- sapply(occ_coord$mean,function(x) {rbinom(nbr_iteration,1,x)},simplify=TRUE)

## store in a 3D array [iteration,site,year]

nbr_site <- length(unique(occ_coord$site))
nbr_year <- length(unique(occ_year$Yearnumeric))

trial_y <- array(NA,dim=c(nbr_iteration,nbr_site,nbr_year))

	for(y in occ_year$Yearnumeric){
		trial_y[,,y] <- trial_iteration[,which(occ_coord.df$YEAR==occ_year$Year[y])]
	}

## shift estimate
full_europe[] <- NA

## For year 2006-2014
M_y <- matrix(NA,nrow=1000,ncol=9)
SD_y <- matrix(NA,nrow=1000,ncol=9)
q25_y <- matrix(NA,nrow=1000,ncol=9)
q75_y <- matrix(NA,nrow=1000,ncol=9)

for (y in occ_year$Yearnumeric){
	per_gricell_obs <- as.matrix(apply(trial_y[,,y],1,function(x) by(x,point_cell[1:nbr_site],FUN=sum)))
	per_gricell_obs[per_gricell_obs>0] <- 1

	M_ <- c()
	SD_ <- c()
	q25_ <- c()
	q75_ <- c()
 	gc.id <- as.numeric(rownames(per_gricell_obs))

	for (it in 1:nbr_iteration){
		full_europe[gc.id] <- per_gricell_obs[,it]
		t.df <- as.data.frame(rasterToPoints(full_europe))
		M_[it] <- mean(t.df$y[t.df$layer==1]/1000)
		SD_[it] <- sd(t.df$y[t.df$layer==1]/1000)
		q25_[it] <- quantile(t.df$y[t.df$layer==1]/1000,0.25)
		q75_[it] <- quantile(t.df$y[t.df$layer==1]/1000,0.75)
	}

	M_y[,y] <- M_
	SD_y[,y] <- SD_
	q25_y[,y] <- q25_
	q75_y[,y] <- q75_
}

save.image(paste0(species_name,".RData"))
}
}
## Estimate linear velocity of northward shift (i.e Km/year) of:

R

library(raster); library(sp); library(rgdal)

## library(rgdal); library(rgeos); library(SDMTools); library(maptools)

if(Sys.info()[4]=="WLD-3VHP992"){
	setwd("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
} else {
	setwd("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/occ_shift_result")
}

shift_table <-data.frame()
l_Rdata <- list.files()

for (i in l_Rdata){

load(i)
print(species_name)

y <- 1
## 1. Core
a <- cbind(c(M_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a[,1]~a[,2])
## hist(a[a[,2]==occ_year$Year[y],1])
shift_a <- coefficients(lm(a[,1]~a[,2]))[2]

## 2. 25th percentile (southern tail)
a25 <- cbind(c(q25_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a25[,1]~a25[,2])
## hist(a25[a25[,2]==occ_year$Year[y],1])
shift_a25 <- coefficients(lm(a25[,1]~a25[,2]))[2]

## 3. 75th percentile (northern front)
a75 <- cbind(c(q75_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a75[,1]~a75[,2])
## hist(a75[a75[,2]==occ_year$Year[y],1])
shift_a75 <- coefficients(lm(a75[,1]~a75[,2]))[2]

## 4. expansion of range (1.67*StandardDeviation)
b <- cbind(c(SD_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(1.67*b[,1]~b[,2])
## hist(1.67*b[b[,2]==occ_year$Year[y],1])
change_core_dim <- coefficients(lm(1.67*b[,1]~b[,2]))[2]


shift_table <- rbind(shift_table,data.frame(species=species_name,core_shift=shift_a,southern_shift=shift_a25,northern_shift=shift_a75,core_range_change=change_core_dim))
}

########################
## this might be better
########################

shift_table_resample  <-data.frame()
for (i in l_Rdata){

load(i)
print(species_name)

## 1. Core
a <- cbind(c(M_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a[,1]~a[,2])
## hist(a[a[,2]==occ_year$Year[y],1])
shift_a <- c()
a <- matrix(a[,1],nrow=nbr_iteration)
time_y <- 1:ncol(a)
for (it in 1:nbr_iteration){
shift_a[it] <- coefficients(lm(a[it,]~time_y))[2]
}; shift_a <- quantile(shift_a,c(0.05,0.5,0.95))

## 2. 25th percentile (southern tail)
a25 <- cbind(c(q25_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a25[,1]~a25[,2])
## hist(a25[a25[,2]==occ_year$Year[y],1])
shift_a25 <- c()
a25 <- matrix(a25[,1],nrow=nbr_iteration)
time_y <- 1:ncol(a25)
for (it in 1:nbr_iteration){
shift_a25[it] <- coefficients(lm(a25[it,]~time_y))[2]
}; shift_a25 <- quantile(shift_a25,c(0.05,0.5,0.95))


## 3. 75th percentile (northern front)
a75 <- cbind(c(q75_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(a75[,1]~a75[,2])
## hist(a75[a75[,2]==occ_year$Year[y],1])
shift_a75 <- c()
a75 <- matrix(a75[,1],nrow=nbr_iteration)
time_y <- 1:ncol(a75)
for (it in 1:nbr_iteration){
shift_a75[it] <-  coefficients(lm(a75[it,]~time_y))[2]
}; shift_a75 <- quantile(shift_a75,c(0.05,0.5,0.95))


## 4. expansion of range (1.67*StandardDeviation)
b <- cbind(c(SD_y),rep(occ_year$Year,rep(nbr_iteration,nbr_year)))
## boxplot(1.67*b[,1]~b[,2])
## hist(1.67*b[b[,2]==occ_year$Year[y],1])
shift_b <- c()
b <- matrix(b[,1],nrow=nbr_iteration)
time_y <- 1:ncol(b)
for (it in 1:nbr_iteration){
shift_b[it] <- coefficients(lm(b[it,]~time_y))[2]
}; change_core_dim <- quantile(shift_b,c(0.05,0.5,0.95))

shift_table_resample <- rbind(shift_table_resample,data.frame(species=species_name,percentile=c("5%","50%","95%"),core_shift=shift_a,southern_shift=shift_a25,northern_shift=shift_a75,core_range_change=change_core_dim))
}

## Northern shift
shift_north <- c(unlist(by(shift_table_resample$core_shift>0,shift_table_resample$species,FUN=sum)))
(shift_north_sp <- names(shift_north)[shift_north==3])
mean(shift_table_resample$core_shift[(as.character(shift_table_resample$species) %in% shift_north_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$core_shift[(as.character(shift_table_resample$species) %in% shift_north_sp) & shift_table_resample$percentile=="50%"])

northern_shift <- c(unlist(by(shift_table_resample$northern_shift>0,shift_table_resample$species,FUN=sum)))
(northern_shift_sp <- names(northern_shift)[northern_shift==3])
mean(shift_table_resample$northern_shift[(as.character(shift_table_resample$species) %in% northern_shift_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$northern_shift[(as.character(shift_table_resample$species) %in% northern_shift_sp) & shift_table_resample$percentile=="50%"])

southern_shift <- c(unlist(by(shift_table_resample$southern_shift>0,shift_table_resample$species,FUN=sum)))
(southern_shift_sp <- names(southern_shift)[southern_shift==3])
mean(shift_table_resample$southern_shift[(as.character(shift_table_resample$species) %in% southern_shift_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$southern_shift[(as.character(shift_table_resample$species) %in% southern_shift_sp) & shift_table_resample$percentile=="50%"])

## Southern shift
shift_north <- c(unlist(by(shift_table_resample$core_shift<0,shift_table_resample$species,FUN=sum)))
(shift_north_sp <- names(shift_north)[shift_north==3])
mean(shift_table_resample$core_shift[(as.character(shift_table_resample$species) %in% shift_north_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$core_shift[(as.character(shift_table_resample$species) %in% shift_north_sp) & shift_table_resample$percentile=="50%"])

northern_shift <- c(unlist(by(shift_table_resample$northern_shift<0,shift_table_resample$species,FUN=sum)))
(northern_shift_sp <- names(northern_shift)[northern_shift==3])
mean(shift_table_resample$northern_shift[(as.character(shift_table_resample$species) %in% northern_shift_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$northern_shift[(as.character(shift_table_resample$species) %in% northern_shift_sp) & shift_table_resample$percentile=="50%"])

southern_shift <- c(unlist(by(shift_table_resample$southern_shift<0,shift_table_resample$species,FUN=sum)))
(southern_shift_sp <- names(southern_shift)[southern_shift==3])
mean(shift_table_resample$southern_shift[(as.character(shift_table_resample$species) %in% southern_shift_sp) & shift_table_resample$percentile=="50%"])
range(shift_table_resample$southern_shift[(as.character(shift_table_resample$species) %in% southern_shift_sp) & shift_table_resample$percentile=="50%"])

## Number of site occupied across latitudinal gradient
library(ggplot2)
library(dplyr)
shift_table <-data.frame()
l_Rdata <- list.files()[grep("RData",list.files())]


load(l_Rdata[24])
print(species_name)

## all site occupied
full_europe[] <- NA
full_europe[as.numeric(rownames(per_gricell_obs))] <- 1
bin_size <- 150000
sample_dist <- as.data.frame(rasterToPoints(full_europe))
sample_dist$bin <- (sample_dist$y-min(sample_dist$y)+bin_size)%/%bin_size

per_gricell_obs <- as.data.frame(apply(trial_y[1:100,,1:4],1,function(x) by(x,point_cell[1:dim(trial_y)[2]],FUN=sum)))
per_gricell_obs[per_gricell_obs>0] <- 1

	full_europe[as.numeric(rownames(per_gricell_obs))] <- per_gricell_obs[,1]
	occ_dist <- as.data.frame(rasterToPoints(full_europe))
	occ_dist$bin <- (occ_dist$y-min(sample_dist$y)+bin_size)%/%bin_size

#####

per_gricell_obs1 <- as.data.frame(apply(trial_y[1:100,,2:5],1,function(x) by(x,point_cell[1:dim(trial_y)[2]],FUN=sum)))
per_gricell_obs1[per_gricell_obs1>0] <- 1

per_gricell_obs2 <- as.data.frame(apply(trial_y[1:100,,6:9],1,function(x) by(x,point_cell[1:dim(trial_y)[2]],FUN=sum)))
per_gricell_obs2[per_gricell_obs2>0] <- 1

## Extinction Colonisation
change_per_grid <- as.matrix(per_gricell_obs2)-as.matrix(per_gricell_obs1)

ch.layer <- matrix(NA,ncol=ncol(change_per_grid),nrow=nrow(change_per_grid))
for(it in 1:ncol(change_per_grid)){
	full_europe[as.numeric(rownames(per_gricell_obs))] <- change_per_grid [,it]
	change_dist <- as.data.frame(rasterToPoints(full_europe))
	ch.layer[,it] <- change_dist$layer
}

change_dist$bin <- (change_dist$y-min(change_dist$y)+bin_size)%/%bin_size

extinction <- as.data.frame(apply(ch.layer,2,function(x) by(x,change_dist$bin,function(y) sum(y==-1))))
colonisation <- as.data.frame(apply(ch.layer,2,function(x) by(x,change_dist$bin,function(y) sum(y==1))))

## stasis present/absent
change_per_grid <- as.matrix(per_gricell_obs2)+as.matrix(per_gricell_obs1)
ch.layer <- matrix(NA,ncol=ncol(change_per_grid),nrow=nrow(change_per_grid))

for(it in 1:ncol(change_per_grid)){
	full_europe[as.numeric(rownames(per_gricell_obs))] <- change_per_grid [,it]
	change_dist <- as.data.frame(rasterToPoints(full_europe))
	ch.layer[,it] <- change_dist$layer
}

change_dist$bin <- (change_dist$y-min(change_dist$y)+bin_size)%/%bin_size

stasis_present <- as.data.frame(apply(ch.layer,2,function(x) by(x,change_dist$bin,function(y) sum(y==2))))
stasis_absent <- as.data.frame(apply(ch.layer,2,function(x) by(x,change_dist$bin,function(y) sum(y==0))))

full_europe[as.numeric(rownames(per_gricell_obs))] <- 1
sampling_dist <- as.data.frame(rasterToPoints(full_europe))
sampling_dist$bin <- (sampling_dist$y-min(sampling_dist$y)+bin_size)%/%bin_size

change_dist <- data.frame(bin=as.numeric(unlist(dimnames(table(sampling_dist$bin[sampling_dist$layer==1])))),latitude=(as.numeric(unlist(dimnames(table(sampling_dist$bin[sampling_dist$layer==1]))))*min(sampling_dist$y))+(bin_size%/%2),count_sampling=as.numeric(table(sampling_dist$bin[sampling_dist$layer==1])))
change_dist$stasis_present <- apply(stasis_present,1,median)
change_dist$extinction <- apply(extinction,1,median)
change_dist$colonisation <- apply(colonisation,1,median)

change_dist$stasis_present_prop <- floor((change_dist$stasis_present+0.0000001)/(change_dist$count_sampling)*100)
change_dist$extinction_prop <- floor((change_dist$extinction+0.0000001)/(change_dist$count_sampling)*100)
change_dist$colonisation_prop <- floor((change_dist$colonisation+0.0000001)/(change_dist$count_sampling)*100)
change_dist$full_sampling_prop <-100

long_data_change <- data.frame(latitude=rep(change_dist$latitude,change_dist$extinction_prop),change=3)
long_data_change <- rbind(long_data_change, data.frame(latitude=rep(change_dist$latitude,change_dist$stasis_present_prop),change=1))
long_data_change <- rbind(long_data_change, data.frame(latitude=rep(change_dist$latitude,change_dist$colonisation_prop),change=2))
long_data_change <- rbind(long_data_change, data.frame(latitude=rep(change_dist$latitude,(100-change_dist$colonisation_prop-change_dist$stasis_present_prop-change_dist$extinction_prop)),change=4))

prop_change <- table(long_data_change$change,long_data_change$latitude)
dev.new()
barplot(prop_change, horiz=TRUE,col=c(rgb(250/255,200/255,0/255),"lightblue3",rgb(255/255,100/255,0/255),"grey90"))

long_data_change_raw <- data.frame(latitude=rep(change_dist$latitude,change_dist$extinction),change=3)
long_data_change_raw <- rbind(long_data_change_raw, data.frame(latitude=rep(change_dist$latitude,change_dist$stasis_present),change=1))
long_data_change_raw <- rbind(long_data_change_raw, data.frame(latitude=rep(change_dist$latitude,change_dist$colonisation),change=2))
long_data_change_raw <- rbind(long_data_change_raw, data.frame(latitude=rep(change_dist$latitude,(change_dist$count_sampling-change_dist$colonisation-change_dist$stasis_present-change_dist$extinction)),change=4))

prop_change_raw <- table(long_data_change_raw$change,long_data_change_raw$latitude)
dev.new()
barplot(prop_change_raw, horiz=TRUE,col=c(rgb(250/255,200/255,0/255),"lightblue3",rgb(255/255,100/255,0/255),"grey90"),legend = c("constant presence", "recent colonization", "recent extinction","absence"))


### Produce a distribution map for period 1 and 2 
library(climateExtract); library(raster); library(sp); library(rgdal); library(rgeos); library(SDMTools); library(maptools)

## GET spatial layers
wd <- getwd()

if(Sys.info()[4]=="WLD-3VHP992"){
	working_d <- "C:/Users/Public/Documents/gis_data"
	setwd(file.path(working_d,"physical_vector_layers/ne_50m_physical/"))
	ocean <- readOGR(".","ne_50m_ocean")
	graticule <- readOGR(".","ne_50m_graticules_10")
	setwd(file.path(working_d,"cultural_vector_layers/ne_50m_cultural/"))
	country <- readOGR(".","ne_50m_admin_0_boundary_lines_land")
	country_p <- readOGR(".","ne_50m_admin_0_sovereignty")
} else {
	working_d <- "/Users/retoschmucki/Documents/GIS layers/natural_earth_vector"
	setwd(file.path(working_d,"50m_physical"))
	ocean <- readOGR(".","ne_50m_ocean")
	setwd("/Users/retoschmucki/Documents/GIS layers/natural_earth_vector/50m_physical/ne_50m_graticules_all/")
	graticule <- readOGR(".","ne_50m_graticules_10")
	setwd(file.path(working_d,"50m_cultural/"))
	country <- readOGR(".","ne_50m_admin_0_boundary_lines_land")
	country_p <- readOGR(".","ne_50m_admin_0_sovereignty")	
}
setwd(wd)

## Crop around European BMS sampling points
CP<-as(extent(-30,50,15,90), "SpatialPolygons")

proj4string(CP) <- CRS("+init=epsg:4326")
country_c <- crop(country,CP, byid=TRUE)
ocean_c <- crop(ocean, CP, byid=TRUE)
graticule_c <- crop(graticule,ocean_c,byid=TRUE)
country_p_c <- crop(country_p,CP,byid=TRUE)

laea.proj <- "+init=epsg:3035"
country_c <- spTransform(country_c,CRS=CRS(laea.proj))
country_p_c <- spTransform(country_p_c,CRS=CRS(laea.proj))
ocean_c <- spTransform(ocean_c,CRS=CRS(laea.proj))
graticule_c <- spTransform(graticule_c,CRS=CRS(laea.proj))

col_ramp <- colorRampPalette(c(rgb(255/255,0/255,0/255),"lightblue4"))(2)

par(mfrow=c(1,2))

full_europe[as.numeric(rownames(per_gricell_obs))] <- rowMeans(per_gricell_obs1)
image(full_europe,axes=T,col=col_ramp,zlim=c(0,1))
plot(country_p_c[country_p_c@data$sovereignt%in% c("United Kingdom","Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
plot(country_c,lwd=0.3,add=T)
polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
## add a legend for the color gradient
## plot(full_europe, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,1),axis.args=list(cex.axis=1,tck=-0.5))


full_europe[as.numeric(rownames(per_gricell_obs))] <- rowMeans(per_gricell_obs2)
image(full_europe,axes=T,col=col_ramp,zlim=c(0,1))
plot(country_p_c[country_p_c@data$sovereignt%in% c("United Kingdom","Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
plot(country_c,lwd=0.3,add=T)
polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
## add a legend for the color gradient
## plot(full_europe, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,1),axis.args=list(cex.axis=1,tck=-0.5))
