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

library(climateExtract)
library(rgdal)
library(rgeos)
library(raster)
library(SDMTools)
library(maptools)
library(sp)

climate_data <- extract_nc_value(2006,2014,local_file=FALSE,clim_variable='mean temp',grid_size=0.25)
monthly_mean <- temporal_mean(climate_data,"monthly")

## full map
spring_summer_months <- as.numeric(monthly_mean$year_month$month) >= 4 & as.numeric(monthly_mean$year_month$month) <= 9
spring_summer_values <- monthly_mean$value_array[,,spring_summer_months]

yearly_spring_summer_meanT <- array(NA,c(length(monthly_mean$longitude),length(monthly_mean$latitude),length(unique(monthly_mean$year_month$year))))

for (y in as.character(unique(monthly_mean$year_month$year))){
	dim3 <- which(as.character(unique(monthly_mean$year_month$year))==y)
	A <- spring_summer_values[,,as.character(monthly_mean$year_month$year)[spring_summer_months]==y]
	yearly_spring_summer_meanT[,,dim3] <- apply(A,c(1,2),mean)
}

species_name <- "Anthocharis_cardamines"

## on MAC /Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/Audusseau_et_al_revised0303 (002)Reto_edit.docx
coord_zone <- read.csv(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIcoord_zone.csv"))
occ_summary <- read.csv(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIoutdaytrimmed.csv"))

## on PC
coord_zone <- read.csv(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIcoord_zone.csv"))
occ_summary <- read.csv(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jags_model_1_",species_name,"_ES_FR_UK_NL_DE_FIoutdaytrimmed.csv"))

coord_zone.sp <- coord_zone
coordinates(coord_zone.sp) <- ~longitude + latitude
proj4string(coord_zone.sp) <- CRS("+init=epsg:3035")
coord_zone.wgs84 <- as.data.frame(spTransform(coord_zone.sp,CRS=CRS("+init=epsg:4326")))
points <- coord_zone.wgs84[,c("site","longitude","latitude")]
names(points) <- c("site_id","longitude","latitude")

point.monthly_mean <- point_grid_extract(monthly_mean,points)

point.monthly_mean.spring_summer <- point.monthly_mean[spring_summer_months,]

value_peryear <- data.frame()

for (per_year in unique(floor(point.monthly_mean.spring_summer$date_extract/100))) {
	val_year <- apply(point.monthly_mean.spring_summer[floor(point.monthly_mean.spring_summer$date_extract/100) == per_year,],2,mean)
	value_peryear <- rbind(value_peryear,val_year)
}

value_peryear <- t(value_peryear)
##names(value_peryear) <- as.character(unique(floor(point.monthly_mean.spring_summer$date_extract/100)))
site_coord.spring_summerT <- cbind(coord_zone,value_peryear[-1,])
names(site_coord.spring_summerT) <- c(names(coord_zone),paste0("y_",as.character(unique(floor(point.monthly_mean.spring_summer$date_extract/100)))))



## produce a map
## set your working directory to containing the physical and cultural folders

## on MAC
working_d <- "/Users/retoschmucki/Documents/GIS layers/natural_earth_vector"
setwd(file.path(working_d,"50m_physical"))
ocean <- readOGR(".","ne_50m_ocean")
setwd("/Users/retoschmucki/Documents/GIS layers/natural_earth_vector/50m_physical/ne_50m_graticules_all/")
graticule <- readOGR(".","ne_50m_graticules_10")
setwd(file.path(working_d,"50m_cultural/"))
country <- readOGR(".","ne_50m_admin_0_boundary_lines_land")
country_p <- readOGR(".","ne_50m_admin_0_sovereignty")

## on PC
working_d <- "C:/Users/Public/Documents/gis_data"
setwd(file.path(working_d,"physical_vector_layers/ne_50m_physical/"))
ocean <- readOGR(".","ne_50m_ocean")
graticule <- readOGR(".","ne_50m_graticules_10")
setwd(file.path(working_d,"cultural_vector_layers/ne_50m_cultural/"))
country <- readOGR(".","ne_50m_admin_0_boundary_lines_land")
country_p <- readOGR(".","ne_50m_admin_0_sovereignty")

## Crop around European BMS sampling points

CP<-as(extent(-30,50,15,90), "SpatialPolygons")

laea.proj <- "+init=epsg:3035"

proj4string(CP) <- CRS("+init=epsg:4326")
country_c <- crop(country,CP,byid=TRUE)
country_p_c <- crop(country_p,CP,byid=TRUE)
ocean_c <- crop(ocean, CP, byid=TRUE)
graticule_c <- crop(graticule,ocean_c,byid=TRUE)

country_c <- spTransform(country_c,CRS=CRS("+init=epsg:3035"))
country_p_c <- spTransform(country_p_c,CRS=CRS("+init=epsg:3035"))
ocean_c <- spTransform(ocean_c,CRS=CRS("+init=epsg:3035"))
graticule_c <- spTransform(graticule_c,CRS=CRS("+init=epsg:3035"))


full_r <- raster(extent(2763687,5733766,1605279,5075336),crs=crs(laea.proj))
res(full_r) <- 50000
col_ramp <- colorRampPalette(c("blue","lightsteelblue3", "yellow", "orange", "red"))(50)
	coordinates(site_coord.spring_summerT) <- ~ longitude+latitude
	proj4string(site_coord.spring_summerT) = CRS(laea.proj)

## var_raster <- rasterize(site_coord.spring_summerT, full_r,1, fun=sum)
## probability of at least one observation per 20 trials 
## var_raster <- rasterize(site_coord.spring_summerT, full_r,1, fun=function (x,na.rm = TRUE) {v <- rep(mean(x),20); v[1:length(x)]<- x; 1-prod(1-v)})

## Temperature
var_raster <- rasterize(site_coord.spring_summerT, full_r,site_coord.spring_summerT$y_2009, fun=mean)

image(var_raster,axes=T,col=col_ramp,zlim=c(0,25))
plot(country_p_c[country_p_c@data$sovereignt%in% c("Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
plot(country_c,lwd=0.3,add=T)
polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
## add a legend for the color gradient
plot(var_raster, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,25),axis.args=list(cex.axis=1,tck=-0.5))

##### HEAD
## species_name <- "Anthocharis_cardamines"
## species_name <- "Apatura_iris"
##########

########################################
## Estimating distribution of occupancy
########################################

## set the OCC prob per site/year
occ_value <- occ_summary[substr(occ_summary$X,1,2)=="z[",]

occ_coord <- data.frame()
for (y in 1:length(unique(occ_coord$YEAR))){
	y1 <-cbind(occ_value[grep(paste0(y,"]$"),occ_value$X, value=FALSE),],coord_zone,YEAR=(min(unique(occ_coord$YEAR))-1+y))
	occ_coord <- rbind(occ_coord,y1)
}	

trial_1000 <- sapply(occ_coord$mean,function(x) {rbinom(1000,1,x)},simplify=TRUE)


## build spatial points to overlay
laea.proj <- "+init=epsg:3035"
wgs84 <- "+init=epsg:4326"

coordinates(occ_coord) <- ~ longitude+latitude
proj4string(occ_coord) = CRS(laea.proj)
## occ_coord_wgs84 <- spTransform(occ_coord,CRS=CRS(wgs84))

## build raster across Europe
full_europe <- raster(extent(2763687,5733766,1605279,5075336),crs=crs(laea.proj))
res(full_europe) <- 50000
full_europe[] <- 1:ncell(full_europe)

point_cell <- extract(full_europe,occ_coord)

occ_coord.df <- as.data.frame(occ_coord)
occ_coord.df$grid_id <- point_cell

length(unique(occ_coord$site))
trial_y <- array(NA,dim=c(1000,length(unique(occ_coord$site)),length(unique(occ_coord$YEAR))))
for(y in 1:9){
trial_y[,,y] <- trial_1000[,which(occ_coord.df$YEAR==y+2005)]
}

## shift estimate
full_europe[] <- NA

## For year 2006-2014
M_y <- matrix(NA,nrow=1000,ncol=9)
SD_y <- matrix(NA,nrow=1000,ncol=9)
q25_y <- matrix(NA,nrow=1000,ncol=9)
q75_y <- matrix(NA,nrow=1000,ncol=9)

for (y in 1:9){
per_gricell_obs <- as.data.frame(apply(trial_y[,,y],1,function(x) by(x,point_cell[1:dim(trial_y)[2]],FUN=sum)))
per_gricell_obs[per_gricell_obs>0] <- 1

M_ <- c()
SD_ <- c()
q25_ <- c()
q75_ <- c()
for (it in 1:dim(per_gricell_obs)[2]){
	full_europe[as.numeric(rownames(per_gricell_obs))] <- per_gricell_obs[,it]
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

a <- cbind(c(M_y),rep(1:9,rep(1000,9)))
boxplot(a[,1]~a[,2])
lm(a[,1]~a[,2])

a25 <- cbind(c(q25_y),rep(1:9,rep(1000,9)))
boxplot(a25[,1]~a25[,2])
lm(a25[,1]~a25[,2])

a75 <- cbind(c(q75_y),rep(1:9,rep(1000,9)))
boxplot(a75[,1]~a75[,2])
summary(lm(a75[,1]~a75[,2]))

b <- cbind(c(SD_y),rep(1:9,rep(1000,9)))
boxplot(b[,1]~b[,2])
lm(b[,1]~b[,2])

## MAP distribution
it=1
for (it in 1:100){
full_europe[] <- NA
full_europe[as.numeric(rownames(per_gricell_obs))] <- per_gricell_obs[,it]

col_ramp <- colorRampPalette(c('orange',"green"))(2)

image(full_europe,axes=T,col=col_ramp,zlim=c(0,1))
plot(country_p_c[country_p_c@data$sovereignt%in% c("Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
plot(country_c,lwd=0.3,add=T)
polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
## add a legend for the color gradient
plot(full_europe, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,1),axis.args=list(cex.axis=1,tck=-0.5))
}

## Compute and map uncertainty (sd)
str(per_gricell_obs)
t <- rowMeans(per_gricell_obs)
sd_prob <- sqrt(t*(1-t))

full_europe[as.numeric(rownames(per_gricell_obs))] <- sd_prob

col_ramp <- colorRampPalette(c("blue","lightsteelblue3", "yellow", "orange", "red"))(25)

image(full_europe,axes=T,col=col_ramp,zlim=c(0,0.5))
plot(country_p_c[country_p_c@data$sovereignt%in% c("Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
plot(country_c,lwd=0.3,add=T)
polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
## add a legend for the color gradient
plot(full_europe, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,0.5),axis.args=list(cex.axis=1,tck=-0.5))


## shift estimate


coord.df <- occ_coord.df[occ_coord.df$YEAR==2006,c("site","latitude","longitude","grid_id")]


t <- data.frame(grid_id=as.numeric(rownames(per_gricell_obs)),occurrence=per_gricell_obs[,1])
t <- merge(coord.df,t,by=c("grid_id"))



test <- apply(trial_y,c(2,3),function(x) mean(x))
str(test)

occ_cell <- matrix(NA,nrow=length(unique(point_cell)),ncol=dim(trial_1000)[1])

for (li in 1:dim(trial_1000)[1]){
	occ_cell[,li] <- as.numeric(by(trial_1000[li,], point_cell, function(x) if(sum(x)==0){0}else{1},simplify=TRUE))
}

occ_cell <- as.data.frame(occ_cell)
occ_cell$cell_id <- as.numeric(names(by(trial_10000[li,], point_cell, function(x) if(sum(x)==0){0}else{1},simplify=TRUE)))

rowMeans(occ_cell[,-10001])


## on MAC
load(paste0("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jagsoutput2.Rdata"))
## on PC
load(paste0("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/",species_name,"/results/jagsoutput2.Rdata"))

mean_occ_p <- apply(out2$BUGSoutput$sims.list$z,c(2,3),mean)




prob_at_least_one.in_many <- function (x,trial=20, na.rm = TRUE) {v <- rep(mean(x),trial); v[1:length(x)]<- x; 1-prod(1-v)}
prob_from_one <- function(x,...) {set.seed(seedset);sample(x,1,replace=TRUE)}

## PER YEAR
nbr_it <- 10000
shift.array <- array(data=NA, dim=c(2,5,nbr_it))

for (z in 1:nbr_it){

density.curves <- list()
y.range <- data.frame()
x.range <- data.frame()
quantile.value <- data.frame()

seedset <- sample(c(500:1000),1)

for(y in c(2:9)){

var_raster2 <- rasterize(site_coord.spring_summerT, full_r,mean_occ_p[,y], fun=prob_from_one)
test.rast2df <- as.data.frame(rasterToPoints(var_raster2))
bin_sim <- sapply(test.rast2df$layer,FUN=function(x){rbinom(1,25,x)})

all.grid <- data.frame(id=names(table(test.rast2df$y)),sample_effort=as.numeric(table(test.rast2df$y)))
northing <- data.frame(id=names(table(rep(test.rast2df$y,bin_sim))),Northing_raw=as.numeric(table(rep(test.rast2df$y,bin_sim))))
merged.northing <- merge(all.grid,northing,by=c("id"))
merged.northing$stdr_effort <- (merged.northing$Northing_raw/merged.northing$sample_effort)

d <- density(rep(as.numeric(as.character(merged.northing$id)),merged.northing$stdr_effort))
y.range <- rbind(y.range,range(d$y))
x.range <- rbind(x.range,range(d$x))
density.curves[[y]] <- d

quantile.value <-  rbind(quantile.value,as.numeric(quantile(rep(as.numeric(as.character(merged.northing$id)),merged.northing$stdr_effort),c(0.1,0.25,0.5,0.75,0.90))))

}

quantile.value$Year <- 1:dim(quantile.value)[1]
names(quantile.value) <- c("p10","p25","p50","p75","p90","Year")
shift.array[,,z] <- matrix(c(as.numeric(coefficients( glm(p10~Year, data = quantile.value))),as.numeric(coefficients( glm(p25~Year, data = quantile.value))),as.numeric(coefficients( glm(p50~Year, data = quantile.value))),
	as.numeric(coefficients( glm(p75~Year, data = quantile.value))),as.numeric(coefficients( glm(p90~Year, data = quantile.value)))),nrow=2,ncol=5)

}

d10 <- density(shift.array[2,1,])
d50 <- density(shift.array[2,2,])
d90 <- density(shift.array[2,3,])

par(mfrow=c(3,1))
plot(d10,col='red',main="d10")
plot(d50,col='blue',main="d50")
plot(d90,col='green',main="d90")



dev.new()
col_ramp <- colorRampPalette(c("blue","green", "magenta"))(9)
plot(density.curves[[2]],xlim=range(x.range),ylim=range(y.range),col=col_ramp[1],lwd=2)

Sys.sleep(2)
for( i in 3:9){
lines(density.curves[[i]],col=col_ramp[i],lwd=2)
Sys.sleep(2)
}

dev.new()
plot(c(2007:2014),quantile.value[,1])
plot(c(2007:2014),quantile.value[,2],col='red')
plot(c(2007:2014),quantile.value[,"p75"],col='red')







# var_raster1 <- rasterize(site_coord.spring_summerT, full_r,mean_occ_p[,y], fun=prob_at_least_one.in_many)
# dev.new()

# image(var_raster2,axes=F,col=col_ramp,zlim=c(0,1))
# plot(country_p_c[country_p_c@data$sovereignt%in% c("Italy","Switzerland","Belgium","Luxembourg","Russia","Denmark","Austria","Poland","Czech Republic","Estonia","Ireland"),],col="white",border=NA,add=T)
# plot(country_c,lwd=0.3,add=T)
# polygon(c(5352923,5352923,5809489,5809489,5352923),c(2024942,3548678,3548678,2024942,2024942),col="white",border=NA)
# plot(ocean_c,col="lightsteelblue1",add=T)
# plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
# ## add a legend for the color gradient
# plot(var_raster, smallplot=c(0.85, 0.88, 0.17, .55),col=col_ramp, legend.only=TRUE,zlim=c(0,1),axis.args=list(cex.axis=1,tck=-0.5))
# title(y)

# hist(rep(test.rast2df$y,bin_sim))
# North_index_m <- mean(rep(test.rast2df$y,bin_sim))
# North_index_sd <- sd(rep(test.rast2df$y,bin_sim))

# East_index_m <- mean(rep(test.rast2df$x,bin_sim))
# East_index_sd <- sd(rep(test.rast2df$x,bin_sim))

# par(mfrow=c(2,1))
# hist(rep(test.rast2df$y,bin_sim),main="Northing")
# hist(rep(test.rast2df$x,bin_sim),main="Easting")
# hist(test.rast2df$y)
# ## clearly this is biased by the distribution of our sampling points across Europe

all.grid <- data.frame(id=names(table(test.rast2df$y)),sample_effort=as.numeric(table(test.rast2df$y)))
northing <- data.frame(id=names(table(rep(test.rast2df$y,bin_sim))),Northing_raw=as.numeric(table(rep(test.rast2df$y,bin_sim))))
merged.northing <- merge(all.grid,northing,by=c("id"))
merged.northing$stdr_effort <- (merged.northing$Northing_raw/merged.northing$sample_effort) ## *merged.northing$id)/sum(merged.northing$id)

sum(merged.northing$stdr_effort*as.numeric(as.character(merged.northing$id)))/sum(merged.northing$stdr_effort)

density_curve <- list()
d <- density(rep(as.numeric(as.character(merged.northing$id)),merged.northing$stdr_effort))

plot(d,xlim=c(1500000,5000000),col='blue')
