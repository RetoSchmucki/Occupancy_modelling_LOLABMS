## =================================
## get the occupancy model results
## Date 17.10.2016
## Author Reto Schmucki
## =================================

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
## on MAC /Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/Audusseau_et_al_revised0303 (002)Reto_edit.docx

coord_zone <- read.csv("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/Araschnia_levana/results/jags_model_1_Araschnia_levana_ES_FR_UK_NL_DE_FIcoord_zone.csv")
## on PC
coord_zone <- read.csv("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/LOLA_BMS/estimate_all_AUGUST/countries_UNI_2006_sd_67/Araschnia_levana/results/jags_model_1_Anthocharis_cardamines_ES_FR_UK_NL_DE_FIcoord_zone.csv ")

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
#var_raster <- rasterize(site_coord.spring_summerT, full_r,1, fun=function (x,na.rm = TRUE) {v <- rep(mean(x),20); v[1:length(x)]<- x; 1-prod(1-v)})

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

species_name <- "Anthocharis_cardamines"
species_name <- "Apatura_iris"
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
