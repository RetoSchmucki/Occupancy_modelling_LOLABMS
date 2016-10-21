R

rm(list=ls())

library(climateExtract); library(raster); library(sp); library(rgdal); library(rgeos); library(SDMTools); library(maptools)

## GET spatial layers

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

setwd("/Users/retoschmucki/")

## GET climate data
climate_data <- extract_nc_value(2006,2014,local_file=FALSE,clim_variable='mean temp',grid_size=0.25)
monthly_mean <- temporal_mean(climate_data,"monthly")

s_m <- seq(from=4,to=(9*12),by=12)
e_m <- seq(from=9,to=(9*12),by=12)
A_S <- array(NA,dim=c(dim(monthly_mean$value_array)[1],dim(monthly_mean$value_array)[2],9))
for (y in 1:9){
A <- monthly_mean$value_array[,,s_m[y]:e_m[y]]
A_S[,,y] <- apply(A,c(1,2),mean)
}

sp.points <- raster("tg_0.25deg_reg_v14.0.nc")
sp.points[] <-apply(t(A_S[,,9]),2, rev) 

## Calculate region summer mean temperature change

n_Tchange <- A_S[,,9] 
n_Tchange[] <- NA 
y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,176:201,y]),year=y))
}
summary(lm(y_temp$temp~y_temp$y))
n_Tchange[,c(176:201)] <- as.numeric(coefficients(lm(y_temp$temp~y_temp$y))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,151:175,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(151:175)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,126:150,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(126:150)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,101:125,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(101:125)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,76:100,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(76:100)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,51:75,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(51:75)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,26:50,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(26:50)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])

y_temp <- data.frame()
for(y in c(1:9)){
y_temp <- rbind(y_temp,data.frame(temp=c(A_S[,1:25,y]),year=y))
}
summary(lm(y_temp$temp~jitter(y_temp$y)))
n_Tchange[,c(1:25)] <- as.numeric(coefficients(lm(y_temp$temp~jitter(y_temp$y)))[2])
sp.points[] <-apply(t(n_Tchange),2, rev) 


## Crop around European BMS sampling points
CP<-as(extent(-30,50,15,90), "SpatialPolygons")

proj4string(CP) <- CRS("+init=epsg:4326")
country_c <- crop(country,CP, byid=TRUE)
ocean_c <- crop(ocean, CP, byid=TRUE)
graticule_c <- crop(graticule,ocean_c,byid=TRUE)
sp.points_c <- crop(sp.points, CP)

col_ramp <- colorRampPalette(c("yellow","red"))(12)
pdf("temp_trend.pdf")
image(sp.points,axes=F,col=col_ramp,zlim=c(-0.1,0.25))
plot(country_c,lwd=0.3,add=T)
plot(ocean_c,col="lightsteelblue1",add=T)
plot(graticule_c, col="gray15",lty=2,lwd=0.5,add=T)
polygon(c(45,45,70,70,45),c(25,80,80,25,25),col="white",border=NA)
plot(sp.points_c, smallplot=c(0.74, 0.76, 0.25, .65),col=col_ramp, legend.only=TRUE,zlim=c(-0.1,0.2),axis.args=list(cex.axis=1,tck=-0.5))

dev.off()



plot(c(A_S[,,1]))
points(1:length(c(A_S[,,2])),c(A_S[,,2]),col='red')
points(1:length(c(A_S[,,2])),c(A_S[,,3]),col='magenta')
points(1:length(c(A_S[,,2])),c(A_S[,,4]),col='blue')

mean(c(A_S[,,1]),na.rm=TRUE)
mean(c(A_S[,,2]),na.rm=TRUE)
mean(c(A_S[,,3]),na.rm=TRUE)
mean(c(A_S[,,4]),na.rm=TRUE)
mean(c(A_S[,,5]),na.rm=TRUE)
mean(c(A_S[,,6]),na.rm=TRUE)
mean(c(A_S[,,7]),na.rm=TRUE)
mean(c(A_S[,,8]),na.rm=TRUE)
mean(c(A_S[,,9]),na.rm=TRUE)

image(A_S[,,1])
dev.new()
image(A_S[,,9])

A_S[is.na(A_S)] <- -999
temp_trend <- apply(A_S,c(1,2),function(x) coefficients(lm(x~c(1:9)))[2])
temp_trend[0] <- NA
temp_trend[temp_trend>10] <- NA
temp_trend[temp_trend<(-10)] <- NA
image(temp_trend)

## build raster across Europe with the appropriate resolution for the grid cell, 50x50Km

gc.resolution <- 50000

full_europe <- raster(extent(2763687,5733766,1605279,5075336),crs=crs(laea.proj))
res(full_europe) <- gc.resolution
proj4string(full_europe) <- CRS("+init=epsg:3035")
full_europe[] <- 1:ncell(full_europe)
gr_cell.points <- as.data.frame(rasterToPoints(full_europe))

names(gr_cell.points) <- c("longitude","latitude","site_id")

coordinates(gr_cell.points) <- ~ longitude+latitude
proj4string(gr_cell.points) = CRS("+init=epsg:3035")
gr_cell.points_wgs84 <- as.data.frame(spTransform(gr_cell.points,CRS=CRS("+init=epsg:4326")))

point.monthly_mean <- point_grid_extract(monthly_mean,gr_cell.points_wgs84)

year <- substr(point.monthly_mean$date_extract,start=1,stop=4)
month <- substr(point.monthly_mean$date_extract,start=5,stop=6)
point.summer_monthly_mean <- point.monthly_mean[as.numeric(month) %in% c(4:9),]

summer.annual.mean <- by(point.summer_monthly_mean[,-1],year[as.numeric(month) %in% c(4:9)],FUN=mean,na.rm=TRUE)



var_raster <- rasterize(site_coord.spring_summerT, full_europe,site_coord.spring_summerT$y_2009, fun=mean)

image(A_S[,,1])
A_S[is.na(A_S)] <- -999
temp_trend <- apply(A_S,c(1,2),function(x) coefficients(lm(x~c(1:9))[2]))

