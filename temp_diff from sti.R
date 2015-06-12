rm(list=ls())
setwd("/home/reto/Documents/LOLA_BMS/Climate_EU")
# setwd("/Users/retoschmucki/Desktop")

# setwd("/home/reto/Desktop/ECAdata")
library(RODBC)
library(ncdf4)
library(raster)
library(rgeos)
library(chron)
library(FNN)
library(tcltk)
library(zoo)
library(classInt) 
library(colorRamps)
library(RColorBrewer)

#######################
# AVERAGE TEMPERATURE
########################

tw.ncdf <- nc_open("tg_0_25deg_reg_1995_2014_v11_0.nc")
lon <- ncvar_get(tw.ncdf,"longitude")
lat <- ncvar_get(tw.ncdf,"latitude")
nlon <- dim(lon)
nlat <- dim(lat)

# set day since in the data
day_since <- ncatt_get( tw.ncdf,"time")$units
timeserie_length <- length(ncvar_get(tw.ncdf,"time"))
fillvalue <- ncatt_get(tw.ncdf, "tg", "_FillValue")
day_vals <- ncvar_get(tw.ncdf,"time")

init_year <- as.numeric(strsplit(unlist(strsplit(gsub('days since ','',day_since),'-',fixed=TRUE)),' ',fixed=TRUE)[[1]][1])
init_month <- as.numeric(strsplit(unlist(strsplit(gsub('days since ','',day_since),'-',fixed=TRUE)),' ',fixed=TRUE)[[2]][1])
init_day <- as.numeric(strsplit(unlist(strsplit(gsub('days since ','',day_since),'-',fixed=TRUE)),' ',fixed=TRUE)[[3]][1])

avg_temp_transect <- data.frame(julianday=day_vals)
avg_temp_transect$day <- month.day.year(avg_temp_transect$julianday, c(month = init_month, day =init_day, year = init_year))$day
avg_temp_transect$month <- month.day.year(avg_temp_transect$julianday, c(month = init_month, day =init_day, year = init_year))$month
avg_temp_transect$year <- month.day.year(avg_temp_transect$julianday, c(month = init_month, day =init_day, year = init_year))$year


first.year <- 2000
last.year <- 2014

firstday<-avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year]
lastday<-avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday)),count=c(nlon,nlat,(lastday-firstday)+1))
tmp.array[tmp.array == fillvalue$value] <- NA


annual.mean <- array(NA,c(nlon,nlat,(last.year-first.year)+1))

# first.toget <-1

for( y in first.year:last.year){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	annual.mean[,,y-(first.year-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}


# build grid climate
lonlat <- expand.grid(lon,lat)
tmp.vec <- as.vector(tmp.array[,,1])
tmp.df01 <- data.frame(cbind(lonlat,tmp.vec))
names(tmp.df01) <- c("lon","lat","tg")
tmp.df_noNA <- tmp.df01[!is.na(tmp.df01$tg),]


# extract transect coordinates
if(Sys.info()[["nodename"]] =='reto-Precision-T1650'){
  uid_pwd <- read.csv("/home/reto/Dropbox/LOLA_BMS/lola_connection.csv",stringsAsFactors = F,header=F) 
} else { 
  uid_pwd <- read.csv("/User/retoschmucki/Dropbox/LOLA_BMS/lola_connection.csv",stringsAsFactors = F,header=F) 
}

ch.lola <- odbcConnect("LOLAdb", uid=uid_pwd[1,2], pwd=uid_pwd[2,2],case="postgresql") 

# retrive transect coord (long-lat)

site_coord_query<- paste0("SELECT DISTINCT
transect_id as site_id,
ST_X(geom) as lon,
ST_Y(geom) as lat
FROM
bms_update.all_bms_transects_coord as p
WHERE
monitoring_type = 1 AND
substring(transect_id from 1 for 2) != \'IL\'")

site_coord <- sqlQuery(ch.lola,site_coord_query,rows_at_time=1)


# get index of closest point in the climate grid
nc_index<-data.frame(site_id=NA,x_index=NA,y_index=NA,x_coord=NA,y_coord=NA)

pb <- tkProgressBar(title = "progress bar", min = 0,max = dim(site_coord)[1], width = 300)

for (i in 1:dim(site_coord)[1]){
  setTkProgressBar(pb, i, label=paste( round(i/dim(site_coord)[1]*100, 0),"% extracted"))
  nnindex<-get.knnx(tmp.df_noNA[,-3],site_coord[i,-1],1)
  nc_index <- rbind(nc_index,data.frame(site_id=site_coord$site_id[i],x_index=which(lon ==tmp.df_noNA[nnindex$nn.index,-3]$lon),y_index=which(lat ==tmp.df_noNA[nnindex$nn.index,-3]$lat),x_coord=site_coord$lon[i],y_coord=site_coord$lat[i]))
}
nc_index <- nc_index[-1,]


# get site per zone per species
site_zone <- read.csv("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEWJUNE/countries_MULTI_2000_sd_1/Aglais_io/jags_model_1_Aglais_io_ES_UK_NL_FIcoord_zone.csv")

z <- "V4"

index_toget <- unique(nc_index[as.character(nc_index$site_id) %in% as.character(site_zone[site_zone[,z]==1,"site"]),c("x_index","y_index")])

year_mean_tempZ <- data.frame()
for (y in 1:dim(annual.mean)[3]){
for (i in 1:length(index_toget$x_index)){
year_mean_tempZ <- rbind(year_mean_tempZ, data.frame(y=c(first.year:last.year)[y],mean_temp=annual.mean[index_toget$x_index[i],index_toget$y_index[i],y]))
}
}

# get species sti
get.sti <- function(ch.sql,species){
  
  transect_TempI_query <- paste("SELECT
		species_name,
		temp_mean,
		temp_sd
		FROM
		climate_data.sti
		WHERE
		species_name = \'",species,"\'
		",sep='')
  
  Species_STI <- sqlQuery(ch.sql,transect_TempI_query ,rows_at_time=1)
  
  return(Species_STI)
}
species <- "Aglais io"

year_mean_tempZ$sti <- get.sti(ch.lola,species)[,"temp_mean"]
year_mean_tempZ$sti_anomali <- year_mean_tempZ$mean_temp - year_mean_tempZ$sti

mean_anomalie <- stats::aggregate(year_mean_tempZ$sti_anomali,by=list(year_mean_tempZ$y),function(x) mean(x,rm.na=T))
names(mean_anomalie) <- c("year","sti_anomalie")

plot(mean_anomalie$year,mean_anomalie$sti_anomalie,ylim=c(-2.5*get.sti(ch.lola,species)[,"temp_sd"],2.5*get.sti(ch.lola,species)[,"temp_sd"]),type="n")
abline(h=0.67*get.sti(ch.lola,species)[,"temp_sd"],col="red",lty=2)
abline(h=0,col="grey")
abline(h=-0.67*get.sti(ch.lola,species)[,"temp_sd"],col="red",lty=2)
points(mean_anomalie$year,mean_anomalie$sti_anomalie,type="l")

##########


mean_climate <- apply(annual.mean,c(1,2),mean,na.rm=T)


lonlat <- expand.grid(lon,lat)
longlat_nona <- lonlat[!is.na(as.vector(mean_climate[,])),]
lonlat_mean_climate <- data.frame(cbind(longlat_nona,as.vector(mean_climate[,])[!is.na(as.vector(mean_climate[,]))]))
names(lonlat_mean_climate) <- c("long","lat","meanclimat")


nclr <- 12
plotclr <- blue2red(nclr) 
class <-  classIntervals(lonlat_mean_climate$meanclimat,n=nclr,style="fixed", fixedBreaks=c(-12,-5,-2,0,2,5,8,10,15,20,26))
colcode <- findColours(class,plotclr)
leg <- leg <-gsub("\\["," ",names(attr(colcode, "table")))
leg <-gsub("\\]"," ",leg)
leg <-gsub("\\)"," ",leg)
leg <-gsub("\\,"," : ",leg)

dev.new()
plot(lonlat_mean_climate$long,lonlat_mean_climate$lat,col=colcode,pch=15,main=paste(first.year,last.year))
legend('bottomright', legend=leg,fill=attr(colcode, "palette"))




# interest period
first.year.p <- 2001
last.year.p <- 2014

firstday.p <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year.p]
lastday.p <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year.p]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday.p)),count=c(nlon,nlat,(lastday.p-firstday.p)+1))
tmp.array[tmp.array == fillvalue$value] <- NA

annual.mean.p <- array(NA,c(nlon,nlat,(last.year.p-first.year.p)+1))

# first.toget <-1

for( y in first.year.p:last.year.p){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	annual.mean.p[,,y-(first.year.p-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday.p+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday.p+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y]-firstday.p+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}

mean_climate.p <- apply(annual.mean.p,c(1,2),mean,na.rm=T)

lonlat <- expand.grid(lon,lat)
longlat_nona <- lonlat[!is.na(as.vector(mean_climate.p[,])),]
lonlat_mean_climate.p <- data.frame(cbind(longlat_nona,as.vector(mean_climate.p[,])[!is.na(as.vector(mean_climate.p[,]))]))
names(lonlat_mean_climate.p) <- c("long","lat","meanclimat")

nclr <- 12
plotclr <- blue2red(nclr) 
class <-  classIntervals(lonlat_mean_climate.p$meanclimat,n=nclr,style="fixed", fixedBreaks=c(-12,-5,-2,0,2,5,8,10,15,20,26))
colcode <- findColours(class,plotclr)
leg <- leg <-gsub("\\["," ",names(attr(colcode, "table")))
leg <-gsub("\\]"," ",leg)
leg <-gsub("\\)"," ",leg)
leg <-gsub("\\,"," : ",leg)

dev.new()
plot(lonlat_mean_climate.p$long,lonlat_mean_climate.p$lat,col=colcode,pch=15,main=paste(first.year.p,last.year.p))
legend('bottomright', legend=leg,fill=attr(colcode, "palette"))

clim_diff <- array(NA,dim=dim(annual.mean.p))

for (y in 1:dim(annual.mean.p)[3]){
clim_diff[,,y] <- annual.mean.p[,,y] - mean_climate
}

mean_tempdiff <- apply(clim_diff,c(1,2),mean,na.rm=T)

mean_tempdiff[mean_tempdiff < -3 | mean_tempdiff > 3] <- NA

lonlat <- expand.grid(lon,lat)
longlat_nona <- lonlat[!is.na(as.vector(mean_tempdiff[,])),]
lonlat_mean_tempdiff <- data.frame(cbind(longlat_nona,as.vector(mean_tempdiff[,])[!is.na(as.vector(mean_tempdiff[,]))]))
names(lonlat_mean_tempdiff) <- c("long","lat","tempdiff")


nclr <- 9
plotclr <- blue2red(nclr) 
class <-  classIntervals(lonlat_mean_tempdiff$tempdiff,n=nclr,style="fixed", fixedBreaks=c(-3.0,-1.0,-0.5,0,0.5,1.0,1.5,3.0))
colcode <- findColours(class,plotclr)
leg <- leg <-gsub("\\["," ",names(attr(colcode, "table")))
leg <-gsub("\\]"," ",leg)
leg <-gsub("\\)"," ",leg)
leg <-gsub("\\,"," : ",leg)

dev.new()
plot(lonlat_mean_tempdiff$long,lonlat_mean_tempdiff$lat,col=colcode,pch=15,main='meantempdiff')
legend('bottomright', legend=leg,fill=attr(colcode, "palette"))


#############################################
#############################################

first.year <- 1971
last.year <- 2000

firstday<-avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year]
lastday<-avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday)),count=c(nlon,nlat,(lastday-firstday)+1))
tmp.array[tmp.array == fillvalue$value] <- NA


summer.mean.c <- array(NA,c(nlon,nlat,(last.year-first.year)+1))

summer.m <- c(4,5,6,7,8,9)
# first.toget <-1

for( y in first.year:last.year){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	summer.mean.c[,,y-(first.year-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}

summer_climate <- apply(summer.mean.c,c(1,2),mean,na.rm=T)



# interest period
first.year.p <- 2001
last.year.p <- 2014

firstday.p <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year.p]
lastday.p <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year.p]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday.p)),count=c(nlon,nlat,(lastday.p-firstday.p)+1))
tmp.array[tmp.array == fillvalue$value] <- NA

summer.mean.p <- array(NA,c(nlon,nlat,(last.year.p-first.year.p)+1))

summer.m <- c(4,5,6,7,8,9)

for( y in first.year.p:last.year.p){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	summer.mean.p[,,y-(first.year.p-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday.p+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday.p+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% summer.m]-firstday.p+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}


clim_diff.s <- array(NA,dim=dim(summer.mean.p))

for (y in 1:dim(summer.mean.p)[3]){
clim_diff.s[,,y] <- summer.mean.p[,,y] - summer_climate
}

mean_tempdiff.s <- apply(clim_diff.s,c(1,2),mean,na.rm=T)

mean_tempdiff.s[mean_tempdiff.s < -3 | mean_tempdiff.s > 3] <- NA

lonlat <- expand.grid(lon,lat)
longlat_nona <- lonlat[!is.na(as.vector(mean_tempdiff.s[,])),]
lonlat_mean_tempdiff.s <- data.frame(cbind(longlat_nona,as.vector(mean_tempdiff.s[,])[!is.na(as.vector(mean_tempdiff.s[,]))]))
names(lonlat_mean_tempdiff.s) <- c("long","lat","tempdiff")


nclr <- 9
plotclr <- blue2red(nclr) 
class <-  classIntervals(lonlat_mean_tempdiff.s$tempdiff,n=nclr,style="fixed", fixedBreaks=c(-3.0,-1.0,-0.5,0,0.5,1.0,1.5,3.0))
colcode <- findColours(class,plotclr)
leg <- leg <-gsub("\\["," ",names(attr(colcode, "table")))
leg <-gsub("\\]"," ",leg)
leg <-gsub("\\)"," ",leg)
leg <-gsub("\\,"," : ",leg)

dev.new()
plot(lonlat_mean_tempdiff.s$long,lonlat_mean_tempdiff.s$lat,col=colcode,pch=15,main='meantempdiff.summer')
legend('bottomright', legend=leg,fill=attr(colcode, "palette"))


#############################
#############################

first.year <- 1971
last.year <- 2000

firstday<-avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year]
lastday<-avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday)),count=c(nlon,nlat,(lastday-firstday)+1))
tmp.array[tmp.array == fillvalue$value] <- NA


winter.mean.c <- array(NA,c(nlon,nlat,(last.year-first.year)+1))

winter.m <- c(1,2,3,10,11,12)
# first.toget <-1

for( y in first.year:last.year){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	winter.mean.c[,,y-(first.year-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}

winter_climate <- apply(winter.mean.c,c(1,2),mean,na.rm=T)



# interest period
first.year.p <- 2001
last.year.p <- 2014

firstday.p <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==first.year.p]
lastday.p <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==last.year.p]

tmp.array <- ncvar_get(tw.ncdf,"tg",start=c(1,1,which(day_vals==firstday.p)),count=c(nlon,nlat,(lastday.p-firstday.p)+1))
tmp.array[tmp.array == fillvalue$value] <- NA

winter.mean.p <- array(NA,c(nlon,nlat,(last.year.p-first.year.p)+1))

winter.m <- c(1,2,3,10,11,12)

for( y in first.year.p:last.year.p){

	# dayone <- avg_temp_transect$julianday[avg_temp_transect$day==1&avg_temp_transect$month==1&avg_temp_transect$year==y]
	# dayend <- avg_temp_transect$julianday[avg_temp_transect$day==31&avg_temp_transect$month==12&avg_temp_transect$year==y]
		
	# last.toget <-first.toget+(dayend-dayone)-1


	winter.mean.p[,,y-(first.year.p-1)] <- apply(tmp.array[,,avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday.p+1],c(1,2),mean, na.rm=T)
	
	print(c(head(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday.p+1,1),tail(avg_temp_transect$julianday[avg_temp_transect$year==y & avg_temp_transect$month %in% winter.m]-firstday.p+1,1)))
	# first.toget <- first.toget+(dayend-dayone)

}


clim_diff.w <- array(NA,dim=dim(winter.mean.p))

for (y in 1:dim(winter.mean.p)[3]){
clim_diff.w[,,y] <- winter.mean.p[,,y] - winter_climate
}

mean_tempdiff.w <- apply(clim_diff.w,c(1,2),mean,na.rm=T)

mean_tempdiff.w[mean_tempdiff.w < -3 | mean_tempdiff.w > 3] <- NA

lonlat <- expand.grid(lon,lat)
longlat_nona <- lonlat[!is.na(as.vector(mean_tempdiff.w[,])),]
lonlat_mean_tempdiff.w <- data.frame(cbind(longlat_nona,as.vector(mean_tempdiff.w[,])[!is.na(as.vector(mean_tempdiff.w[,]))]))
names(lonlat_mean_tempdiff.w) <- c("long","lat","tempdiff")


nclr <- 9
plotclr <- blue2red(nclr) 
class <-  classIntervals(lonlat_mean_tempdiff.w$tempdiff,n=nclr,style="fixed", fixedBreaks=c(-3.0,-1.0,-0.5,0,0.5,1.0,1.5,3.0))
colcode <- findColours(class,plotclr)
leg <- leg <-gsub("\\["," ",names(attr(colcode, "table")))
leg <-gsub("\\]"," ",leg)
leg <-gsub("\\)"," ",leg)
leg <-gsub("\\,"," : ",leg)

dev.new()
plot(lonlat_mean_tempdiff.w$long,lonlat_mean_tempdiff.w$lat,col=colcode,pch=15,main='meantempdiff.winter')
legend('bottomright', legend=leg,fill=attr(colcode, "palette"))





