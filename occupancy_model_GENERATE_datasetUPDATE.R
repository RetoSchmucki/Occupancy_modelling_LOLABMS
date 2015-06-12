# AUTHOR: RETO SCHMUCKI
# UPDATED: 28/01/2015
# ==========================
# Script to extract species count data 

# This script generate a .RData file with all object needed to run an occupancy model with 
# no covariate across the full European range where BMS have data for the period 2000-2012. 
 
#***************************************
# MODULE 1
# extract information from the database
# and structure the entry data file with
# presence/absence and covariates
#***************************************

# UNIVOLTINE
# species_list <- c('Anthocharis cardamines','Apatura iris','Aphantopus hyperantus',"Araschnia levana",'Argynnis aglaja','Argynnis paphia','Callophrys rubi','Carterocephalus palaemon','Favonius quercus','Gonepteryx rhamni','Hipparchia semele','Maniola jurtina','Nymphalis polychloros','Thymelicus lineola','Coenonympha tullia','Hesperia comma','Polyommatus coridon','Pyronia tithonus','Thymelicus sylvestris')
# voltinism <- "UNI"

# MULTIVOLTINE
 species_list <- c('Aglais io','Aglais urticae','Boloria selene','Celastrina argiolus','Colias crocea','Limenitis camilla','Lycaena phlaeas','Ochlodes sylvanus','Papilio machaon','Pararge aegeria','Pieris napi','Pieris rapae','Plebejus argus','Polygonia c-album','Polyommatus icarus')
 voltinism <- "MULTI"

library(RODBC)
library(rgdal)
library(lattice)

if(Sys.info()[["nodename"]] =='reto-Precision-T1650'){
uid_pwd <- read.csv("/home/reto/Dropbox/LOLA_BMS/lola_connection.csv",stringsAsFactors = F,header=F) 
} else { 
uid_pwd <- read.csv("/User/retoschmucki/Dropbox/LOLA_BMS/lola_connection.csv",stringsAsFactors = F,header=F) 
}

ch.lola <- odbcConnect("LOLAdb", uid=uid_pwd[1,2], pwd=uid_pwd[2,2],case="postgresql") 

# LOCATE flight curve file
if (Sys.info()[['nodename']] == 'reto-Precision-T1650') {
flight_curve_folder <- "/home/reto/Documents/LOLA_BMS/flight_curves_update2014/"
bd  <- "/usr/bin/jags"
OCC_wd <- "/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEWJUNE"
} else {
if (Sys.info()[['nodename']] == 'reto-OptiPlex-745') {
flight_curve_folder <- "/home/reto/Documents/LOLA_BMS/flight_period_curves/"
bd  <- "/usr/bin/jags"
OCC_wd <- "/home/reto/Documents/LOLA_BMS/Occupancy_modelling"
} else {
flight_curve_folder <- "/Users/retoschmucki/Documents/LOLA_project/Occupancy Model project/all_flight_period_curves/"
bd  <- "/usr/local/bin/JAGS"
OCC_wd <- "/Users/retoschmucki/Documents/LOLA_project/Occupancy Model Project/Occupancy_modellingLOLA_aixUNI"}} # on my MAC


###########################################

#MCMC parameters
ni <- 15000
nb <- 10000
nt <- 10

# select model
# =====================================

modelnbr <- 1

jags_model <- paste("OCCmodel_",modelnbr,"_jagscript_tidev_zone.txt",sep="")

country_of_flight <- noquote(c("\'ES\',\'FR\',\'UK\',\'NL\',\'DE\',\'FI\'")) # LINE 240

start.year <- 2000
end.year <- 2014
min_nbr_year <- 3

GRIDsize <- 100000
nbr_siteperGRID <- 30

slice_sd <- 0.67 #core slice extend mesured in Tidev SD 

year_toget <- noquote(c(paste(as.character(c(start.year:end.year)),collapse=",")))

if (start.year==2000) country <- noquote(c("\'ES\',\'UK\',\'NL\',\'FI\'"))
if (start.year==2006) country <- noquote(c("\'ES\',\'FR\',\'UK\',\'NL\',\'DE\',\'FI\'"))

folder <- paste("countries",voltinism,start.year,"sd",gsub("0.","",slice_sd),sep="_")

source("/home/reto/Dropbox/LOLA_BMS/Occupancy Modelling Project/R-Scripts/Occupancy_modelling_LOLABMS/OCC_model_dataGENERATION_functionUPDATE.R")

# build a grid with bms transect distribution

site_visit <- get.site_visit(ch.sql=ch.lola,country=country,year_toget=year_toget)

site_visit_sub <- site_subset.min_year(dataset=site_visit,visit_year=site_visit$visit_year,site_id=site_visit$transect_id,min_nbr_year=min_nbr_year)

laea.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
site_coord <- site_visit_sub[,c("transect_id","lat","long")]

site_subset_spatial <- site_subset.sp_grid(site_coord=site_coord,GRIDsize=GRIDsize,nbr_siteperGRID=nbr_siteperGRID,projectionCRS=laea.proj)

for (occsp in species_list){

# occsp <- species_list[3]
# set parameter for occupancy model
# ==============================================================

species <- occsp

site_with_obs <- get.site_with_obs(ch.sql=ch.lola,species=species,country=country,year_toget=year_toget)

#map of sites where the species was observed
plot(site_subset_spatial$long,site_subset_spatial$lat,col="blue")
points(site_subset_spatial$long[as.character(site_subset_spatial$transect_id) %in% as.character(site_with_obs$transect_id)],site_subset_spatial$lat[as.character(site_subset_spatial$transect_id) %in% as.character(site_with_obs$transect_id)],col="red")

# Retrieve and organized data
# =============================================================

year_sampled <- unique(site_visit$visit_year[as.character(site_visit$transect_id) %in% as.character(site_subset_spatial$transect_id)])
year_sampled <- year_sampled[year_sampled >= 1980 & !is.na(year_sampled)]
year_sampled <- year_sampled[order(year_sampled)]

site_visit_sub <- site_visit[as.character(site_visit$transect_id) %in% as.character(site_subset_spatial$transect_id),]

# 1.2 Generate a sampling date dataframe with Julian day
sampling_date <- data.frame(climate_zone=site_visit_sub$climate_zone,site_year=paste(site_visit_sub$transect_id,site_visit_sub$visit_year,sep='_'),date=paste(site_visit_sub$visit_day,site_visit_sub$visit_month,site_visit_sub$visit_year,sep='/'),year=site_visit_sub$visit_year)
sampling_date$julian_day <- strptime(sampling_date$date, "%d/%m/%Y")$yday+1

# =========================================================================
# TRIM SAMPLING DATE FOR THE SPECIES CLOSURE WINDOW DEFINED FROM THE FLIGHT
# PERIOD USING THE GAM DEVELOPED IN DENNIS ET AL. 2013 
# The flight period is computed a the scale of the bioclimatic region for
# each year
# =========================================================================

# extract the flight period curve computed at the bioclimatic region scale
sp_pheno_data <- read.csv(paste(flight_curve_folder,gsub(" ","_",species,fixed=TRUE),"_",gsub(",","_",gsub("\'","",as.character(country_of_flight))),".csv",sep=""), header=T, stringsAsFactors=F)
colnames(sp_pheno_data) <- c("CLIMATIC_ZONE","SPECIES","YEAR","WEEK","DAYNO","DAYNO_adj","NM")

window_obs <- get.obs_window(sp_pheno_data=sp_pheno_data,prct_cut=0.10)

# SPLIT sampling site and year per level of information available about the phenology; => 1.pheno for that year available, 2.pheno for that year missing, 3. pheno for that zone missing
zone_with_pheno <- unique(paste(window_obs$climate_zone,window_obs$Year))
sampling_date_withpheno_year <- sampling_date[paste(sampling_date$climate_zone,sampling_date$year) %in% zone_with_pheno,]
sampling_date_withpheno_noyear <- sampling_date[(paste(sampling_date$climate_zone,0) %in% zone_with_pheno) & !(paste(sampling_date$climate_zone,sampling_date$year) %in% zone_with_pheno),]
sampling_date_withoutpheno <- sampling_date[!sampling_date$climate_zone %in% unique(window_obs$climate_zone),]

# CASE 1. with precise phenology info
if(dim(sampling_date_withpheno_year)[1]>0){
	sampling_date_withpheno_year$NM <- NA
	match_index <- match(paste(sampling_date_withpheno_year$climate_zone,sampling_date_withpheno_year$year,sampling_date_withpheno_year$julian_day,sep='_'),paste(window_obs$climate_zone,window_obs$Year,window_obs$Dayno,sep='_'))
	sampling_date_withpheno_year$NM <- window_obs$NM[match_index]}

# CASE 2. with phenology info averaged across years (specific year missing)
if(dim(sampling_date_withpheno_noyear)[1]>0){
	sampling_date_withpheno_noyear$NM <- NA
	match_index <- match(paste(sampling_date_withpheno_noyear$climate_zone,0,sampling_date_withpheno_noyear$julian_day,sep='_'),paste(window_obs$climate_zone,0,window_obs$Dayno,sep='_'))
	sampling_date_withpheno_noyear$NM <- window_obs$NM[match_index]}

# CASE 3. no phenology available for the region, use the average in the nearest bioclimatic region
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step 1. identify closest zone
if(dim(sampling_date_withoutpheno)[1]>0){
	sampled_zone <- data.frame(zonename=unique(as.character(sampling_date_withoutpheno$climate_zone[!is.na(sampling_date_withoutpheno$climate_zone)])),zonenumber=match(substr(unique(as.character(sampling_date_withoutpheno$climate_zone[!is.na(sampling_date_withoutpheno$climate_zone)])),1,1),toupper(letters)))
	pheno_zone_available <- data.frame(zonename=unique(as.character(window_obs$climate_zone)),zonenumber=match(substr(unique(as.character(window_obs$climate_zone)),1,1),toupper(letters)))
	sampled_zone$closest_z <- NA

	for (i in 1:length(sampled_zone$zonenumber)){
		sampled_zone$closest_z[i] <- as.character(pheno_zone_available$zonename)[which(abs(sampled_zone$zonenumber[i]-pheno_zone_available$zonenumber)==min(abs(sampled_zone$zonenumber[i]-pheno_zone_available$zonenumber)))]}

# step 2. merge info from the closest zone
	match_index <- match(as.character(sampling_date_withoutpheno$climate_zone),sampled_zone$zonename)
	sampling_date_withoutpheno$nearest_zone <- NA
	sampling_date_withoutpheno$nearest_zone <- sampled_zone$closest_z[match_index]
	sampling_date_withoutpheno$NM <- NA
	match_index <- match(paste(sampling_date_withoutpheno$nearest_zone,0,sampling_date_withoutpheno$julian_day,sep='_'),paste(window_obs$climate_zone,0,window_obs$Dayno,sep='_'))
	sampling_date_withoutpheno$NM <- window_obs$NM[match_index]
	sampling_date_withoutpheno <- sampling_date_withoutpheno[,c("climate_zone","site_year","date","year","julian_day","NM")]}
#======
# END

# RECOMBINE all transects
sampling_date_withphenowindow <- rbind(sampling_date_withpheno_year,sampling_date_withpheno_noyear,sampling_date_withoutpheno)

# CLIP date within the observation window (NA)
sampling_date_withphenowindow <- sampling_date_withphenowindow[!is.na(sampling_date_withphenowindow$NM),]

# ORDER per transect_id and year
sampling_date_withphenowindow$transect_id <- substr(as.character(sampling_date_withphenowindow$site_year),1,nchar(as.character(sampling_date_withphenowindow$site_year))-5)
sampling_date_withphenowindow <- sampling_date_withphenowindow[order(sampling_date_withphenowindow$transect_id,sampling_date_withphenowindow$year,sampling_date_withphenowindow$julian_day),]

# ====================================================================
# SITES AND DATES CORRESPONDING TO THE CLOSURE WINDOW FOR SPECIES XXX
# ====================================================================

# 1.3 Generate a list of sampling date per site_year and extract the maximum length, but less than 27 weeks
nbr_samp_peryear <- split(sampling_date_withphenowindow$julian_day,sampling_date_withphenowindow$site_year)
ex_transect <- names(unlist(lapply(nbr_samp_peryear,length)))[as.numeric(unlist(lapply(nbr_samp_peryear,length)))>27]
sampling_date <- sampling_date_withphenowindow[!sampling_date_withphenowindow$site_year %in% ex_transect,c("climate_zone","site_year","date","year","julian_day","NM","transect_id")]

# MAXIMUM number of count observed and list of site sampled
max.nbr_count <- max(unlist(lapply(split(sampling_date$julian_day,sampling_date$site_year),length)))
site_sampled <- unique(sampling_date$transect_id)


# ==================================================================
# MODULE 2
# Impute the observed values in the Count
# and the Julian day matrix
# ==================================================================


species_visit <- get.species_count(ch.sql=ch.lola,country=country,species=species,year_toget=year_toget)

transect_length <- get.transect_length(ch.sql=ch.lola,country=country)

transect_TempI <- get.transect_tempindex(ch.sql=ch.lola,country=country)

Species_STI <- get.sti(ch.sql=ch.lola,species=species)

# COMPUTE temperature deviation for each transect
transect_TempI_dev <- data.frame(transect_id=transect_TempI$transect_id,TIDEV=(Species_STI$temp_mean-transect_TempI$manntmp)/Species_STI$temp_sd,coorx=transect_TempI$long,coory=transect_TempI$lat)

# COMPUTE sampling date with observation
species_sampling_date <- data.frame(site_year=paste(species_visit$transect_id,species_visit$visit_year,sep='_'),date=paste(species_visit$visit_day,species_visit$visit_month,species_visit$visit_year,sep='/'),year=species_visit$visit_year,climate_zone=species_visit$climate_zone)
species_sampling_date$julian_day <- strptime(species_sampling_date$date, "%d/%m/%Y")$yday+1
species_sampling_date$count <- species_visit$individual_count

# TRIM around the observation window with sampling_date 
species_sampling_date <- species_sampling_date[paste(as.character(species_sampling_date$site_year),species_sampling_date$julian_day,sep='_') %in% paste(as.character(sampling_date$site_year),sampling_date$julian_day,sep='_'),]

# ADD ZEROS in sampling visit when no individual where observed
zero_obs_date <- sampling_date[which(!paste(sampling_date$site_year,sampling_date$julian_day,sep="_") %in% paste(species_sampling_date$site_year,species_sampling_date$julian_day,sep="_")),c("site_year","date","year","climate_zone","julian_day")]
zero_obs_date$count <- 0
species_sampling_date <- rbind(species_sampling_date,zero_obs_date)

# ORDER species_sampling_date and merge NM
species_sampling_date <- species_sampling_date[order(as.character(species_sampling_date$site_year),species_sampling_date$julian_day),]
species_sampling_date <- merge(species_sampling_date,sampling_date[,c("site_year","julian_day","NM")],by=c("site_year","julian_day"))


#########################################
# BUILD EMPTY OBJECT TO RECEIVE THE DATA
#########################################

# GENERATE an empty[NA] count per visit matrix
count_mat <- data.frame(matrix(NA,ncol=max.nbr_count,nrow=(length(year_sampled)*length(site_sampled))))
names(count_mat) <- paste('Count',1:max.nbr_count,sep='')

# GENERATE an empty[NA] julian days matrix COVARIATE
day_mat <- data.frame(matrix(NA,ncol=max.nbr_count,nrow=(length(year_sampled)*length(site_sampled))))
names(day_mat) <- paste('Day',1:max.nbr_count,sep='')

# GENERATE an empty[NA] NM matrix COVARIATE
nm_mat <- data.frame(matrix(NA,ncol=max.nbr_count,nrow=(length(year_sampled)*length(site_sampled))))
names(nm_mat) <- paste('NM',1:max.nbr_count,sep='')

# GENERATE an empty[NA] transect length matrix COVARIATE
tl_mat <- data.frame(matrix(NA,ncol=1,nrow=(length(year_sampled)*length(site_sampled))))
names(tl_mat) <- paste('TL',1,sep='')

# GENERATE an empty[NA] transect latitude-longitude matrix COVARIATE
coord_mat <- data.frame(matrix(NA,ncol=2,nrow=(length(year_sampled)*length(site_sampled))))
names(coord_mat) <- paste(c('LAT','LONG'),1,sep='')

# GENERATE an empty[NA] transect temperature index deviance matrix COVARIATE
TempIdev_mat <- data.frame(matrix(NA,ncol=1,nrow=(length(year_sampled)*length(site_sampled))))
names(TempIdev_mat) <- paste('TIDEV',1,sep='')

# GENERATE a site_year matrix
site.year_mat <- data.frame(Site=rep(site_sampled,rep(length(year_sampled),length(site_sampled))),Year=year_sampled, site_year=paste(rep(site_sampled,rep(length(year_sampled),length(site_sampled))),year_sampled,sep='_'))


# ASSOCIATE the list of days with count and NM
julian_day.list <- split(species_sampling_date$julian_day,as.character(species_sampling_date$site_year))
species_count.list <- split(species_sampling_date$count,as.character(species_sampling_date$site_year))
NM.list <- split(species_sampling_date$NM,as.character(species_sampling_date$site_year))

# FEED the julian day covariate
max.len <- max(sapply(julian_day.list, length))
corrected.list <- lapply(julian_day.list, function(x) {c(x, rep(NA, max.len - length(x)))})
julian_mat <- do.call(rbind, corrected.list)
julian_mat <- julian_mat[row.names(julian_mat) %in% as.character(site.year_mat$site_year),]
row_index <- match(row.names(julian_mat),as.character(site.year_mat$site_year))
day_mat[row_index,]<-julian_mat

# FEED the observed count and convert to presence absence
max.len <- max(sapply(species_count.list, length))
corrected.list <- lapply(species_count.list, function(x) {c(x, rep(NA, max.len - length(x)))})
obs_mat <- do.call(rbind, corrected.list)
obs_mat <- obs_mat[row.names(obs_mat) %in% as.character(site.year_mat$site_year),]
row_index <- match(row.names(obs_mat),as.character(site.year_mat$site_year))
count_mat[row_index,]<-obs_mat
count_mat[count_mat > 0] <- 1

# FEED the NM (flight curve)
max.len <- max(sapply(NM.list, length))
corrected.list <- lapply(NM.list, function(x) {c(x, rep(NA, max.len - length(x)))})
NM_mat <- do.call(rbind, corrected.list)
NM_mat <- NM_mat[row.names(NM_mat) %in% as.character(site.year_mat$site_year),]
row_index <- match(row.names(NM_mat),as.character(site.year_mat$site_year))
nm_mat[row_index,]<-NM_mat

# FEED transect length covariate
row_index <- match(as.character(site.year_mat$Site),as.character(transect_length$transect_id))
site_length <- transect_length$transect_length[row_index]
tl_mat[,1] <- site_length

# FEED latitude covariate
row_index <- match(as.character(site.year_mat$Site),as.character(transect_length$transect_id))
site_lat <- transect_length$lat[row_index]
coord_mat[,1] <- site_lat

# FEED longitude covariate
row_index <- match(as.character(site.year_mat$Site),as.character(transect_length$transect_id))
site_long <- transect_length$long[row_index]
coord_mat[,2] <- site_long

# FEED STI covariate
row_index <- match(as.character(site.year_mat$Site),as.character(transect_TempI_dev$transect_id))
site_TIDEV <- transect_TempI_dev$TIDEV[row_index]
TempIdev_mat[,1] <- site_TIDEV 

# combine data in one dataframe
entry_mat <- cbind(site.year_mat[,-3],count_mat,day_mat,nm_mat,tl_mat,coord_mat,TempIdev_mat)
entry_mat <- entry_mat[!is.na(entry_mat$LAT1),]
entry_mat <- entry_mat[substr(as.character(entry_mat$Site),1,2)%in%unique(substr(as.character(site_with_obs$transect_id),1,2)),]
entry_mat <- entry_mat[order(entry_mat$Year,entry_mat$Site),]

country_with_obs <- unique(substr(as.character(site_with_obs$transect_id),1,2))

# *******************************
# MODULE 3
# Define working directory and
# generate files names
# *******************************
 
setwd(OCC_wd)

if (!file.exists(file.path(OCC_wd,folder,gsub(" ","_",species)))){
dir.create(file.path(OCC_wd,folder,gsub(" ","_",species)),recursive = TRUE)
} else {cat(paste("file",file.path(OCC_wd,folder),"exist! \n"))}
pad <- file.path(OCC_wd,folder,gsub(" ","_",species))

setwd(pad)

#####################################################
# Data management 
#####################################################

# Read data

site_lat_mat <- data.frame(site=as.character(entry_mat[entry_mat$Year == start.year,"Site"]),latitude=entry_mat[entry_mat$Year == start.year,"LAT1"],longitude=entry_mat[entry_mat$Year == start.year,"LONG1"],tidev=entry_mat[entry_mat$Year == start.year,"TIDEV1"])

tidev_mat_zone <- matrix(0,nrow=length(site_lat_mat$site),ncol=4)

tidev_mat_zone[site_lat_mat$tidev < (-1*slice_sd),1] <- 1
tidev_mat_zone[site_lat_mat$tidev >= (-1*slice_sd) & site_lat_mat$tidev < 0 ,2] <- 1
tidev_mat_zone[site_lat_mat$tidev >= 0 & site_lat_mat$tidev <= slice_sd ,3] <- 1
tidev_mat_zone[site_lat_mat$tidev > slice_sd,4] <- 1

tidev_matdf <-  as.data.frame(tidev_mat_zone)
tidev_matdf <- cbind(site_lat_mat,tidev_matdf)
tidev_matdf$suitable <- 0
tidev_matdf$suitable[as.character(tidev_matdf$site) %in% as.character(site_with_obs$transect_id)] <- 1

sum_col <- colSums(is.na(entry_mat))
nr <- min(which(sum_col[grep("Count", names(sum_col))] == dim(entry_mat)[1])) # nr of visits
if(is.infinite(nr)){ nr <- length(grep("Count", names(entry_mat)))}

site <- entry_mat[,"Site"]
usite <- unique(site)
(nsite <- length(usite))
nvisit <- nr
uyear <- unique(entry_mat[,"Year"])
minyear  <- min(uyear)
nyear <- length(uyear)
time_interval <- diff(uyear)[1]

(nbr_site.per_zone <- data.frame(Z1=colSums(tidev_mat_zone)[1],Z2=colSums(tidev_mat_zone)[2],Z3=colSums(tidev_mat_zone)[3],Z4=colSums(tidev_mat_zone)[4]))
(nbr_site.per_zone_suitable <- data.frame(Z1=sum(tidev_matdf$suitable==1&tidev_matdf$V1==1),Z2=sum(tidev_matdf$suitable==1&tidev_matdf$V2==1),Z3=sum(tidev_matdf$suitable==1&tidev_matdf$V3==1),Z4=sum(tidev_matdf$suitable==1&tidev_matdf$V4==1)))

# Write lists of years and sites for species under consideration, for post-processing
# =====================================================================================

Yearnumeric <- 1:nyear
Year <- uyear
temp1 <- cbind(Yearnumeric, Year)

#####################################################################################
# Convert calender years into year 1, 2 etc. required for calculation of trends
# The sequence of years need to be complete and regular
#####################################################################################

oyear <- (uyear-(minyear-1))
year <- (uyear-(minyear-1))
sumX  <- sum(oyear)
sumX2 <- sum(oyear*oyear)
maxyear <- max(oyear)


# Fill matrices with data
M <- array(NA,dim=c(nsite,nvisit,nyear))
NM   <-array(NA,dim=c(nsite,nvisit,nyear))

# Covariate columns in inputfile (one for each visit)
firstcount <- which(names(entry_mat)== "Count1") 
lastcount <- which(names(entry_mat)== paste("Count", nr, sep=""))

firstNM <- which(names(entry_mat)== "NM1") 
lastNM <- which(names(entry_mat)== paste("NM", nr, sep=""))


for (t in 1:nyear) {
   kp <- entry_mat[,"Year"]==((minyear-time_interval)+(t*time_interval)) 	

dcount <- entry_mat[kp, firstcount:lastcount] 	# observations in M 
   M[,,t] <- as.matrix(dcount)

dnm<-entry_mat[kp, firstNM:lastNM]			# covariate NM from GAM
	NM[,,t]<-as.matrix(dnm)
}

# CONSTANT covariates
TL   <-array(entry_mat[1:nsite,"TL1"],dim=c(nsite,1,1))
TIDEV <- array(entry_mat[1:nsite,"TIDEV1"],dim=c(nsite,1,1))
SUITABLE <- matrix(tidev_matdf$suitable,nrow=nsite,ncol=1)

nbr_suit_zone <- c(sum(tidev_matdf$V1 * tidev_matdf$suitable),sum(tidev_matdf$V2 * tidev_matdf$suitable),sum(tidev_matdf$V3 * tidev_matdf$suitable),sum(tidev_matdf$V4 * tidev_matdf$suitable))

#######################################
# Standardize all covariates
#######################################

cat("\n*** Values used for covariate standardisation *** \n")

mn1 <-mean(NM,na.rm=TRUE)
sd1 <-sqrt(var(NM[1:length(NM)],na.rm=TRUE))
cat("Mean and sd for standardising NM:", mn1, sd1, "\n")
NM <- (NM-mn1)/sd1

mn2 <-mean(TL,na.rm=TRUE)
sd2 <-sqrt(var(TL[1:length(TL)],na.rm=TRUE))
cat("Mean and sd for standardising TL:", mn2, sd2, "\n")
TL <- (TL-mn2)/sd2

mn3 <-mean(TIDEV,na.rm=TRUE)
sd3 <-sqrt(var(TIDEV[1:length(TIDEV)],na.rm=TRUE))
cat("Mean and sd for standardising TIDEV:", mn3, sd3, "\n")
TIDEV <- (TIDEV-mn3)/sd3

# Impute means for missing values
NM[is.na(NM)] <- 0
TL[is.na(TL)] <- 0
TIDEV[is.na(TIDEV)] <- 0

covarinfo <- cbind(mn1,sd1,mn2,sd2,mn3,sd3) # table

# DEFINE JOB NAME (FILES NAMES)

# jobnaam <- gsub(" ","_",paste("jags_model",modelnbr,species_visit$fauna_europea_species[1],gsub(",","_",gsub("\'","",as.character(country))),sep='_'))
jobname <- gsub(" ","_",paste("jags_model",modelnbr,species_visit$fauna_europea_species[1],gsub(",","_",gsub("\'","",as.character(country))),sep='_'))

jobnamein   <-  jobname                               
jobnameout  <-  paste(jobname, "outdaytrimmed.csv",sep="")     	# standard output
jobnameoutc <-  paste(jobname, "outc.csv",sep="")    		# additional info
jobnameoutj <-  paste(jobname, "year.csv", sep="")   		# list of years
jobnameoutp <-  paste(jobname, "plot.csv", sep="")   		# list of sites
jobnameoutcoord <-  paste(jobname, "coord.csv", sep="")   	# list of sites coordinates
jobnameouts <-  paste(jobname, "sims.csv", sep="")   		# simulations in vertical format
jobnamezonecoord <-  paste(jobname, "coord_zone.csv", sep="")   		# simulations in vertical format

#####################################
# Define the MODEL WITHOUT COVARIATE
#####################################

#if(start.year == 2006){
#	if(Sys.info()[["nodename"]] =='reto-Precision-T1650'){
#	source('/home/reto/Dropbox/LOLA_BMS/Occupancy Modelling Project/R-Scripts/OCC_model_jagscript_tidev_zoneMEAN.r')
#	} else {
#	source("/Users/retoschmucki/Dropbox/LOLA_BMS/Occupancy Modelling Project/R-Scripts/OCC_model_jagscript_tidev_zoneMEAN.r")
#	}

#} else {
#	if(Sys.info()[["nodename"]] =='reto-Precision-T1650'){
	source('/home/reto/Dropbox/LOLA_BMS/Occupancy Modelling Project/R-Scripts/Occupancy_modelling_LOLABMS/randomPhiGam_OCCmodel.R')
#	} else {
#	source("/Users/retoschmucki/Dropbox/LOLA_BMS/Occupancy Modelling Project/R-Scripts/Occupancy_modelling_LOLABMS/randomPhiGam_OCCmodel.R")
#	}
#}
# Read priors and iteration settings 
psi_min	     <- 0
psi_max	     <- 1
bnm_min	     <- -10
bnm_max	     <- 10
btl_min	     <- -10
btl_max	     <- 10
bphi_min     <- -10
bphi_max     <- 10
bgam_min     <- -10
bgam_max     <- 10

write.csv(temp1, file = jobnameoutj, row.names = FALSE)
write.csv(site_lat_mat,file=jobnameoutcoord,row.names=FALSE)
write.csv(covarinfo, file = jobnameoutc, row.names=FALSE)
write.csv(tidev_matdf,file= jobnamezonecoord , row.names=FALSE)

save("M","nsite","nvisit","nyear",
"sumX","sumX2","maxyear","NM","TL","TIDEV","tidev_mat_zone","SUITABLE","nbr_suit_zone","psi_min","psi_max","bnm_min","bnm_max","btl_min","btl_max","bphi_min","bphi_max","bgam_min","bgam_max","species","country_of_flight","country","start.year","end.year","year_toget","min_nbr_year","GRIDsize","nbr_siteperGRID","slice_sd","country_with_obs",
"jobnamein","jobnameout","jobnameoutc","jobnameoutj","jobnameoutp","jobnameoutcoord","jobnameouts","jobnamezonecoord",file=paste(species,".RData",sep=""))

}
