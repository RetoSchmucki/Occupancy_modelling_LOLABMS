

# 1.
# FUNCTION to retrieve site visits from the LOLA_BMS database
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.site_visit <- function(ch.sql,country,year_toget){

site_visit_query<- paste("SELECT DISTINCT
all_species_count_visit.transect_id,
visit_year,
visit_month,
visit_day,
ST_Y(ST_Transform(geom,3035)) as LAT,
ST_X(ST_Transform(geom,3035)) as LONG,
all_bms_transects_coord.eco_zone_name as climate_zone
FROM
bms_update.all_species_count_visit
LEFT JOIN bms_update.all_bms_transects_coord ON all_species_count_visit.transect_id = all_bms_transects_coord.transect_id
WHERE
substring(all_species_count_visit.transect_id from 1 for 2) IN (",country,") AND
visit_year in (",year_toget,") AND
all_bms_transects_coord.altitude < 3000 AND
all_species_count_visit.monitoring_type = 1
ORDER BY
transect_id,
visit_year,
visit_month,
visit_day",sep="")

site_visit <- sqlQuery(ch.sql,site_visit_query,rows_at_time=1)

return(site_visit)
}
#========
# END


# 2.
# FUNCTION to subset site with a minimum number of monitoring year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

site_subset.min_year <- function(dataset,visit_year,site_id,min_nbr_year){

	site_nbr_year <- aggregate(site_visit$visit_year,by=list(site_visit$transect_id),function(x) sum(!duplicated(x)))
	names(site_nbr_year)<- c("transect_id","nbr_monitor_year")
	site_withmin_nbr_year <- site_nbr_year[site_nbr_year$nbr_monitor_year >= min_nbr_year,]
	site_visit_sub <- site_visit[as.character(site_visit$transect_id) %in% as.character(site_withmin_nbr_year$transect_id),]

return(site_visit_sub)
}
#======
# END


# 4.
###########
# FUNCTION to retrieve site with species count
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.site_with_obs <- function(ch.sql,species,country,year_toget){

species_site_query <- paste("SELECT
distinct(all_species_count_visit.transect_id)
FROM
bms_update.all_species_count_visit
WHERE
substring(all_species_count_visit.transect_id from 1 for 2) IN (",country,") AND
fauna_europea_species = \'",as.name(species),"\' AND
visit_year in (",year_toget,")
ORDER BY
all_species_count_visit.transect_id",sep='')

site_with_obs <- sqlQuery(ch.sql,species_site_query,rows_at_time=1)

return(site_with_obs)
}

#=======
# END


# 5.
# FUNCTION to subset sites to homogenize spatial distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

site_subset.sp_grid <- function(site_coord,GRIDsize,nbr_siteperGRID,projectionCRS){

	require(rgdal)
	
	# build spatial point object
	poi <- site_coord
	poi <- poi[!duplicated(poi),]
	poi.sp <- poi
	coordinates(poi.sp) <- ~ long+lat
	proj4string(poi.sp) = CRS(projectionCRS)

	# build a grid based on the points
	bb <- bbox(poi.sp)
	cs <- c(GRIDsize,GRIDsize)
	cc <- bb[,1]+(cs/2)
	cd <- ceiling(diff(t(bb))/cs)
	grd <- GridTopology(cellcentre.offset=cc,cellsize=cs,cells.dim=cd)
	sp_grd <- SpatialGridDataFrame(grd,data=data.frame(id=1:prod(cd)),proj4string=CRS(projectionCRS))

	# get point per gridcell
	point.grid <- cbind(poi,over(poi.sp,sp_grd))
	nbr.pointpergrid <- aggregate(as.character(point.grid$transect_id),by=list(point.grid$id),function(x) length(x))
	
	# random subset for gridcell having more than the maximum per cell specified [nbr_siteperGRID]
	site.tokeep <- as.character(point.grid$transect_id)[!point.grid$id %in% nbr.pointpergrid$Group.1[nbr.pointpergrid$x > nbr_siteperGRID]] 
		for (cell in nbr.pointpergrid$Group.1[nbr.pointpergrid$x > nbr_siteperGRID]){
		site.tokeep <- c(site.tokeep,sample(as.character(point.grid$transect_id)[point.grid$id == cell],nbr_siteperGRID,replace=F))}


	poi.sub <- site_coord[as.character(site_coord$transect_id) %in% site.tokeep,c("transect_id","lat","long")]
	poi.sub <- poi.sub[!duplicated(poi.sub),]

return(poi.sub)}

#====================
# END


# 6.
# FUNCTION to get the observation window for each year and climate zone 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.obs_window <- function(sp_pheno_data,prct_cut=0.10){

	# add the mean phenology for each zone with at least one phenology curve
	avg_nm_cz <- aggregate(sp_pheno_data$NM,by=list(sp_pheno_data$CLIMATIC_ZONE, sp_pheno_data$DAYNO),function(x) mean(x,na.rm=T))
	names(avg_nm_cz) <- c("CLIMATIC_ZONE","DAYNO","NM")
	avg_nm_cz$YEAR <- 0
	avg_nm_cz$SPECIES <- species

	# add this to the original phenology data
	sp_pheno_data2 <- rbind(sp_pheno_data[,c("CLIMATIC_ZONE","SPECIES","YEAR","DAYNO","NM")],avg_nm_cz[,c("CLIMATIC_ZONE","SPECIES","YEAR","DAYNO","NM")])

window_obs <- data.frame()

	for (y in unique(sp_pheno_data2$YEAR)){
		pheno_data_year <- sp_pheno_data2[sp_pheno_data2$YEAR == y,]
		zonelist <- unique(pheno_data_year$CLIMATIC_ZONE)
		for(cz in unique(zonelist)){
			pheno_data_year_zone <- pheno_data_year[pheno_data_year$CLIMATIC_ZONE == cz,]
			day_window <- data.frame(climate_zone=cz,Year=y,Dayno=pheno_data_year_zone$DAYNO[pheno_data_year_zone$NM >= prct_cut*max(pheno_data_year_zone$NM)],NM=pheno_data_year_zone$NM[pheno_data_year_zone$NM >= prct_cut*max(pheno_data_year_zone$NM)])
			window_obs <- rbind(window_obs,day_window) 
		} #cz
	} #y

return(window_obs)
}
#======
# END


# 7.
# FUNCTION to retrieve species count from the LOLA-BMS database
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.species_count <- function(ch.sql,country,species,year_toget){

	species_visit_query <- paste("SELECT DISTINCT
		fauna_europea_species,
		all_species_count_visit.transect_id,
		visit_year,
		visit_month,
		visit_day,
		individual_count,
		all_bms_transects_coord.eco_zone_name as climate_zone
		FROM
		bms_update.all_species_count_visit
		LEFT JOIN bms_update.all_bms_transects_coord ON all_species_count_visit.transect_id = all_bms_transects_coord.transect_id
		WHERE
		substring(all_species_count_visit.transect_id from 1 for 2) IN (",country,") AND
		fauna_europea_species = \'",as.name(species),"\' AND
		visit_year in (",year_toget,")
		ORDER BY
		all_species_count_visit.transect_id,
		visit_year,
		visit_month,
		visit_day",sep='')

	species_visit<- sqlQuery(ch.sql,species_visit_query,rows_at_time=1)

return(species_visit)
}

#=======
# END


# 8.
# FUNCTION to retrieve transect length
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.transect_length <- function(ch.sql,country){

	transect_length_query <- paste("SELECT
		transect_id,
		transect_length,
		ST_Y(ST_Transform(geom,3035)) as LAT,
		ST_X(ST_Transform(geom,3035)) as LONG
		FROM
		bms_update.all_bms_transects_coord
		WHERE
		substring(transect_id from 1 for 2) IN (",country,")
		ORDER BY
		transect_id",sep='')

	transect_length <- sqlQuery(ch.sql,transect_length_query,rows_at_time=1)
	transect_length$lat <- as.numeric(as.character(transect_length$lat))
	transect_length$long <- as.numeric(as.character(transect_length$long))

return(transect_length)
}
#======
# END


# 9.
# FUNCTION to retrieve transect temperature index
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.transect_tempindex <- function(ch.sql,country){

	transect_TempI_query <- paste("SELECT
		p.transect_id,
		g.cgrsname,
		manntmp,
		ST_Y(ST_Transform(p.geom,3035)) as lat,
		ST_X(ST_Transform(p.geom,3035)) as long
		FROM
		bms_update.all_bms_transects_coord as p,
		climate_data.cgrs_latlon_grid as g,
		climate_data.cgrs_sti 
		WHERE
		ST_intersects(p.geom,g.geom) AND
		cgrs_sti.cgrsname = g.cgrsname ",sep='')

	transect_TempI <- sqlQuery(ch.sql,transect_TempI_query ,rows_at_time=1)

return(transect_TempI)
}
#======
# END

# 10.
# FUNCTION to retrieve species temperature index (sti) and sd of the sti
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#======
# END

