library(arm)

voltinisme <- "UNI"
#voltinisme <- "MULTI"

sp_list <- list.files(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/",sep=""))
sp_list <- sp_list[!sp_list %in% c("launch_multiple_uni.sh", "sub.sh", "sub.sh~", "reto.sh", "launch_jags_OCMtidev_zoneLOCAL.R", "launch_multiple_uni2.sh", "launch_multiple_multi2.sh", "launch_multiple_multi.sh")]


# sink(file = 'fooLONG.txt', append = T)

est_df <- data.frame()

for (species1 in sp_list) {
  
  site_coord <- read.csv(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/", 
                               species1, "/results/jags_model_1_", species1, "_ES_UK_NL_FIcoord_zone.csv", sep = ""), header = T, stringsAsFactors = F)
  load(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/", 
             species1, "/", gsub("_", " ", species1), ".RData", sep = ""))
  load(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/", 
             species1, "/results/jagsoutput.Rdata", sep = ""))
  

  occupancy_zone <- data.frame()
  
  for (z in 1:4) {
    
    quant_y <- data.frame()
    
    for (y in 2:13) {
      quant <- (out$BUGSoutput$sims.list$psi.fs_z[, z, y] - mean(out$BUGSoutput$sims.list$psi.fs_z[,z, -1]))/sd(out$BUGSoutput$sims.list$psi.fs_z[, z, -1])
      quant_y <- rbind(quant_y, quant)
    }
    
    if (sum(quant_y) == 0 | sum(!is.na(quant_y)) == 0) {
      est_df <- rbind(est_df, data.frame(species = species1, z = z, mean_est = 0, lci = 0, uci = 0))
    } else {
      est_dist <- round(quantile(apply(quant_y, 2, function(x) lm(x ~ c(1:12) - 1)$coefficients), c(0.025, 0.5, 0.975)), 4)
      est_df <- rbind(est_df, data.frame(species = species1, z = z, mean_est = est_dist[2], lci = est_dist[1], uci = est_dist[3]))
    }
  }  # zone
}  # species

# multipanel plot

est_df <- merge(est_df, data.frame(species = as.character(est_df$species[est_df$z == 4])[order(est_df$mean[est_df$z == 4])], ord = 1:length(est_df$species[est_df$z == 
                4])), by = c("species"))
est_df <- est_df[order(est_df$ord, est_df$z), ]

par(mfrow = c(5, 1))

for (z in 4:1) {
  
  if (z == 1) {
    par(mar = c(2, 2, 1, 0))
  } else {
    par(mar = c(2, 2, 1, 0))
  }
  
  plot(1:length(est_df$mean[est_df$z == z]),est_df$mean[est_df$z == z], ylim = c(min(est_df$lci[est_df$z == z]), max(est_df$uci[est_df$z == z])), 
       xaxt = "n", xlab = "", ylab = "estimate", cex = 0.6)
  
  if (z == 1) 
    axis(1, at = 1:length(est_df$mean[est_df$z == z]), labels = gsub("_", " ", as.character(est_df$species[est_df$z == z])), tick = FALSE, 
         las = 2, cex = 0.8)
  
  segments(y0 = est_df$lci[est_df$z == z], x0 = c(1:length(est_df$mean[est_df$z == z])), y1 = est_df$uci[est_df$z == z], x1 = c(1:length(est_df$mean[est_df$z == z])))
  abline(h = 0, lty = 2, col = "grey")
}

# plot(est_df$mean[est_df$z == 2], 1:length(est_df$mean[est_df$z == 2]), xlim = c(min(est_df$lci[est_df$z == 2]), max(est_df$uci[est_df$z == 2])))
# segments(x0 = est_df$lci[est_df$z == 2], y0 = c(1:length(est_df$mean[est_df$z == 1])), x1 = est_df$uci[est_df$z == 1], y1 = c(1:length(est_df$mean[est_df$z == 2])))
# abline(v = 0, lty = 2, col = "grey")
# 
# plot(occupancy_zone1$year, occupancy_zone1$mean, ylim = c(min(occupancy_zone1[, -1]), max(occupancy_zone1[, -1])))
# abline(h = 0, col = "grey", lty = 2)
# points(occupancy_zone1$year, occupancy_zone1$lci, col = "red", lty = 2, type = "l")
# points(occupancy_zone1$year, occupancy_zone1$uci, col = "red", lty = 2, type = "l")
