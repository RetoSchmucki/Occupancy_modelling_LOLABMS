library(arm)

voltinisme <- "UNI"
#voltinisme <- "MULTI"

sp_list <- list.files(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/",sep=""))
sp_list <- sp_list[!sp_list %in% c("launch_multiple_uni.sh", "sub.sh", "sub.sh~", "reto.sh", "launch_jags_OCMtidev_zoneLOCAL.R", "launch_multiple_uni2.sh", "launch_multiple_multi2.sh", "launch_multiple_multi.sh")]

beta_phi <- data.frame()

# for (species1 in sp_list) {
  
  species1 <- sp_list[2]
  
  load(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/", 
             species1, "/", gsub("_", " ", species1), ".RData", sep = ""))
  load(paste("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEW/ResultOCC/temporal_trend_gam_phi_estimates/countries_",voltinisme,"_2000_sd_67/", 
             species1, "/results/jagsoutput.Rdata", sep = ""))
  
  
round(quantile(apply(out$BUGSoutput$sims.list$beta.gam1, 1, function(x) lm(x ~ c(1:12) - 1)$coefficients),c(0.025, 0.5, 0.975)), 4)

  
data1 <- (cbind(1:12,c(t(out$BUGSoutput$sims.list$beta.gam1))))
  
plot(data1[,1],data1[,2])

summary(lm(data1[,2]~data1[,1]-1))



  