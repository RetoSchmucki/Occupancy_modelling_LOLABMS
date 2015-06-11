#################################################################################
# Model 5 with Temperature index deviation
#
#	logit(gam) <- alpha[t] + beta*TIDEV
#	logit(phi) <- alpha[t] + beta*TIDEV
#
#	========================================================================
#
#	Intercept vary with time.
#
#################################################################################

# Save BUGS description of the model to working directory
sink("OCCmodel_1_jagscript_tidev_zone.txt")
cat("
model {

# Priors 
# state model priors

psi ~ dunif(psi_min,psi_max)

# observation model priors 
for (t in 1:nyear) {
	alpha.p[t] ~ dnorm(mu.lp, tau.lp)	# p random year intercept
} #t

mu.lp ~ dnorm(0, 0.01)
tau.lp <- 1 / (sd.lp * sd.lp)
sd.lp ~ dunif(0, 5)

bnm.p ~ dunif(bnm_min, bnm_max)
btl.p ~ dunif(btl_min, btl_max)

#############
# PHI priors
#############

# beta.phi1 model priors
#========================
for (t in 2:nyear) {
	beta.phi1[t-1] <- muphi1 + betaphi1*y[t-1] + epsilonphi1[t-1]
  epsilonphi1[t-1] ~ dnorm(0,tauphi1)
  } #t

muphi1 ~ dnorm(0, 0.01)
betaphi1 ~ dnorm(0,0.01) I(-10,10)
sigmaphi1 ~ dunif(0,10)
tauphi1 <- pow(sigmaphi1,-2)
sigmaphi1_2 <- pow(sigmaphi1,2)


# beta.phi2 model priors
#========================
for (t in 2:nyear) {
	beta.phi2[t-1] <- muphi2 + betaphi2*y[t-1] + epsilonphi2[t-1]
  epsilonphi2[t-1] ~ dnorm(0,tauphi2)
  } #t

muphi2 ~ dnorm(0, 0.01)
betaphi2 ~ dnorm(0,0.01) I(-10,10)
sigmaphi2 ~ dunif(0,10)
tauphi2 <- pow(sigmaphi2,-2)
sigmaphi2_2 <- pow(sigmaphi2,2)


# beta.phi3 model priors
#========================
for (t in 2:nyear) {
	beta.phi3[t-1] <- muphi3 + betaphi3*y[t-1] + epsilonphi3[t-1]
  epsilonphi3[t-1] ~ dnorm(0,tauphi3)
  } #t

muphi3 ~ dnorm(0, 0.01)
betaphi3 ~ dnorm(0,0.01) I(-10,10)
sigmaphi3 ~ dunif(0,10)
tauphi3 <- pow(sigmaphi3,-2)
sigmaphi3_2 <- pow(sigmaphi3,2)


# beta.phi4 model priors
#========================
for (t in 2:nyear) {
	beta.phi4[t-1] <- muphi4 + betaphi4*y[t-1] + epsilonphi4[t-1]
  epsilonphi4[t-1] ~ dnorm(0,tauphi4)
  } #t

muphi4 ~ dnorm(0, 0.01)
betaphi4 ~ dnorm(0,0.01) I(-10,10)
sigmaphi4 ~ dunif(0,10)
tauphi4 <- pow(sigmaphi4,-2)
sigmaphi4_2 <- pow(sigmaphi4,2)

#############
# GAM priors
#############

# beta.gam1 model priors
#========================
for (t in 2:nyear) {
beta.gam1[t-1] <- mugam1 + betagam1*y[t-1] + epsilongam1[t-1]
epsilongam1[t-1] ~ dnorm(0,taugam1)
} #t

mugam1 ~ dnorm(0, 0.01)
betagam1 ~ dnorm(0,0.01) I(-10,10)
sigmagam1 ~ dunif(0,10)
taugam1 <- pow(sigmagam1,-2)
sigmagam1_2 <- pow(sigmagam1,2)


# beta.gam2 model priors
#========================
for (t in 2:nyear) {
beta.gam2[t-1] <- mugam2 + betagam2*y[t-1] + epsilongam2[t-1]
epsilongam2[t-1] ~ dnorm(0,taugam2)
} #t

mugam2 ~ dnorm(0, 0.01)
betagam2 ~ dnorm(0,0.01) I(-10,10)
sigmagam2 ~ dunif(0,10)
taugam2 <- pow(sigmagam2,-2)
sigmagam2_2 <- pow(sigmagam2,2)


# beta.gam3 model priors
#========================
for (t in 2:nyear) {
beta.gam3[t-1] <- mugam3 + betagam3*y[t-1] + epsilongam3[t-1]
epsilongam3[t-1] ~ dnorm(0,taugam3)
} #t

mugam3 ~ dnorm(0, 0.01)
betagam3 ~ dnorm(0,0.01) I(-10,10)
sigmagam3 ~ dunif(0,10)
taugam3 <- pow(sigmagam3,-2)
sigmagam3_2 <- pow(sigmagam3,2)


# beta.gam4 model priors
#========================
for (t in 2:nyear) {
beta.gam4[t-1] <- mugam4 + betagam4*y[t-1] + epsilongam4[t-1]
epsilongam4[t-1] ~ dnorm(0,taugam4)
} #t

mugam4 ~ dnorm(0, 0.01)
betagam4 ~ dnorm(0,0.01) I(-10,10)
sigmagam4 ~ dunif(0,10)
taugam4 <- pow(sigmagam4,-2)
sigmagam4_2 <- pow(sigmagam4,2)



#====================================================



# State model
#==============

for (i in 1:nsite){
	z[i,1] ~ dbern(psi)	# occupancy (binomial)

	for (t in 2:nyear){
		muZ[i,t] <- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gam[i,t-1]
		
		logit(phi[i,t-1]) <- beta.phi1[t-1]*tidev_mat_zone[i,1] + beta.phi2[t-1]*tidev_mat_zone[i,2] + beta.phi3[t-1]*tidev_mat_zone[i,3] + beta.phi4[t-1]*tidev_mat_zone[i,4]
		
		logit(gam[i,t-1]) <- beta.gam1[t-1]*tidev_mat_zone[i,1] + beta.gam2[t-1]*tidev_mat_zone[i,2] + beta.gam3[t-1]*tidev_mat_zone[i,3] + beta.gam4[t-1]*tidev_mat_zone[i,4]
		
		z[i,t] ~ dbern(muZ[i,t])
	} #t

} #i


# Observation model 
#====================

for (t in 1:nyear){
	for(i in 1:nsite){
		for(j in 1:nvisit) {
		Py[i,j,t]<- z[i,t]*p[i,j,t]

		logit(p[i,j,t]) <- alpha.p[t] + bnm.p*NM[i,j,t] + btl.p*TL[i,1,1]	# NM:flight period, TL:transect length

		M[i,j,t] ~ dbern(Py[i,j,t])

		Presi[i,j,t] <- abs(M[i,j,t]-p[i,j,t])

		y.new[i,j,t] ~ dbern(Py[i,j,t])

		Presi.new[i,j,t] <- abs(y.new[i,j,t]-p[i,j,t])

		} #j
	}#i
}#t

# Derived parameters state model
# =======================================================

# Finite sample occupancy
for (t in 1:nyear) {
	psi.fs[t] <- sum(z[1:nsite,t])/nsite
}

# Finite sample occupancy ZONE
for (zo in 1:4){
	for (t in 1:nyear) {
	psi.fs_z[zo,t] <- sum(z[1:nsite,t] * tidev_mat_zone[1:nsite,zo]) / (sum(tidev_mat_zone[1:nsite,zo]) + 0.00000001)
	}
}

for (zo in 1:4){
	for (t in 1:nyear) {
	psi.suit_z[zo,t] <- sum(z[1:nsite,t] * tidev_mat_zone[1:nsite,zo] * SUITABLE[1:nsite,1]) / (sum(tidev_mat_zone[1:nsite,zo] * SUITABLE[1:nsite,1]) + 0.00000001)
	}
}

# Finite sample persistence ZONE suitable long

for (zo in 1:4){
	for (t in 1:nyear-1) {
		for (i in 1:nsite){
		OC_T1phi[zo,i,t] <- z[i,t] * tidev_mat_zone[i,zo] * SUITABLE[i,1]
		OC_T2phi[zo,i,t] <- z[i,t+1] * tidev_mat_zone[i,zo] * SUITABLE[i,1]
		}
	phi.suit_z[zo,t] <- sum((OC_T1phi[zo,1:nsite,t] + OC_T2phi[zo,1:nsite,t]) == 2) / (sum(OC_T1phi[zo,1:nsite,t]) + 0.00000001)
	}
}

# ==========================================================
# Finite sample colonization ZONE suitable long
for (zo in 1:4){
	for (t in 1:nyear-1) {
		for (i in 1:nsite){
		OC_T1gam[zo,i,t] <- z[i,t] * tidev_mat_zone[i,zo] * SUITABLE[i,1]
		OC_T2gam[zo,i,t] <- z[i,t+1] * tidev_mat_zone[i,zo] * SUITABLE[i,1]
		UC_T1gam[zo,i,t] <- (z[i,t]==0) * tidev_mat_zone[i,zo] * SUITABLE[i,1]
		}
	gam.suit_z[zo,t] <- sum(OC_T1gam[zo,1:nsite,t] < OC_T2gam[zo,1:nsite,t]) / (sum(UC_T1gam[zo,1:nsite,t]) + 0.00000001)
	}
}
# ==========================================================

# end of model formulation
}
",fill=TRUE)
sink()
