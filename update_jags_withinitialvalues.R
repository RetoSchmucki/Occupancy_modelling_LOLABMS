# setwd("/home/reto/Documents/LOLA_BMS/Occupancy Model project/Occupancy_modellingNEWJUNE/Anthocharis_cardamines")
library(R2jags)
load(list.files()[grep(".RData",list.files())])
load("results/jagsoutput.Rdata")

modelnbr <- 1
ni <- 10000
nb <- 5000
nt <- 10
nch <- 3

jags_model <- paste0("OCCmodel_",modelnbr,"_jagscript_tidev_zone.txt")

analysis.par <- data.frame(parmeters=c("species","country_of_flight","country","country_with_obs","year_toget","min_nbr_year","nbr_per_gridcell","ni","nb","nt","jags_model"),values=c(species,country_of_flight,country,paste(country_with_obs,collapse=","),year_toget,min_nbr_year,nbr_siteperGRID,ni,nb,nt,jags_model))

y <- c(1:nyear)

data <- list("M", "y", "nsite", "nvisit", "nyear", "sumX", "sumX2", "maxyear","NM","TL","psi_min","psi_max","bnm_min","bnm_max","btl_min","btl_max","bphi_min","bphi_max","bgam_min","bgam_max","tidev_mat_zone","SUITABLE")


# (2) Starting values
warn <- -1
zst <- apply(M,c(1,3),max,na.rm=TRUE)  	# inits per site and year; c(1,3) means: inits for dimension 1-3 = site, visit and year
warn <- 0
zst[zst==(-Inf)] <- 1 			# inits for -infinite: 1
zst[zst==(Inf)]  <- 1			# inits for infinite : 1
pst <- rep(0.5, nyear) 

inits2 <- function() {list (z=zst, 
alpha.p=apply(out$BUGSoutput$sims.list$alpha.p,2,function(x) rnorm(1,mean(x),sd(x))),
muphi1=apply(out$BUGSoutput$sims.list$muphi1,2,function(x) rnorm(1,mean(x),sd(x))),
muphi2=apply(out$BUGSoutput$sims.list$muphi2,2,function(x) rnorm(1,mean(x),sd(x))),
muphi3=apply(out$BUGSoutput$sims.list$muphi3,2,function(x) rnorm(1,mean(x),sd(x))),
muphi4=apply(out$BUGSoutput$sims.list$muphi4,2,function(x) rnorm(1,mean(x),sd(x))),
mugam1=apply(out$BUGSoutput$sims.list$mugam1,2,function(x) rnorm(1,mean(x),sd(x))),
mugam2=apply(out$BUGSoutput$sims.list$mugam2,2,function(x) rnorm(1,mean(x),sd(x))),
mugam3=apply(out$BUGSoutput$sims.list$mugam3,2,function(x) rnorm(1,mean(x),sd(x))),
mugam4=apply(out$BUGSoutput$sims.list$mugam4,2,function(x) rnorm(1,mean(x),sd(x))),
betaphi1=apply(out$BUGSoutput$sims.list$betaphi1,2,function(x) rnorm(1,mean(x),sd(x))),
betaphi2=apply(out$BUGSoutput$sims.list$betaphi2,2,function(x) rnorm(1,mean(x),sd(x))),
betaphi3=apply(out$BUGSoutput$sims.list$betaphi3,2,function(x) rnorm(1,mean(x),sd(x))),
betaphi4=apply(out$BUGSoutput$sims.list$betaphi4,2,function(x) rnorm(1,mean(x),sd(x))),
betagam1=apply(out$BUGSoutput$sims.list$betagam1,2,function(x) rnorm(1,mean(x),sd(x))),
betagam2=apply(out$BUGSoutput$sims.list$betagam2,2,function(x) rnorm(1,mean(x),sd(x))),
betagam3=apply(out$BUGSoutput$sims.list$betagam3,2,function(x) rnorm(1,mean(x),sd(x))),
betagam4=apply(out$BUGSoutput$sims.list$betagam4,2,function(x) rnorm(1,mean(x),sd(x))),
sigmaphi1=apply(out$BUGSoutput$sims.list$sigmaphi1,2,function(x) runif(1,min(x),min(max(x),5))),
sigmaphi2=apply(out$BUGSoutput$sims.list$sigmaphi2,2,function(x) runif(1,min(x),min(max(x),5))),
sigmaphi3=apply(out$BUGSoutput$sims.list$sigmaphi3,2,function(x) runif(1,min(x),min(max(x),5))),
sigmaphi4=apply(out$BUGSoutput$sims.list$sigmaphi4,2,function(x) runif(1,min(x),min(max(x),5))),
sigmagam1=apply(out$BUGSoutput$sims.list$sigmagam1,2,function(x) runif(1,min(x),min(max(x),5))),
sigmagam2=apply(out$BUGSoutput$sims.list$sigmagam2,2,function(x) runif(1,min(x),min(max(x),5))),
sigmagam3=apply(out$BUGSoutput$sims.list$sigmagam3,2,function(x) runif(1,min(x),min(max(x),5))),
sigmagam4=apply(out$BUGSoutput$sims.list$sigmagam4,2,function(x) runif(1,min(x),min(max(x),5))),
bnm.p=apply(out$BUGSoutput$sims.list$bnm.p,2,function(x) rnorm(1,mean(x),sd(x))),
btl.p=apply(out$BUGSoutput$sims.list$btl.p,2,function(x) rnorm(1,mean(x),sd(x))))}


parameters <- c("z", "psi.fs", "psi.fs_z", "psi.suit_z", "phi.suit_z", "gam.suit_z",
                "bnm.p", "btl.p", "alpha.p",
                "muphi1", "muphi2", "muphi3", "muphi4",
                "sigmaphi1_2", "sigmaphi2_2", "sigmaphi3_2", "sigmaphi4_2",
                "betaphi1", "betaphi2", "betaphi3", "betaphi4",
                "beta.phi1", "beta.phi2", "beta.phi3", "beta.phi4",
                "mugam1", "mugam2", "mugam3", "mugam4",
                "sigmagam1_2", "sigmagam2_2", "sigmagam3_2", "sigmagam4_2",
                "betagam1", "betagam2", "betagam3", "betagam4",
                "beta.gam1", "beta.gam2", "beta.gam3", "beta.gam4")

library(R2jags)

jags_fit <- function (data, inits, parameters.to.save, model.file = "model.bug", 
                      n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
                      n.thin = max(1, floor((n.iter - n.burnin)/1000)), DIC = TRUE, 
                      working.directory = NULL, jags.seed = 123, refresh = n.iter/50, 
                      progress.bar = "text", digits = 4, RNGname = c("Wichmann-Hill", 
                                                                     "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")) 
{
  if (!is.null(working.directory)) {
    working.directory <- path.expand(working.directory)
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
  }
  else {
    savedWD <- getwd()
    working.directory <- savedWD
  }
  if (is.character(data) && length(data) == 1 && regexpr("\\.txt$", 
                                                         data) > 0) {
    if (all(basename(data) == data)) {
      fn2 <- file.path(working.directory, data)
      if (normalizePath(fn2) != normalizePath(data)) {
        try(file.copy(fn2, data, overwrite = TRUE))
      }
    }
    if (!file.exists(data)) {
      stop("File", data, "does not exist")
    }
    if (file.info(data)["size"] == 0) {
      stop("Empty data file ", data)
    }
    e <- new.env()
    eval(parse(data), e)
    data <- as.list(e)
  }
  else if (is.character(data) || (is.list(data) && all(sapply(data, 
                                                              is.character)))) {
    dlist <- lapply(as.list(data), get, envir = parent.frame(1))
    names(dlist) <- unlist(data)
    data <- dlist
  }
  else if (!is.list(data)) {
    stop("data must be a character vector of object names, a list of object names, or a list of objects")
  }
  if (is.function(model.file)) {
    temp <- tempfile("model")
    temp <- if (is.R() || .Platform$OS.type != "windows") {
      paste(temp, "txt", sep = ".")
    }
    else {
      gsub("\\.tmp$", ".txt", temp)
    }
    write.model(model.file, con = temp, digits = digits)
    model.file <- gsub("\\\\", "/", temp)
    if (!is.R()) 
      on.exit(file.remove(model.file), add = TRUE)
  }
  if (DIC) {
    parameters.to.save <- c(parameters.to.save, "deviance")
    load.module("dic", quiet = TRUE)
  }
  if (n.burnin > 0) {
    n.adapt <- n.burnin
  }
  else {
    n.adapt <- 100
  }
  if (!missing(inits) && !is.function(inits) && !is.null(inits) && 
      (length(inits) != n.chains)) {
    stop("Number of initialized chains (length(inits)) != n.chains")
  }
  RNGname <- match.arg(RNGname)
  if (RNGname %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                     "Super-Duper", "Mersenne-Twister")) {
    RNGname <- paste("base::", RNGname, sep = "")
  }
  else {
    stop("The name of the RNG is not correctly provided!")
  }
  load.module("glm")
  init.values <- vector("list", n.chains)
  if (missing(inits)) {
    for (i in 1:n.chains) {
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- abs(.Random.seed[i + 
                                                       1])
    }
  }
  else if (is.null(inits)) {
    for (i in 1:n.chains) {
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- abs(.Random.seed[i + 
                                                       1])
    }
  }
  else if (is.function(inits)) {
    if (any(names(formals(inits)) == "chain")) {
      for (i in 1:n.chains) {
        init.values[[i]] <- inits(chain = i)
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- abs(.Random.seed[i + 
                                                         1])
      }
    }
    else {
      for (i in 1:n.chains) {
        init.values[[i]] <- inits()
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- abs(.Random.seed[i + 
                                                         1])
      }
    }
  }
  else {
    if (!is.list(inits)) {
      stop("Invalid inits")
    }
    if (length(inits) != n.chains) {
      stop("Number of initialized chains (length(inits)) != n.chains")
    }
    for (i in 1:n.chains) {
      init.values[[i]] <- inits[[i]]
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- abs(.Random.seed[i + 
                                                       1])
    }
  }
  m <- jags.model(model.file, data = data, inits = init.values, 
                  n.chains = n.chains, n.adapt = 0)
  adapt(m, n.iter = n.adapt, by = refresh, progress.bar = progress.bar, 
        end.adaptation = TRUE)
  samples <- coda.samples(model = m, variable.names = parameters.to.save, 
                          n.iter = (n.iter - n.burnin), thin = n.thin, by = refresh, 
                          progress.bar = progress.bar)
  fit <- mcmc2bugs(samples, model.file = model.file, program = "jags", 
                   DIC = DIC, DICOutput = NULL, n.iter = n.iter, n.burnin = n.burnin, 
                   n.thin = n.thin)[c("sims.array","sims.list","summary","pD","DIC")]
  out <- list(model = m, BUGSoutput = fit, parameters.to.save = parameters.to.save, 
              model.file = model.file, n.iter = n.iter, DIC = DIC)
  class(out) <- "rjags"
  return(out)
}
environment(jags_fit) <- environment(jags)

set.seed(1234)

out2 <- jags_fit(data, inits2, parameters, jags_model, n.chains=nch, n.iter = ni, n.thin=nt, n.burnin=nb, DIC=TRUE, working.directory=NULL)
