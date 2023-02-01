## load packages ----
library(dplyr)
library(nimble)
library(parallel)

## bring in data ----
elk <- read.csv(
    "indata/elk_share.csv"
    )  %>% 
    dplyr::mutate(
        ID = as.factor(ID),
        ts = as.POSIXct(ts)
    ) 
elk.split <- split(elk, elk$ID) ## split data frame into list
elk.parturition <- read.csv(
    "indata/elkvit_share.csv"
    )  %>% 
    dplyr::mutate(
        ID = as.factor(ID),
        TruePart = as.POSIXct(TruePart)
    ) 
deer <- read.csv(
    "indata/deer_share.csv"
    )  %>% 
    dplyr::mutate(
        ID = as.factor(ID),
        ts = as.POSIXct(ts)
    )
deer.split <- split(deer, deer$ID)
deer.parturition <- read.csv(
    "indata/deervit_share.csv"
    )  %>% 
    dplyr::mutate(
        ID = as.factor(ID),
        TruePart = as.POSIXct(TruePart)
    ) 

## function for lcpm ----
lcmp.nimble <- function(data,vit.data){
  require(nimble)
  require(dplyr)
  niter <-  200000
  nburnin <- 100000
  
  ## data for model
  locations <- matrix(c(data$easting.centered,data$northing.centered),ncol = 2)
  
  ## model
  model.lcpm <- nimbleCode({
    # Priors for mean (attractant)
    chi[1:2]<-c(0,0)
    mu[1,1:2]~dmnorm(chi[1:2],cov=Sigma[1:2,1:2])
    mu[2,1:2]~dmnorm(chi[1:2],cov=Sigma[1:2,1:2])
    
    ## priors for Delta
    Delta[1,1]~dunif(0,100000)
    Delta[2,2]~dunif(0,100000)
    Delta[2,1]<-0
    Delta[1,2]<-Delta[2,1]
    
    rho.AR[1,1]~dunif(0.1,0.8) ## single rho.AR
    M[1,2]<-0
    M[2,1]<-0
    M[1,1]<-rho.AR[1,1]
    M[2,2]<-rho.AR[1,1]
    I.M[1:2,1:2]<-I[1:2,1:2]-M[1:2,1:2]
    ## z indexing
    prob[1:T]<-rep(1/T,T)
    tau[1,1]~dcat(prob[1:T])
    tau[2,1] <- 0
    ## spike and slab prior
    ind[1] ~ dbern(gamma[1])
    gamma[1] ~ dbeta(0.5,0.5)
    ## change-point based on spike and slab
    cp[1] <- tau[ind[1] + 1,1]
    z[1,1]<-1
    for(j in 2:T){
      z[j,1] <- (j >= cp[1]) + 1
    }
    ##likelihood
    y[1,1:2]~dmnorm(mu[z[1,1],1:2],cov=Delta[1:2,1:2])
    psi[1,1:2] <- mu[z[1,1],1:2]
    for(i in 2:T){
      psi[i,1:2]<-(y[i-1,1:2]%*%M[1:2,1:2]+mu[z[i,1],1:2]%*%I.M[1:2,1:2])[1,1:2]
      y[i,1:2]~dmnorm(psi[i,1:2],cov=Delta[1:2,1:2])
    }
  })
  
  data <- list(
    ## this is the data
    y = locations,
    I = diag(2),
    Sigma = diag(5000,2,2)
  )
  const <- list(
    T = nrow(data)
  )
  inits <- list(mu = cbind(rnorm(2,locations[1,1],0.01),
                            rnorm(2,locations[1,2],0.01)),
                Delta = diag(5,2,2),
                rho.AR = diag(0.5,1),
                tau = matrix(1,nrow=2,ncol=1),ind = 1,
                gamma = 0.5)
  params <- c("mu","Delta",
              "rho.AR","tau",
              "cp")
  
  results <- nimbleMCMC(
    code = model.lcpm,
    constants = const,
    data = data,
    inits = inits,
    monitors = params,
    niter = niter,
    nburnin = nburnin
  )
  ## find the true part
  vit.date <- vit.data %>% 
    dplyr::filter(ID == unique(data$ID))
  
  ## find the true part ts
  tmp.truepart<-which(data$ts==vit.date$TruePart)
  ## create data frame
  df.lcpm <- data.frame(results[,],ID=id) %>% 
    dplyr::mutate(taudiff = cp.1. - tmp.truepart,
                  model = "lcpm")
  
  return(df.lcpm)
}

## run for elk ----
## create cluster to run in parallel
this_cluster <- makeCluster(6)

parallel::clusterExport(this_cluster,"lcpm.nimble")
system.time({
    lcpm.out.elk <- parLapply(
        cl = this_cluster,
        X = elk.split,
        fun = lcpm.nimble,
        vit.data = elk.parturition
    )
})
stopCluster(this_cluster)
elk.lcpm <- do.call(rbind,lapply(lcpm.out.elk,as.data.frame))
date_print <- format(Sys.Date(),"%Y%m%d")
write.csv(elk.lcpm,paste("outdata/lcpm_elk_",date_print,".csv",sep=""),row.names=F)

## run for deer ----

## create cluster to run in parallel
this_cluster <- makeCluster(6)

parallel::clusterExport(this_cluster,"lcpm.nimble")
system.time({
    lcpm.out.deer <- parLapply(
        cl = this_cluster,
        X = deer.split,
        fun = lcpm.nimble,
        vit.data = deer.parturition
    )
})
stopCluster(this_cluster)
deer.lcpm <- do.call(rbind,lapply(lcpm.out.deer,as.data.frame))
date_print <- format(Sys.Date(),"%Y%m%d")
write.csv(deer.lcpm,paste("outdata/lcpm_deer_",date_print,".csv",sep=""),row.names=F)

