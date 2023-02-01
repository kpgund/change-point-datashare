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

## function mmcpm - elk ----
mmcpm.nimble <- function(data,vit.data){
  require(nimble)
  require(dplyr)
  
  niter <-  200000
  nburnin <-  100000
  
  data <- data[-1,] ## taking out the first and the last rows with NAs
  data <- data[-nrow(data),]  
  
  T <- nrow(data) ## T is length of data
  
  ## create data
  ones <- numeric(T)+1 ## for the model
  steps <- data$sl/1000 ## make sl km 
  steps <- ifelse(steps==0,0.01,steps) ## steps cannot be 0, replace with very small number
  turns <- data$ta
  
  ## model
  modfile.mmcpm <- nimbleCode({
    # priors on movement parameters
    # step lengths
    a[1] ~ dgamma(0.001,0.001)
    a[2] ~ dgamma(0.001,0.001)
    b[1] ~ dunif(0, 8)
    b[2] ~ dunif(0, 8)
    
    # turning angles 
    mu[1] ~ dunif(-1, 6)
    mu[2] ~ dunif(-1, 6)
    rho[1] ~ dunif(0.0,1.0)
    rho[2] ~ dunif(0.0,1.0)
    
    pi <- 3.14159265359
    prob[1:T]<-rep(1/T,T)
    tau[1,1] ~ dcat(prob[1:T])
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
    
    #Likelihood  
    for (i in 1:T){
      # likelihood for steplengths
      steps[i] ~ dweib(shape=a[z[i,1]],scale=b[z[i,1]])
      # likelihood for turns
      ones[i] ~ dbern( wc[i] )
      wc[i] <- (1/(2*pi)*(1-pow(rho[z[i,1]],2))/(1+pow(rho[z[i,1]],2)-2*rho[z[i,1]]*cos(turns[i]-mu[z[i,1]])))/100
      turns[i] ~ dunif(-3.14159265359, 3.14159265359)
    }
  })
  
  inits <- list(a = c(0.5,0.5), b = c(0.5,0.5), 
                mu = c(0.5,0.5), rho = c(0.5,0.5),
                tau = matrix(1,nrow=2,ncol=1),ind = 1,
                gamma = 0.5)
  dat <- list(
    ## this is the data
    steps = steps,
    turns = turns,
    ones = ones
  )
  consts <- list(
    ## constants
    T = T
  )
  params <- c("a","b",
              "mu","rho",
              "tau","cp")
  
  results <- nimbleMCMC(
    code = modfile.mmcpm,
    constants = consts,
    data = dat,
    inits = inits,
    monitors = params,
    niter = niter,
    nburnin = nburnin
  )
  
  ## find the true part
  vit.date <- vit.data %>% 
    dplyr::filter(ID == unique(data$ID))
  
  ## find the true part ts
  tmp.truepart <- which(data$ts==vit.date$TruePart)
  ## create data frame
  df.mmcpm <- data.frame(results[,],ID=unique(data$ID)) %>% 
    dplyr::mutate(
        taudiff = cp.1. - tmp.truepart,
        model = "mmcpm"
    )
  
  return(df.mmcpm)
}
## run for elk ----

## create cluster to run in parallel
this_cluster <- makeCluster(6)

parallel::clusterExport(this_cluster,"mmcpm.nimble")
system.time({
    mmcpm.out.elk <- parLapply(
        cl = this_cluster,
        X = elk.split,
        fun = mmcpm.nimble,
        vit.data = elk.parturition
    )
})
stopCluster(this_cluster)
elk.mmcpm <- do.call(rbind,lapply(mmcpm.out.elk,as.data.frame))
date_print <- format(Sys.Date(),"%Y%m%d")
write.csv(elk.mmcpm,paste("outdata/mmcpm_elk_",date_print,".csv",sep=""),row.names=F)

## function mmcpm - deer ----
mmcpm.exp.nimble <- function(data,vit.data){
  require(nimble)
  require(dplyr)
  
  niter <-  200000
  nburnin <-  100000
  
  data <- data[-1,] ## taking out the first and the last rows with NAs
  data <- data[-nrow(data),]  
  
  T<-nrow(data) ## T is length of data
  
  ## create data
  ones<-numeric(T)+1 ## for the model
  steps<-data$sl/1000 ## make sl km 
  steps<- ifelse(steps==0,0.01,steps) ## steps cannot be 0, replace with very small number
  turns<-data$ta
  
  ## exponential model file
  modfile.exp.mmcpm <- nimbleCode( {
    # priors on movement parameters
    # step lengths
    lambda[1] ~ dgamma(0.001,0.001)
    lambda[2] ~ dgamma(0.001,0.001)
    
    # turning angles 
    mu[1] ~ dunif(-1, 6)
    mu[2] ~ dunif(-1, 6)
    rho[1] ~ dunif(0.0,1.0)
    rho[2] ~ dunif(0.0,1.0)
    
    pi <- 3.14159265359
    prob[1:T]<-rep(1/T,T)
    tau[1,1] ~ dcat(prob[1:T])
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
    
    #Likelihood  
    for (i in 1:T){
      # likelihood for steplengths
      steps[i] ~ dexp(lambda[z[i,1]]) # a is the shape parameter
      # likelihood for turns
      ones[i] ~ dbern( wc[i] )
      wc[i] <- (1/(2*pi)*(1-pow(rho[z[i,1]],2))/(1+pow(rho[z[i,1]],2)-2*rho[z[i,1]]*cos(turns[i]-mu[z[i,1]])))/100
      turns[i] ~ dunif(-3.14159265359, 3.14159265359)
    }
  })
  
  ## inital values
  inits <- list(lambda = c(0.5,0.5),
                mu = c(0.5,0.5), rho = c(0.5,0.5),
                tau = matrix(1,nrow=2,ncol=1),ind = 1,
                gamma = 0.5)
  ## data
  dat <- list(
    ## this is the data
    steps = steps,
    turns = turns,
    ones = ones
  )
  consts <- list(
    ## constants
    T = T
  )
  params <- c("lambda",
              "mu","rho",
              "tau","cp")
  
  results <- nimbleMCMC(
    code = modfile.exp.mmcpm,
    constants = consts,
    data = dat,
    inits = inits,
    monitors = params,
    niter = niter,
    nburnin = nburnin
  )
  
    ## find the true part
  vit.date <- vit.data %>% 
    dplyr::filter(ID == unique(data$ID))
  
  ## find the true part ts
  tmp.truepart <- which(data$ts==vit.date$TruePart)
  ## create data frame
  df.mmcpm <- data.frame(results[,],ID=unique(data$ID)) %>% 
    dplyr::mutate(
        taudiff = cp.1. - tmp.truepart,
        model = "mmcpm"
    )
  
  return(df.mmcpm)
}
## run for deer ----

## create cluster to run in parallel
this_cluster <- makeCluster(6)

parallel::clusterExport(this_cluster,"mmcpm.nimble")
system.time({
    mmcpm.out.deer <- parLapply(
        cl = this_cluster,
        X = deer.split,
        fun = mmcpm.exp.nimble,
        vit.data = deer.parturition
    )
})
stopCluster(this_cluster)
deer.mmcpm <- do.call(rbind,lapply(mmcpm.out.deer,as.data.frame))
date_print <- format(Sys.Date(),"%Y%m%d")
write.csv(elk.mmcpm,paste("outdata/mmcpm_deer_",date_print,".csv",sep=""),row.names=F)
