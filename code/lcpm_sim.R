## running the LCPM on varying post-event durations, fix intervals, and 
library(nimble)
library(dplyr)
library(parallel)

## bring in the data 
lcpm.sim <- readRDS(
    "indata/lcpmSimData.RData"
)
## nested list
## first nest is HMD (Duratin), HMF (fix interval), HBC (magnitude of change)
## next nest is the different scenarios for each
## third nest is for each simulated datasets

## mmcpm function ----
lcpm_sim_func <- function(seed,data){
  require(nimble)
  require(dplyr)
  
    ## gcp model file
  modfile.lcpm <- nimbleCode({
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
  
  data.lcpm <- list(
    y = data[1:nrow(data),1:2],
    I = diag(2),
    Sigma = diag(5000,2,2)
  )
  consts.lcpm <- list(
    T = nrow(data)
  )
  inits.lcpm <- list(
    mu = cbind(
        rnorm(2,data[1,1],0.01),
        rnorm(2,data[1,2],0.01)
    ),
    Delta = diag(5,2,2),
    rho.AR = diag(0.5,1),
    tau = matrix(1,nrow=2,ncol=1),
    ind = 1,
    gamma = 0.5
    )
  
  params.lcpm <- c("mu","Delta",
                  "rho.AR","tau",
                  "cp")
  
  ## number of simulations
  n.sim <- 1
  ## number of data sets to simulation
  n.data <- 50
  ## number of iterations
  n.iter <- n.iter
  ## number of burnin
  n.burn <- n.burn
  
  ## run in 
  results <- nimbleMCMC(
    code = modfile.lcpm,
    constants = consts.lcpm,
    data = data.lcpm,
    inits = inits.lcpm,
    monitors = params.lcpm,
    niter = n.iter,
    nburnin = n.burn,
    setSeed = seed
  )
  
    ## find sim info and detail to find true part
  siminfo <- unique(data$siminfo)
  if(siminfo == ("HMF")){
    siminfoDetail <- unique(data$siminfoDetail)
    if(siminfoDetail == "15min"){
        true.part <- 192    
    }
    if(siminfoDetail == "30min"){
        true.part <- 97
    }
    if(siminfoDetail == "60min"){
        true.part <- 49
    }
    }else{
        true.part <- 192
    }

    ## return data frame of results 
  data.frame(results[,]) %>% 
    dplyr::mutate(
        cpdiff = cp.1. - true.part,
        datasim = rep(unique(data$datasim)),
        simInfo = rep(unique(data$siminfo)),
        simInfoDetail = rep(unique(data$siminfoDetail)),
        hmd = rep(unique(data$howMuchData)),
        hmf = rep(unique(data$howManyFixes))) 
}

## apply function in loop for all scenarios ----
for(q in 1:3){ ## for each question (duration, fix rate, magnitude)
    tmp.q <- lcpm.sim[[q]]
    scenario <- length(tmp.q)
    for(s in 1:scenario){ ## for each scenario in each question
        dat <- tmp.q[[s]]
        dat.split <- split(
            dat,
            dat$datasim
        )
        cl <- makeCluster(6)

        out.l <- parLapply(
            cl = cl,
            X = dat.split, 
            fun = lcpm_sim_func,
            seed = 10120,
            n.iter = 1000, ## change to higher iterations for inference
            n.burn = 500
        )

        out.df <- do.call(
            rbind,
            lapply(
                out.l,
                as.data.frame
            )
        ) ## create master data.frame for all data lengths and simulations (i x k)

        date_print <- format(Sys.Date(),"%Y%m%d")
        siminfo <- unique(dat.split[[1]]$siminfo)
        siminfoDetail <- unique(dat.split[[1]]$siminfoDetail)
        write.csv(
            out.df,
            paste("outdata/lcpm_",siminfo,"_",siminfoDetail,"_",date_print,".csv",sep=""),
            row.names = FALSE
        ) ## save master data.frame

        stopCluster(cl)
    }
}