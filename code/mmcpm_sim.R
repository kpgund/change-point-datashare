## running the MMCPM on varying post-event durations, fix intervals, and 
library(nimble)
library(dplyr)
library(parallel)

## bring in the data 
mmcpm.sim <- readRDS(
    "mmcpmSimData.RData"
)
## nested list
## first nest is HMD (Duratin), HMF (fix interval), HBC (magnitude of change)
## next nest is the different scenarios for each
## third nest is for each simulated datasets

## mmcpm function ----
mmcpm_sim_func <- function(seed,data,n.iter,n.burn){
  require(nimble)
  require(dplyr)
  
  ## model
    mmcpm.modfile <-nimbleCode({
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
  
  mmcpm.inits <- list(a = c(0.5,0.5), b = c(0.5,0.5), 
                    mu = c(0.5,0.5), rho = c(0.5,0.5),
                    tau = matrix(1,nrow=2,ncol=1),ind = 1,
                    gamma = 0.5)
  data.mmcpm <- list(
    ## this is the data
    steps = data$steps,
    turns = data$turns,
    ones = data$ones
  )
  consts.mmcpm <- list(
    ## constants
    T = nrow(data)
  )
  
  ## number of simulations
  n.sim <- 1
  ## number of data sets to simulation
  n.data <- 50
  ## number of iterations
  n.iter <- n.iter
  ## number of burnin
  n.burn <- n.burn
  
  params.mmcpm <- c("a","b",
                  "mu","rho",
                  "tau","cp")
  ## run in 
  results <- nimbleMCMC(
    code =  mmcpm.modfile,
    constants = consts.mmcpm,
    data = data.mmcpm,
    inits = mmcpm.inits,
    monitors = params.mmcpm,
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

  ## make df
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
    tmp.q <- mmcpm.sim[[q]]
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
            fun = mmcpm_sim_func,
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
            paste("outdata/mmcpm_",siminfo,"_",siminfoDetail,"_",date_print,".csv",sep=""),
            row.names = FALSE
        ) ## save master data.frame

        stopCluster(cl)
    }
}