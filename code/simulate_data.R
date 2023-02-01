######
##
## creating dataframe and writing CSV files for simulations
##
######

library(lubridate)
library(dplyr)
library(circular)
library(mvtnorm)

## read in MMCPM data ----
## read in the latest categorize_outputs
mmcpm.sim.params <- read.csv(
    "indata/summary_posteriormean_mmcpm.csv"
)

## create output for mmcpmSimData
mmcpmSimData <- list()

## HBC - MMCPM ---- 
## set seed
set.seed(10120) ## NEED TO RUN THIS FIRST VERY IMPORTANT 
## number of simulations
n.sim <- 1
## number of data sets to simulation
n.data <- 50

HBC.data.final <- list()
HBC.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## simulate the data
  T<-385 ##(gives 48 hours pre/post at 15 min fixes)
  ## create time stamp
  ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)

  # estimated tau in elk
  a<-list(matrix(c(mmcpm.sim.params$sl_shape_before,
                   mmcpm.sim.params$sl_shape_before),nrow=1), ## no diff sl
          matrix(c(mmcpm.sim.params$sl_shape_before,
                   mmcpm.sim.params$sl_shape_after),nrow=1), ## diff in sl
          matrix(c(mmcpm.sim.params$sl_shape_before,
                   mmcpm.sim.params$sl_shape_before),nrow=1), ## no diff sl
          matrix(c(mmcpm.sim.params$sl_shape_before,
                   mmcpm.sim.params$sl_shape_after),nrow=1)) ## diff in sl
  b<-matrix(c(mmcpm.sim.params$sl_scale_before,
              mmcpm.sim.params$sl_scale_before),nrow=1) ## b is not differing so 
  ## turning angle
  mu<-list(matrix(c(mmcpm.sim.params$ta_mu_before,
                    mmcpm.sim.params$ta_mu_before),nrow=1), ## no diff ta
           matrix(c(mmcpm.sim.params$ta_mu_before,
                    mmcpm.sim.params$ta_mu_before),nrow=1), ## no diff ta
           matrix(c(mmcpm.sim.params$ta_mu_before,
                    mmcpm.sim.params$ta_mu_after),nrow=1), ## diff ta
           matrix(c(mmcpm.sim.params$ta_mu_before,
                    mmcpm.sim.params$ta_mu_after),nrow=1)) ## diff ta
  rho<-list(matrix(c(mmcpm.sim.params$ta_rho_before,
                     mmcpm.sim.params$ta_rho_before),nrow=1), ## no diff ta
            matrix(c(mmcpm.sim.params$ta_rho_before
                     ,mmcpm.sim.params$ta_rho_before),nrow=1), ## no diff ta
            matrix(c(mmcpm.sim.params$ta_rho_before,
                     mmcpm.sim.params$ta_rho_after),nrow=1), ## diff ta
            matrix(c(mmcpm.sim.params$ta_rho_before,
                     mmcpm.sim.params$ta_rho_after),nrow=1)) ## diff ta
  
  
  opt <- c("no_diff","diff_sl","diff_ta","diff_both")
  tau<-192
  
  for(k in 1:length(opt)){
    ## turning angles - wrapped cauchy (package)
    turns<- rep(NA,T)
    steps<-rep(NA,T)
    z<-rep(0,T)
    ## create steps and turns based on parameters
    for(i in 1:T){
      z[i]<-ifelse(i<tau,0,1)
      z.idx<-z[i]+1
      steps[i] <- rweibull(1, shape = a[[k]][1, z.idx], scale = b[, z.idx])
      turns[i] <- rwrappedcauchy(1, mu = circular(mu[[k]][1,z.idx],zero=pi,units="radians"), 
                                 rho = rho[[k]][1,z.idx])
    }
    ## ones
    ones<-numeric(T)+1
    ## create the dataset
    sim.data.final <- data.frame(ones = ones,
                                 steps = steps,
                                 turns = turns,
                                 state = z,
                                 ts = ts,
                                 datasim = j,
                                 siminfo = "HBC",
                                 siminfoDetail = opt[k],
                                 hmd = "48hr",
                                 hmf = "15min",
                                 hbc = opt[k])
    
    HBC.data.final[[k]] <- sim.data.final
  }
  HBC.data[[j]] <- HBC.data.final
}

## put length of data into own data frames
hbc.both <- NULL;hbc.none <- NULL;hbc.sl <- NULL;hbc.ta <- NULL
for(i in 1:length(HBC.data)){
  tmp.l <- HBC.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hbc <- unique(tmp$hbc)
    if(hbc == "no_diff"){
      hbc.none <- rbind(hbc.none,
                        tmp)
    }
    if(hbc == "diff_sl"){
      hbc.sl <- rbind(hbc.sl,
                      tmp)
    }
    if(hbc == "diff_ta"){
      hbc.ta <- rbind(hbc.ta,
                      tmp)
    }
    if(hbc == "diff_both"){
      hbc.both <- rbind(hbc.both,
                        tmp)
    }
  }
}

hbc.list <- list(hbc.none,hbc.sl,
                 hbc.ta,hbc.both)
## put data into final container
mmcpmSimData[[1]] <- hbc.list

## HMF - MMCPM ----
## set seed
set.seed(10120)
HMF.data.final <- list()
HMF.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## simulate the data
  T<-385 ##(gives 48 hours pre/post at 15 min fixes)
  
  ## turning angles - wrapped cauchy (package)
  turns<- rep(NA,T)
  steps<-rep(NA,T)
  z<-rep(0,T)
  ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)
  ts[192] ## as of right now, parturition is my birth day and time
  
  # estimated tau in elk
  a<-matrix(c(mmcpm.sim.params$sl_shape_before,
              mmcpm.sim.params$sl_shape_after),nrow=1) ## need to change to mean
  b<-matrix(c(mmcpm.sim.params$sl_scale_before,
              mmcpm.sim.params$sl_scale_after),nrow=1) ## need to change to sd
  mu<-matrix(c(mmcpm.sim.params$ta_mu_before,
               mmcpm.sim.params$ta_mu_after),nrow=1)
  rho<-matrix(c(mmcpm.sim.params$ta_rho_before,
                mmcpm.sim.params$ta_rho_after),nrow=1) # must be an interval between 0 and 1
  tau<-192
  
  ## create steps and turns based on parameters
  for(i in 1:T){
    z[i]<-ifelse(i<tau,0,1)
    z.idx<-z[i]+1
    steps[i] <- rweibull(1, shape = a[1, z.idx], scale = b[1, z.idx])
    turns[i] <- turns[i] <- rwrappedcauchy(1, mu = circular(mu[1,z.idx],zero=pi,units="radians"), 
                                           rho = rho[1,z.idx])
  }
## truncate and save as needed in a loop
  how.much.data <- 24 ## in hours after change
  how.much.data.quotes <- c("24hr")  ## for hours after change
  how.many.fixes <- c(1,2,4)
  how.many.fixes.quotes <- c("15min","30min","60min")
  ## ones
  ones<-numeric(T)+1
  ## create the dataset
  sim.data.final <- data.frame(ones = ones,
                               steps = steps,
                               turns = turns,
                               state = z,
                               ts = ts,
                               datasim = j,
                               siminfo = "hmf")
  
  ## take the full data set and truncate "fixes"
  for(k in 1:length(how.many.fixes)){
    if(k == 1){
      sim.data.fix <- sim.data.final %>% 
        dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data))) %>% 
        dplyr::mutate(
            siminfoDetail = how.many.fixes.quotes[k],
            hmd = how.much.data.quotes,
            hmf = how.many.fixes.quotes[k]
        )
      # write.csv(sim.data.fix,paste("simulations/sim_dataset_HMD_HMF_BCP",how.much.data.quotes[k],
      #                            how.many.fixes.quotes[k],".csv",sep="_"),row.names=FALSE)
      
    }else{
      sim.data.fix <- sim.data.final %>% 
        dplyr::filter(row_number() %% how.many.fixes[k] == 1) %>%
        dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data))) %>% 
        dplyr::mutate(
            siminfoDetail = how.many.fixes.quotes[k],
            hmd = how.much.data.quotes,
            hmf = how.many.fixes.quotes[k]
        )
      # write.csv(sim.data.fix,paste("simulations/sim_dataset_HMD_HMF_BCP",how.much.data.quotes[k],
      #                              how.many.fixes.quotes[k],".csv",sep="_"),row.names=FALSE)
      # assign(paste("sim_dataset_",how.much.data.quotes[k],"_15min",sep=""),sim.data.15min)
    }
    HMF.data.final[[k]] <- sim.data.fix
  }
  HMF.data[[j]] <- HMF.data.final
}

## put length of data into own data frames
hmf.15 <- NULL;hmf.30 <- NULL;hmf.60 <- NULL
for(i in 1:length(HMF.data)){
  tmp.l <- HMF.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hmf <- unique(tmp$hmf)
    if(hmf == "15min"){
      hmf.15 <- rbind(hmf.15,
                      tmp)
    }
    if(hmf == "30min"){
      hmf.30 <- rbind(hmf.30,
                      tmp)
    }
    if(hmf == "60min"){
      hmf.60 <- rbind(hmf.60,
                      tmp)
    }
  }
}

hmf.list <- list(hmf.15,hmf.30,
                 hmf.60)
## put data into final container
mmcpmSimData[[2]] <- hmf.list

## HMD - MMCPM ----
## set seed
set.seed(10120)
HMD.data.final <- list()
HMD.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## simulate the data
  T<-385 ##(gives 48 hours pre/post at 15 min fixes)
  
  ## turning angles - wrapped cauchy (package)
  turns<- rep(NA,T)
  steps<-rep(NA,T)
  z<-rep(0,T)
  ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)
  ts[192] ## as of right now, parturition is my birth day and time
  
  # estimated tau in elk
  a<-matrix(c(mmcpm.sim.params$sl_shape_before,
              mmcpm.sim.params$sl_shape_after),nrow=1) ## need to change to mean
  b<-matrix(c(mmcpm.sim.params$sl_scale_before,
              mmcpm.sim.params$sl_scale_after),nrow=1) ## need to change to sd
  mu<-matrix(c(mmcpm.sim.params$ta_mu_before,
               mmcpm.sim.params$ta_mu_after),nrow=1)
  rho<-matrix(c(mmcpm.sim.params$ta_rho_before,
                mmcpm.sim.params$ta_rho_after),nrow=1) # must be an interval between 0 and 1
  tau<-192
  
  ## create steps and turns based on parameters
  for(i in 1:T){
    z[i]<-ifelse(i<tau,0,1)
    z.idx<-z[i]+1
    steps[i] <- rweibull(1, shape = a[1, z.idx], scale = b[1, z.idx])
    turns[i] <- rwrappedcauchy(1, mu = circular(mu[1,z.idx],zero=pi,units="radians"), rho = rho[1,z.idx])
  }
  
  ## create data frame of data
  
  ## ones
  ones<-numeric(T)+1
  
  sim.data.final <- data.frame(ones = ones,
                               steps = steps,
                               turns = turns,
                               state = z,
                               ts = ts,
                               datasim = j,
                               siminfo = "hmd")
  
  ## truncate and save as needed in a loop
  how.much.data <- c(48,24,12,6,3) ## in hours
  how.much.data.quotes <- c("48hr","24hr","12hr","6hr","3hr") 
  
  ## take the full data set and truncate
  for(k in 1:length(how.much.data)){
    sim.data.15min <- sim.data.final %>% 
      dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data[k])+minutes(15))) %>% 
      dplyr::mutate(
        siminfoDetail = how.much.data.quotes[k],
        hmd = how.much.data.quotes[k],
        hmf = "15min")
    # write.csv(sim.data.15min,paste("simulations/sim_dataset_HMD_",how.much.data.quotes[k],".csv",sep=""),row.names=FALSE)
    
    HMD.data.final[[k]] <- sim.data.15min
  }
  HMD.data[[j]] <- HMD.data.final
}
## put length of data into own data frames
hmd.48 <- NULL;hmd.24 <- NULL;hmd.12 <- NULL;hmd.06 <- NULL;hmd.03 <- NULL
for(i in 1:length(HMD.data)){
  tmp.l <- HMD.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hmd <- unique(tmp$hmd)
    if(hmd == "48hr"){
      hmd.48 <- rbind(hmd.48,
                      tmp)
    }
    if(hmd == "24hr"){
      hmd.24 <- rbind(hmd.24,
                      tmp)
    }
    if(hmd == "12hr"){
      hmd.12 <- rbind(hmd.12,
                      tmp)
    }
    if(hmd == "6hr"){
      hmd.06 <- rbind(hmd.06,
                      tmp)
    }
    if(hmd == "3hr"){
      hmd.03 <- rbind(hmd.03,
                      tmp)
    }
  }
}

hmd.list <- list(hmd.48,hmd.24,
                 hmd.12,hmd.06,hmd.03)
## put data into final container
mmcpmSimData[[3]] <- hmd.list

## save mmcpm simulated data as .RDS ----
saveRDS(
    mmcpmSimData,
    "indata/mmcpmSimData.RData"
)


## read in data LCPM ----
## read in the latest categorize_outputs
max.cent.diff <- read.csv(
    "indata/summary_posteriormean_lcpm.csv"
)
max.cent.N.var <- 19674.51 ## for variance 

## create output for mmcpmSimData
lcpmSimData <- list()

## HBC- LCPM ----
## set seed
set.seed(10120) ## NEED TO RUN THIS FIRST VERY IMPORTANT 
## number of simulations
n.sim <- 1
## number of data sets to simulation
n.data <- 50

## HBC for the loop
diff.for.HBC <- c(1,0.75,0.5,0.25,0)
diff.for.HBC.quotes <- c("100","75","50","25","0")


mu <- matrix(c(max.cent.diff$Ecent_mean_before, ## mean
               max.cent.diff$Ecent_mean_after,
               max.cent.diff$Ncent_mean_before,
               max.cent.diff$Ncent_mean_after),nrow=2,byrow=F)
## figuring out triangle stuff
new.points <- approx(mu,n=5)
new.df<- do.call(cbind,lapply(new.points,as.data.frame))
colnames(new.df) <- c("x","y")
new.df$dist <- as.factor(c("100","75","50","25","0"))


HBC.data.final <- list()
HBC.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## create 4 subsets of the data
  for(k in 1:(length(diff.for.HBC)-1)){
    T<-385 ##(gives 48 hours pre/post at 15 min fixes)
    
    mu.1 <- matrix(c(new.df[5,1],new.df[5,2]),ncol=2)
    mu.2 <- matrix(c(new.df[k,1],new.df[k,2]),ncol=2)
    mu <- rbind(mu.1,mu.2)
    # mean.sig <- c(mean(gcp.sim.params$Ecent_var),mean(gcp.sim.params$Ncent_var))
    Delta<-matrix(c(max.cent.N.var,0,0,max.cent.N.var),nrow=2,byrow=2) ## variance
    z<-rep(0,T)
    z.idx <- z[1]+1
    y<-matrix(0,nrow=T,ncol=2)
    y[1,]<-rmvnorm(1,mu[z.idx,],Delta)
    rho<-(.75)
    M<-rho*diag(2)
    ## create the time stamp sequence (just to help with part date (still my bday))
    ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)
    tau<-192 ## change point is halfway
    
    ## loop to simulate data points
    for (i in 2:T){
      z[i]<-ifelse(i<tau,0,1)
      z.idx<-z[i]+1
      psi<-t(M%*%y[i-1,])+t((diag(2)-M)%*%(mu[z.idx,]))
      y[i,]<-rmvnorm(1,psi,Delta)
    }
    
    ## create sim data frame
    sim.data.final <- data.frame(x = y[,1],
                                 y = y[,2],
                                 state = z,
                                 ts = ts,
                                 datasim = j,
                                 siminfo = "hbc",
                                 siminfoDetail = diff.for.HBC.quotes[k],
                                 hmd = "48hr",
                                 hmf = "15min",
                                 hbc = diff.for.HBC.quotes[k])
    
    ## put the df into list
    HBC.data.final[[k]] <- sim.data.final
  }
  ## put that list into another list
  HBC.data[[j]] <- HBC.data.final
}

## put length of data into own data frames
hbc.100 <- NULL;hbc.75 <- NULL;hbc.50 <- NULL;hbc.25 <- NULL
for(i in 1:length(HBC.data)){
  tmp.l <- HBC.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hbc <- unique(tmp$siminfoDetail)
    if(hbc == "100"){
      hbc.100 <- rbind(hbc.100,
                        tmp)
    }
    if(hbc == "75"){
      hbc.75 <- rbind(hbc.75,
                      tmp)
    }
    if(hbc == "50"){
      hbc.50 <- rbind(hbc.50,
                      tmp)
    }
    if(hbc == "25"){
      hbc.25 <- rbind(hbc.25,
                        tmp)
    }
  }
}

hbc.list <- list(hbc.100,hbc.75,
                 hbc.50,hbc.25)

lcpmSimData[[1]] <- hbc.list

## HMF - LCPM ----
## set seed
set.seed(10120)
HMF.data.final <- list()
HMF.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## simulate the data
  T<-385 ##(gives 48 hours pre/post at 15 min fixes)
  
  mu <- matrix(c(max.cent.diff$Ecent_mean_before,
                 max.cent.diff$Ecent_mean_after,
                 max.cent.diff$Ncent_mean_before,
                 max.cent.diff$Ncent_mean_after),nrow=2)
  # mean.sig <- c(mean(gcp.sim.params$Ecent_var),mean(gcp.sim.params$Ncent_var))
  Delta<-matrix(c(max.cent.N.var,0,0,max.cent.N.var),nrow=2,byrow=2) ## variance
  z<-rep(0,T)
  z.idx <- z[1]+1
  y<-matrix(0,nrow=T,ncol=2)
  y[1,]<-rmvnorm(1,mu[z.idx,],Delta)
  #
  rho<-matrix(c(.8,0,0,.8),byrow=T,nrow=2)
  M<-rho*diag(2)
  tau<-192
  z<-rep(0,T)
  ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)
  
  ## create locations based on parameters
  for (i in 2:T){
    z[i]<-ifelse(i<tau,0,1)
    z.idx<-z[i]+1
    psi<-t(M%*%y[i-1,])+t((diag(2)-M)%*%(mu[z.idx,]))
    y[i,]<-rmvnorm(1,psi,Delta)
  }
  
  ## create data frame of data
  sim.data.final <- data.frame(long = y[,1],
                               lat = y[,2],
                               state = z,
                               ts = ts,
                               datasim = j,
                               siminfo = "hmf")
  
  ## truncate and save as needed in a loop
  how.much.data <- 24 ## in hours after change
  how.much.data.quotes <- c("24hr")  ## for hours after change
  how.many.fixes <- c(1,2,4)
  how.many.fixes.quotes <- c("15min","30min","60min")
  
  ## take the full data set and truncate "fixes"
  for(k in 1:length(how.many.fixes)){
    if(k == 1){
      sim.data.fix <- sim.data.final %>% 
        dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data))) %>% 
        dplyr::mutate(
            siminfoDetail = how.many.fixes.quotes[k],
            hmd = how.much.data.quotes,
            hmf = how.many.fixes.quotes[k])
      # write.csv(sim.data.fix,paste("simulations/sim_dataset_HMD_HMF_BCP",how.much.data.quotes[k],
      #                            how.many.fixes.quotes[k],".csv",sep="_"),row.names=FALSE)
      
    }else{
      sim.data.fix <- sim.data.final %>% 
        dplyr::filter(row_number() %% how.many.fixes[k] == 1) %>%
        dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data))) %>% 
        dplyr::mutate(
            siminfoDetail = how.many.fixes.quotes[k],
            hmd = how.much.data.quotes,
            hmf = how.many.fixes.quotes[k])
      # write.csv(sim.data.fix,paste("simulations/sim_dataset_HMD_HMF_BCP",how.much.data.quotes[k],
      #                              how.many.fixes.quotes[k],".csv",sep="_"),row.names=FALSE)
      # assign(paste("sim_dataset_",how.much.data.quotes[k],"_15min",sep=""),sim.data.15min)
    }
    HMF.data.final[[k]] <- sim.data.fix
  }
  HMF.data[[j]] <- HMF.data.final
}

## put length of data into own data frames
hmf.15 <- NULL;hmf.30 <- NULL;hmf.60 <- NULL
for(i in 1:length(HMF.data)){
  tmp.l <- HMF.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hmf <- unique(tmp$siminfoDetail)
    if(hmf == "15min"){
      hmf.15 <- rbind(hmf.15,
                      tmp)
    }
    if(hmf == "30min"){
      hmf.30 <- rbind(hmf.30,
                      tmp)
    }
    if(hmf == "60min"){
      hmf.60 <- rbind(hmf.60,
                      tmp)
    }
  }
}

hmf.list <- list(hmf.15,hmf.30,
                 hmf.60)
lcpmSimData[[2]] <- hmf.list

## HMD - LCPM ----
# set seed
set.seed(10120)


HMD.data.final <- list()
HMD.data <- list()
## create 3 data sets to be run by 3 functions
for(j in 1:n.data){
  ## simulate the data
  T<-385 ##(gives 48 hours pre/post at 15 min fixes)
  
  mu <- matrix(c(max.cent.diff$Ecent_mean_before,
                 max.cent.diff$Ecent_mean_after,
                 max.cent.diff$Ncent_mean_before,
                 max.cent.diff$Ncent_mean_after),nrow=2)
  # mean.sig <- c(mean(gcp.sim.params$Ecent_var),mean(gcp.sim.params$Ncent_var))
  Delta<-matrix(c(max.cent.N.var,0,0,max.cent.N.var),nrow=2,byrow=2) ## variance
  z<-rep(0,T)
  z.idx <- z[1]+1
  y<-matrix(0,nrow=T,ncol=2)
  y[1,]<-rmvnorm(1,mu[z.idx,],Delta)
  #
  rho<-0.8
  M<-rho*diag(2)
  tau<-192
  z<-rep(0,T)
  ts <- seq(ymd_hm("1994-05-17 12:28"),by="15 min",length.out=T)
  
  ## create locations based on parameters
  for (i in 2:T){
    z[i]<-ifelse(i<tau,1,2)
    z.idx<-z[i]
    psi<-t(M%*%y[i-1,])+t((diag(2)-M)%*%(mu[z.idx,]))
    y[i,]<-rmvnorm(1,psi,Delta)
  }
  
  ## create data frame of data
  sim.data.final <- data.frame(x = y[,1],
                               y = y[,2],
                               state = z,
                               ts = ts,
                               datasim = j,
                               siminfo = "hmd")
  
  ## truncate and save as needed in a loop
  how.much.data <- c(48,24,12,6,3) ## in hours
  how.much.data.quotes <- c("48hr","24hr","12hr","6hr","3hr") 
  
  ## take the full data set and truncate
  for(k in 1:length(how.much.data)){
    sim.data.15min <- sim.data.final %>% 
      dplyr::filter(ts <= (ymd_hm("1994-05-19 12:13")+hours(how.much.data[k])+minutes(15))) %>% 
      dplyr::mutate(
        siminfoDetail = how.much.data.quotes[k],
        hmd = how.much.data.quotes[k],
        hmf = "15min")
    # write.csv(sim.data.15min,paste("simulations/sim_dataset_HMD_HMF_BCP_",how.much.data.quotes[k],"_15min.csv",sep=""),row.names=FALSE)
    # assign(paste("sim_dataset_",how.much.data.quotes[k],"_15min",sep=""),sim.data.15min)
    
    HMD.data.final[[k]] <- sim.data.15min
  }
  HMD.data[[j]] <- HMD.data.final
}

## put length of data into own data frames
hmd.48 <- NULL;hmd.24 <- NULL;hmd.12 <- NULL;hmd.06 <- NULL;hmd.03 <- NULL
for(i in 1:length(HMD.data)){
  tmp.l <- HMD.data[[i]]
  for(j in 1:length(tmp.l)){
    tmp <- tmp.l[[j]]
    hmd <- unique(tmp$siminfoDetail)
    if(hmd == "48hr"){
      hmd.48 <- rbind(hmd.48,
                      tmp)
    }
    if(hmd == "24hr"){
      hmd.24 <- rbind(hmd.24,
                      tmp)
    }
    if(hmd == "12hr"){
      hmd.12 <- rbind(hmd.12,
                      tmp)
    }
    if(hmd == "6hr"){
      hmd.06 <- rbind(hmd.06,
                      tmp)
    }
    if(hmd == "3hr"){
      hmd.03 <- rbind(hmd.03,
                      tmp)
    }
  }
}

hmd.list <- list(hmd.48,hmd.24,
                 hmd.12,hmd.06,hmd.03)
lcpmSimData[[3]] <- hmd.list

## save lcpm simulated data as .RDS ----
saveRDS(
    lcpmSimData,
    "indata/lcpmSimData.RData"
)
