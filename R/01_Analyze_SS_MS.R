#Program: Sims 07.22.22.R
#Developed by: BES 
#Purpose: Simulation study for single v multi-span bridging estimators
#Created: 07.22.22, last updated 01.31.23


t1 <- Sys.time()

library(geex)
library(data.table)

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
subdir_name <- as.character(args[1])

user_home_dir <- " " # enter path for top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")

# load SS/MS estimation functions 
source('00_Estimators_01.31.23.R')

sim <- Sys.getenv("SLURM_ARRAY_TASK_ID") 

setwd(dir_path)
alldat<-fread(file = paste("DGM_",subdir_name,".csv",sep=""))

dat<-alldat[alldat$iter==sim,]
dat$miss<-as.numeric(is.na(dat$cd4wk8obs))

#estimate IPTWs
dat320<-dat[dat$S==1,]
dat388<-dat[dat$S==2,]

dat320$prob3<-0.5
dat320$Wa<-ifelse(dat320$A==2,(dat320$prob3)^(-1),(1-dat320$prob3)^(-1))

dat388$prob3<-0.5
dat388$Wa<-ifelse(dat388$A==2,(dat388$prob3)^(-1),(1-dat388$prob3)^(-1))

trtall<-rbind(dat320,dat388)
trtall<-trtall[,c('part.id','Wa')]

dat<-merge(dat,trtall,by="part.id")

# estimate IPMWs
dat320.1<-dat[dat$S==1 & dat$A==1,]
dat320.2<-dat[dat$S==1 & dat$A==2,]
dat388.2<-dat[dat$S==2 & dat$A==2,]
dat388.3<-dat[dat$S==2 & dat$A==3,]

missmod.320.1 <- glm((miss==0) ~ as.factor(IDU),family = binomial(), data = dat320.1)
dat320.1$probnmiss<-predict(missmod.320.1, dat320.1, type="response")
dat320.1$Wm<-(dat320.1$probnmiss)^(-1)

missmod.320.2 <- glm((miss==0) ~ as.factor(IDU),family = binomial(), data = dat320.2)
dat320.2$probnmiss<-predict(missmod.320.2, dat320.2, type="response")
dat320.2$Wm<-(dat320.2$probnmiss)^(-1)

missmod.388.2 <- glm((miss==0) ~ as.factor(IDU),family = binomial(), data = dat388.2)
dat388.2$probnmiss<-predict(missmod.388.2, dat388.2, type="response")
dat388.2$Wm<-(dat388.2$probnmiss)^(-1)

missmod.388.3 <- glm((miss==0) ~ as.factor(IDU),family = binomial(), data = dat388.3)
dat388.3$probnmiss<-predict(missmod.388.3, dat388.3, type="response")
dat388.3$Wm<-(dat388.3$probnmiss)^(-1)

miss.all<-rbind(dat320.1,dat320.2,dat388.2,dat388.3)
miss.all<-miss.all[,c('part.id','Wm')]

dat<-merge(dat,miss.all,by="part.id")

#Estimate sampling weights
ssdat<-dat[dat$A==1|dat$A==3,]
transp.model.ss <- glm((S==1) ~ as.factor(IDU) + blcd4, family = binomial(), data = ssdat)
ssdat$prob320.ss<-predict(transp.model.ss, ssdat, type="response")
ssdat$Wtr.ss<-ifelse(ssdat$S==1,(1-ssdat$prob320.ss)/ssdat$prob320.ss,1)

transp.model.ms <- glm((S==1) ~ as.factor(IDU) + blcd4, family = binomial(), data = dat)
dat$prob320.ms<-predict(transp.model.ms, dat, type="response")
dat$Wtr.ms<-ifelse(dat$S==1,(1-dat$prob320.ms)/dat$prob320.ms,1)

#calculate overall weights
dat$WTFIN<-dat$Wa*dat$Wtr.ms*dat$Wm
ssdat$WTFIN<-ssdat$Wa*ssdat$Wtr.ss*ssdat$Wm

#derive variables needed in geex
dat$dual<-ifelse(dat$A==1,1,0)
dat$triple<-ifelse(dat$A==2,1,0)
dat$triple.loc<-ifelse(dat$A==2 & dat$S==2,1,0)
dat$triple.dis<-ifelse(dat$A==2 & dat$S==1,1,0)
dat$quad<-ifelse(dat$A==3,1,0)
dat$ACTG320<-ifelse(dat$S==1,1,0)
dat$sstransp<-ifelse(dat$A==2,0,1)
dat$nmissout<-1-dat$miss

#### estimate the treatment effect based on four proposed estimators
#naive - continuous outcome
naive.3<-with(dat, mean(cd4wk8obs[A == '3' & miss=='0']))
naive.1<-with(dat, mean(cd4wk8obs[A == '1' & miss=='0']))
naive.ss<-with(dat, mean(cd4wk8obs[A == '3' & miss=='0']))-with(dat, mean(cd4wk8obs[A == '1' & miss=='0'])) 

#naive - binary outcome
naive.3.bin<-with(dat, mean(cd4wk8obs_bin[A == '3' & miss=='0']))
naive.1.bin<-with(dat, mean(cd4wk8obs_bin[A == '1' & miss=='0']))
naive.ss_bin<-with(dat, mean(cd4wk8obs_bin[A == '3' & miss=='0']))-with(dat, mean(cd4wk8obs_bin[A == '1' & miss=='0'])) 

t2 <- Sys.time()
print(round(t2 - t1, 4))
SS_res_naive<-SS_naive_est(data=dat,
                           estfun=estfun_naive_SS, 
                           EL.IV=naive.3,
                           ED.IV=naive.1,
                           RD.IV=naive.ss
)

t3 <- Sys.time()
print(round(t3 - t2, 4))
SS_res_naive_bin<-SS_naive_est(data=dat,
                               estfun=estfun_naive_SS_bin, 
                               EL.IV=naive.3.bin,
                               ED.IV=naive.1.bin,
                               RD.IV=naive.ss_bin
)
colnames(SS_res_naive_bin)<-c("ATE.naive.SS.bin","seATE.naive.SS.bin")

naive.2.l<-with(dat, mean(cd4wk8obs[A == '2' & S=='2' & miss=='0']))
naive.2.d<-with(dat, mean(cd4wk8obs[A == '2' & S=='1' & miss=='0']))
naive.ms<-with(dat, mean(cd4wk8obs[A == '3' & miss=='0']))-with(dat, mean(cd4wk8obs[A == '2' & S=='2' & miss=='0']))+with(dat, mean(cd4wk8obs[A == '2' & S=='1' & miss=='0']))-with(dat, mean(cd4wk8obs[A == '1' & miss=='0'])) 
naive.ms.diag<-naive.2.l-naive.2.d

naive.2.l.bin<-with(dat, mean(cd4wk8obs_bin[A == '2' & S=='2' & miss=='0']))
naive.2.d.bin<-with(dat, mean(cd4wk8obs_bin[A == '2' & S=='1' & miss=='0']))
naive.ms_bin<-with(dat, mean(cd4wk8obs_bin[A == '3' & miss=='0']))-with(dat, mean(cd4wk8obs_bin[A == '2' & S=='2' & miss=='0']))+with(dat, mean(cd4wk8obs_bin[A == '2' & S=='1' & miss=='0']))-with(dat, mean(cd4wk8obs_bin[A == '1' & miss=='0'])) 
naive.ms.diag_bin<-naive.2.l.bin-naive.2.d.bin

t4 <- Sys.time()
print(round(t4 - t3, 4))
MS_res_naive<-MS_naive_est(data=dat,
                           estfun=estfun_naive_MS, 
                           EL3.IV=naive.3,
                           EL2.IV=naive.2.l,
                           ED2.IV=naive.2.d,
                           ED1.IV=naive.1,
                           RD.IV=naive.ms,
                           Diag.IV=naive.ms.diag
)
colnames(MS_res_naive)<-c("ATE.naive.MS","seATE.naive.MS","Naive.MS.Diag","seNaive.MS.Diag")

t5 <- Sys.time()
print(round(t5 - t4, 4))
MS_res_naive_bin<-MS_naive_est(data=dat,
                               estfun=estfun_naive_MS_bin, 
                               EL3.IV=naive.3.bin,
                               EL2.IV=naive.2.l.bin,
                               ED2.IV=naive.2.d.bin,
                               ED1.IV=naive.1.bin,
                               RD.IV=naive.ms_bin,
                               Diag.IV=naive.ms.diag_bin
)
colnames(MS_res_naive_bin)<-c("ATE.naive.MS.bin","seATE.naive.MS.bin","Naive.MS.Diag.bin","seNaive.MS.Diag.bin")

### SS estimator
dat.obs<-dat[dat$miss==0,]
ssdat.obs<-ssdat[ssdat$miss==0,]
ss.Ehat_L3<-sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==3)*ssdat.obs$cd4wk8obs)/sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==3))
ss.Ehat_D1<-sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==1)*ssdat.obs$cd4wk8obs)/sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==1))
ATE_ss.quad1<-ss.Ehat_L3-ss.Ehat_D1

ss.Ehat_L3_bin<-sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==3)*ssdat.obs$cd4wk8obs_bin)/sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==3))
ss.Ehat_D1_bin<-sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==1)*ssdat.obs$cd4wk8obs_bin)/sum(ssdat.obs$WTFIN*as.numeric(ssdat.obs$A==1))
ATE_ss.quad1_bin<-ss.Ehat_L3_bin-ss.Ehat_D1_bin

#continuous outcome
t6 <- Sys.time()
print(round(t6 - t5, 4))
SS_res<-SS_est(data=dat,
               estfun=estfun_SS, 
               tr_formula=ACTG320 ~ sstransp + sstransp:(IDU + blcd4) -1, 
               md_formula=nmissout ~ dual + dual:(IDU)-1, 
               ml_formula=nmissout ~ quad + quad:(IDU)-1,
               EL.IV=ss.Ehat_L3,
               ED.IV=ss.Ehat_D1,
               RD.IV=ATE_ss.quad1
)

#binary outcome
t7 <- Sys.time()
print(round(t7 - t6, 4))
SS_res_bin<-SS_est(data=dat,
                   estfun=estfun_SS_bin, 
                   tr_formula=ACTG320 ~ sstransp + sstransp:(IDU + blcd4) -1, 
                   md_formula=nmissout ~ dual + dual:(IDU)-1, 
                   ml_formula=nmissout ~ quad + quad:(IDU)-1,
                   EL.IV=ss.Ehat_L3_bin,
                   ED.IV=ss.Ehat_D1_bin,
                   RD.IV=ATE_ss.quad1_bin
)
colnames(SS_res_bin)<-c("ATE.SS.bin","seATE.SS.bin")


### MS estimator
ms.Ehat_L3<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==3)*as.numeric(dat.obs$S==2)*dat.obs$cd4wk8obs)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==3)*as.numeric(dat.obs$S==2))
ms.Ehat_L2<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==2)*dat.obs$cd4wk8obs)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==2))
ms.Ehat_D2<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==1)*dat.obs$cd4wk8obs)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==1))
ms.Ehat_D1<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==1)*as.numeric(dat.obs$S==1)*dat.obs$cd4wk8obs)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==1)*as.numeric(dat.obs$S==1))

ATE_ms.quad1<-ms.Ehat_L3-ms.Ehat_L2+ms.Ehat_D2-ms.Ehat_D1
Diagnostic<-ms.Ehat_L2-ms.Ehat_D2

ms.Ehat_L3_bin<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==3)*as.numeric(dat.obs$S==2)*dat.obs$cd4wk8obs_bin)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==3)*as.numeric(dat.obs$S==2))
ms.Ehat_L2_bin<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==2)*dat.obs$cd4wk8obs_bin)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==2))
ms.Ehat_D2_bin<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==1)*dat.obs$cd4wk8obs_bin)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==2)*as.numeric(dat.obs$S==1))
ms.Ehat_D1_bin<-sum(dat.obs$WTFIN*as.numeric(dat.obs$A==1)*as.numeric(dat.obs$S==1)*dat.obs$cd4wk8obs_bin)/sum(dat.obs$WTFIN*as.numeric(dat.obs$A==1)*as.numeric(dat.obs$S==1))

ATE_ms.quad1_bin<-ms.Ehat_L3_bin-ms.Ehat_L2_bin+ms.Ehat_D2_bin-ms.Ehat_D1_bin
Diagnostic_bin<-ms.Ehat_L2_bin-ms.Ehat_D2_bin

#continuous outcome
t8 <- Sys.time()
print(round(t8 - t7, 4))
MS_res<-MS_est(data=dat,
               estfun=estfun_MS, 
               tr_formula=ACTG320 ~ IDU + blcd4, 
               md1_formula=nmissout ~ dual + dual:(IDU) -1, 
               md2_formula=nmissout ~ triple.dis + triple.dis:(IDU) -1, 
               ml2_formula=nmissout ~ triple.loc + triple.loc:(IDU) -1,               
               ml3_formula=nmissout ~ quad + quad:(IDU)-1,
               EL3.IV=ms.Ehat_L3,
               EL2.IV=ms.Ehat_L2,
               ED2.IV=ms.Ehat_D2,
               ED1.IV=ms.Ehat_D1,
               RD.IV=ATE_ms.quad1,
               Diag.IV=Diagnostic
)
colnames(MS_res)<-c("ATE.MS","seATE.MS","MS.Diag","seMS.Diag")

#binary outcome
t9 <- Sys.time()
print(round(t9 - t8, 4))
MS_res_bin<-MS_est(data=dat,
                   estfun=estfun_MS_bin, 
                   tr_formula=ACTG320 ~ IDU + blcd4, 
                   md1_formula=nmissout ~ dual + dual:(IDU) -1, 
                   md2_formula=nmissout ~ triple.dis + triple.dis:(IDU) -1, 
                   ml2_formula=nmissout ~ triple.loc + triple.loc:(IDU) -1,               
                   ml3_formula=nmissout ~ quad + quad:(IDU)-1,
                   EL3.IV=ms.Ehat_L3_bin,
                   EL2.IV=ms.Ehat_L2_bin,
                   ED2.IV=ms.Ehat_D2_bin,
                   ED1.IV=ms.Ehat_D1_bin,
                   RD.IV=ATE_ms.quad1_bin,
                   Diag.IV=Diagnostic_bin
)
colnames(MS_res_bin)<-c("ATE.MS.bin","seATE.MS.bin","MS.Diag.bin","seMS.Diag.bin")

t10 <- Sys.time()
print(round(t10 - t9, 4))
print(round(t10 - t1, 4))

#combine results and output
res<-cbind(SS_res_naive, SS_res_naive_bin, SS_res, SS_res_bin, MS_res_naive, MS_res_naive_bin, MS_res, MS_res_bin)

output_filename <- paste(dir_path,"/results_", sim, ".csv", sep="")
write.csv(res, output_filename)


