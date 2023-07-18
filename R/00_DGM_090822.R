#Program: 00_DGM_090722.R
#Developed by: BES
#Purpose: Simulation study for single v multi-span bridging estimators
#Created: 07.22.22, last updated 09.08.22


user_home_dir <- " " # enter path for top-level directory for the project


## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
num.sims <- as.numeric(args[1]) #number of simulations
Scen <- as.numeric(args[2]) #scenario (1-5)
num.part.NF <- as.numeric(args[3]) #number of participants, non-focal population
num.part.F <- as.numeric(args[4]) #number of participants, target population
subdir_name <- as.character(args[5]) #directory where output will be stored

DGM<-rep(NA, num.sims)

library(boot)

set.seed(09082022)

for (k in 1:num.sims){

  #assign trial indicator and randomize treatment
  num.part<-num.part.NF+num.part.F #total number of participants
  part.id<-seq(1:num.part)
  dat<-as.data.frame(part.id)
  dat$S<-ifelse(dat$part.id<=num.part.NF,1,2) #sampling indicator
  dat$rand<-rbinom(num.part,1,0.5)
  dat$A<-ifelse(dat$S==1,1+dat$rand,2+dat$rand) #treatment level

  #generate measured confounders - IDU and BL CD4, UC, potential outcomes (scenario-specific), and missing indicator for each scenario
 if(Scen==1) { 
    dat$IDUprob<-0.25
    dat$IDU<-rbinom(num.part, 1, dat$IDUprob) 
    dat$mu1<-(175-10*dat$IDU)
    dat$blcd4pre<-rnorm(num.part, mean = (dat$mu1), sd = 30)
    dat$blcd4<-ifelse(dat$blcd4pre<0,0,dat$blcd4pre)
    dat$mu2_1<-50-5*dat$IDU+1.1*dat$blcd4
    dat$mu2_2<-80-5*dat$IDU+1.1*dat$blcd4
    dat$mu2_3<-110-5*dat$IDU+1.1*dat$blcd4
    dat$mprob<-0.15
}

if(Scen==2) {
  dat$IDUprob<-ifelse(dat$S==1,0.5,0.2)
  dat$IDU<-rbinom(num.part, 1, dat$IDUprob) 
  dat$mu1<-(175+10*(dat$S==2)-20*dat$IDU)
  dat$blcd4pre<-rnorm(num.part, mean = (dat$mu1), sd = 30)
  dat$blcd4<-ifelse(dat$blcd4pre<0,0,dat$blcd4pre)
  dat$mu2_1<-35-80*dat$IDU+1.0*dat$blcd4 
  dat$mu2_2<-30-10*dat$IDU+1.1*dat$blcd4
  dat$mu2_3<-40+20*dat$IDU+1.2*dat$blcd4
  dat$mprob<-ifelse(dat$S==1,inv.logit(-2.0+0.5*dat$IDU),inv.logit(-2.1+0.5*dat$IDU))
}

if(Scen==3) {
  dat$IDUprob<-ifelse(dat$S==1,0.5,0.2)
  dat$IDU<-rbinom(num.part, 1, dat$IDUprob) 
  dat$mu1<-(175+10*(dat$S==2)-20*dat$IDU)
  dat$blcd4pre<-rnorm(num.part, mean = (dat$mu1), sd = 30)
  dat$blcd4<-ifelse(dat$blcd4pre<0,0,dat$blcd4pre)
  dat$mu2_1<-35-80*dat$IDU+1.0*dat$blcd4 
  dat$mu2_2<-ifelse(dat$S==1,45-10*dat$IDU+1.1*dat$blcd4,40+10*dat$IDU+1.0*dat$blcd4)
  dat$mu2_3<-40+20*dat$IDU+1.2*dat$blcd4
  dat$mprob<-ifelse(dat$S==1,inv.logit(-2.0+0.5*dat$IDU),inv.logit(-2.1+0.5*dat$IDU))
}

if(Scen==4) {
  dat$IDUprob<-ifelse(dat$S==1,0.5,0.2)
  dat$IDU<-rbinom(num.part, 1, dat$IDUprob) 
  dat$mu1<-(175+10*(dat$S==2)-20*dat$IDU)
  dat$blcd4pre<-rnorm(num.part, mean = (dat$mu1), sd = 30)
  dat$blcd4<-ifelse(dat$blcd4pre<0,0,dat$blcd4pre)
  dat$mu2_1<-45-80*dat$IDU+1.0*dat$blcd4-30*(dat$S==1)
  dat$mu2_2<-30-10*dat$IDU+1.1*dat$blcd4
  dat$mu2_3<-30+20*dat$IDU+1.2*dat$blcd4+20*(dat$S==2)
  dat$mprob<-ifelse(dat$S==1,inv.logit(-2.0+0.5*dat$IDU),inv.logit(-2.1+0.5*dat$IDU))
}

if(Scen==5) {
  dat$IDUprob<-ifelse(dat$S==1,0.5,0.2)
  dat$IDU<-rbinom(num.part, 1, dat$IDUprob) 
  dat$mu1<-(175+10*(dat$S==2)-20*dat$IDU)
  dat$blcd4pre<-rnorm(num.part, mean = (dat$mu1), sd = 30)
  dat$blcd4<-ifelse(dat$blcd4pre<0,0,dat$blcd4pre)
  dat$mu2_1<-45-80*dat$IDU+1.0*dat$blcd4-30*(dat$S==1)
  dat$mu2_2<-30-10*dat$IDU+1.1*dat$blcd4+10*(dat$S==1)
  dat$mu2_3<-30+20*dat$IDU+1.2*dat$blcd4+20*(dat$S==2)
  dat$mprob<-ifelse(dat$S==1,inv.logit(-2.0+0.5*dat$IDU),inv.logit(-2.1+0.5*dat$IDU))
}
  #generate potential outcomes, continuous and binary
  
  dat$cd4wk8pre_1<-rnorm(num.part, mean = (dat$mu2_1), sd = 20)
  dat$cd4wk8pre_2<-rnorm(num.part, mean = (dat$mu2_2), sd = 20)
  dat$cd4wk8pre_3<-rnorm(num.part, mean = (dat$mu2_3), sd = 20)
  
  dat$cd4wk8_1<-ifelse(dat$cd4wk8pre_1<0,0,dat$cd4wk8pre_1)
  dat$cd4wk8_2<-ifelse(dat$cd4wk8pre_2<0,0,dat$cd4wk8pre_2)
  dat$cd4wk8_3<-ifelse(dat$cd4wk8pre_3<0,0,dat$cd4wk8pre_3)
  
  dat$cd4wk8_1_bin<-ifelse(dat$cd4wk8_1>250,1,0)
  dat$cd4wk8_2_bin<-ifelse(dat$cd4wk8_2>250,1,0)
  dat$cd4wk8_3_bin<-ifelse(dat$cd4wk8_3>250,1,0)

  #generate missing data indicator
  dat$miss<-rbinom(num.part,1,dat$mprob)

  #generate observed outcome - dependent on assigned treatment, missing data indicator, and potential outcomes
  dat$cd4wk8obs<-ifelse(dat$miss==0,dat$cd4wk8_1*(dat$A==1)+dat$cd4wk8_2*(dat$A==2)+dat$cd4wk8_3*(dat$A==3),NA)
  dat$cd4wk8obs_bin<-ifelse(dat$miss==0,dat$cd4wk8_1_bin*(dat$A==1)+dat$cd4wk8_2_bin*(dat$A==2)+dat$cd4wk8_3_bin*(dat$A==3),NA)
  
  iter<-dat[,c('part.id','S','A','IDU','blcd4','cd4wk8obs','cd4wk8obs_bin')]
  iter$iter<-k

  DGM<-rbind(DGM,iter)
  
}

  ### create a directory for these sims
  dir_path <- paste(user_home_dir,subdir_name,"/",sep="")
  dir.create(dir_path)
  setwd(dir_path)
  DGM2<-as.data.frame(DGM[-1,])
  write.csv(DGM2, file = paste("DGM_",subdir_name,".csv",sep=""))



