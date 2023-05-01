#Program: Sims 07.22.22.R
#Developed by: BES
#Purpose: Simulation study for single v multi-span bridging estimators, determine truth for each scenario empirically
#Created: 07.22.22, last updated 09.07.22

library(boot)

#specify number of participants in each population and the number of scenarios
num.part.NF<-20000000
num.part.F<-20000000
num.Scen<-5 

set.seed(09082022)

truth.res<-rep(NA, num.Scen)

for (Scen in 1:num.Scen){

#assign trial indicator and randomize treatment
num.part<-num.part.NF+num.part.F
part.id<-seq(1:num.part)
dat<-as.data.frame(part.id)
dat$S<-ifelse(dat$part.id<=num.part.NF,1,2)
dat$rand<-rbinom(num.part,1,0.5)
dat$A<-ifelse(dat$S==1,1+dat$rand,2+dat$rand)

#generate measured confounders - IDU and BL CD4, UC, and potential outcomes (scenario-specific) for each scenario
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
  
dat$cd4wk8pre_1<-rnorm(num.part, mean = (dat$mu2_1), sd = 20)
dat$cd4wk8pre_2<-rnorm(num.part, mean = (dat$mu2_2), sd = 20)
dat$cd4wk8pre_3<-rnorm(num.part, mean = (dat$mu2_3), sd = 20)
  
dat$cd4wk8_1<-ifelse(dat$cd4wk8pre_1<0,0,dat$cd4wk8pre_1)
dat$cd4wk8_2<-ifelse(dat$cd4wk8pre_2<0,0,dat$cd4wk8pre_2)
dat$cd4wk8_3<-ifelse(dat$cd4wk8pre_3<0,0,dat$cd4wk8pre_3)

dat$cd4wk8_1_bin<-ifelse(dat$cd4wk8_1>250,1,0)
dat$cd4wk8_2_bin<-ifelse(dat$cd4wk8_2>250,1,0)
dat$cd4wk8_3_bin<-ifelse(dat$cd4wk8_3>250,1,0)

truth<-with(dat, mean(cd4wk8_3[S=='2']))-with(dat, mean(cd4wk8_1[S=='2'])) 
truth_bin<-with(dat, mean(cd4wk8_3_bin[S=='2']))-with(dat, mean(cd4wk8_1_bin[S=='2'])) 
  
#combine results
res<-cbind(Scen, truth, truth_bin)
  
truth.res<-rbind(truth.res,res)
  
rm(list= ls()[!(ls() %in% c('num.Scen','truth.res','num.part.NF','num.part.F'))])

}


truth.res2<-as.data.frame(truth.res[-1,])

write.csv(truth.res2, "Sim_truth_091222.csv")

