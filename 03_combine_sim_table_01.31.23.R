#Program: Combine Sim results and output table and Figures
#Developed by: BES
#Created: 07.28.22, last updated 01.31.23

#bring in true value for each scenario
truth<-as.data.frame(read.csv("Sim_truth_091222.csv"))

#bring in results for each scenario and sample size considered
Sc1_pr<-as.data.frame(read.csv("results_N1_1000_N2_400_SC1.csv"))
Sc2_pr<-as.data.frame(read.csv("results_N1_1000_N2_400_SC2.csv"))
Sc3_pr<-as.data.frame(read.csv("results_N1_1000_N2_400_SC3.csv"))
Sc4_pr<-as.data.frame(read.csv("results_N1_1000_N2_400_SC4.csv"))
Sc5_pr<-as.data.frame(read.csv("results_N1_1000_N2_400_SC5.csv"))

Sc1_alt<-as.data.frame(read.csv("results_N1_1000_N2_2000_SC1.csv"))
Sc2_alt<-as.data.frame(read.csv("results_N1_1000_N2_2000_SC2.csv"))
Sc3_alt<-as.data.frame(read.csv("results_N1_1000_N2_2000_SC3.csv"))
Sc4_alt<-as.data.frame(read.csv("results_N1_1000_N2_2000_SC4.csv"))
Sc5_alt<-as.data.frame(read.csv("results_N1_1000_N2_2000_SC5.csv"))

Sc1_alt2<-as.data.frame(read.csv("results_N1_400_N2_1000_SC1.csv"))
Sc2_alt2<-as.data.frame(read.csv("results_N1_400_N2_1000_SC2.csv"))
Sc3_alt2<-as.data.frame(read.csv("results_N1_400_N2_1000_SC3.csv"))
Sc4_alt2<-as.data.frame(read.csv("results_N1_400_N2_1000_SC4.csv"))
Sc5_alt2<-as.data.frame(read.csv("results_N1_400_N2_1000_SC5.csv"))


#calculate summary measures for continuous outcome (bias, coverage, ASE, ESE, SER, diagnostic characteristics)
calc_bias <- function(ATEest,seATEest,truth,scen,est,diag,sediag){
  ATEbias<-mean(ATEest, na.rm = TRUE)-truth
  ATEpctbias<-mean((ATEest-truth)/truth)*100
  ll<-ATEest-1.96*seATEest
  ul<-ATEest+1.96*seATEest
  cov <-100*mean((ll <= truth & ul >= truth))
  ASE<-mean(na.omit(seATEest))
  ESE<-sd(na.omit(ATEest))
  SER<-ASE/ESE
  RMSE<-sqrt(ATEbias^2+ESE^2)
  lldiag<-diag-1.96*sediag
  uldiag<-diag+1.96*sediag
  meandiag<-mean(diag, na.rm = TRUE)
  diagzero<-100*mean((lldiag <= 0 & uldiag >= 0))
  results<-c(scen, "&", round(ATEbias,digits=1),"&", round(ASE, digits=2),"&", round(ESE, digits=2),"&", round(SER, digits=2),"&",round(RMSE, digits=1),"&", round(cov, digits=0),"&", round(meandiag,digits=1),"&", round(diagzero, digits=0),'\\')
}

#summarize sim results for all 4 estimators and combine
prep_data <- function(dat,Scen){
  naive.MS<-calc_bias(dat$ATE.naive.MS,dat$seATE.naive.MS,truth$truth[Scen],Scen,'Naive MS',dat$Naive.MS.Diag,dat$seNaive.MS.Diag)
  naive.SS<-calc_bias(dat$ATE.naive.SS,dat$seATE.naive.SS,truth$truth[Scen],Scen,'Naive SS',dat$Naive.MS.Diag,dat$seNaive.MS.Diag)
  MS<-      calc_bias(dat$ATE.MS,dat$seATE.MS,truth$truth[Scen],Scen,'Bridging MS',dat$MS.Diag,dat$seMS.Diag)
  SS<-      calc_bias(dat$ATE.SS,dat$seATE.SS,truth$truth[Scen],Scen,'Bridging SS',dat$MS.Diag,dat$seMS.Diag)
  combine<-rbind(naive.MS,naive.SS,MS,SS)
}

#calculate summary measures for binary outcome (bias, coverage, ASE, ESE, SER, diagnostic characteristics)
calc_bias_bin <- function(ATEest,seATEest,truth,scen,est,diag,sediag){
  ATEbias<-mean(100*ATEest, na.rm = TRUE)-100*truth
  ATEpctbias<-mean((100*ATEest-100*truth)/(100*truth))*100
  ll<-ATEest-1.96*seATEest
  ul<-ATEest+1.96*seATEest
  cov <-100*mean((ll <= truth & ul >= truth))
  ASE<-mean(na.omit(100*seATEest))
  ESE<-sd(na.omit(100*ATEest))
  SER<-ASE/ESE
  RMSE<-sqrt(ATEbias^2+ESE^2)
  lldiag<-diag-1.96*sediag
  uldiag<-diag+1.96*sediag
  diagzero<-100*mean((lldiag <= 0 & uldiag >= 0))
  meandiag<-mean(100*diag, na.rm = TRUE)
  results<-c(scen, "&", round(ATEbias,digits=1),"&", round(ASE, digits=2),"&", round(ESE, digits=2),"&", round(SER, digits=2),"&",round(RMSE, digits=1),"&", round(cov, digits=0),"&", round(meandiag,digits=1),"&",round(diagzero, digits=0),'\\')
}

#summarize sim results for all 4 estimators and combine
prep_data_bin <- function(dat,Scen){
  naive.MS<-calc_bias_bin(dat$ATE.naive.MS.bin,dat$seATE.naive.MS.bin,truth$truth_bin[Scen],Scen,'Naive MS',dat$Naive.MS.Diag.bin,dat$seNaive.MS.Diag.bin)
  naive.SS<-calc_bias_bin(dat$ATE.naive.SS.bin,dat$seATE.naive.SS.bin,truth$truth_bin[Scen],Scen,'Naive SS',dat$Naive.MS.Diag.bin,dat$seNaive.MS.Diag.bin)
  MS<-      calc_bias_bin(dat$ATE.MS.bin,dat$seATE.MS.bin,truth$truth_bin[Scen],Scen,'Bridging MS',dat$MS.Diag.bin,dat$seMS.Diag.bin)
  SS<-      calc_bias_bin(dat$ATE.SS.bin,dat$seATE.SS.bin,truth$truth_bin[Scen],Scen,'Bridging SS',dat$MS.Diag.bin,dat$seMS.Diag.bin)
  combine<-rbind(naive.MS,naive.SS,MS,SS)
}

#run the above functions for each sample size and scenario considered, then output summary tables
Sc1<-prep_data(Sc1_pr,1)
Sc2<-prep_data(Sc2_pr,2)
Sc3<-prep_data(Sc3_pr,3)
Sc4<-prep_data(Sc4_pr,4)
Sc5<-prep_data(Sc5_pr,5)
Prim_all_cont<-rbind(Sc1,Sc2,Sc3,Sc4,Sc5)
write.csv(Prim_all_cont, "SimSummaryPrim_cont_013123.csv")

Sc1.bin<-prep_data_bin(Sc1_pr,1)
Sc2.bin<-prep_data_bin(Sc2_pr,2)
Sc3.bin<-prep_data_bin(Sc3_pr,3)
Sc4.bin<-prep_data_bin(Sc4_pr,4)
Sc5.bin<-prep_data_bin(Sc5_pr,5)
Prim_all_bin<-rbind(Sc1.bin,Sc2.bin,Sc3.bin,Sc4.bin,Sc5.bin)
write.csv(Prim_all_bin, "SimSummaryPrim_bin_013123.csv")

Sc1.s<-prep_data(Sc1_alt,1)
Sc2.s<-prep_data(Sc2_alt,2)
Sc3.s<-prep_data(Sc3_alt,3)
Sc4.s<-prep_data(Sc4_alt,4)
Sc5.s<-prep_data(Sc5_alt,5)
Sec_all_cont<-rbind(Sc1.s,Sc2.s,Sc3.s,Sc4.s,Sc5.s)
write.csv(Sec_all_cont, "SimSummarySec_013123.csv")

Sc1.bin.s<-prep_data_bin(Sc1_alt,1)
Sc2.bin.s<-prep_data_bin(Sc2_alt,2)
Sc3.bin.s<-prep_data_bin(Sc3_alt,3)
Sc4.bin.s<-prep_data_bin(Sc4_alt,4)
Sc5.bin.s<-prep_data_bin(Sc5_alt,5)
Sec_all_bin<-rbind(Sc1.bin.s,Sc2.bin.s,Sc3.bin.s,Sc4.bin.s,Sc5.bin.s)
write.csv(Sec_all_bin, "SimSummarySec_bin_013123.csv")

Sc1.s2<-prep_data(Sc1_alt2,1)
Sc2.s2<-prep_data(Sc2_alt2,2)
Sc3.s2<-prep_data(Sc3_alt2,3)
Sc4.s2<-prep_data(Sc4_alt2,4)
Sc5.s2<-prep_data(Sc5_alt2,5)
Sec2_all_cont<-rbind(Sc1.s2,Sc2.s2,Sc3.s2,Sc4.s2,Sc5.s2)
write.csv(Sec2_all_cont, "SimSummarySec2_013123.csv")

Sc1.bin.s2<-prep_data_bin(Sc1_alt2,1)
Sc2.bin.s2<-prep_data_bin(Sc2_alt2,2)
Sc3.bin.s2<-prep_data_bin(Sc3_alt2,3)
Sc4.bin.s2<-prep_data_bin(Sc4_alt2,4)
Sc5.bin.s2<-prep_data_bin(Sc5_alt2,5)
Sec2_all_bin<-rbind(Sc1.bin.s2,Sc2.bin.s2,Sc3.bin.s2,Sc4.bin.s2,Sc5.bin.s2)
write.csv(Sec2_all_bin, "SimSummarySec2_bin_013123.csv")

