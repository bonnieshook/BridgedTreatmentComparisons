#Program: 00_Estimators_05.13.22.R
#Developed by: BES
#Purpose: This program contains functions for fusion M-estimation
#Created: 05.13.22
#Last Updated: 05.13.22

library(numDeriv)
library(boot)
library(geex)

### single-span estimating equations
geex_SS <- function(data, estfun, tr_formula, md_formula, ml_formula, EL.IV, ED.IV, RD.IV){
  
  t_model  <- glm(tr_formula, data = data, family =binomial)
  md_model <- glm(md_formula, data = data, family =binomial)
  ml_model <- glm(ml_formula, data = data, family =binomial)
  models <- list(t=t_model, md=md_model, ml=ml_model)
  
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    #root_control = setup_root_control(start = coef(t_model),coef(md_model),coef(ml_model),EL.IV,ED.IV,RD.IV)),
    roots = c(coef(t_model),coef(md_model),coef(ml_model),EL.IV,ED.IV,RD.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
}

#continuous outcome
estfun_SS <- function(data,models){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$miss==1,0,data$cd4wk8obs)
  
  Xt <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$t))
  T_pos  <- (1):(ncol(Xt))
  T_scores  <- grab_psiFUN(models$t, data)
  
  Xmd <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md))
  Md_pos  <- (ncol(Xt)+1):(ncol(Xt)+ncol(Xmd))
  Md_scores  <- grab_psiFUN(models$md, data)

  Xml <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml))
  Ml_pos  <- (ncol(Xt)+ncol(Xmd)+1):(ncol(Xt)+ncol(Xmd)+ncol(Xml))
  Ml_scores  <- grab_psiFUN(models$ml, data)
  
  function(theta){
    p<-length(theta)
    SLprob <- 0.5 #388 prob triple
    SDprob <- 0.5 #320 prob triple
    SW <-(1-S)*(Triple*SLprob^(-1)+(1-Triple)*(1-SLprob)^(-1))+S*(Triple*SDprob^(-1)+(1-Triple)*(1-SDprob)^(-1)) 
    TRprob<-plogis(Xt %*% theta[T_pos])
    TRW <-S*Dual*(1-TRprob)/TRprob+(1-S)*1+S*(1-Dual)*1
    MDprob<-plogis(Xmd %*% theta[Md_pos])
    MLprob<-plogis(Xml %*% theta[Ml_pos])
    MW<-S*Dual*(MDprob^(-1))+(1-S)*Quad*(MLprob^(-1))+Triple*1
    WTFIN<-SW*TRW*MW
    CML.num<-WTFIN*(1-S)*Quad*NotMiss*Y
    CML.den<-WTFIN*(1-S)*Quad*NotMiss
    CMD.num<-WTFIN*S*Dual*NotMiss*Y
    CMD.den<-WTFIN*S*Dual*NotMiss 
    
    # SS estimators
    c(T_scores(theta[T_pos]),
      Md_scores(theta[Md_pos]),
      Ml_scores(theta[Ml_pos]),
      CML.num-CML.den*theta[p-2],
      CMD.num-CMD.den*theta[p-1],
      theta[p-2]-theta[p-1]-theta[p])
  }
}

#binary outcome
estfun_SS_bin <- function(data,models){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$miss==1,0,data$cd4wk8obs_bin)
  
  Xt <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$t))
  T_pos  <- (1):(ncol(Xt))
  T_scores  <- grab_psiFUN(models$t, data)
  
  Xmd <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md))
  Md_pos  <- (ncol(Xt)+1):(ncol(Xt)+ncol(Xmd))
  Md_scores  <- grab_psiFUN(models$md, data)
  
  Xml <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml))
  Ml_pos  <- (ncol(Xt)+ncol(Xmd)+1):(ncol(Xt)+ncol(Xmd)+ncol(Xml))
  Ml_scores  <- grab_psiFUN(models$ml, data)
  
  function(theta){
    p<-length(theta)
    SLprob <- 0.5 #388 prob triple
    SDprob <- 0.5 #320 prob triple
    SW <-(1-S)*(Triple*SLprob^(-1)+(1-Triple)*(1-SLprob)^(-1))+S*(Triple*SDprob^(-1)+(1-Triple)*(1-SDprob)^(-1)) 
    TRprob<-plogis(Xt %*% theta[T_pos])
    TRW <-S*Dual*(1-TRprob)/TRprob+(1-S)*1+S*(1-Dual)*1
    MDprob<-plogis(Xmd %*% theta[Md_pos])
    MLprob<-plogis(Xml %*% theta[Ml_pos])
    MW<-S*Dual*(MDprob^(-1))+(1-S)*Quad*(MLprob^(-1))+Triple*1
    WTFIN<-SW*TRW*MW
    CML.num<-WTFIN*(1-S)*Quad*NotMiss*Y
    CML.den<-WTFIN*(1-S)*Quad*NotMiss
    CMD.num<-WTFIN*S*Dual*NotMiss*Y
    CMD.den<-WTFIN*S*Dual*NotMiss 
    
    # SS estimators
    c(T_scores(theta[T_pos]),
      Md_scores(theta[Md_pos]),
      Ml_scores(theta[Ml_pos]),
      CML.num-CML.den*theta[p-2],
      CMD.num-CMD.den*theta[p-1],
      theta[p-2]-theta[p-1]-theta[p])
  }
}


#format output from geex
SS_est <- function(data,estfun,tr_formula,md_formula,ml_formula,EL.IV,ED.IV,RD.IV) {
  
  geex_results_ss <- geex_SS(data,estfun,tr_formula,md_formula,ml_formula,EL.IV,ED.IV,RD.IV)
  
  ATE.SS<-geex_results_ss@estimates[length(geex_results_ss@estimates)]
  seATE.SS <- sqrt(geex_results_ss@vcov[length(geex_results_ss@estimates),length(geex_results_ss@estimates)])
  ATE_ests.SS<-cbind(ATE.SS,seATE.SS)
  
  return(ATE_ests.SS)
}


### multi-span estimating equations
geex_MS <- function(data, estfun, tr_formula, md1_formula, md2_formula, ml2_formula, ml3_formula, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV){
  
  t_model  <- glm(tr_formula, data = data, family =binomial)
  md1_model <- glm(md1_formula, data = data, family =binomial)
  md2_model <- glm(md2_formula, data = data, family =binomial)
  ml2_model <- glm(ml2_formula, data = data, family =binomial)
  ml3_model <- glm(ml3_formula, data = data, family =binomial)
  models <- list(t=t_model, md1=md1_model,md2=md2_model, ml2=ml2_model, ml3=ml3_model)
  
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    #root_control = setup_root_control(start = c(coef(t_model),coef(md1_model),coef(md2_model),coef(ml2_model),coef(ml3_model),EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV)),
    roots = c(coef(t_model),coef(md1_model),coef(md2_model),coef(ml2_model),coef(ml3_model),EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
}

#continuous outcome
estfun_MS <- function(data,models){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$nmissout==1,data$cd4wk8obs,0)
  
  Xt <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$t))
  T_pos  <- (1):(ncol(Xt))
  T_scores  <- grab_psiFUN(models$t, data)
  
  Xmd1 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md1))
  Md1_pos  <- (ncol(Xt)+1):(ncol(Xt)+ncol(Xmd1))
  Md1_scores  <- grab_psiFUN(models$md1, data)

  Xmd2 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md2))
  Md2_pos  <- (ncol(Xt)+ncol(Xmd1)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2))
  Md2_scores  <- grab_psiFUN(models$md2, data)
  
  Xml2 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml2))
  Ml2_pos  <- (ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2))
  Ml2_scores  <- grab_psiFUN(models$ml2, data)
 
  Xml3 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml3))
  Ml3_pos  <- (ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2)+ncol(Xml3))
  Ml3_scores  <- grab_psiFUN(models$ml3, data)
   
  function(theta){
    p<-length(theta)
    SLprob <- 0.5 #388 prob triple
    SDprob <- 0.5 #320 prob triple
    SW <-(1-S)*(Triple*SLprob^(-1)+(1-Triple)*(1-SLprob)^(-1))+S*(Triple*SDprob^(-1)+(1-Triple)*(1-SDprob)^(-1)) 
    TRprob<-plogis(Xt %*% theta[T_pos])
    TRW <-S*(1-TRprob)/TRprob+(1-S)*1
    MD1prob<-plogis(Xmd1 %*% theta[Md1_pos])
    MD2prob<-plogis(Xmd2 %*% theta[Md2_pos])
    ML2prob<-plogis(Xml2 %*% theta[Ml2_pos])
    ML3prob<-plogis(Xml3 %*% theta[Ml3_pos])
    MW<-S*Dual*(MD1prob^(-1))+S*Triple*(MD2prob^(-1))+(1-S)*Triple*(ML2prob^(-1))+(1-S)*Quad*(ML3prob^(-1))
    WTFIN<-SW*TRW*MW
    CML3.num<-WTFIN*(1-S)*Quad*NotMiss*Y
    CML3.den<-WTFIN*(1-S)*Quad*NotMiss
    CML2.num<-WTFIN*(1-S)*Triple*NotMiss*Y
    CML2.den<-WTFIN*(1-S)*Triple*NotMiss
    CMD2.num<-WTFIN*S*Triple*NotMiss*Y
    CMD2.den<-WTFIN*S*Triple*NotMiss 
    CMD1.num<-WTFIN*S*Dual*NotMiss*Y
    CMD1.den<-WTFIN*S*Dual*NotMiss     
    # MS estimators
    c(T_scores(theta[T_pos]),
      Md1_scores(theta[Md1_pos]),
      Md2_scores(theta[Md2_pos]),
      Ml2_scores(theta[Ml2_pos]),
      Ml3_scores(theta[Ml3_pos]),
      CML3.num-CML3.den*theta[p-5],
      CML2.num-CML2.den*theta[p-4],
      CMD2.num-CMD2.den*theta[p-3],
      CMD1.num-CMD1.den*theta[p-2],
      theta[p-5]-theta[p-4]+theta[p-3]-theta[p-2]-theta[p-1],
      theta[p-4]-theta[p-3]-theta[p])
  }
}

#binary outcome
estfun_MS_bin <- function(data,models){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$nmissout==1,data$cd4wk8obs_bin,0)
  
  Xt <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$t))
  T_pos  <- (1):(ncol(Xt))
  T_scores  <- grab_psiFUN(models$t, data)
  
  Xmd1 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md1))
  Md1_pos  <- (ncol(Xt)+1):(ncol(Xt)+ncol(Xmd1))
  Md1_scores  <- grab_psiFUN(models$md1, data)
  
  Xmd2 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$md2))
  Md2_pos  <- (ncol(Xt)+ncol(Xmd1)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2))
  Md2_scores  <- grab_psiFUN(models$md2, data)
  
  Xml2 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml2))
  Ml2_pos  <- (ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2))
  Ml2_scores  <- grab_psiFUN(models$ml2, data)
  
  Xml3 <- grab_design_matrix(data = data,rhs_formula = grab_fixed_formula(models$ml3))
  Ml3_pos  <- (ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2)+1):(ncol(Xt)+ncol(Xmd1)+ncol(Xmd2)+ncol(Xml2)+ncol(Xml3))
  Ml3_scores  <- grab_psiFUN(models$ml3, data)
  
  function(theta){
    p<-length(theta)
    SLprob <- 0.5 #388 prob triple
    SDprob <- 0.5 #320 prob triple
    SW <-(1-S)*(Triple*SLprob^(-1)+(1-Triple)*(1-SLprob)^(-1))+S*(Triple*SDprob^(-1)+(1-Triple)*(1-SDprob)^(-1)) 
    TRprob<-plogis(Xt %*% theta[T_pos])
    TRW <-S*(1-TRprob)/TRprob+(1-S)*1
    MD1prob<-plogis(Xmd1 %*% theta[Md1_pos])
    MD2prob<-plogis(Xmd2 %*% theta[Md2_pos])
    ML2prob<-plogis(Xml2 %*% theta[Ml2_pos])
    ML3prob<-plogis(Xml3 %*% theta[Ml3_pos])
    MW<-S*Dual*(MD1prob^(-1))+S*Triple*(MD2prob^(-1))+(1-S)*Triple*(ML2prob^(-1))+(1-S)*Quad*(ML3prob^(-1))
    WTFIN<-SW*TRW*MW
    CML3.num<-WTFIN*(1-S)*Quad*NotMiss*Y
    CML3.den<-WTFIN*(1-S)*Quad*NotMiss
    CML2.num<-WTFIN*(1-S)*Triple*NotMiss*Y
    CML2.den<-WTFIN*(1-S)*Triple*NotMiss
    CMD2.num<-WTFIN*S*Triple*NotMiss*Y
    CMD2.den<-WTFIN*S*Triple*NotMiss 
    CMD1.num<-WTFIN*S*Dual*NotMiss*Y
    CMD1.den<-WTFIN*S*Dual*NotMiss     
    # MS estimators
    c(T_scores(theta[T_pos]),
      Md1_scores(theta[Md1_pos]),
      Md2_scores(theta[Md2_pos]),
      Ml2_scores(theta[Ml2_pos]),
      Ml3_scores(theta[Ml3_pos]),
      CML3.num-CML3.den*theta[p-5],
      CML2.num-CML2.den*theta[p-4],
      CMD2.num-CMD2.den*theta[p-3],
      CMD1.num-CMD1.den*theta[p-2],
      theta[p-5]-theta[p-4]+theta[p-3]-theta[p-2]-theta[p-1],
      theta[p-4]-theta[p-3]-theta[p])
  }
}

#format output from geex
MS_est <- function(data, estfun, tr_formula, md1_formula, md2_formula, ml2_formula, ml3_formula, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV) {
  
  geex_results_ms <- geex_MS(data, estfun, tr_formula, md1_formula, md2_formula, ml2_formula, ml3_formula, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV)
  
  ATE.MS<-geex_results_ms@estimates[length(geex_results_ms@estimates)-1]
  seATE.MS <- sqrt(geex_results_ms@vcov[length(geex_results_ms@estimates)-1,length(geex_results_ms@estimates)-1])
  ATE_ests.MS<-cbind()
  
  Diag<-geex_results_ms@estimates[length(geex_results_ms@estimates)]
  seDiag <- sqrt(geex_results_ms@vcov[length(geex_results_ms@estimates),length(geex_results_ms@estimates)])

  MS_ests<-cbind(ATE.MS,seATE.MS,Diag,seDiag)
  
  return(MS_ests)
}



### naive single-span estimating equations
geex_naive_SS <- function(data, estfun, EL.IV, ED.IV, RD.IV){
  
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    #root_control = setup_root_control(start = c(EL.IV,ED.IV,RD.IV))
    roots = c(EL.IV,ED.IV,RD.IV),
    compute_roots = FALSE
    )
}

#continuous outcome
estfun_naive_SS <- function(data){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$miss==1,0,data$cd4wk8obs)
 
  function(theta){
    p<-length(theta)
    CML.num<-(1-S)*Quad*NotMiss*Y
    CML.den<-(1-S)*Quad*NotMiss
    CMD.num<-S*Dual*NotMiss*Y
    CMD.den<-S*Dual*NotMiss 
    
    # SS estimators
    c(CML.num-CML.den*theta[p-2],
      CMD.num-CMD.den*theta[p-1],
      theta[p-2]-theta[p-1]-theta[p])
  }
}

#binary outcome
estfun_naive_SS_bin <- function(data){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$miss==1,0,data$cd4wk8obs_bin)
  
  function(theta){
    p<-length(theta)
    CML.num<-(1-S)*Quad*NotMiss*Y
    CML.den<-(1-S)*Quad*NotMiss
    CMD.num<-S*Dual*NotMiss*Y
    CMD.den<-S*Dual*NotMiss 
    
    # SS estimators
    c(CML.num-CML.den*theta[p-2],
      CMD.num-CMD.den*theta[p-1],
      theta[p-2]-theta[p-1]-theta[p])
  }
}

#format output from geex
SS_naive_est <- function(data,estfun,EL.IV,ED.IV,RD.IV) {
  
  geex_results_ss <- geex_naive_SS(data,estfun,EL.IV,ED.IV,RD.IV)
  
  ATE.naive.SS<-geex_results_ss@estimates[length(geex_results_ss@estimates)]
  seATE.naive.SS <- sqrt(geex_results_ss@vcov[length(geex_results_ss@estimates),length(geex_results_ss@estimates)])
  ATE_ests.naive.SS<-cbind(ATE.naive.SS,seATE.naive.SS)
  
  return(ATE_ests.naive.SS)
}



### naive multi-span estimating equations
geex_naive_MS <- function(data, estfun, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV){
  
  m_estimate(
    estFUN = estfun, 
    data   = data, 
    #root_control = setup_root_control(start = c(EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV))
    roots = c(EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV),
    compute_roots = FALSE,
    )
}

#continuous outcome
estfun_naive_MS <- function(data){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$nmissout==1,data$cd4wk8obs,0)
 
  function(theta){
    p<-length(theta)
    CML3.num<-(1-S)*Quad*NotMiss*Y
    CML3.den<-(1-S)*Quad*NotMiss
    CML2.num<-(1-S)*Triple*NotMiss*Y
    CML2.den<-(1-S)*Triple*NotMiss
    CMD2.num<-S*Triple*NotMiss*Y
    CMD2.den<-S*Triple*NotMiss 
    CMD1.num<-S*Dual*NotMiss*Y
    CMD1.den<-S*Dual*NotMiss     
    # MS estimators
    c(CML3.num-CML3.den*theta[p-5],
      CML2.num-CML2.den*theta[p-4],
      CMD2.num-CMD2.den*theta[p-3],
      CMD1.num-CMD1.den*theta[p-2],
      theta[p-5]-theta[p-4]+theta[p-3]-theta[p-2]-theta[p-1],
      theta[p-4]-theta[p-3]-theta[p])
  }
}

#binary outcome
estfun_naive_MS_bin <- function(data){
  
  S<-data$ACTG320
  Dual<-data$A==1
  Triple<-data$A==2
  Quad<-data$A==3
  NotMiss<-data$nmissout
  Y<-ifelse(data$nmissout==1,data$cd4wk8obs_bin,0)
  
  function(theta){
    p<-length(theta)
    CML3.num<-(1-S)*Quad*NotMiss*Y
    CML3.den<-(1-S)*Quad*NotMiss
    CML2.num<-(1-S)*Triple*NotMiss*Y
    CML2.den<-(1-S)*Triple*NotMiss
    CMD2.num<-S*Triple*NotMiss*Y
    CMD2.den<-S*Triple*NotMiss 
    CMD1.num<-S*Dual*NotMiss*Y
    CMD1.den<-S*Dual*NotMiss     
    # MS estimators
    c(CML3.num-CML3.den*theta[p-5],
      CML2.num-CML2.den*theta[p-4],
      CMD2.num-CMD2.den*theta[p-3],
      CMD1.num-CMD1.den*theta[p-2],
      theta[p-5]-theta[p-4]+theta[p-3]-theta[p-2]-theta[p-1],
      theta[p-4]-theta[p-3]-theta[p])
  }
}


#format output from geex
MS_naive_est <- function(data, estfun, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV) {
  
  geex_results_ms <- geex_naive_MS(data, estfun, EL3.IV,EL2.IV, ED2.IV, ED1.IV, RD.IV, Diag.IV)
  
  ATE.MS<-geex_results_ms@estimates[length(geex_results_ms@estimates)-1]
  seATE.MS <- sqrt(geex_results_ms@vcov[length(geex_results_ms@estimates)-1,length(geex_results_ms@estimates)-1])
  
  Diag<-geex_results_ms@estimates[length(geex_results_ms@estimates)]
  seDiag <- sqrt(geex_results_ms@vcov[length(geex_results_ms@estimates),length(geex_results_ms@estimates)])
  
  MS_ests<-cbind(ATE.MS,seATE.MS,Diag,seDiag)
  
  return(MS_ests)
}