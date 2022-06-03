#' Get bias estimate
#'
#' @description Calculation of performance measure of estimation: Bias
#' as the difference between the estimated and actual values
#'
#' @param vec Estimated value from simulation
#' @param theta True value
#'
#' @return Bias of effect estimation
#' @export
#'
bias<-function(vec, theta){
  return(mean(vec)-theta)
}



#' Get several performance measures of (effect estimation) simulation studies
#'
#' @description Calculation of several performance measure for the estimates
#' obtained from the simulation repeats
#'
#' @param vecEst vector of estimated effect values from simulation repeats
#' @param vecSe vector of estimated standard error from simulation repeats
#' @param vec.p vector of significance (probability of being sampled from null hypothesis) from simulation repeats
#' @param theta true value for effect
#' @param alpha Type I error rate = significance level for testing the null hypothesis
#'
#' @return Vector of Performance measures:
#' 1. Average of the estimates of the parameter of interest
#' 2. Empirical standard error as an assessment of the estimation uncertainty
#' 3. Bias = mean deviation of the estimates from the true value of the parameter of interest is an indicator of accuracy
#' 4. Coverage of a confidence interval is a measurement that can be used to control the Type I error
#' 5. Empirical power = the proportion of simulation samples in which the H0 of no effect is rejected at a significance level when H0 is false
#' @export
#'
performanceMeas<-function(vecEst, vecSe, vec.p, theta, alpha=0.05){

  z<-stats::qnorm(p=1-alpha/2)
  CI<-data.frame(CI_low=vecEst-z*vecSe, CI_up=vecEst+z*vecSe)

  res<-c(
    AvEst=mean(vecEst),
    EmpSE=stats::sd(vecEst),
    Bias=bias(vecEst,theta),
    Coverage=mean(as.numeric((CI$CI_low<=theta)&(theta<=CI$CI_up))),
    Power=mean(as.numeric(vec.p<=alpha))
  )

  return(res)

}


#' Data loss simulation for given study data
#'
#' @description Simulate a loss of data points within a full study data obtained
#' within I clusters and observed over K time points
#'
#' @param data.all given data of a study (can be simulated or obtained)
#' @param I Number of clusters within the study
#' @param K Number of time points within the study
#' @param B Type of cluster loss: 0 = No cluster loss, 1 = Cluster missing at random, 2 = Cluster is missing at beginning, 3 = Cluster is missing at end of the study
#' @param C Number of clusters which getting lost
#' @param D Number of individual loss (chosen randomly from all individuals within the study)
#'
#' @return Matrix of data which is remained after data loss
#' @export
#'
Data.loss.SWD<-function(data.all, I,K, B="0", C=0, D=0){

  ## data deletion depends on:  B,C, D      ##

  if(C==0){data.delete<-data.all
  }else{

    #number of cluster loss
    clusternr<-sample(1:I, size=C, replace=FALSE)

    #timepoint of cluster loss
    switch(B,
           #B0: no Cluster missing
           "0"={timepoint<-0},
           #B1: Cluster missing at random
           "1"={timepoint<-sample(2:K, size=C, replace=TRUE)  },
           #B2: Cluster is missing at beginning
           "2"={timepoint<-sample((1:round(K/4))+1, size=C, replace=TRUE)  },
           #B3: Cluster is missing at end
           "3"={timepoint<-sample((K-round(K/4)+1):K, size=C, replace=TRUE)  }
    )

    dc<-NULL
    for(i in 1:C){

      which.tp<-timepoint[i]:K
      dc<-rbind(dc, cbind(rep(clusternr[i],length(which.tp)), timepoint[i]:K))
    }
    colnames(dc)<-c("cluster", "measurement")

    # given vector of  timepoints and corresponding vector of cluster
    # collect rows whcih has to delete
    collect.all<-NULL
    for(i in 1:dim(dc)[1]){
      collect<-unlist(sapply(1:dim(data.all)[1], function(j){

        if(identical(unname(unlist(data.all[j,c("cluster", "measurement")]),force=TRUE),unname(unlist(dc[i,]),force=TRUE))){return(j)}
      }))
      collect.all<-c(collect.all,collect )
    }
    # now delete rows in data
    data.delete<-data.all[-collect.all, ]

  }
  if(D!=0){#individiual loss

    loss.ind<-sample(1:dim(data.delete)[1], size=D, replace=FALSE)
    data.delete<-data.delete[-loss.ind, ]
  }

  return(data.delete)
}


#' Get performance measures from effect estimation within all
#' repeats of the simulation
#'
#' @description Simulation for ??????
#' Repeats anzSim times the following steps
#'   1. Determining the design matrix regarding the chosen design (Design, Number of cluster, time points, individuals)
#'      knowing Fidelity pattern and possible data loss
#'   2. Sample data of the given study design and expected data loss using the package "samplingDataCRT"
#'   3. Estimation of Effects using linear mixed model estimation
#'
#' @param anzSim Number of simulation repeats
#' @param type Study design type = "cross-sec" for cross-sectional or  "long" for longitudinal
#' @param sigma.1 Within variability or error variance (Hughes&Hussey sigma)
#' @param sigma.2  Between individual variance (only used when the individuals are followed over time: longitudinal study)
#' @param sigma.3 Between clusters variability (Hughes&Hussey tau)
#' @param K Number of time points (measurement)
#' @param I Number of cluster
#' @param J Number of Individuals (=individuals per cluster)
#' @param mu.0 Baseline mean within the model specification
#' @param theta Intervention effect
#' @param betas Time trend could be included
#' @param X Design matrix for building linear model of given ooptimal study
#' @param X.A Design matrix for building linear model of sampling real data given real setting
#' @param B.cond Condition for type of cluster loss: 0 = No cluster loss, 1 = Cluster missing at random, 2 = Cluster is missing at beginning, 3 = Cluster is missing at end of the study
#' @param C.cond Condition for number of clusters which getting lost
#' @param D.cond Condition for Number of individual loss (chosen randomly from all individuals within the study)
#'
#' @seealso  Reference for linear model formulation of stepped wedge cluster randomized trials
#' by Hussey and Hughes "Design and analysis of stepped wedge cluster randomized trials",
#' Contemporary Clinical Trials, 2007 \url{https://www.sciencedirect.com/science/article/pii/S1551714406000632}
#' as well as performance measurements for simulation studies
#' by Burton et al. “The design of simulation studies in medical statistics”,
#' Stat Med, 2006 \url{https://pubmed.ncbi.nlm.nih.gov/16947139/}
#'
#' @return Vector of Estimates summaru=ies and performance measures (see function \code{performanceMeas}):
#' 1. Average of the estimates of the parameter of interest
#' 2. Empirical standard error as an assessment of the estimation uncertainty
#' 3. Bias = mean deviation of the estimates from the true value of the parameter of interest is an indicator of accuracy
#' 4. Coverage of a confidence interval is a measurement that can be used to control the Type I error
#' 5. Empirical power = the proportion of simulation samples in which the H0 of no effect is rejected at a significance level when H0 is false
#'
#' @export
simulation<-function(anzSim,type, sigma.1,sigma.2=NULL,sigma.3,K,J,I,mu.0,theta,betas,
                        X, X.A, B.cond="0",C.cond=0, D.cond=0){

  # Desingmatrix Daten laut Studeinedesign
  D<-samplingDataCRT::completeDataDesignMatrix(J, X)
  # Desingmatrix Daten f?r reale Daten
  if(!is.null(X.A)){A<-samplingDataCRT::completeDataDesignMatrix(J, X.A)
  }else{A<-NULL}


  parameters<-c(mu.0, betas, theta)

  #Covarianzmatrix
  if(type=="cross-sec"){
    V<-samplingDataCRT::CovMat.Design(K=K, J=J, I=I, sigma.1=sigma.1, sigma.3=sigma.3)
  }
  if(type=="long"){
    V<-samplingDataCRT::CovMat.Design(K=K, J=J, I=I,sigma.1, sigma.2, sigma.3)
  }


  res.all<-NULL
  for(s in 1:anzSim){ ##Repeats simulation

    # sample Data
    # sample I cluster from distribution of cluster, each ave the same mean vector
    sample.data<-samplingDataCRT::sampleData(type, K,J,I, D, A, V, parameters )
    #Loss of Data
    data.delete<-Data.loss.SWD(sample.data, I,K,B.cond,C.cond, D.cond)
    # print(xtabs(~cluster+measurement, data=data.delete))

    ####    Estimation by linear mixed model ##
    # timepoint 1 <- 0
    data.delete$measurement<-as.factor((as.numeric(data.delete$measurement)-1))
    #cross-sec
    if(type=="cross-sec"){

      #lm.res<-lmer(val~ intervention + measurement+(1|cluster)+ (1|subject), data=sample.data)
      lm.res<-lme4::lmer(val~ intervention + measurement+(1|cluster), data=data.delete)
      lm.0<-lme4::lmer(val~ measurement+(1|cluster), data=data.delete)
    }

    #longitudinal
    if(type=="long"){

      lm.res<-lme4::lmer(val~ intervention + measurement+(1|cluster)+ (1|subject), data=data.delete)
      #summary(lm.res)
      lm.0<-lme4::lmer(val~ measurement+(1|cluster)+ (1|subject), data=data.delete)
      #summary(lm.res)
    }


    #random effects
    res.ran<-as.data.frame(summary(lm.res)$varcor)[,5]^2
    names(res.ran)<-as.data.frame(summary(lm.res)$varcor)[,1]
    #fixed efects
    res.fix<-lme4::fixef(lm.res)[-1]
    #SE of fixed effect estimates
    SEs<-stats::coef(summary(lm.res))[-1,"Std. Error"]
    names(SEs)<-paste("SE.",names(SEs), sep="")
    #anova for intervention
    anova.p<-stats::anova(lm.res,lm.0)[8][2,1]

    res.all<-rbind(res.all, c(res.ran,res.fix, SEs, p.intervention=anova.p))
  }

  #Summary of all repeats
  summEst<-colMeans(res.all)
  names(summEst)<-paste(names(summEst), "Mean",".")
  perf.intervention<-performanceMeas(res.all[,"intervention"], res.all[,"SE.intervention"], res.all[,"p.intervention"], theta)
  names(perf.intervention)<-paste(names(perf.intervention),"Intervention",sep=".")
  #c(performanceMeas(est["intervention",], est["SE.intervention",], est["p.intervention",], theta))
  return(c(summEst, perf.intervention))
}


#' Implementation of parallel matrix
#'
#' @description Given a design of a hypothetical cluster parallel randomized study
#' and their pattern of fidelity over time
#' conduct the function the corresponding design matrix
#' which can be used for linear regression modeling
#'
#' @param nC Number of cluster
#' @param nT Number of time points
#' @param nSw Number of cluster with control condition 0 (nC-nSw is then the number of cluster with treatment condition 1)
#' @param pattern Fidelity pattern (deviation from 1)  for each time point
#'
#' @return Design matrix of size nT x nC for a parrallel cluster randomized trial with
#' values of 0 for control condition and values between 0 and 1 for degree of implementation of treatment (fidelity)
#' @export
#'
implemMatrix.parallel<-function (nC, nT, nSw, pattern){
  if (length(pattern)  == nT) {
    ma <- rbind(t(replicate(nSw, rep(0, nT))), t(replicate(nC -nSw, pattern)))
    return(ma)
  }
  else {
    stop("The length of the pattern must be the same as the number of timepoints.")
  }
}



#' Implementation of crossover matrix
#'
#' @description Given a design of a hypothetical cluster randomized stepped wedge study and their pattern of fidelity over time
#' conduct the function the corresponding design matrix
#' which can be used for linear regression modeling
#'
#' @param nC Number of cluster
#' @param nT Number of time points
#' @param swP Number of timepoints the for clusters be at control
#' @param nSw Number of cluster with switches from control condition 0 to treatment condition 1 at each time point
#' @param pattern Fidelity pattern (deviation from 1)  for each time point
#'
#' @return Design matrix of size nT x nC for a cross-over cluster randomized trial with
#' values of 0 for control condition and values between 0 and 1 for degree of implementation of treatment (fidelity)' @return
#' @export
#
implemMatrix.crossover<-function (nC, nT, nSw,swP=NULL, pattern){
  if (is.null(swP)) {
    swP <- ceiling(nT/2)
  }
  if (length(pattern)  == swP) {

    ma<- rbind(t(replicate(nSw, c(rep(0, swP),  pattern))),
               t(replicate(nC - nSw, c(pattern, rep(0, nT - swP))))
               )
    return(ma)
  }
  else {
    stop("time points of intervention must be the same as the length of the pattern")
  }
}

