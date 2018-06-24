#Npat is the number of patient
#Ncov is the number of covariates here, which is fixed at 4 for now(including AUC^2)
#Comp_prob_0 and Comp_prob_1 are the probability for different CR
#CR is complete remission status
#mus and sds are the means and sds for AUC
#sigma is the sd for the survival time
#beta is the linear coefficients for lifetime mean
#censormean is the "mean" for censor time
#prob is the probability used when simulate life time

#We have four covaraites for generating the truth: Age, CR, AUC and AUC^2

#method=1: two mixture of lognormals
#method=2: two mixture of weibulls

#mixture equals to 1 means 2 component
#mixture equals to 0 means only one component, i.e. no mixture


#' Simulate data
#' @param Npat A positive integer indicating the number of patients to be simulated.
#' 
#' @param method A binary value (either 0 or 1). If method is set to 0, the survival times are simulated from a Weibull distribution.
#' If method is set to 1, the survival times are simulated from a Lognormal distribution. Default value is 1. 
#' 
#' @param mixture A binary value (either 0 or 1). If method is set to 0, the survival times are simulated from a mixture of two 
#' distributions. If method is set to 1, the survival times are simulated from single distribution. Default value is 1. 
#' 
#' @param censor A binary value (either 0 or 1). If censor is set to 1, survival times of patients will be censored (at roughly 25%).
#' If censor is set to 0, the survival times of patients will not be censored.
#'  Default value is 0. 
#' 
#' @return Returns a list composed of the patients' generated attributes. The covariates are Age, AUC, and CR.
#' 
#' @examples 
#' #Simulate Data for 10 patients
#' Npat<- 10 
#' data<-simulate_data(Npat)
#' 
#' @export
simulate_data <- function(Npat, method=1, mixture=1, censor=0)
{  
  Comp_prob=c(0.55,0.35,0.1)
  mus=c(4.5,5.5,7.25);sds=sqrt(c(1,1,1)); sigma = 0.4; 
  shape = 1.5; censcale = 650 #parameters used for weibull distribution
  censormean = log(650)
  
  if(mixture) prob = 0.4 else
    prob = 0
  
  Data <- NULL
  Data$Npat <- Npat
  #Simulate Age from discrete uniform
  # Age <- Data$Age <- sample(c(16:65), Npat, replace=TRUE)
  
  #Simulate CR
  CR <- Data$CR <- rbinom(Npat,1,0.5)
  
  #Simulate AUC
  # AUC <- numeric(Npat)
  # for (i in 1:Npat)
  # {
  #     components <- sample(1:length(mus),prob=Comp_prob,size=Npat,replace=TRUE)
  #     AUC[i]=rnorm(1,mean=mus[components],sd=sds[components])
  # } 
  # Data$AUC <- AUC
  True_Age<-c(59, 46 ,52 ,56 ,37 ,54 ,53 ,27 ,39 ,60 ,49 ,58 ,57 ,35 ,48 ,46 ,50 ,36 ,49 ,36 ,65 ,16 ,53 ,50 ,27 ,55 ,37 ,55 ,13 ,57 ,51,
              36 ,34 ,59 ,40 ,50 ,56 ,53 ,53 ,51 ,31 ,33 ,44 ,61 ,51 ,38 ,57 ,64 ,55 ,50 ,36 ,52 ,42 ,28 ,24 ,53 ,56 ,60 ,58 ,55 ,42 ,61,
              43 ,54 ,36 ,46 ,61 ,52 ,23 ,55 ,53 ,55 ,61 ,24 ,59 ,62 ,42 ,52 ,32 ,62 ,37 ,58 ,22 ,44 ,55 ,63 ,56 ,55 ,50 ,48 ,32 ,52 ,65,
              38 ,13 ,28 ,48 ,35 ,53 ,52 ,24 ,32 ,65 ,51 ,31 ,50 ,47 ,44 ,39 ,21 ,20 ,48 ,27 ,45 ,22 ,51 ,54 ,62 ,49 ,23 ,35 ,57 ,42 ,45,
              54 ,44 ,58 ,58 ,44 ,53 ,38 ,43 ,61 ,35 ,45 ,43 ,26 ,44 ,32 ,59 ,60 ,41 ,56 ,57 ,36 ,50 ,53 ,31 ,27 ,27 ,29)
  remove<-c(13,16)
  True_Age<-True_Age[! True_Age %in% remove]
  add<-seq(25,35)
  True_Age<-c(True_Age,add,add)
  True_AUC<-c(4.218909 , 2.621726 ,4.274692 ,4.600881 ,5.114714 ,5.506069 ,5.232825 ,5.309927 ,6.875603 ,5.417502
              ,4.669992 ,4.845335 ,5.816077 ,7.507137 ,6.371383 ,5.777296 ,4.390200 ,5.258493 ,7.324009 ,5.463360
              ,5.626651 ,4.588791 ,5.902899 ,5.548625 ,4.630019 ,5.959815 ,4.310302 ,5.466724 ,6.012101 ,5.240228
              ,5.327913 ,6.161089 ,5.204301 ,3.914655 ,4.885487 ,5.923654 ,5.158748 ,4.760329 ,4.412944 ,6.405455
              ,4.118528 ,5.240039 ,3.556275 ,5.763631 ,7.208160 ,4.344777 ,4.269055 ,4.392512 ,3.609402 ,7.228775
              ,4.449121 ,5.367043 ,5.869129 ,5.898398 ,4.231518 ,4.743104 ,5.169906 ,6.275440 ,7.359592 ,6.656494
              ,5.595894 ,5.889592 ,3.812498 ,4.158808 ,5.286502 ,4.547545 ,4.040851 ,6.152578 ,5.703770 ,4.567299
              ,5.949434 ,5.927006 ,6.182740 ,4.504007 ,3.541518 ,4.355357 ,5.444196 ,5.474398 ,5.804139 ,4.888809
              ,5.412979 ,5.409267 ,4.373623 ,4.640287 ,6.705499 ,5.671171 ,4.638472 ,4.277929 ,4.535522 ,4.714513
              ,5.485711 ,4.285346 ,6.070886 ,7.058000 ,6.770000 ,6.909000 ,5.369000 ,4.266000 ,4.995000 ,5.290000
              ,4.668000 ,3.449000 ,4.694000 ,3.976000 ,3.435000 ,3.535000 ,4.908000 ,3.936000 ,3.881000 ,3.995000
              ,4.122000 ,4.472000 ,5.830000 ,5.594000 ,4.956000 ,5.797000 ,4.151000 ,4.127000 ,4.013000 ,5.744000
              ,4.158000 ,5.615000 ,4.943000 ,5.092000 ,4.181000 ,4.613000 ,3.912000 ,5.077000 ,6.111000 ,5.733000
              ,3.598000 ,5.764000 ,4.906000 ,3.928000 ,6.550000 ,4.599000 ,5.043000 ,7.151000 ,3.385000 ,4.416000
              ,5.340000 ,3.621000 ,4.615000 ,5.255000 ,4.966000 ,5.219000 ,4.168000 ,4.045000 ,5.169000 ,8.259000
              ,2.931000)
  Age <- Data$Age <- sample(True_Age,Npat,replace=TRUE)
  AUC<- Data$AUC <-sample(True_AUC,Npat,replace=TRUE)
  #Simulate
  lifetimes<-numeric(Npat)
  Ncov = 7
  dpat <- matrix(0, Npat, Ncov)
  dpat[,1]=scale(Age)
  dpat[,2] = AUC
  dpat[,3] = CR
  dpat[,4] = AUC^2
  dpat[,5] = AUC*scale(Age)
  dpat[,6] = AUC*CR
  # dpat[,7] = log(1+exp(50+(AUC*(scale(Age)-0.3)^3*2.1*abs((0.5-CR)))))
  dpat[,7] = AUC*scale(Age)*abs(CR)
  if(method ==1)
  {
    beta1<-numeric(Ncov+1)
    #For log exp as 7th covariate 
    # beta2<-c(4.9,-0.10, 0.7,  0.1, -0.084,-0.07, 0.27,-0.05)
    # beta2<-c(5.5,-0.1, 0.7,  0.1, -0.09,-0.07, 0.15,-0.05)
    
    #Low impact 3 way inter
    # beta2<-c(3.5,-0.10, 0.7,  0.1, -0.075,-0.07, 0.06,-0.1)
    #3 way inter 
    # beta2<-c(5,-0.10, 0.7,  0.1, -0.075,-0.07, 0.06,-0.3)
    #3 way inter a
    # beta2<-c(3.5,-0.10, 0.7,  0.1, -0.075,-0.07, 0.5,-0.3)
    beta2<-c(4,-0.10, 0.7,  0.3, -0.1,-0.068, 0.2,-0.18)
    
    beta1<-beta2
    
  } else
  {
    beta2<-c(5, -0.05, 0.05,  0.15, -0.15)
    beta1<-c(3.1, -0.5,0.45,0.4, -0.2)
  }
  
  lifetimes <- rep(0, Npat)
  
  if(method==1){
    Data$truth='lognorm'
    censtimes<-numeric(Npat)
    for (i in 1:Npat){
      pat_cov<- c(1,dpat[i,])
      if (runif(1)<prob)
      {
        meanlog<- sum(beta1*pat_cov)
        lifetimes[i]<-rlnorm(1,meanlog=meanlog,sdlog=sigma)
      }
      else{
        meanlog<- sum(beta2*pat_cov)
        lifetimes[i]<-rlnorm(1,meanlog=meanlog,sdlog=sigma)
      }
      censtimes[i]<-rlnorm(1,meanlog=meanlog+0.41,sdlog=sigma)
    }
    # censtimes <-rlnorm(Npat,meanlog=censormean,sdlog=sigma)
    
  }else{#method==2
    Data$truth='weib'
    for (i in 1:Npat){
      pat_cov<- c(1,dpat[i,])
      u<-runif(1)
      if (u<prob){
        scale<- exp(sum(beta1*pat_cov))
        lifetimes[i]<-rweibull(1,shape=shape,scale=scale)
      }
      else{scale<- exp(sum(beta2*pat_cov))
      lifetimes[i]<-rweibull(1,shape=shape,scale=scale)
      }
    }
    censtimes <-rweibull(Npat,shape=shape,scale=censcale)
  }
  
  
  if (censor==0)
  {
    life_cens = lifetimes
    death = rep(1, Npat)
  } else
  {
    life_cens<-pmin(lifetimes,censtimes)
    death <- (lifetimes < censtimes) #death = 1 means observe death
  }
  
  Data$OS <- life_cens
  Data$death <- death
  Data$beta1<-beta1
  Data$beta2<-beta2
  Data$prob<-prob
  Data$sigma<-sigma
  Data$dpat<-dpat
  # plot(density(log(Data$OS)))
  # sum(Data$death)
  return(Data)
}

# bad<-which(Data$OS==min(Data$OS))
# c(1,dpat[bad,])*beta2
