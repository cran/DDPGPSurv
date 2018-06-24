#' Posterior mean survival estimation using Gibbs MCMC Outputs
#' @param mcmc A list-The MCMC output from mcmc_Gibbs.
#' @param new_pat A matrix consisting of the covariates of new sample points. Each row represents one sample point.
#' It can also be vector consisting of the covariates of a new sample point. 
#' @param if_plot A logical variable indicating whether the mean survival estimation should be plotted.
#' 1 means estimation should be plotted; 0 means otherwise. Defualt value is 0, i.e. not plotted.
#' @param quantiles A two-dimensional vector determining the quantiles for the confidence bounds to be plotted. The 
#' first value is the lower quantile and the second value is the upper quantile. Default is (0.025,0.975)
#' @param cov_col An integer defining which covariate column in new_pat to plot on the x-axis for the mean survival figure. 
#' 
#' @examples 
#' ########################################
#' #Source dependent packages
#' library(MASS)
#' library(mc2d)
#' library(mvnfast)
#' library(survival)
#' #Simulate Data using built-in data generation
#' ##
#' seed=1
#' set.seed(seed)
#' Npat=5
#' data <- simulate_data(Npat)
#' 
#' ########################################
#' #Run MCMC 
#' ########################################
#' #Inputs for mcmc 
#' response <- log(data$OS)
#' covariate <- cbind(scale(data$Age),data$AUC,data$CR)
#' censor_status <- data$death
#' mcmc_settings<-NULL
#' mcmc_settings$nskip<-10
#' mcmc_settings$nburn<-50
#' mcmc_settings$ndisplay<-100
#' mcmc_settings$nsave<-20
#' mcmc_settings$sigma_jump<-c(0,4,2.5,4,2.6)
#' ###################
#' #Run MCMC function
#' mcmc_Gibbs<-mcmc_DDPGP(response,covariate,censor_status,mcmc_settings)
#' ########################################
#' #Plotting Mean Survival Estimation 
#' ########################################
#' #Parameters in Plotting Mean Survival Estimation 
#' example_CR <- 1
#' example_Age <- 1
#' range_AUC <- seq(2.6, 7, 0.1)
#' new_pat_1<-cbind(example_Age,range_AUC,example_CR)
#' if_plot=1
#' DPGP_mean<-DDPGP_meansurvival(mcmc_Gibbs,new_pat_1,if_plot,cov_col=2)
#' 
#' @return A list composed of the estimation for the mean survival across the values for the covariate specified.
#' The mean_survival is the average of the mean survival calculuation for all iterations, and lower_quant and upper_quant
#' is that for the quantiles specified. optimal is the optimal covariate value (which is specified in cov_col).
#' meansurvival_all is the mean survival for all iterations.

#' @export
#'


# mcmc<-mcmc_Gibbs_out_real; example_Age<-0; range_AUC <- seq(2.6, 7, 0.1); example_CR<-1
# new_pat<-cbind(example_Age,range_AUC,example_CR); quantiles=c(0.025,0.975); cov_col=2; if_plot=1
DDPGP_meansurvival<- function(mcmc,new_pat,if_plot = 0,quantiles=c(0.025,0.975),cov_col=1)
{
  covariate<-mcmc$covariate
  Npat<-mcmc$Npat
  iter_saved = length(mcmc$H)
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each sample point
  d[,1] = 1
  d[,2:Ncov] = covariate

  cov<-d
  
  #If the new_pat is a vetor, change it to a matrix
  if(length(dim(new_pat))==0){
    new_pat <- t(new_pat)
  }
  
  #Quantile Functions
  quantile_lower<-function(x){
    quantile(x,probs=quantiles[1])
  }
  quantile_upper<-function(x){
    quantile(x,probs=quantiles[2])
  }
  #build a matrix to record the density estimation for new sample points.
  N_new <- dim(new_pat)[1]
  Survival_lower<-Survival_upper<-   meansurvival_AUC <- matrix(NA, nrow = iter_saved, ncol = N_new)
  
  tmp <- matrix(0, (iter_saved), length(range))
  for(id_pat in 1:N_new){
    if (id_pat%%5==0){
      cat("Patient Number", id_pat,'of',N_new,'complete \n')
    }
    dpat=c(1,new_pat[id_pat,])
    for (iter in 1:iter_saved)
    {
      M=mcmc$M[iter];r=mcmc$r[,iter];
      pat_cov=dpat;
      beta_0=mcmc$beta_0
      H=mcmc$H[iter]
      betah=matrix(t(mcmc$betah[1:H,,iter]), ncol=Ncov);
      muh=mcmc$muh[1:H,,iter]; sigma2 = mcmc$sigma2[iter];ns=mcmc$ns[,iter]
      wh<-(ns/Npat)[1:H]
      theta<-mcmc$theta[,iter]
      sigma2_cov<-mcmc$sigma2_cov[iter]
      covariance<-update_cov(Ncov, Npat, cov,theta,sigma2_cov)
      inv_covariance = chol2inv(chol(covariance))
      nj <- table(r)
      
      # new_cov<-rbind(pat_cov,cov)
      # new_covar<-update_cov(Ncov, Npat+1,new_cov,theta,sigma2_cov)
      # c1<-as.vector(new_covar[2:(Npat+1),1])

      coeff = (t(cov)-pat_cov)/(theta)
      c1 = exp(-colSums(coeff^2))*sigma2_cov^2
      

      # cluster_var <- as.numeric(sigma2 + sigma2_cov^2+0.01 - c1%*%inv_covariance%*%c1)
      cluster_mean <- t(betah%*%pat_cov)+c1%*%inv_covariance%*%t(muh-t(cov%*%t(betah)))

      #cluster_mean[which(cluster_mean>max(response))] = max(response)
      # meansurvival_AUC[iter,id_pat]= sum(wh * exp(cluster_mean+0.5*cluster_var)) + M/(Npat+M)*exp(beta_0%*%pat_cov)
      meansurvival_AUC[iter,id_pat]= sum(wh * cluster_mean)
      # meansurvival_AUC[iter,id_pat]= sum(wh * cluster_mean+0.5*cluster_var)
    }
    
    
  }
  meansurvival_AUC_mean<-apply(meansurvival_AUC, 2, mean)
  Survival_lower=apply((meansurvival_AUC),2,quantile_lower)
  Survival_upper=apply((meansurvival_AUC),2,quantile_upper)
  id_max_meansurvival =  which(meansurvival_AUC_mean==max(meansurvival_AUC_mean))
  cov_star<-new_pat[,cov_col]
  if(if_plot){
    # plot(cov_star,meansurvival_AUC_mean,col='black',type="l", xlab='Covariate',
    #      main = paste0("Mean Survival Plot"))
    # legend("topleft",paste("Covariate Max=", cov_star[id_max_meansurvival],sep=""))
    # lines(c(cov_star[id_max_meansurvival],cov_star[id_max_meansurvival]),c(min(Survival_lower),max(Survival_upper)))
    # lines(cov_star,Survival_upper, col='red')
    # lines(cov_star,Survival_lower, col='red')
    plot(cov_star,meansurvival_AUC_mean,col='black',type="l", xlab='Covariate',
         main = paste0("Mean Survival Plot"),ylab="Log of Mean Survival")
    legend("topleft",paste("Covariate Max=", cov_star[id_max_meansurvival],sep=""))
    lines(c(cov_star[id_max_meansurvival],cov_star[id_max_meansurvival]),c(min(Survival_lower),max(Survival_upper)))
    lines(cov_star,Survival_upper, col='red')
    lines(cov_star,Survival_lower, col='red')
    
  }
  out <- NULL
  out$mean_survival<-meansurvival_AUC_mean
  out$lower_quant <- Survival_lower
  out$upper_quant <- Survival_upper
  out$optimal<-cov_star[id_max_meansurvival]
  out$meansurvival_all<-meansurvival_AUC
  return(out)
}


# true_meansurvival <- function(data,new_pat, if_plot = 0.5, cov_col=1,col='black',perc_shade=60)
# { 
#   prob<-data$prob
#   sigma<-data$sigma
#   shape<-data$shape
#   #If the new_pat is a vetor, change it to a matrix
#   if(length(dim(new_pat))==0){
#     new_pat <- t(new_pat)
#   }
#   #build a matrix to record the density for new sample points.
#   N_new <- dim(new_pat)[1]
#   meansurvival_truth <- numeric(N_new)
#   beta1<-data$beta1
#   beta2<-data$beta2
#   if (data$truth=='lognorm'){
#     for(id_pat in 1:N_new){
#       meanlog1<- sum(beta1*new_pat[id_pat,])
#       meanlog2<- sum(beta2*new_pat[id_pat,])
#       meansurvival_truth[id_pat]=(prob)*(meanlog1)+(1-prob)*(meanlog2)
#     }
#   }else if (data$truth=='weib'){
#     for(id_pat in 1:N_new){
#       scale1<- exp(sum(beta1*new_pat[id_pat,]))
#       scale2<- exp(sum(beta2*new_pat[id_pat,]))
#       meanweib1<-log(scale1)+digamma(1)/shape
#       meanweib2<-log(scale2)+digamma(1)/shape
#       meansurvival_truth[id_pat]=(prob)*meanweib1+(1-prob)*meanweib2
#     }
#   }
#   id_max_meansurvival_truth =  which(meansurvival_truth==max(meansurvival_truth))
#   # mycol <- t_col(col, perc = perc_shade)
#   if(if_plot==1){
#     plot(new_pat[,cov_col],meansurvival_truth,col=col,type='l',main='True Survival Function')
#     legend("topright",paste("True Optimal AUC=", range_AUC[id_max_meansurvival_truth],sep=""))
#     lines(c(range_AUC[id_max_meansurvival_truth],range_AUC[id_max_meansurvival_truth]),c(0,max(meansurvival_truth)),col=col)
#     
#   }
#   if(if_plot==0.5){
#     lines(new_pat[,cov_col],meansurvival_truth,col=col)
#     # legend("topright",paste("True Optimal AUC=", range_AUC[id_max_meansurvival_truth],sep=""),box.col=col,bg=col)
#     legend("topright",paste("True Optimal AUC=", range_AUC[id_max_meansurvival_truth],sep=""))
#     lines(c(range_AUC[id_max_meansurvival_truth],range_AUC[id_max_meansurvival_truth]),c(0,max(meansurvival_truth)),col=col)
#     
#   }
#   return(meansurvival_truth)
# }
