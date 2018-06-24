#' Gibbs Sampling MCMC function
#'
#' @param response An n-dimensional vector consisting of the verall survival time of the patients in the dataset (n is the number of patients).
#' @param covariate An nxp dimensional matrix consisting of the relevant covariates for each patient (p is the number of covariates)
#' @param censor_status An n-dimensional vector with the censoring status of the response (0 for censored, 1 for uncensored)
#' @param mcmc_settings A list for MCMC setup. nskip is the number of iterations to skip between saved iterations. 
#' nburn is the number of total iterations to burn. ndisplay is number of iterations per which the display message will 
#' appear. nsave is the total number of iterations to save. The total number of iterations (including those not 
#' saved) will be mcmc_settings$nburn+mcmc_settings$nskip*mcmc_settings$nsave. mcmc_settings$sigma_jump 
#' is a (p+2)-dimensionial vector, where p is the number of 
#' covariates. These are the jump sizes in the M-H algorithm. The first value in this vector
#' refers to the jump size of the 1's covariate column (which is typically 0). The next 
#' p values are the jump sizes of the covariates in the order of the covariate vector. The
#' last value in this vector is the jump size for the variance term (sigma_0) in the 
#' covariance matrix. See example. 
#' @param lambda1 A hyper-parameter in prior for the M (M~Gamma(lambda1, lambda2)). Default value is 1.
#' @param lambda2 A hyper-parameter in prior for the M (M~Gamma(lambda1, lambda2)). Default value is 1.
#' @param delta1 A hyper-parameter in prior for the sigma2 (1/(sigma2^2)~Gamma(delta1, delta2)). Default value is 4.
#' @param delta2 A hyper-parameter in prior for the sigma2 (1/(sigma2^2)~Gamma(delta1, delta2)). Default value is 3.
#' @param mu0 A hyper-parameter for the mean in the normal prior for the Metropolis Hastings algorithm for the scale parameter (l) in the covariance matrix. 
#' Default value is 0.
#' @param tau2 A hyper-parameter for the variance in the normal prior for the Metropolis Hastings algorithm for the scale parameter (l) in the covariance matrix. 
#' Default value is 10.
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
#' #Plotting Density/Survival/Hazard Estimation 
#' ########################################
# #Parameters in Plotting Estimation for Functions
# range=seq(2,8,1)
# example_AUC <- 5
# example_CR <- 1
# example_Age <- 1
# new_pat<-cbind(example_Age,example_AUC,example_CR)
# if_plot=1
# ###################
# #Plot DDP-GP Density Esimation
# a=DDPGP_Dens(mcmc_Gibbs,new_pat,range, if_plot)
# 
# #Plot DDP-GP Survival Esimation
# b=DDPGP_Surv(mcmc_Gibbs,new_pat,range, if_plot)
# 
# #Plot DDP-GP Hazard Esimation
# c=DDPGP_Haz(mcmc_Gibbs,new_pat,range, if_plot)
# ########################################
# #Plotting Mean Survival Estimation
# ########################################
# #Parameters in Plotting Mean Survival Estimation
# range_AUC <- seq(2.6, 7, 0.1)
# new_pat_1<-cbind(example_Age,range_AUC,example_CR)
# if_plot=1
# DPGP_mean<-DDPGP_meansurvival(mcmc_Gibbs,new_pat_1,if_plot,cov_col=2)
#' 
#' @return A list composed of the mcmc outputs. These outputs can be inputted into the other functions in this package.
#'
#' @export
mcmc_DDPGP<-function(response,covariate,censor_status,mcmc_settings, lambda1 = 1, lambda2 = 1,
                                delta1 = 4, delta2 = 3, mu0 = 0, tau2 = 10){
  mcmc<-NULL
  # if (hyperparam==0){
  #   mcmc<-mcmc_Gibbs(response,covariate,censor_status,mcmc_settings,lambda1 = 1, lambda2 = 1,
  #                     delta1 = 4, delta2 = 3)
  # } else if (hyperparam==1){
    mcmc<-mcmc_Gibbs_Hyperparam(response,covariate,censor_status,mcmc_settings,lambda1 = lambda1, lambda2 = lambda2,
                      delta1 = delta1, delta2 = delta2, mu0 =mu0, tau2 = tau2)
  # }
  return(mcmc)
}


  
#' Posterior density function estimation using MCMC results
#' @param mcmc A list-The MCMC output from mcmc_Gibbs or FiniteDP
#' @param new_pat A matrix consisting of the covariates of new sample points. Each row represents one sample point. 
#' It can also be vector consisting of the covariates of a new sample point. 
#' @param range A vector consisting of the values to be evaluated in the plot for density estimation. 
#' @param if_plot A logical variable indicating whether the density estimation should be plotted. 
#' 1 means estimation should be plotted; 0 means otherwise. Defualt value is 0, i.e. not plotted. 
#' @param quantiles A two-dimensional vector determining the quantiles for the confidence bounds to be plotted. The 
#' first value is the lower quantile and the second value is the upper quantile. Default is (0.025,0.975)
#' @param color A string determining the color of the plot. 
#' 
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
#' #Plotting Density Estimation 
#' ########################################
#' #Parameters in Plotting Estimation for Functions 
#' range=seq(2,8,1)
#' example_AUC <- 5
#' example_CR <- 1
#' example_Age <- 1
#' new_pat<-cbind(example_Age,example_AUC,example_CR)
#' if_plot=1
#' ###################
#' #Plot DDP-GP Density Esimation
#' a=DDPGP_Dens(mcmc_Gibbs,new_pat,range, if_plot)
#' 
#' @return A list composed of the estimation for the density across the given range. The Density_estimation is the 
#' mean across all iterations, and the Density_lower_quant and Density_upper_quant are the density estimations for the
#' specified quantiles.
#' Each row represents the estimation for each new sample point.
#' 
#' @export
DDPGP_Dens <- function(mcmc,new_pat,range, if_plot = 0, quantiles=c(0.025,0.975),color='green'){
  dens<-NULL
    dens<-DDPGP_density_Hyperparam(mcmc,new_pat,range, if_plot, quantiles,color)
  return(dens)
}




#' Posterior survival function estimation using Gibbs MCMC results
#' @param mcmc A list-The MCMC output from mcmc_Gibbs 
#' @param new_pat A matrix consisting of the covariates of new sample points. Each row represents one sample point. 
#' It can also be vector consisting of the covariates of a new sample point. 
#' @param range A vector consisting of the values to be evaluated in the plot for survival estimation. 
#' @param if_plot A logical variable indicating whether the survival estimation should be plotted. 
#' 1 means estimation should be plotted; 0 means otherwise. Defualt value is 0, i.e. not plotted. 
#' @param quantiles A two-dimensional vector determining the quantiles for the confidence bounds to be plotted. The 
#' first value is the lower quantile and the second value is the upper quantile. Default is (0.025,0.975)
#' @param color A string determining the color of the plot. 
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
#' #Plotting Survival Estimation 
#' ########################################
#' #Parameters in Plotting Estimation for Functions 
#' range=seq(2,8,1)
#' example_AUC <- 5
#' example_CR <- 1
#' example_Age <- 1
#' new_pat<-cbind(example_Age,example_AUC,example_CR)
#' if_plot=1
#' ###################
#' #Plot DDP-GP Survival Esimation
#' b=DDPGP_Surv(mcmc_Gibbs,new_pat,range, if_plot)
#' 
#' @return A list composed of the estimation for the survival across the given range. The Survival_estimation is the 
#' mean across all iterations, and the Survival_lower_quant and Survival_upper_quant are the survival function estimations for the
#' specified quantiles.
#' 
#' @export
DDPGP_Surv <- function(mcmc,new_pat,range, if_plot = 0,quantiles=c(0.025,0.975),color='green'){
  surv<-NULL
    surv<-DDPGP_survival_Hyperparam(mcmc,new_pat,range, if_plot, quantiles,color)
  return(surv)
}





#' Posterior hazard function estimation using Gibbs MCMC results
#' @param mcmc A list-The MCMC output from mcmc_Gibbs 
#' @param new_pat A matrix consisting of the covariates of new sample points. Each row represents one sample point. 
#' It can also be vector consisting of the covariates of a new sample point. 
#' @param range A vector consisting of the values to be evaluated in the plot for hazard estimation. 
#' @param if_plot A logical variable indicating whether the hazard estimation should be plotted. 
#' 1 means estimation should be plotted; 0 means otherwise. Defualt value is 0, i.e. not plotted. 
#' @param quantiles A two-dimensional vector determining the quantiles for the confidence bounds to be plotted. The 
#' first value is the lower quantile and the second value is the upper quantile. Default is (0.025,0.975)
#' @param color A string determining the color of the plot. 
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
#' #Plotting Hazard Estimation 
#' ########################################
#' #Parameters in Plotting Estimation for Functions 
#' range=seq(2,8,1)
#' example_AUC <- 5
#' example_CR <- 1
#' example_Age <- 1
#' new_pat<-cbind(example_Age,example_AUC,example_CR)
#' if_plot=1
#' ###################
#' #Plot DDP-GP Hazard Esimation
#' c=DDPGP_Haz(mcmc_Gibbs,new_pat,range, if_plot)
#'
#' @return A list composed of the estimation for the hazard across the given range. The Hazard_estimation is the 
#' mean across all iterations, and the Hazard_lower_quant and Hazard_upper_quant are the hazard function estimations for the
#' specified quantiles.
#' 
#' @export
DDPGP_Haz <- function(mcmc,new_pat,range, if_plot = 0,quantiles=c(0.025,0.975),color='green'){
  hazard<-NULL
    hazard<-DDPGP_hazard_Hyperparam(mcmc,new_pat,range, if_plot, quantiles,color)
  return(hazard)
}