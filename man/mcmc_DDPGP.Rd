% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main_Functions.R
\name{mcmc_DDPGP}
\alias{mcmc_DDPGP}
\title{Gibbs Sampling MCMC function}
\usage{
mcmc_DDPGP(response, covariate, censor_status, mcmc_settings, lambda1 = 1,
  lambda2 = 1, delta1 = 4, delta2 = 3, mu0 = 0, tau2 = 10)
}
\arguments{
\item{response}{An n-dimensional vector consisting of the verall survival time of the patients in the dataset (n is the number of patients).}

\item{covariate}{An nxp dimensional matrix consisting of the relevant covariates for each patient (p is the number of covariates)}

\item{censor_status}{An n-dimensional vector with the censoring status of the response (0 for censored, 1 for uncensored)}

\item{mcmc_settings}{A list for MCMC setup. nskip is the number of iterations to skip between saved iterations. 
nburn is the number of total iterations to burn. ndisplay is number of iterations per which the display message will 
appear. nsave is the total number of iterations to save. The total number of iterations (including those not 
saved) will be mcmc_settings$nburn+mcmc_settings$nskip*mcmc_settings$nsave. mcmc_settings$sigma_jump 
is a (p+2)-dimensionial vector, where p is the number of 
covariates. These are the jump sizes in the M-H algorithm. The first value in this vector
refers to the jump size of the 1's covariate column (which is typically 0). The next 
p values are the jump sizes of the covariates in the order of the covariate vector. The
last value in this vector is the jump size for the variance term (sigma_0) in the 
covariance matrix. See example.}

\item{lambda1}{A hyper-parameter in prior for the M (M~Gamma(lambda1, lambda2)). Default value is 1.}

\item{lambda2}{A hyper-parameter in prior for the M (M~Gamma(lambda1, lambda2)). Default value is 1.}

\item{delta1}{A hyper-parameter in prior for the sigma2 (1/(sigma2^2)~Gamma(delta1, delta2)). Default value is 4.}

\item{delta2}{A hyper-parameter in prior for the sigma2 (1/(sigma2^2)~Gamma(delta1, delta2)). Default value is 3.}

\item{mu0}{A hyper-parameter for the mean in the normal prior for the Metropolis Hastings algorithm for the scale parameter (l) in the covariance matrix. 
Default value is 0.}

\item{tau2}{A hyper-parameter for the variance in the normal prior for the Metropolis Hastings algorithm for the scale parameter (l) in the covariance matrix. 
Default value is 10.}
}
\value{
A list composed of the mcmc outputs. These outputs can be inputted into the other functions in this package.
}
\description{
Gibbs Sampling MCMC function
}
\examples{
########################################
#Source dependent packages
library(MASS)
library(mc2d)
library(mvnfast)
library(survival)
#Simulate Data using built-in data generation
##
seed=1
set.seed(seed)
Npat=5
data <- simulate_data(Npat)

########################################
#Run MCMC 
########################################
#Inputs for mcmc 
response <- log(data$OS)
covariate <- cbind(scale(data$Age),data$AUC,data$CR)
censor_status <- data$death
mcmc_settings<-NULL
mcmc_settings$nskip<-10
mcmc_settings$nburn<-50
mcmc_settings$ndisplay<-100
mcmc_settings$nsave<-20
mcmc_settings$sigma_jump<-c(0,4,2.5,4,2.6)
###################
#Run MCMC function
mcmc_Gibbs<-mcmc_DDPGP(response,covariate,censor_status,mcmc_settings)
########################################
#Plotting Density/Survival/Hazard Estimation 
########################################

}
