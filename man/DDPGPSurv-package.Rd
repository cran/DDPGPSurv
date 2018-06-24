\name{DDPGPSurv-package}
\alias{DDPGPSurv-package}
\alias{DDPGPSurv}
\docType{package}
\title{DDP-GP Survival Analysis Package}
\description{
\packageDescription{DDPGPSurv}
}
\details{
The Dependent Dirichlet Process-Gaussian Process model is summarized below. The Dependent Dirichlet Process has a 
Gaussian Process as its base measure.
\deqn{y_i|\mu_h(x_i),\sigma^2=\sum_{h=1}^{\infty}w_{h}N(\mu_{h}(x_i),\sigma^2)}
\deqn{\mu_h(x)\sim GP(x\beta_h,C(\cdot,\cdot))}
\deqn{C(x_{i},x_{j})=\sigma_0^2exp\{-\Vert \frac{x_{i}-x_{j}}{l^2}\Vert^{2}\}+I_{\{i=j\}}J^{2}}




Here, the mean of the Gaussian Process is modeled using linear regression with regression coefficients \eqn{\beta_h} and J is a small diagonal perturbation to the squared exponential covariance function (set to 0.01).

To complete the model specification, independent hyperpriors are assumed:
\deqn{\beta_{h}\sim N(\beta_{0},\Sigma_{0})}
\deqn{\sigma^{-2}\sim Gamma(\delta_{1},\delta_{2})}
We use a hyperparamter \eqn{M} to control the weights (\eqn{w_{h}}) in the model:
\deqn{M\sim Gamma(\lambda_{1},\lambda_{2})}
}
\author{
\packageAuthor{DDPGPSurv}

Maintainer: \packageMaintainer{DDPGPSurv}
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
#Parameters in Plotting Estimation for Functions 
range=seq(2,8,1)
example_AUC <- 5
example_CR <- 1
example_Age <- 0
new_pat<-cbind(example_Age,example_AUC,example_CR)
if_plot=1
###################
#Plot DDP-GP Density Esimation
a=DDPGP_Dens(mcmc_Gibbs,new_pat,range, if_plot)

#Plot DDP-GP Survival Esimation
b=DDPGP_Surv(mcmc_Gibbs,new_pat,range, if_plot)

#Plot DDP-GP Hazard Esimation
c=DDPGP_Haz(mcmc_Gibbs,new_pat,range, if_plot)
########################################
#Plotting Mean Survival Estimation 
########################################
#Parameters in Plotting Mean Survival Estimation 
range_AUC <- seq(2.6, 7, 0.1)
new_pat_1<-cbind(example_Age,range_AUC,example_CR)
if_plot=1
DDPGP_mean<-DDPGP_meansurvival(mcmc_Gibbs,new_pat_1,if_plot,cov_col=2)
}