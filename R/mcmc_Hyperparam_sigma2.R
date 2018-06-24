
update_cov <- function(Ncov, Npat, cov, theta,sigma2_cov)
{
  cova = matrix(0, Npat, Npat)
  for ( i in 1:Ncov)
  {
    cova = cova + outer(cov[,i],cov[,i],'-')^2/(theta[i]^2)
  }
  return(sigma2_cov^2*exp(-cova) + 0.01*diag(Npat))
}

#nlm(L,c(1,1,1),muh=muh, Npat, cov, r=r, Ncov)


L <- function(theta, muh, Npat, cov, Ncov,betah,r,H,sigma2_cov,mu0,tau2)
{
  if (is.vector(muh)){
    muh<-matrix(muh,nrow=1)
  }
  if (is.vector(betah)){
    betah<-matrix(betah,nrow=1)
  }
  covariance = update_cov(Ncov,Npat,cov,theta,sigma2_cov)
  sum1 = 0
  sum2 = 0
  for (h in 1:H)
  {
    if (any(r==h)){
      Ah = which(r==h)
      mu=cov[Ah,]%*%betah[h,]
      sum1 = sum1 + t(muh[h,Ah]-mu)%*%chol2inv(chol(covariance[Ah, Ah]))%*%(muh[h,Ah]-mu)
      if (length(Ah)==1) sum2 = sum2 + log(covariance[Ah,Ah]) else
        sum2 = sum2 + determinant(covariance[Ah, Ah], logarithm=TRUE)$modulus
    }}
  return(-1/2*sum2-1/2*sum1- sum((theta-mu0)^2/(2*tau2)))
  # return(-1/2*sum2-1/2*sum1- sum((theta-mu0)^2/(2*tau2))  - ((sigma2_cov-mu0)^2/(2*tau2)) )
  
  #Change the prior for sigma0 to normal 
  #return(-1/2*sum2-1/2*sum1)
  
}

update_theta_comp<- function(theta,muh,Npat,cov,Ncov,betah,i,sigma_jump,r,H,sigma2_cov,mu0,tau2){
  prop_theta<-theta
  prop_theta[i]<- rnorm(1,theta[i],sigma_jump[i])
  # if (prop_theta[i]<0){
  #   prop_theta[i]=-prop_theta[i]
  # }
  loglike.curr<- L(theta,muh,Npat,cov,Ncov,betah,r,H,sigma2_cov,mu0,tau2)
  loglike.prop<- L(prop_theta,muh,Npat,cov,Ncov,betah,r,H,sigma2_cov,mu0,tau2)
  ratio<- exp(loglike.prop-loglike.curr)
  if (runif(1)<ratio){
    theta_updated<-prop_theta
    acc=1
  } else{
    theta_updated<-theta
    acc=0
  }
  return(list(theta_updated=theta_updated,acc=acc))
}

# theta=mcmc$theta[,iter-1]
# muh=mcmc$muh[1:mcmc$H[iter],,iter]
# betah=mcmc$betah[1:mcmc$H[iter],,iter]
# r=mcmc$r[,iter]
# H=mcmc$H[iter]
update_theta <- function(theta,muh,Npat,cov,Ncov,betah,sigma_jump,r,H,sigma2_cov,mu0,tau2)
{
  accept<-numeric(Ncov)
  for(i in 1:Ncov){
    out<-update_theta_comp(theta,muh,Npat,cov,Ncov,betah,i,sigma_jump,r,H,sigma2_cov,mu0,tau2)
    theta<-out$theta_updated
    accept[i]<-out$acc
  }
  return(list(theta=theta,accept=accept))
}

#just one theta for all covariates
update_theta_single <- function(theta,muh,Npat,cov,Ncov,betah,sigma_jump,r,H,sigma2_cov,mu0,tau2)
{
  prop_theta<- rep(rnorm(1,theta[1],sigma_jump[1]),Ncov)
  # if (prop_theta[i]<0){
  #   prop_theta[i]=-prop_theta[i]
  # }
  loglike.curr<- L(theta,muh,Npat,cov,Ncov,betah,r,H,sigma2_cov,mu0,tau2)
  loglike.prop<- L(prop_theta,muh,Npat,cov,Ncov,betah,r,H,sigma2_cov,mu0,tau2)
  ratio<- exp(loglike.prop-loglike.curr)
  if (runif(1)<ratio){
    theta_updated<-prop_theta
    acc=1
  } else{
    theta_updated<-theta
    acc=0
  }

  accept<-rep(acc,Ncov)
  return(list(theta=theta_updated,accept=accept))
}

update_sigma2_cov<- function(theta,muh,Npat,cov,Ncov,betah,sigma_jump,r,H,sigma2_cov,mu0,tau2){
  prop_sigma2<-rnorm(1,sigma2_cov,sigma_jump[Ncov+1])
  # if (prop_sigma2<0){
  #   prop_sigma2=-prop_sigma2
  # }
  loglike.curr<- L(theta,muh,Npat,cov,Ncov,betah,r,H,sigma2_cov,mu0,tau2)
  loglike.prop<- L(theta,muh,Npat,cov,Ncov,betah,r,H,prop_sigma2,mu0,tau2)
  ratio<- exp(loglike.prop-loglike.curr)
  if (runif(1)<ratio){
    sigma2_updated<-prop_sigma2
    acc=1
  } else{
    sigma2_updated<-sigma2_cov
    acc=0
  }
  return(list(sigma2_cov=sigma2_updated,accept=acc))
}

mcmc_Gibbs_Hyperparam_singletheta<-function(response,covariate,censor_status,mcmc_settings,lambda1 = 1, lambda2 = 1,
                                delta1 = 4, delta2 = 3,mu0=0, tau2=100)
{
  Niter=mcmc_settings$nsave*mcmc_settings$nskip+mcmc_settings$nburn
  burn.in=mcmc_settings$nburn
  skip<-mcmc_settings$nskip
  sigma_jump<-mcmc_settings$sigma_jump
  Npat <- length(response) #Number of patients
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each patient
  d[,1] = 1
  d[,2:Ncov] = covariate
  #Define covariance matrix for GP
  cov=d
  #I used Finite DP to approximate DP. Of course, we can change it later
  #prior for sigma2 let mean=1, var=0.5
  # assign('lambda1', lambda1, envir = .GlobalEnv)
  # assign('lambda2', lambda2, envir = .GlobalEnv)
  # assign('delta1', delta1, envir = .GlobalEnv)
  # assign('delta2', delta2, envir = .GlobalEnv)
  # assign('mu0', mu0, envir = .GlobalEnv)
  # assign('tau2', tau2, envir = .GlobalEnv)
  beta_0 <- rep(0,Ncov) #Prior for betah
  
  ########################################################################
  #mcmc draws for each variable
  mcmc <- NULL
  mcmc$M <- rep(NA, Niter)
  mcmc$H <- rep(NA,Niter)
  mcmc$sigma2 <- rep(NA, Niter)
  mcmc$muh <- array(NA, c(Npat, Npat, Niter))
  mcmc$betah <- array(NA, c(Npat, Ncov, Niter))
  mcmc$r <- array(NA, c(Npat, Niter))
  mcmc$ns <- array(NA, c(Npat, Niter))
  #Here mcmc$imputey is for the censoring setup. If censoring exsits, then we need
  #to impute y. But for now, imputey=response
  mcmc$imputey <- matrix(NA, Npat, Niter)
  mcmc$theta<-matrix(NA,Ncov,Niter)
  mcmc$sigma2_cov<-numeric(Niter)
  #####################
  set.seed(1)
  ##Initialize variables
  H<- 5
  initial <- init(response, Npat,H)
  mcmc$M[1] = initial$M
  mcmc$H[1] = H
  mcmc$sigma2[1] = 10
  mcmc$muh[1:H,,1] = initial$muh
  mcmc$r[,1] = initial$r
  for(j in 1:H){mcmc$ns[j,1]=sum(mcmc$r[,1]==j)}
  #assign the initial values using a survival regression first
  lfit = survreg(Surv(exp(response),censor_status)~d,dist="weibull")
  mcmc$betah[1:H,,1] = matrix(rep(lfit$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc$imputey[,1] = response
  mcmc$theta[,1]=rep(1,Ncov)
  theta_accept<-numeric(Ncov)
  mcmc$sigma2_cov[1]<-1
  sigma2_accept<-0
  
  ###################################################################################
  #run MCMC Gibbs sampling
  ptm <- proc.time()
  
  for (iter in 2:Niter)
  {
    # print(iter)
    if (iter%%mcmc_settings$ndisplay==0){
      cat("Iteration", iter,'of',Niter,'complete \n')
      # print(mcmc$theta[,iter-1])
      # print(sigma2_accept)
      # print(theta_accept)
    }
    
    covariance<-update_cov(Ncov, Npat, cov,mcmc$theta[,iter-1],mcmc$sigma2_cov[iter-1])
    inv_covariance = chol2inv(chol((covariance)))
    #lambda1 = 1; lambda2 = 1;delta1 = 4; delta2 = 3
    tmp = update_r(mcmc$sigma2[iter-1],
                   matrix(mcmc$muh[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1]),
                   Npat, censor_status, response,
                   mcmc$r[,iter-1],
                   covariance,cov,
                   mcmc$H[iter-1],
                   mcmc$M[iter-1],
                   mcmc$ns[1:mcmc$H[iter-1],iter-1],beta_0,Ncov)
    mcmc$r[,iter] = tmp$r
    mcmc$imputey[,iter] = tmp$impute_y
    mcmc$H[iter] = tmp$H
    mcmc$muh[1:mcmc$H[iter],,iter] = tmp$muh
    
    # sigma2<-mcmc$sigma2[iter-1]
    # muh<-matrix(mcmc$muh[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1])
    # betah<-matrix(mcmc$betah[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1])
    # r<-mcmc$r[,iter-1]
    # H<- mcmc$H[iter-1]
    # M<- mcmc$M[iter-1]
    # ns<-mcmc$ns[1:mcmc$H[iter-1],iter-1]
    
    mcmc$betah[1:mcmc$H[iter],,iter] = update_betah(Ncov,mcmc$muh[1:mcmc$H[iter],,iter],
                                                    cov, covariance,
                                                    mcmc$r[,iter], inv_covariance,
                                                    mcmc$H[iter],beta_0)
    # print(mcmc$betah[1:mcmc$H[iter],,iter])
    # betah_save<-mcmc$betah[1:mcmc$H[iter],,iter]
    mcmc$muh[1:mcmc$H[iter],,iter] = update_muh(mcmc$r[,iter],
                                                mcmc$sigma2[iter-1],
                                                mcmc$imputey[,iter],
                                                Npat, covariance, cov,
                                                mcmc$betah[1:mcmc$H[iter],,iter],
                                                inv_covariance,
                                                mcmc$H[iter])
    # print(mcmc$muh[1:mcmc$H[iter],,iter])
    
    mcmc$ns[1:mcmc$H[iter],iter] <- update_ns(mcmc$r[,iter],mcmc$H[iter])
    mcmc$sigma2[iter] = update_sigma2(mcmc$r[,iter],
                                      mcmc$imputey[,iter],
                                      mcmc$muh[1:mcmc$H[iter],,iter], Npat,
                                      mcmc$H[iter], delta1, delta2)
    mcmc$M[iter] = update_M(mcmc$M[iter-1], mcmc$H[iter],lambda1,lambda2,Npat)
    out<- update_theta_single(mcmc$theta[,iter-1],mcmc$muh[1:mcmc$H[iter],,iter],Npat,cov,Ncov,mcmc$betah[1:mcmc$H[iter],,iter],sigma_jump,mcmc$r[,iter],mcmc$H[iter],mcmc$sigma2_cov[iter-1],mu0,tau2)
    mcmc$theta[,iter]<-out$theta
    out_sigma<-update_sigma2_cov(mcmc$theta[,iter],mcmc$muh[1:mcmc$H[iter],,iter],Npat,cov,Ncov,mcmc$betah[1:mcmc$H[iter],,iter],sigma_jump,mcmc$r[,iter],mcmc$H[iter],mcmc$sigma2_cov[iter-1],mu0,tau2)
    mcmc$sigma2_cov[iter]<-out_sigma$sigma2_cov
    
    if (iter>burn.in){
      theta_accept=theta_accept+out$accept
      sigma2_accept=sigma2_accept+out_sigma$accept
    }
    
    # print(mcmc$theta[,iter])
  }
  id<-seq(burn.in+skip,Niter,skip)
  mcmc1<-NULL
  mcmc1$sigma2=mcmc$sigma2[id]
  mcmc1$betah=mcmc$betah[,,id]
  mcmc1$muh=mcmc$muh[,,id]
  mcmc1$ns=mcmc$ns[,id]
  mcmc1$r=mcmc$r[,id]
  mcmc1$imputey=mcmc$imputey[,id]
  mcmc1$H=mcmc$H[id]
  mcmc1$M=mcmc$M[id]
  mcmc1$covariate=covariate
  mcmc1$Npat=Npat
  mcmc1$beta_0=beta_0
  mcmc1$theta<-mcmc$theta[,id]
  mcmc1$theta_accept<-theta_accept
  mcmc1$sigma2_cov<-mcmc$sigma2_cov[id]
  mcmc1$sigma2_accept<-sigma2_accept
  
  mcmc1$Finite=FALSE  
  mcmc1$Hyper=TRUE
  
  # mcmc1$Niter=Niter
  # mcmc1$burn.in=burn.in
  return(mcmc1)
}



mcmc_Gibbs_Hyperparam<-function(response,covariate,censor_status,mcmc_settings,lambda1 = 1, lambda2 = 1,
                                delta1 = 4, delta2 = 3,mu0=0, tau2=100)
{
  Niter=mcmc_settings$nsave*mcmc_settings$nskip+mcmc_settings$nburn
  burn.in=mcmc_settings$nburn
  skip<-mcmc_settings$nskip
  sigma_jump<-mcmc_settings$sigma_jump
  Npat <- length(response) #Number of patients
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each patient
  d[,1] = 1
  d[,2:Ncov] = covariate
  #Define covariance matrix for GP
  cov=d
  #I used Finite DP to approximate DP. Of course, we can change it later
  #prior for sigma2 let mean=1, var=0.5
  # assign('lambda1', lambda1, envir = .GlobalEnv)
  # assign('lambda2', lambda2, envir = .GlobalEnv)
  # assign('delta1', delta1, envir = .GlobalEnv)
  # assign('delta2', delta2, envir = .GlobalEnv)
  # assign('mu0', mu0, envir = .GlobalEnv)
  # assign('tau2', tau2, envir = .GlobalEnv)
  beta_0 <- rep(0,Ncov) #Prior for betah
  
  ########################################################################
  #mcmc draws for each variable
  mcmc <- NULL
  mcmc$M <- rep(NA, Niter)
  mcmc$H <- rep(NA,Niter)
  mcmc$sigma2 <- rep(NA, Niter)
  mcmc$muh <- array(NA, c(Npat, Npat, Niter))
  mcmc$betah <- array(NA, c(Npat, Ncov, Niter))
  mcmc$r <- array(NA, c(Npat, Niter))
  mcmc$ns <- array(NA, c(Npat, Niter))
  #Here mcmc$imputey is for the censoring setup. If censoring exsits, then we need
  #to impute y. But for now, imputey=response
  mcmc$imputey <- matrix(NA, Npat, Niter)
  mcmc$theta<-matrix(NA,Ncov,Niter)
  mcmc$sigma2_cov<-numeric(Niter)
  #####################
  set.seed(1)
  ##Initialize variables
  H<- 5
  initial <- init(response, Npat,H)
  mcmc$M[1] = initial$M
  mcmc$H[1] = H
  mcmc$sigma2[1] = 10
  mcmc$muh[1:H,,1] = initial$muh
  mcmc$r[,1] = initial$r
  for(j in 1:H){mcmc$ns[j,1]=sum(mcmc$r[,1]==j)}
  #assign the initial values using a survival regression first
  lfit = survreg(Surv(exp(response),censor_status)~d,dist="weibull")
  mcmc$betah[1:H,,1] = matrix(rep(lfit$coefficients[-2], H),H,Ncov, byrow=T)
  mcmc$imputey[,1] = response
  mcmc$theta[,1]=rep(1,Ncov)
  theta_accept<-numeric(Ncov)
  mcmc$sigma2_cov[1]<-1
  sigma2_accept<-0
  
  ###################################################################################
  #run MCMC Gibbs sampling
  ptm <- proc.time()
  
  for (iter in 2:Niter)
  {
    # print(iter)
    if (iter%%mcmc_settings$ndisplay==0){
      cat("Iteration", iter,'of',Niter,'complete \n')
      # print(mcmc$theta[,iter-1])
      # print(sigma2_accept)
      # print(theta_accept)
    }
    
    covariance<-update_cov(Ncov, Npat, cov,mcmc$theta[,iter-1],mcmc$sigma2_cov[iter-1])
    inv_covariance = chol2inv(chol((covariance)))
    #lambda1 = 1; lambda2 = 1;delta1 = 4; delta2 = 3
    tmp = update_r(mcmc$sigma2[iter-1],
                        matrix(mcmc$muh[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1]),
                        Npat, censor_status, response,
                        mcmc$r[,iter-1],
                        covariance,cov,
                        mcmc$H[iter-1],
                        mcmc$M[iter-1],
                        mcmc$ns[1:mcmc$H[iter-1],iter-1],beta_0,Ncov)
    mcmc$r[,iter] = tmp$r
    mcmc$imputey[,iter] = tmp$impute_y
    mcmc$H[iter] = tmp$H
    mcmc$muh[1:mcmc$H[iter],,iter] = tmp$muh

    # sigma2<-mcmc$sigma2[iter-1]
    # muh<-matrix(mcmc$muh[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1])
    # betah<-matrix(mcmc$betah[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1])
    # r<-mcmc$r[,iter-1]
    # H<- mcmc$H[iter-1]
    # M<- mcmc$M[iter-1]
    # ns<-mcmc$ns[1:mcmc$H[iter-1],iter-1]

    mcmc$betah[1:mcmc$H[iter],,iter] = update_betah(Ncov,mcmc$muh[1:mcmc$H[iter],,iter],
                                                    cov, covariance,
                                                    mcmc$r[,iter], inv_covariance,
                                                    mcmc$H[iter],beta_0)
    # print(mcmc$betah[1:mcmc$H[iter],,iter])
    # betah_save<-mcmc$betah[1:mcmc$H[iter],,iter]
    mcmc$muh[1:mcmc$H[iter],,iter] = update_muh(mcmc$r[,iter],
                                                mcmc$sigma2[iter-1],
                                                mcmc$imputey[,iter],
                                                Npat, covariance, cov,
                                                mcmc$betah[1:mcmc$H[iter],,iter],
                                                inv_covariance,
                                                mcmc$H[iter])
    # print(mcmc$muh[1:mcmc$H[iter],,iter])

    mcmc$ns[1:mcmc$H[iter],iter] <- update_ns(mcmc$r[,iter],mcmc$H[iter])
    mcmc$sigma2[iter] = update_sigma2(mcmc$r[,iter],
                                      mcmc$imputey[,iter],
                                      mcmc$muh[1:mcmc$H[iter],,iter], Npat,
                                      mcmc$H[iter], delta1, delta2)
    mcmc$M[iter] = update_M(mcmc$M[iter-1], mcmc$H[iter],lambda1,lambda2,Npat)
    out<- update_theta(mcmc$theta[,iter-1],mcmc$muh[1:mcmc$H[iter],,iter],Npat,cov,Ncov,mcmc$betah[1:mcmc$H[iter],,iter],sigma_jump,mcmc$r[,iter],mcmc$H[iter],mcmc$sigma2_cov[iter-1],mu0,tau2)
    mcmc$theta[,iter]<-out$theta
    out_sigma<-update_sigma2_cov(mcmc$theta[,iter],mcmc$muh[1:mcmc$H[iter],,iter],Npat,cov,Ncov,mcmc$betah[1:mcmc$H[iter],,iter],sigma_jump,mcmc$r[,iter],mcmc$H[iter],mcmc$sigma2_cov[iter-1],mu0,tau2)
    mcmc$sigma2_cov[iter]<-out_sigma$sigma2_cov
    
    if (iter>burn.in){
      theta_accept=theta_accept+out$accept
      sigma2_accept=sigma2_accept+out_sigma$accept
    }
    
    # print(mcmc$theta[,iter])
  }
  id<-seq(burn.in+skip,Niter,skip)
  mcmc1<-NULL
  mcmc1$sigma2=mcmc$sigma2[id]
  mcmc1$betah=mcmc$betah[,,id]
  mcmc1$muh=mcmc$muh[,,id]
  mcmc1$ns=mcmc$ns[,id]
  mcmc1$r=mcmc$r[,id]
  mcmc1$imputey=mcmc$imputey[,id]
  mcmc1$H=mcmc$H[id]
  mcmc1$M=mcmc$M[id]
  mcmc1$covariate=covariate
  mcmc1$Npat=Npat
  mcmc1$beta_0=beta_0
  mcmc1$theta<-mcmc$theta[,id]
  mcmc1$theta_accept<-theta_accept
  mcmc1$sigma2_cov<-mcmc$sigma2_cov[id]
  mcmc1$sigma2_accept<-sigma2_accept
  
  mcmc1$Finite=FALSE  
  mcmc1$Hyper=TRUE
  
  # mcmc1$Niter=Niter
  # mcmc1$burn.in=burn.in
  return(mcmc1)
}

