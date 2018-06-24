
DDPGP_density_survival_Hyperparam <- function(mcmc,new_pat,range, if_plot = 0, quantiles=c(0.025,0.975),color='green')
{
  beta_0<-mcmc$beta_0
  covariate<-mcmc$covariate
  Npat<-mcmc$Npat
  iter_saved = length(mcmc$H)
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each sample point
  d[,1] = 1
  d[,2:Ncov] = covariate
  
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
  Survival_lower<-Survival_upper<- Survival_estimation<- Density_lower<-Density_upper<-Density_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  
  tmp_s<- tmp <- matrix(0, (iter_saved), length(range))
  for(id_pat in 1:N_new){
    dpat=c(1,new_pat[id_pat,])
    for (iter in 1:iter_saved)
    {
      M=mcmc$M[iter];r=mcmc$r[,iter];
      pat_cov=dpat;wh=mcmc$wh[,iter];
      muh=mcmc$muh[,,iter]; betah=mcmc$betah[,,iter];H=mcmc$H[iter]
      sigma2 = mcmc$sigma2[iter]
      theta<-mcmc$theta[,iter]
      sigma2_cov<-mcmc$sigma2_cov[iter]
      cov<-d
      covariance<-update_cov(Ncov, Npat, cov,theta,sigma2_cov)
      inv_covariance = chol2inv(chol(covariance))
      nj <- table(r)
      # ths <- as.numeric(names(nj)) # unique values
      # coeff = (t(cov)-pat_cov)
      # c1 = exp(-colSums(coeff^2*theta^2))
      # 
      new_cov<-rbind(pat_cov,cov)
      new_covar<-update_cov(Ncov, Npat+1,new_cov,theta,sigma2_cov)
      c1<-as.vector(new_covar[2:(Npat+1),1])
      
      ft=0
      Ft=0
      for (h in 1:H)
      {
        mean=betah[h,]%*%pat_cov+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,])
        ft <- ft + nj[h]/(Npat+M)*dnorm((range), mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
        Ft <- Ft + nj[h]/(Npat+M)*pnorm(range, mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
        
      }
      ft<- ft + M/(Npat+M)*dnorm((range),mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      Ft<- Ft + M/(Npat+M)*pnorm(range,mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      tmp[iter,] = ft   
      tmp_s[iter,] = 1-Ft

    }

    Density_estimation[id_pat,] = apply(tmp, 2, mean)
    Density_lower[id_pat,]=apply(tmp,2,quantile_lower)
    Density_upper[id_pat,]=apply(tmp,2,quantile_upper)
    Survival_estimation[id_pat,] = apply(tmp_s, 2, mean)
    Survival_lower[id_pat,]=apply(tmp_s,2,quantile_lower)
    Survival_upper[id_pat,]=apply(tmp_s,2,quantile_upper)
    
    if(if_plot==1){
      plot(range,Density_upper[id_pat,],type="l", lty = 2,col=color,
           main = paste0("Density estimation for sample point ", as.character(id_pat)),xlab='Log of Survival Time',ylab='Density Function')
      lines(range,Density_lower[id_pat,],lty = 2, col = color)
      lines(range,Density_estimation[id_pat,], col = color)
    }else if (if_plot==0.5){
      lines(range,Density_estimation[id_pat,], col = color)
      lines(range,Density_lower[id_pat,],lty = 2, col = color)
      lines(range,Density_upper[id_pat,],lty = 2, col = color)
    }
  }
  output<-NULL
  output$Density_estimation<-Density_estimation
  output$Density_lower_quant<-Density_lower
  output$Density_upper_quant<-Density_upper
  output$Survival_estimation<-Survival_estimation
  output$Survival_lower_quant<-Survival_lower
  output$Survival_upper_quant<-Survival_upper
  return(output)
}


DDPGP_density_Hyperparam <- function(mcmc,new_pat,range, if_plot = 0, quantiles=c(0.025,0.975),color='green')
{
  beta_0<-mcmc$beta_0
  covariate<-mcmc$covariate
  Npat<-mcmc$Npat
  iter_saved = length(mcmc$H)
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each sample point
  d[,1] = 1
  d[,2:Ncov] = covariate
  
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
  Density_lower<-Density_upper<-Density_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  
  tmp <- matrix(0, (iter_saved), length(range))
  for(id_pat in 1:N_new){
    dpat=c(1,new_pat[id_pat,])
    for (iter in 1:iter_saved)
    {
      M=mcmc$M[iter];r=mcmc$r[,iter];
      pat_cov=dpat;wh=mcmc$wh[,iter];
      muh=mcmc$muh[,,iter]; betah=mcmc$betah[,,iter];H=mcmc$H[iter]
      sigma2 = mcmc$sigma2[iter]
      theta<-mcmc$theta[,iter]
      sigma2_cov<-mcmc$sigma2_cov[iter]
      cov<-d
      covariance<-update_cov(Ncov, Npat, cov,theta,sigma2_cov)
      inv_covariance = chol2inv(chol(covariance))
      nj <- table(r)
      # ths <- as.numeric(names(nj)) # unique values
      # coeff = (t(cov)-pat_cov)
      # c1 = exp(-colSums(coeff^2*theta^2))
      # 
      new_cov<-rbind(pat_cov,cov)
      new_covar<-update_cov(Ncov, Npat+1,new_cov,theta,sigma2_cov)
      c1<-as.vector(new_covar[2:(Npat+1),1])
      
      ft=0
      for (h in 1:H)
      {
        mean=betah[h,]%*%pat_cov+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,])
        ft <- ft + nj[h]/(Npat+M)*dnorm((range), mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
      }
      ft<- ft + M/(Npat+M)*dnorm((range),mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      tmp[iter,] = ft
    }
    
    Density_estimation[id_pat,] = apply(tmp, 2, mean)
    Density_lower[id_pat,]=apply(tmp,2,quantile_lower)
    Density_upper[id_pat,]=apply(tmp,2,quantile_upper)
    if(if_plot==1){
      plot(range,Density_upper[id_pat,],type="l", lty = 2,col=color, ylim=c(min(Density_lower),max(Density_upper)),
           main = paste0("Density estimation for sample point ", as.character(id_pat)),xlab='Log of Survival Time',ylab='Density Function')
      lines(range,Density_lower[id_pat,],lty = 2, col = color)
      lines(range,Density_estimation[id_pat,], col = color)
    }else if (if_plot==0.5){
      lines(range,Density_estimation[id_pat,], col = color)
      lines(range,Density_lower[id_pat,],lty = 2, col = color)
      lines(range,Density_upper[id_pat,],lty = 2, col = color)
    }
  }
  output<-NULL
  output$Density_estimation<-Density_estimation
  output$Density_lower_quant<-Density_lower
  output$Density_upper_quant<-Density_upper
  return(output)
}





DDPGP_survival_Hyperparam <- function(mcmc,new_pat,range, if_plot = 0,quantiles=c(0.025,0.975),color='green')
{
  beta_0<-mcmc$beta_0
  covariate<-mcmc$covariate
  Npat<-mcmc$Npat
  iter_saved = length(mcmc$H)
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each sample point
  d[,1] = 1
  d[,2:Ncov] = covariate
  
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
  Survival_lower<-Survival_upper<- Survival_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  
  tmp <- matrix(0, (iter_saved), length(range))
  for(id_pat in 1:N_new){
    dpat=c(1,new_pat[id_pat,])
    for (iter in 1:iter_saved)
    {
      M=mcmc$M[iter];r=mcmc$r[,iter];
      pat_cov=dpat;wh=mcmc$wh[,iter];
      muh=mcmc$muh[,,iter]; betah=mcmc$betah[,,iter];H=mcmc$H[iter]
      sigma2 = mcmc$sigma2[iter]
      theta<-mcmc$theta[,iter]
      sigma2_cov<-mcmc$sigma2_cov[iter]
      cov<-d
      covariance<-update_cov(Ncov, Npat, cov,theta,sigma2_cov)
      inv_covariance = chol2inv(chol(covariance))
      nj <- table(r)
      # ths <- as.numeric(names(nj)) # unique values
      # coeff = (t(cov)-pat_cov)
      # c1 = exp(-colSums(coeff^2*theta^2))
      # 
      new_cov<-rbind(pat_cov,cov)
      new_covar<-update_cov(Ncov, Npat+1,new_cov,theta,sigma2_cov)
      c1<-as.vector(new_covar[2:(Npat+1),1])
      ft=0
      for (h in 1:H)
      {
        mean=betah[h,]%*%pat_cov+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,])
        ft <- ft + nj[h]/(Npat+M)*pnorm(range, mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
      }
      ft<- ft + M/(Npat+M)*pnorm(range,mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      tmp[iter,] = 1-ft
    }
    Survival_estimation[id_pat,] = apply(tmp, 2, mean)
    Survival_lower[id_pat,]=apply(tmp,2,quantile_lower)
    Survival_upper[id_pat,]=apply(tmp,2,quantile_upper)
    if(if_plot==1){
      plot(range,Survival_upper[id_pat,],type="l", lty = 2, col=color, ylim=c(min(Survival_lower),max(Survival_upper)),
           main = paste0("Survival estimation for sample point ", as.character(id_pat)),xlab='Log of Survival Time',ylab='Survival Function')
      lines(range,Survival_lower[id_pat,],lty = 2, col = color)
      lines(range,Survival_estimation[id_pat,], col = color)
    }else if (if_plot==0.5){
      lines(range,Survival_estimation[id_pat,], col = color)
      lines(range,Survival_lower[id_pat,],lty = 2, col = color)
      lines(range,Survival_upper[id_pat,],lty = 2, col = color)
    }
  }
  output<-NULL
  output$Survival_estimation<-Survival_estimation
  output$Survival_lower_quant<-Survival_lower
  output$Survival_upper_quant<-Survival_upper
  return(output)
}



DDPGP_hazard_Hyperparam <- function(mcmc,new_pat,range, if_plot = 0,quantiles=c(0.025,0.975),color='green')
{
  beta_0<-mcmc$beta_0
  covariate<-mcmc$covariate
  iter_saved = length(mcmc$H)
  Npat<-mcmc$Npat
  if (length(covariate)==Npat) Ncov <- 2 else
    Ncov <- dim(covariate)[2]+1
  d <- matrix(0, Npat, Ncov)  #Matrix of all covariates of each sample point
  d[,1] = 1
  d[,2:Ncov] = covariate
  
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
  Hazard_lower<-Hazard_upper<-Hazard_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  
  tmp <- matrix(0, (iter_saved), length(range))
  for(id_pat in 1:N_new){
    dpat=c(1,new_pat[id_pat,])
    for (iter in 1:iter_saved)
    {
      M=mcmc$M[iter];r=mcmc$r[,iter];
      pat_cov=dpat;wh=mcmc$wh[,iter];
      muh=mcmc$muh[,,iter]; betah=mcmc$betah[,,iter];H=mcmc$H[iter]
      sigma2 = mcmc$sigma2[iter]
      theta<-mcmc$theta[,iter]
      sigma2_cov<-mcmc$sigma2_cov[iter]
      cov<-d
      covariance<-update_cov(Ncov, Npat, cov,theta,sigma2_cov)
      inv_covariance = chol2inv(chol(covariance))
      nj <- table(r)
      # ths <- as.numeric(names(nj)) # unique values
      # coeff = (t(cov)-pat_cov)
      # c1 = exp(-colSums(coeff^2*theta^2))
      # 
      new_cov<-rbind(pat_cov,cov)
      new_covar<-update_cov(Ncov, Npat+1,new_cov,theta,sigma2_cov)
      c1<-as.vector(new_covar[2:(Npat+1),1])
      ft=0
      St=0
      for (h in 1:H)
      {
        mean=betah[h,]%*%pat_cov+c1%*%inv_covariance%*%(muh[h,]-cov%*%betah[h,])
        ft <- ft + nj[h]/(Npat+M)*dnorm(range, mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
        St <- St + nj[h]/(Npat+M)*pnorm(range, mean=mean, sd=sqrt(sigma2+sigma2_cov^2+0.01-c1%*%inv_covariance%*%c1))
        
      }
      ft<- ft + M/(Npat+M)*dnorm((range),mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      St<- St + M/(Npat+M)*pnorm((range),mean=beta_0%*%pat_cov, sd=sqrt(sigma2+sigma2_cov^2+0.01+pat_cov%*%diag(Ncov)%*%pat_cov))
      
      tmp[iter,] = ft/(1-St)
    }
    Hazard_estimation[id_pat,] = apply(tmp, 2, mean)
    Hazard_lower[id_pat,]=apply(tmp,2,quantile_lower)
    Hazard_upper[id_pat,]=apply(tmp,2,quantile_upper)
    if(if_plot==1){
      plot(range,Hazard_upper[id_pat,],type="l", lty = 2, col=color,ylim=c(min(Hazard_lower),max(Hazard_upper)),
           main = paste0("Hazard estimation for sample point ", as.character(id_pat)),xlab='Log of Survival Time',ylab='Hazard Function')
      lines(range,Hazard_lower[id_pat,],lty = 2, col = color)
      lines(range,Hazard_estimation[id_pat,], col = color)
    }else if (if_plot==0.5){
      lines(range,Hazard_estimation[id_pat,], col = color)
      lines(range,Hazard_lower[id_pat,],lty = 2, col = color)
      lines(range,Hazard_upper[id_pat,],lty = 2, col = color)
    }
  }
  output<-NULL
  output$Hazard_estimation<-Hazard_estimation
  output$Hazard_lower_quant<-Hazard_lower
  output$Hazard_upper_quant<-Hazard_upper
  return(output)
}


true_density <- function(data,new_pat,range, if_plot = 0.5,col='black')
{ 
  prob<-data$prob
  #If the new_pat is a vetor, change it to a matrix
  if(length(dim(new_pat))==0){
    new_pat <- t(new_pat)
  }
  sigma<-data$sigma
  shape<-data$shape
  #build a matrix to record the density for new sample points.
  N_new <- dim(new_pat)[1]
  Density_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  beta1<-data$beta1
  beta2<-data$beta2
  for(id_pat in 1:N_new){
    if (data$truth=='lognorm'){
    meanlog1<- sum(beta1*new_pat[id_pat,])
    meanlog2<- sum(beta2*new_pat[id_pat,])
    y=(prob)*dnorm(range , meanlog1, sigma)+(1-prob)*dnorm(range , meanlog2, sigma)
    # y=(prob)*dlnorm(exp(range) , meanlog1, sigma)+(1-prob)*dlnorm(exp(range) , meanlog2, sigma)
    
    }else if (data$truth=='weib'){
      scale1<- exp(sum(beta1*new_pat[id_pat,]))
      scale2<-exp(sum(beta2*new_pat[id_pat,]))
      y<-prob*(exp(range)*dweibull(exp(range),shape=shape,scale=scale1))+(1-prob)*(exp(range)*dweibull(exp(range),shape=shape,scale=scale2))
    }
    if(if_plot==1){
      plot(range,y,col=col,type='l',main='True Density Function')
    }
    if(if_plot==0.5){
      lines(range,y,col=col)
    }
    Density_estimation[id_pat,] <- y
  }
  return(Density_estimation)
}

true_survival <- function(data,new_pat,range, if_plot = 0.5,col='black')
{ 
  prob<-data$prob
  #If the new_pat is a vetor, change it to a matrix
  if(length(dim(new_pat))==0){
    new_pat <- t(new_pat)
  }
  sigma<-data$sigma
  shape<-data$shape
  #build a matrix to record the density for new sample points.
  N_new <- dim(new_pat)[1]
  Survival_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  beta1<-data$beta1
  beta2<-data$beta2
  for(id_pat in 1:N_new){
    if (data$truth=='lognorm'){
      meanlog1<- sum(beta1*new_pat[id_pat,])
      meanlog2<- sum(beta2*new_pat[id_pat,])
      y=(prob)*(1-pnorm(range , meanlog1, sigma))+(1-prob)*(1-pnorm(range , meanlog2, sigma))
    }else if (data$truth=='weib'){
      scale1<- exp(sum(beta1*new_pat[id_pat,]))
      scale2<-exp(sum(beta2*new_pat[id_pat,]))
      #log_range <- log(range + 1)
      y<-prob*(1-pweibull(exp(range),shape=shape,scale=scale1))+(1-prob)*(1-pweibull(exp(range),shape=shape,scale=scale2))
    }
    if(if_plot==1){
      plot(range,y,col=col,type='l',lwd = 3,main='True Survival Function vs Estimations')
    }
    if(if_plot==0.5){
      lines(range,y,col=col)
    }
    Survival_estimation[id_pat,] <- y
  }
  return(Survival_estimation)
}


true_hazard <- function(data,new_pat,range, if_plot = 0.5,col='black')
{ 
  prob<-data$prob
  #If the new_pat is a vetor, change it to a matrix
  if(length(dim(new_pat))==0){
    new_pat <- t(new_pat)
  }
  sigma<-data$sigma
  shape<-data$shape
  #build a matrix to record the density for new sample points.
  N_new <- dim(new_pat)[1]
  Hazard_estimation <- matrix(NA, nrow = N_new, ncol = length(range))
  beta1<-data$beta1
  beta2<-data$beta2
  for(id_pat in 1:N_new){
    if (data$truth=='lognorm'){
      meanlog1<- sum(beta1*new_pat[id_pat,])
      meanlog2<- sum(beta2*new_pat[id_pat,])
      y=(prob)*dnorm(range , meanlog1, sigma)+(1-prob)*dnorm(range , meanlog2, sigma)
      z=(prob)*(1-pnorm(range , meanlog1, sigma))+(1-prob)*(1-pnorm(range , meanlog2, sigma))
      
      plot_y=(prob)*dnorm(range , meanlog1, sigma)/(1-pnorm(range , meanlog1, sigma))+(1-prob)*dnorm(range , meanlog2, sigma)/(1-pnorm(range , meanlog2, sigma))
    }else if (data$truth=='weib'){
      scale1<- exp(sum(beta1*new_pat[id_pat,]))
      scale2<-exp(sum(beta2*new_pat[id_pat,]))
      y<-prob*(dweibull(range,shape=shape,scale=scale1))+(1-prob)*(dweibull(range,shape=shape,scale=scale2))
      z<-prob*(1-pweibull(range,shape=shape,scale=scale1))+(1-prob)*(1-pweibull(range,shape=shape,scale=scale2))
      
    }
    plot_y<-y/z
    if(if_plot==1){
      plot(range,plot_y,col=col,type='l',main='True Density Function')
    }
    if(if_plot==0.5){
      lines(range,plot_y,col=col)
    }
    Hazard_estimation[id_pat,] <- plot_y
  }
  return(Hazard_estimation)
}
