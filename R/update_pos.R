init <- function(response, Npat,H)
{
	hc = hclust(dist(response)^2, "cen")
	r = cutree(hc, k=H)
    wh1 <- table(r)/Npat
    idx <- order(wh1,decreasing=T)
    wh <- wh1[idx]
    mh = sapply(split(response,r),mean)
    muh = matrix(mh,H, Npat)
    M = 1
    sigma2 = 0.1
    return(list(muh=muh, wh=wh, M=M,sigma2=sigma2, r=r))
}

update_r<- function(sigma2, muh, Npat, censor_status, response, r, covariance, cov,
                   H, M , ns,beta_0,Ncov)
{
  impute_y = response
  #generated for new cluster
  #set.seed(1)
  A <- cov%*%t(cov)+covariance
  A_inv<-chol2inv(chol(A))
  muh_temp <- muh
  
  test_r=update_r_within_testing(sigma2, Npat, censor_status, response, r, covariance, cov,
                                 H, M , ns,impute_y,
                                 muh_temp,beta_0,Ncov,sample,rtrunc,rmvn,A_inv,A)
  if (dim(test_r$muh)[1]>test_r$H){
    test_r$muh<-test_r$muh[1:test_r$H,]
  }
  return(list(r=test_r$r,impute_y=test_r$impute_y, H = test_r$H, muh = test_r$muh))
}

#sigma2=mcmc$sigma2[iter-1];muh=matrix(mcmc$muh[1:mcmc$H[iter-1],,iter-1],nrow=mcmc$H[iter-1])
#r=mcmc$r[,iter-1];H=mcmc$H[iter-1];M=mcmc$M[iter-1];ns=mcmc$ns[1:mcmc$H[iter-1],iter-1]
# update_r_test<- function(sigma2, muh, Npat, censor_status, response, r, covariance, cov,
#                     H, M, ns)
# {
#   A <- cov%*%t(cov)+covariance
#   A_inv<-chol2inv(chol(A))
#   impute_y = response
#   beta0 <<- rep(0,Ncov)
#   for (i in 1:Npat)
#   {
#     ns[r[i]] = ns[r[i]]-1
#     if (ns[r[i]]==0)
#     {
#       H = H-1
#       muh = muh[-r[i],]
#       ns = ns[-r[i]]
#       for (id in 1:Npat)
#       {
#         if (r[id]>r[i]) r[id] = r[id] - 1
#       }
#     }
#     mu = muh[,i]
#     if (censor_status[i]==1)
#     {
#       ph = dnorm(response[i], m=mu, sd=sqrt(sigma2))
#       ph[H+1] <- dnorm(response[i], m= 0, sd=sqrt(sigma2+A[i,i]))
#       ph <- ph*c(ns,M)
#       r[i] = sample(1:(H+1), 1, prob=ph)
#     }
#     if (censor_status[i]==0)
#     {
#       ph = 1-pnorm(response[i], m=mu, sd=sqrt(sigma2))
#       ph[H+1] <- 1-pnorm(response[i], m= 0, sd=sqrt(sigma2+A[i,i]))
#       ph <- ph*c(ns,M)
#       r[i] = sample(1:(H+1), 1, prob=ph)
#       if (r[i] == (H+1)) impute_y[i] = rtrunc(rnorm, 1, linf=response[i], lsup=Inf, mean = 0, sd = sqrt(sigma2+A[i,i])) else 
#         impute_y[i] = rtrunc(rnorm, 1, linf=response[i], lsup=Inf, mean = mu[r[i]], sd = sqrt(sigma2)) 
#     }
#     if (r[i]==(H+1))
#     {
#       H = H+1
#       ns[H] <- 1
#       t2 = A_inv
#       t2[i,i] = t2[i,i] + 1/sigma2
#       t2_inv = chol2inv(chol(t2))
#       t1 = t2_inv[,i]*impute_y[i]
#       muh_candidate <- rmvn(1, t1, t2_inv)
#       muh = rbind (muh, muh_candidate)
#     } else
#     {
#       ns[r[i]] = ns[r[i]] + 1
#     }
#   }
#   return(list(r=r,impute_y=impute_y, H = H, muh = muh))
# }




update_muh <- function(r, sigma2, response, Npat, covariance, cov, betah, inv_covariance, H)
{
  muh <- matrix(0,H, Npat)     # initialize
  for(h in 1:H){
    # some data assigned to h-th pointmass
    Sh <- which(r==h)
    nh <- length(Sh)
    W = matrix(0, nh, Npat)
    for (i in 1:nh)
    {
      W[i, Sh[i]] = 1
    }
    #W_ji=1 if patient i is in the jth cluster
    var = chol2inv(chol(inv_covariance+1/sigma2*t(W)%*%W))
    if (H==1){
      sel_betah<-betah
    } else{
      sel_betah<-betah[h,]
    }
    mu = var%*%(t(W)%*%response[Sh]/sigma2 + inv_covariance%*%cov%*%sel_betah)
    muh[h,] <- c(rmvn(1, mu, var))
    # no data assinged to h-th pointmass# sample from base measure
    #mu <- cov%*%betah[h,]
    #muh[h,] <- mvrnorm(1, mu, covariance)
  }
  
  return(muh)
}

update_betah <- function(Ncov, muh, cov, covariance, r, inv_covariance, H,beta_0)
{
  betah <- matrix(0, H, Ncov)     # initialize
  var = solve(t(cov)%*%inv_covariance%*%cov+diag(Ncov))
  for(h in 1:H){     # some data assigned to h-th pointmass
    if (H==1){
      sel_muh<-muh
    } else{
      sel_muh<-muh[h,]
    }
    mu = var%*%(t(cov)%*%inv_covariance%*%sel_muh + diag(Ncov)%*%beta_0)
    betah[h,] <- c(rmvn(1, mu, var))
    # no data assinged to h-th pointmass     # sample from base measure
    #mu <- rep(0, Ncov)
    #betah[h,] <- mvrnorm(1, mu, diag(Ncov))}
  }
  return(betah)
}


update_sigma2 <- function(r, response, muh, Npat, H, delta1, delta2)
{
  tmp = rep(0, Npat)
  for (h in 1:H)
  {
    Sh = which(r==h)
    if (H==1){
      tmp[Sh]=muh[Sh]
    } else{
      tmp[Sh] = muh[h,Sh]
    }
  }
  inverse_sigma2 <- rgamma(1, delta1+Npat/2, delta2+sum((response-tmp)^2)/2)
  sigma2 <- 1/inverse_sigma2
  return(sigma2)
}

#r<-mcmc$r[,iter]; M=mcmc$M[iter-1]
# update_wh_and_M <- function(r, M)
# {
#   ## returns: wh
#   vh <- rep(0,H)  # initialize
#   wh <- rep(0,H)
#   V <-  1         # record prod_{g<h} (1-vh_h)
#   for(h in 1:(H-1)){
#     Ah <- which(r==h)
#     Bh <- which(r>h)
#     vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
#     wh[h] <- vh[h]*V
#     V <- V*(1-vh[h])
#   }
#   vh[H] <- 1.0
#   wh[H] <- V
#   #M = rgamma(1, H+lambda1-1, lambda2-sum(log(1-vh[1:(H-1)])) )
#   tmp1 = log(1-vh[1:(H-1)])
#   M <- rgamma(1, lambda1+H-1, lambda2-sum( ifelse(tmp1<(-708.3964), -708.3964, tmp1)  )  )
#   return(list(wh=wh, M=M))
# }

# update_wh <- function(r, M,H)
# {
#   ## returns: wh
#   vh <- rep(0,H)  # initialize
#   wh <- rep(0,H)
#   V <-  1         # record prod_{g<h} (1-vh_h)
#   for(h in 1:(H-1)){
#     Ah <- which(r==h)
#     Bh <- which(r>h)
#     vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
#     wh[h] <- vh[h]*V
#     V <- V*(1-vh[h])
#   }
#   vh[H] <- 1.0
#   wh[H] <- V
#   return(wh)
# }

update_M <- function(Mold,H,lambda1,lambda2,Npat)
{
	m = rbeta(1, (Mold+1), Npat)
	ratio = (lambda1+H-1)/(Npat*(lambda2-log(m)))
	prop = ratio/(1+ratio)
	u = runif(1)
	if (u<prop){Mnew = rgamma(1, lambda1+H, lambda2-log(m))} else{
		Mnew = rgamma(1, lambda1+H-1, lambda2-log(m))}
	return(Mnew)
}

#ns=number of patients in each cluster
update_ns<-function(r,H){
  ns <- rep(NA,H)
  for(h in 1:H){
    ns[h] <- sum(r==h)
  }
  return(ns)
}
