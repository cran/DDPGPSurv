#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R.h>
#include <iostream>
#include <assert.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends("mvnfast")]]
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void good(double a){
  Rcout << "Good " << a <<std::endl; // for debugging purposes
}

// [[Rcpp::export]]
void showD(double b, double a){
  Rcout << "Value " <<a <<" is " << std::endl << b <<std::endl;  //for debugging
}

// [[Rcpp::export]]
void showRowvec(arma::rowvec b){
  Rcout << "Rowvec is " << std::endl << b <<std::endl;  //for debugging
}


// [[Rcpp::export]]
void showMatrix(arma::mat b){
  Rcout << "Matrix is " << std::endl << b <<std::endl; //for debugging
}

// [[Rcpp::export]]
void stab(arma::rowvec b){
}

// [[Rcpp::export]]
List update_r_within_testing(double sigma2,
                             int Npat, arma::rowvec censor_status, arma::rowvec response,
                             arma::rowvec r, arma::mat covariance, arma::mat cov,
                             int H, double M , arma::rowvec ns, arma::rowvec impute_y,
                             arma::mat muh_temp,arma::rowvec beta_0,int Ncov,Function samp,
                             Function rtrunc,Function rmvn,arma::mat A_inv,arma::mat A){
  List ans;
  //Rcpp::Environment package_env("package:mvnfast");
  //Rcpp::Function rmvn=package_env["rmvn"];
  for(int i=0; i<Npat;i++){
    ns[r[i]-1]=ns[r[i]-1]-1;
    if(ns[r[i]-1]==0){
      H= H-1;
      muh_temp.shed_row(r[i]-1);
      ns.shed_col(r[i]-1);
      for(int id=0; id<Npat;id++){
        if(r[id]>r[i]){
          r[id]= r[id]-1;
        }
      }
    }
    arma::colvec mu= muh_temp.col(i);
    arma::rowvec ph(H+1);
    if (censor_status[i]==1)
    {
      for( int k=0; k<H; k++) {
        ph[k] = R::dnorm(response[i], mu[k], sqrt(sigma2),0)*ns[k];
        // showD(ph[k],1);
      }
      ph[H]= R::dnorm(response[i], sum(beta_0%cov.row(i)), sqrt(sigma2+covariance(i,i)+as_scalar(cov.row(i)*arma::eye(Ncov,Ncov)*cov.row(i).t())),0)*M;
      // good(1);
      ph=ph/sum(ph);
      r[i]= as<double>(samp(arma::linspace<arma::uvec>(1,H+1,H+1), 1,0, ph));
    }
    if (censor_status[i] == 0)
    {
      for( int k=0; k<H; k++) {
        ph[k] =(1-R::pnorm(response[i],mu[k],sqrt(sigma2),1,0))*ns[k];
      }
      ph[H]= (1-R::pnorm(response[i],sum(beta_0%cov.row(i)), sqrt(sigma2+covariance(i,i)+as_scalar(cov.row(i)*arma::eye(Ncov,Ncov)*cov.row(i).t())),1,0))*M;
      //stab(ph);
      // good(2);
      // showRowvec(ph);
      ph=ph/sum(ph);
      r[i]= as<double>(samp(arma::linspace<arma::uvec>(1,H+1,H+1), 1,0, ph));
      if (r[i]==(H+1)){
        impute_y[i] = as<double>(rtrunc("rnorm", 1,response[i], arma::datum::inf,0, sqrt(sigma2+A(i,i))));
      } else {
      impute_y[i] = as<double>(rtrunc("rnorm", 1,response[i], arma::datum::inf, mu[r[i]-1], sqrt(sigma2)));
      }
      // if ((response[i]-mu[r[i]-1])>3*sqrt(sigma2)){
      //   impute_y[i]=response[i];
      // }
      // if ((response[i]-mu[r[i]-1])<=3*sqrt(sigma2)){
      //   impute_y[i] = as<double>(rtrunc("rnorm", 1,response[i], arma::datum::inf, mu[r[i]-1], sqrt(sigma2)));
      //   
      // }
      // if (impute_y[i]>pow(10,99)){
      //   impute_y[i]=response[i];
      // }
      //showRowvec(impute_y);
      
    }
    if(r[i]==(H+1)){
      H= H+1;
      ns.resize(ns.size()+1);
      ns[H-1]=1;
      arma::mat t1= (A_inv*cov*beta_0.t());
      t1(i)=t1(i)+response[i]/sigma2;
      arma::mat t2=A_inv;
      t2(i,i)=t2(i,i)+1/sigma2;
      arma::mat t2_chol_inv=chol(t2).i();
      // t2=t2.i();
      t2=t2_chol_inv*t2_chol_inv.t();
      arma::mat mu_temp_c= t2*t1;
      arma::colvec mu_temp=mu_temp_c.col(0);
      arma::rowvec muh_candidate= as<arma::rowvec>(rmvn(1,mu_temp,t2));
      muh_temp= join_cols(muh_temp, muh_candidate);

    }
    else{
      ns[r[i]-1]=ns[r[i]-1]+1;
    }
  }
  // good(3);
  // showRowvec(impute_y);
  ans["ns"] = ns;
  ans["r"]=r;
  ans["impute_y"]=impute_y;
  ans["H"]=H;
  ans["muh"]=muh_temp;
  return(ans);
}


// // [[Rcpp::export]]
// double mean_surviv(SEXP r, arma::mat cov, arma::rowvec pat_cov,int H,arma::mat betah, List mcmc,
//                    arma::mat inv_covariance, arma::mat muh, int Npat, double M,int Ncov, double sigma2,
//                    arma::rowvec beta_0,int nAUC, int Niter, int burn_in,int lag, arma::rowvec range_AUC){
//   Environment myEnv = Environment::global_env();
//   Function mySum = myEnv["mcmc_extract"];
//   arma::mat dpat(nAUC,Ncov,arma::fill::ones);
//   dpat.row(1)=range_AUC;
//   arma::mat meansurvival_AUC((Niter-burn_in)/lag,nAUC);
//   for(int k=0; k<nAUC;k++){
//   good(k);
//     for(int i=0; i<H;i++){
//       int iter = burn_in + i*lag;
//       arma::rowvec pat_cov=dpat.row(k);
//       List mcmc_it=mcmc_extract(mcmc,iter);
//
//     }
//   }
//   return(2);
// }


// [[Rcpp::export]]
arma::rowvec fmeansurviv(arma::rowvec x, SEXP r, arma::mat cov, arma::rowvec pat_cov,int H,arma::mat betah,
                         arma::mat inv_covariance, arma::mat muh, int Npat, double M,int Ncov, double sigma2,
                         arma::rowvec beta_0){
  
  IntegerVector nj=table(as<NumericVector>(r));
  arma::rowvec a=as<arma::rowvec>(nj)/2;
  arma::mat coeff=cov.each_row()-pat_cov;
  arma::rowvec c1=exp(-sum(square(coeff.t()),0));
  arma::rowvec Ft;
  for(int h=0; h<H;h++){
    double norm_mean=as_scalar(betah.row(h)*pat_cov.t()+c1*inv_covariance*(muh.row(h)-cov*betah.row(h)));
    arma::rowvec pvec;
    for(int j=0; j<Npat; j++) {
      pvec[j]=R::pnorm(log(x[j]),norm_mean,sqrt(sigma2),1,0);
    }
    Ft=Ft+nj[h]/(Npat+M)*pvec;
  }
  arma::rowvec pvec1;
  for(int l=0; l<Npat; l++) {
    pvec1[l]=R::pnorm(log(x[l]),as_scalar(beta_0*pat_cov.t()),sqrt(sigma2+as_scalar(pat_cov*arma::eye(Ncov,Ncov)*pat_cov.t())),1,0);
  }
  Ft=Ft+M/(Npat+M)*pvec1;
  arma::rowvec St=1-Ft;
  return(St);
}




// [[Rcpp::export]]
IntegerVector RcppTable(SEXP x) {
  switch (TYPEOF(x)) {
  case INTSXP: return table(as<IntegerVector>(x));
  case REALSXP: return table(as<NumericVector>(x));
  case STRSXP: return table(as<CharacterVector>(x));
  case LGLSXP: return table(as<LogicalVector>(x));
  default: {
    stop("untested SEXP type");
    return R_NilValue;
  }
  }
}
