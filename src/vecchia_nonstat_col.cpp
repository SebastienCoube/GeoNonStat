#include <RcppArmadillo.h>
#include <R.h>
#include <iostream>
#include <vector>
#include <cmath>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

//#include <boost/math/special_functions/bessel.hpp>
//#include <boost/math/special_functions/gamma.hpp>

//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::cube nonstat_vecchia_Linv_col(
    arma::mat log_range,
    arma::mat locs,
    arma::mat NNarray, 
    int num_threads,
    bool compute_derivative){
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  // data dimensions
  int n = locs.n_cols;
  int m = NNarray.n_rows;
  arma::cube result(m, 1 + log_range.n_rows * m * compute_derivative, n);
  result(0,0,0) = 1;
  
  
  // pre-compute expensive matrix exponentials
  arma::mat exprange(log_range.n_rows + log_range.n_rows * log_range.n_rows * compute_derivative, n);
  // loop over every observation  
#pragma omp parallel for
  for(int i=0; i<n; i++){
    arma::mat lograngemat(2,2);
    arma::mat rangemat(2,2);
    arma::vec log_range_col = log_range.col(i); 
    
    
    lograngemat(0,0) = log_range_col(0);
    lograngemat(1,1) = log_range_col(1);
    lograngemat(1,0) = log_range_col(2)/sqrt(2);
    lograngemat(0,1) = lograngemat(1,0);
    rangemat = expmat_sym(lograngemat);
    
    exprange(0, i) = rangemat(0,0); 
    exprange(1, i) = rangemat(1,1); 
    exprange(2, i) = rangemat(0,1); 
    if(compute_derivative){
      lograngemat(0,0) = lograngemat(0,0) + .001;
      rangemat = expmat_sym(lograngemat);
      exprange(3, i) = rangemat(0,0); 
      exprange(4, i) = rangemat(1,1); 
      exprange(5, i) = rangemat(0,1); 
      lograngemat(0,0) = lograngemat(0,0) - .001;
      lograngemat(1,1) = lograngemat(1,1) + .001;
      rangemat = expmat_sym(lograngemat);
      exprange(6, i) = rangemat(0,0); 
      exprange(7, i) = rangemat(1,1); 
      exprange(8, i) = rangemat(0,1); 
      lograngemat(1,1) = lograngemat(1,1) - .001;
      lograngemat(0,1) = lograngemat(0,1) + .001;
      lograngemat(1,0) = lograngemat(0,1);
      rangemat = expmat_sym(lograngemat);
      exprange(9, i) = rangemat(0,0); 
      exprange(10, i) = rangemat(1,1); 
      exprange(11, i) = rangemat(0,1); 
    }
  }

#pragma omp parallel for
  for(int i=1; i<n; i++){
    // Fill subsets of locations, range, and pre-compute expensive matrix exponentials
    int bsize = std::min(i+1,m);
    arma::mat locsub(2, bsize);
    arma::mat rangesub(3 + 9 * compute_derivative, bsize);
    for(int j=0; j<bsize; j++){
      int idx = NNarray(j,i);
      // subsets of locations
      locsub.col(j) = locs.col(idx - 1); 
      // subsets of range
      rangesub.col(j) = exprange.col(idx - 1); 
    }
    // initialize recipient for result
    arma::mat out(m, 1 + log_range.n_rows * m * compute_derivative);
    
    // computing covariance 
    double mahala_dist;
    arma::mat  sigma11(bsize-1, bsize-1);
    arma::vec  sigma12(bsize-1);
    arma::mat hybrid_range_inv(2,2) ;
    // anisotropic case
    if(log_range.n_rows == 3){
      for(int j1=1; j1<bsize; j1++){
        //filling sigma11
        for(int j2=1; j2<j1; j2++){
          
          hybrid_range_inv(0,0) = (rangesub(0, j1) + rangesub(0, j2))*.5;
          hybrid_range_inv(1,0) = (rangesub(2, j1) + rangesub(2, j2))*.5;
          hybrid_range_inv(0,1) = hybrid_range_inv(1,0);
          hybrid_range_inv(1,1) = (rangesub(1, j1) + rangesub(1, j2))*.5;
          hybrid_range_inv = inv_sympd(hybrid_range_inv);
          double mahala_dist = pow(
            pow( (locsub(0, j1)-locsub(0, j2)),2) * hybrid_range_inv(0,0)
            +pow((locsub(1, j1)-locsub(1, j2)),2) * hybrid_range_inv(1,1)
            + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * hybrid_range_inv(1,0), 
            .5
          );
          sigma11 (j1-1, j2-1) = 
            pow(
              (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                (rangesub(0, j2)*rangesub(1, j2) - rangesub(2, j2)*rangesub(2, j2)), 
                .25) * //determinants
                pow(arma::det(hybrid_range_inv), .5) * 
                exp(-mahala_dist) * (1+ mahala_dist);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
          
        }
        sigma11 (j1-1, j1-1) = 1.0001;
        //filling sigma12
        hybrid_range_inv(0,0) = (rangesub(0, j1) + rangesub(0, 0))*.5;
        hybrid_range_inv(1,0) = (rangesub(2, j1) + rangesub(2, 0))*.5;
        hybrid_range_inv(0,1) = hybrid_range_inv(1,0);
        hybrid_range_inv(1,1) = (rangesub(1, j1) + rangesub(1, 0))*.5;
        hybrid_range_inv = inv_sympd(hybrid_range_inv);
        mahala_dist = pow(
          pow( (locsub(0, j1)-locsub(0, 0)),2) * hybrid_range_inv(0,0)
          +pow((locsub(1, j1)-locsub(1, 0)),2) * hybrid_range_inv(1,1)
          + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * hybrid_range_inv(1,0), 
          .5
        );
        sigma12 (j1-1) = 
          pow(
            (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
              (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0)), 
              .25) * //determinants
              pow(arma::det(hybrid_range_inv), .5) * 
              exp(-mahala_dist) * (1+ mahala_dist);
      }
    }
    // solving sigma11
    arma::mat agmis11  = arma::inv_sympd(sigma11) ;
    // computing a vector used everyvhere
    arma::mat salt  = agmis11 * sigma12 ;
    //computing Vecchia approx itself
    double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
    double pow_invcondsd_3  = pow(inverse_cond_sd, 3);
    out(0, 0) = inverse_cond_sd ;
    for(int j=1; j<bsize; j++){
      out(j,0) = - salt(j-1) * inverse_cond_sd;
    }
    
    // computing derivative of Vecchia approx 
    if(compute_derivative){
      // derivative of sigma 11 wrt to parent range parameters
      arma::mat dsigma11(bsize-1, bsize-1);
      // derivative of sigma 12 wrt to parent range parameters
      arma::vec dpsigma12(bsize-1);
      // derivative of sigma 12 wrt to child range parameters
      arma::vec dcsigma12(bsize-1);
      // result of multiplications
      arma::vec salt_dsigma11(bsize-1);
      arma::vec the_rest(bsize-1);
      arma::vec salt_dsigma11_agmis11(bsize-1);
      for(int d_idx=0; d_idx<3; d_idx++){
        // computing derivative of covariance
        // anisotropic case
        if(log_range.n_rows == 3){
          for(int j1=1; j1<bsize; j1++){
            //filling dsigma11
            for(int j2=1; j2<bsize; j2++){
              if(j1!=j2){
                hybrid_range_inv(0,0) = (rangesub((d_idx+1) * 3 + 0, j1) + rangesub((d_idx+1) * 3 + 0, j2))*.5;
                hybrid_range_inv(1,1) = (rangesub((d_idx+1) * 3 + 1, j1) + rangesub((d_idx+1) * 3 + 1, j2))*.5;
                hybrid_range_inv(1,0) = (rangesub((d_idx+1) * 3 + 2, j1) + rangesub((d_idx+1) * 3 + 2, j2))*.5;
                hybrid_range_inv(0,1) = hybrid_range_inv(1,0);
                hybrid_range_inv = inv_sympd(hybrid_range_inv);
                mahala_dist = pow(
                  pow( (locsub(0, j1)-locsub(0, j2)),2) * hybrid_range_inv(0,0)
                  +pow((locsub(1, j1)-locsub(1, j2)),2) * hybrid_range_inv(1,1)
                  + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * hybrid_range_inv(1,0), 
                  .5
                );
                dsigma11 (j1-1, j2-1) = 
                  pow(
                    (rangesub(d_idx * 3 + 0, j1)*rangesub(d_idx * 3 + 1, j1) - rangesub(d_idx * 3 + 2, j1)*rangesub(d_idx * 3 + 2, j1)) * 
                      (rangesub(d_idx * 3 + 0, j2)*rangesub(d_idx * 3 + 1, j2) - rangesub(d_idx * 3 + 2, j2)*rangesub(d_idx * 3 + 2, j2)), 
                      .25) * //determinants
                      pow(arma::det(hybrid_range_inv), .5) * 
                      exp(-mahala_dist) * (1+ mahala_dist);
              }
            }
            dsigma11 (j1-1, j1-1) = 1.0001;
            //filling dpsigma12
            hybrid_range_inv(0,0) = (rangesub((d_idx + 1) * 3 + 0, j1) + rangesub(0, 0))*.5;
            hybrid_range_inv(1,1) = (rangesub((d_idx + 1) * 3 + 1, j1) + rangesub(1, 0))*.5;
            hybrid_range_inv(1,0) = (rangesub((d_idx + 1) * 3 + 2, j1) + rangesub(2, 0))*.5;
            hybrid_range_inv(0,1) = hybrid_range_inv(1,0);
            hybrid_range_inv = inv_sympd(hybrid_range_inv);
            mahala_dist = pow(
              pow( (locsub(0, j1)-locsub(0, 0)),2) * hybrid_range_inv(0,0)
              +pow((locsub(1, j1)-locsub(1, 0)),2) * hybrid_range_inv(1,1)
              + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * hybrid_range_inv(1,0), 
              .5
            );
            dpsigma12 (j1-1) = 
              pow(
                (rangesub((d_idx + 1) * 3 + 0, j1)*rangesub((d_idx + 1) * 3 + 1, j1) - rangesub((d_idx + 1) * 3 + 2, j1)*rangesub((d_idx + 1) * 3 + 2, j1)) * 
                  (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0)), 
                  .25) * //determinants
                  pow(arma::det(hybrid_range_inv), .5) * 
                  exp(-mahala_dist) * (1+ mahala_dist);
            //filling dcsigma12
            hybrid_range_inv(0,0) = (rangesub(0, j1) + rangesub((d_idx + 1) * 3 + 0, 0))*.5;
            hybrid_range_inv(1,1) = (rangesub(1, j1) + rangesub((d_idx + 1) * 3 + 1, 0))*.5;
            hybrid_range_inv(1,0) = (rangesub(2, j1) + rangesub((d_idx + 1) * 3 + 2, 0))*.5;
            hybrid_range_inv(0,1) = hybrid_range_inv(1,0);
            hybrid_range_inv = inv_sympd(hybrid_range_inv);
            mahala_dist = pow(
              pow( (locsub(0, j1)-locsub(0, 0)),2) * hybrid_range_inv(0,0)
              +pow((locsub(1, j1)-locsub(1, 0)),2) * hybrid_range_inv(1,1)
              + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * hybrid_range_inv(1,0), 
              .5
            );
            dcsigma12 (j1-1) = 
              pow(
                (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                  (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0)), 
                  .25) * //determinants
                  pow(arma::det(hybrid_range_inv), .5) * 
                  exp(-mahala_dist) * (1+ mahala_dist);
          }
        }
        dsigma11 = (dsigma11 - sigma11)*1000;
        dpsigma12 = (dpsigma12 - sigma12)*1000;
        dcsigma12 = (dcsigma12 - sigma12)*1000;
        
        // computing gradient of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        out(0, 1 + m*d_idx) = 
          //a = 0, see article's annex
          arma::sum(salt % dcsigma12) * pow_invcondsd_3; // b
        //c = 0
        //rest of the coefficients
        the_rest = 
          - agmis11 *  dcsigma12 * inverse_cond_sd // d
          //e = 0
          - salt * out(0, 1 + m*d_idx); //f
          for(int j2 = 1; j2<bsize; j2++)
          {
            out(j2, 1 + m*d_idx) = the_rest(j2-1);
          }
          
          // case parents' ranges are differentiated 
          for(int j1 = 1; j1<bsize;  j1++)
          {
            // computing stuff         
            //                     |     .    | 
            //                     |     .    | 
            //                     | .........|     
            //        (---------)  |     .    |
            
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)  += salt(j1-1)*dsigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1)+= salt(j2  )*dsigma11(j1-1, j2); //  vertical bar of the cross
            }
            salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            out(0, 1 + j1+ m*d_idx) = 
              // a = 0
              dpsigma12(j1-1) * salt(j1-1)  * pow_invcondsd_3 // b 
              -arma::sum(salt_dsigma11 % salt)     * pow_invcondsd_3/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              out(j2, 1+ j1 + m*d_idx) = - agmis11(j2-1, j1-1) *  dpsigma12(j1-1) * out(0, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * out(0, 0) // e
              - salt(j2-1)                                    * out(0, 1 + j1 + m*d_idx) ; //f
            }
          }
      }
    }
    result.slice(i) = out;
  }
  
  return(result);
}


//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]

arma::vec derivative_sandwich 
(
    arma::cube derivative, 
    arma::vec left_vector, 
    arma::vec right_vector, 
    Rcpp::IntegerMatrix NNarray 
)
{
  int  n = NNarray.nrow() ;
  int  m = NNarray.ncol() ;
  arma::vec res(n);
  res.fill(0);
  // looping over the row idx of the derivative of tilde R
  // # pragma omp parallel
  for(int row_idx = 0; row_idx<n; row_idx++){
    int bsize = std::min(row_idx+1,m);
    {
      // looping over the column idx of the derivative of tilde R (column idx retrieved through NNarray)
      for(int col_idx = 0; col_idx < bsize; col_idx++)
      {
        // looping over the index of the covariance parameter wrt which the derivative has been computed  (retrieved through NNarray)
        for(int covparm_idx = 0; covparm_idx < bsize; covparm_idx++)
        {
          //std::cout << NNarray(row_idx, 0) ; 
          //std::cout << "\n" ; 
          res(NNarray(row_idx, covparm_idx)-1)+= derivative(row_idx, covparm_idx, col_idx)
          * left_vector(row_idx) 
          * right_vector(NNarray(row_idx, col_idx)-1)
          ; 
        }
      }
    }
  }
  return(res);
}

// function used in sufficient ll gradient computation. Returns the gradient of the sum of diagonal terms.
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec log_determinant_derivative 
(
    arma::cube derivative, 
    arma::mat compressed_sparse_chol, 
    Rcpp::IntegerMatrix NNarray
)
{
  int  n = NNarray.nrow() ;
  int  m = NNarray.ncol() ;
  arma::vec res(n);
  res.fill(0);
  // looping over the row idx of the derivative of tilde R
  for(int row_idx = 0; row_idx<n; row_idx++){
    int bsize = std::min(row_idx+1,m);
    {
      // looping over the index of the covariance parameter wrt which the derivative has been computed  (retrieved through NNarray)
      for(int covparm_idx = 0; covparm_idx < bsize; covparm_idx++)
      {
        //std::cout << NNarray(row_idx, 0) ; 
        //std::cout << "\n" ; 
        res(NNarray(row_idx, covparm_idx)-1)+= derivative(row_idx, covparm_idx, 0)/compressed_sparse_chol(row_idx, 0); 
      }
    }
  }
  return(res);
}


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

Rcpp::NumericMatrix nonstat_covmat(
    arma::mat log_range,
    std::string covfun_name,
    arma::mat locs){
  
  // data dimensions
  int n = locs.n_rows;
  int dim = locs.n_cols;
  // objects of interest
  Rcpp::NumericMatrix Covmat(n,n);
  
  
  //case of isotropic range parameters
  if(covfun_name == "nonstationary_exponential_isotropic")
  {
    arma::mat sqeucdist(n, n);
    arma::vec range = exp(log_range) ; 
    // loop over every observation    
    for(int i=0; i<n; i++){
      for(int j=0; j<=i; j++){
        sqeucdist(i, j) = 0  ; 
        for(int k=0;k<dim;k++){ 
          sqeucdist(i, j) +=  pow((locs(i, k)-locs(j, k)), 2) ; 
        }
        sqeucdist(i, i) = 0  ; 
      }
    }
    // compute squexp covariance 
    for(int i=0; i<n; i++){
      for(int j=0; j<=i; j++){
        Covmat (i, j) = 
          pow(range(i) *    range(j)    , dim*.25  ) * 
          pow(range(i)*.5 + range(j)*.5 , dim*(-.5)) * 
          exp(- pow(sqeucdist(i, j)/(range(i)*.5 + range(j)*.5) , .5));
        Covmat (j, i) = Covmat (i, j) ; 
      }
    }
  }
  
  //case of anisotropic range parameters
  if(covfun_name == "nonstationary_exponential_anisotropic")
  {
    if(log_range.n_cols==3)
    {
      // range matrices computed from log_range
      arma::cube range_matrices(2,2,n); 
      for(int i=0; i<n; i++){
        arma::mat logmat(2, 2) ; 
        logmat(0, 0) = log_range(i, 0);
        logmat(1, 1) = log_range(i, 1);
        logmat(1, 0) = log_range(i, 2)/pow(2, .5);
        logmat(0, 1) = log_range(i, 2)/pow(2, .5);
        range_matrices.slice(i) = expmat_sym(logmat) ; 
      }
      for(int i=0; i<n; i++){
        Rcpp::checkUserInterrupt();
        for(int j=0; j<n; j++){
          arma::mat locdiff = locs.row(i) - locs.row(j);
          arma::mat hybrid_range = range_matrices.slice(i)*.5 + range_matrices.slice(j)*.5 ;
          arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          Covmat (i, j) = 
            pow(arma::det(range_matrices.slice(i)) * arma::det(range_matrices.slice(j)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0,0), .5));
          Covmat (j, i) = Covmat (i, j) ; 
        }
      }
    }
  }
  return(Covmat);
}


