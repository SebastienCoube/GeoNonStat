#define ARMA_NO_DEBUG
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
//[[Rcpp::export]]

double maternThingy_(double mahala_dist, bool smoothness15){
  double res;
  if (smoothness15){res = (exp(- mahala_dist) * (1 + mahala_dist));}else{res = (exp(- mahala_dist));}
  return(res);
}

double det_sym22(double a, double d, double b){
  return(a*d - b*b);
}

// double res[3] 
void inv_sym22(double det, double a, double d, double b, double res[3]){
  res[0] = d / det;
  res[1] = a / det;
  res[2] = -b / det;
}

void DetInvSym22(double a, double d, double b, double res[4]){
  res[3] = 1 / (a*d - b*b);
  res[0] = d  *  res[3];
  res[1] = a  *  res[3];
  res[2] = -b * res[3];
}


using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
void vecchia_(
    arma::mat log_range,
    arma::mat locs,
    arma::mat NNarray, 
    int num_threads,
    int start_idx,
    double smoothness,
    bool compute_derivative, 
    arma::cube &result){
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  bool smoothness15 = false; 
  if (smoothness == 1.5)smoothness15 = true;
  double diag_term = 1.001;
  double diff_term = 0.0001;
  double sqrt2 = sqrt(2); 
  // initializing stuff
  int n = locs.n_cols;
  int m = NNarray.n_rows;
  //arma::cube result(m, 1 + log_range.n_rows * m * compute_derivative, n);
  result(0,0,0) = 1;
  
  // pre-compute expensive matrix exponentials in case of aniso 
  int exprange_nrow = 1;
  int exprange_ncol = n;
  if(log_range.n_rows == 3){
    exprange_nrow = 3 + 9*compute_derivative;
  }
  arma::mat exprange(exprange_nrow, exprange_ncol);
  if(log_range.n_rows == 1){
    // loop over every observation  
    exprange = exp(log_range);
  }
  if(log_range.n_rows == 3){
    // loop over every observation  
    
#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; i++){
      arma::mat lograngemat(2,2);
      arma::mat rangemat(2,2);
      arma::vec log_range_col = log_range.col(i); 
      
      lograngemat(0,0) = log_range_col(0);
      lograngemat(1,1) = log_range_col(1);
      lograngemat(1,0) = log_range_col(2)/sqrt2;
      lograngemat(0,1) = lograngemat(1,0);
      rangemat = expmat_sym(lograngemat);
      
      exprange(0, i) = rangemat(0,0); 
      exprange(1, i) = rangemat(1,1); 
      exprange(2, i) = rangemat(0,1); 
      if(compute_derivative){
        lograngemat(0,0) = lograngemat(0,0) + diff_term;
        rangemat = expmat_sym(lograngemat);
        exprange(3, i) = rangemat(0,0); 
        exprange(4, i) = rangemat(1,1); 
        exprange(5, i) = rangemat(0,1); 
        lograngemat(0,0) = lograngemat(0,0) - diff_term;
        lograngemat(1,1) = lograngemat(1,1) + diff_term;
        rangemat = expmat_sym(lograngemat);
        exprange(6, i) = rangemat(0,0); 
        exprange(7, i) = rangemat(1,1); 
        exprange(8, i) = rangemat(0,1); 
        lograngemat(1,1) = lograngemat(1,1) - diff_term;
        lograngemat(0,1) = lograngemat(0,1) + diff_term/sqrt2;
        lograngemat(1,0) = lograngemat(0,1);
        rangemat = expmat_sym(lograngemat);
        exprange(9, i)  = rangemat(0,0); 
        exprange(10, i) = rangemat(1,1); 
        exprange(11, i) = rangemat(0,1); 
      }
    }
    
  }
  
  // Computing Vecchia approx
#pragma omp parallel 
  {
  arma::mat out(m, 1 + log_range.n_rows * m * compute_derivative);
#pragma omp for
  for(int i=start_idx; i<n; i++){
    // Fill subsets of locations, range
    int bsize = std::min(i+1,m);
    arma::mat locsub(2, bsize);
    arma::mat rangesub(exprange_nrow, bsize);
    for(int j=0; j<bsize; j++){
      int idx = NNarray(j,i);
      // subsets of locations
      locsub.col(j) = locs.col(idx - 1); 
      // subsets of range
      rangesub.col(j) = exprange.col(idx - 1); 
    }
    // initialize recipient for result
    // computing covariance 
    double mahala_dist;
    double det_inv_hybrid_range[4];
    arma::mat  sigma11(bsize-1, bsize-1);
    arma::vec  sigma12(bsize-1);
    // anisotropic case
    if(log_range.n_rows == 3){
      for(int j1=1; j1<bsize; j1++){
        //filling sigma11
        for(int j2=1; j2<j1; j2++){
          DetInvSym22((rangesub(0, j1) + rangesub(0, j2))*.5,  
                      (rangesub(1, j1) + rangesub(1, j2))*.5,  
                      (rangesub(2, j1) + rangesub(2, j2))*.5, 
                      det_inv_hybrid_range);
          
          mahala_dist = sqrt(
            (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) * det_inv_hybrid_range[0]
          + (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[1]
          + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[2]
          );
          sigma11 (j1-1, j2-1) = 
            sqrt(sqrt(
                (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                  (rangesub(0, j2)*rangesub(1, j2) - rangesub(2, j2)*rangesub(2, j2))) * //determinants
                  (det_inv_hybrid_range[3])) * 
                  maternThingy_(mahala_dist, smoothness15);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
          
        }
        sigma11 (j1-1, j1-1) = diag_term;
        //filling sigma12
        DetInvSym22((rangesub(0, j1) + rangesub(0, 0))*.5, 
                    (rangesub(1, j1) + rangesub(1, 0))*.5, 
                    (rangesub(2, j1) + rangesub(2, 0))*.5, det_inv_hybrid_range);
        mahala_dist = sqrt(
          (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) * det_inv_hybrid_range[0]
        +(locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[1]
        + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
        );
        sigma12 (j1-1) = 
          sqrt(sqrt(
              (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0))) * //determinants
                (det_inv_hybrid_range[3])) * 
                maternThingy_(mahala_dist, smoothness15);
      }
    }
    // isotropic case
    if(log_range.n_rows == 1){
      for(int j1=1; j1<bsize; j1++){
        //filling sigma11
        for(int j2=1; j2<j1; j2++){
          mahala_dist = sqrt( 
            ( 
                (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) + 
                  (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2))
            )
          /(rangesub(0, j1)*.5 + rangesub(0, j2)*.5));
          sigma11 (j1-1, j2-1) = 
          sqrt(sqrt(rangesub(0, j1) *    rangesub(0, j2))/ 
          (rangesub(0, j1)*.5 + rangesub(0, j2)*.5)) * 
          maternThingy_(mahala_dist, smoothness15);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
        }
        sigma11 (j1-1, j1-1) = diag_term;
        mahala_dist = sqrt( 
          ( (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) )
          /(rangesub(0, j1)*.5 + rangesub(0, 0)*.5));
        sigma12 (j1-1) = 
        sqrt(sqrt(rangesub(0, j1) *    rangesub(0, 0)) / 
        (rangesub(0, j1)*.5 + rangesub(0, 0)*.5)) * 
        maternThingy_(mahala_dist, smoothness15);
      }
    }
    
    // solving sigma11
    arma::mat agmis11  = arma::inv_sympd(sigma11) ;
    // computing a vector used everyvhere
    arma::mat salt  = agmis11 * sigma12 ;
    //computing Vecchia approx itself
    double inverse_cond_sd = 1/sqrt(diag_term- arma::dot(salt, sigma12));
    double pow_invcondsd_3  = inverse_cond_sd*inverse_cond_sd*inverse_cond_sd;
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
      //for(int d_idx=0; d_idx<log_range.n_rows; d_idx++){
      for(int d_idx=0; d_idx<log_range.n_rows; d_idx++){
        
        // computing derivative of covariance
        // column j1 of dsigma11 contains d(sigma11(,j1))/d(range(j1))
        
        // anisotropic case
        if(log_range.n_rows == 3){
          for(int j1=1; j1<bsize; j1++){
            //filling dsigma11
            for(int j2=1; j2<bsize; j2++){
              DetInvSym22(
                (rangesub((d_idx + 1) * 3 + 0, j1) + rangesub(0, j2))*.5,
                (rangesub((d_idx + 1) * 3 + 1, j1) + rangesub(1, j2))*.5,
                (rangesub((d_idx + 1) * 3 + 2, j1) + rangesub(2, j2))*.5,
                det_inv_hybrid_range);
              double mahala_dist = sqrt(
                ((locsub(0, j1)-locsub(0, j2))) * ((locsub(0, j1)-locsub(0, j2))) * det_inv_hybrid_range[0]
              +((locsub(1, j1)-locsub(1, j2))) * ((locsub(1, j1)-locsub(1, j2))) * det_inv_hybrid_range[1]
              + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[2]
              );
              dsigma11 (j1-1, j2-1) = 
                sqrt(sqrt(
                    (rangesub((d_idx + 1) * 3 + 0, j1)*rangesub((d_idx + 1) * 3 + 1, j1) - rangesub((d_idx + 1) * 3 + 2, j1)*rangesub((d_idx + 1) * 3 + 2, j1)) * 
                      (rangesub(0, j2)*rangesub(1, j2) - rangesub(2, j2)*rangesub(2, j2))) * //determinants
                      (det_inv_hybrid_range[3])) * 
                      maternThingy_(mahala_dist, smoothness15);
            }
            dsigma11 (j1-1, j1-1) = diag_term;
            //filling dpsigma12
            DetInvSym22(
              (rangesub((d_idx + 1) * 3 + 0, j1) + rangesub(0, 0))*.5,  
              (rangesub((d_idx + 1) * 3 + 1, j1) + rangesub(1, 0))*.5, 
              (rangesub((d_idx + 1) * 3 + 2, j1) + rangesub(2, 0))*.5, 
              det_inv_hybrid_range);
            mahala_dist = sqrt(
              (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) * det_inv_hybrid_range[0]
            +(locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[1]
            + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
            );
            dpsigma12 (j1-1) = 
              sqrt(sqrt(
                  (rangesub((d_idx + 1) * 3 + 0, j1)*rangesub((d_idx + 1) * 3 + 1, j1) - rangesub((d_idx + 1) * 3 + 2, j1)*rangesub((d_idx + 1) * 3 + 2, j1)) * 
                    (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0))) * //determinants
                    (det_inv_hybrid_range[3])) * 
                    maternThingy_(mahala_dist, smoothness15);
            //filling dcsigma12
            DetInvSym22((rangesub(0, j1) + rangesub((d_idx + 1) * 3 + 0, 0))*.5, 
                        (rangesub(1, j1) + rangesub((d_idx + 1) * 3 + 1, 0))*.5, 
                        (rangesub(2, j1) + rangesub((d_idx + 1) * 3 + 2, 0))*.5,
                        det_inv_hybrid_range);
            mahala_dist = sqrt(
              ((locsub(0, j1)-locsub(0, 0))) * ((locsub(0, j1)-locsub(0, 0))) * det_inv_hybrid_range[0]
            + ((locsub(1, j1)-locsub(1, 0))) * ((locsub(1, j1)-locsub(1, 0))) * det_inv_hybrid_range[1]
            + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
            );
            dcsigma12 (j1-1) = 
              sqrt(sqrt(
                  (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                    (rangesub((d_idx + 1) * 3 + 0, 0 )*rangesub((d_idx + 1) * 3 + 1, 0 ) - 
                    rangesub((d_idx + 1) * 3 + 2, 0 )*rangesub((d_idx + 1) * 3 + 2, 0 ))) * //determinants
                    (det_inv_hybrid_range[3])) * 
                    maternThingy_(mahala_dist, smoothness15);
            
          }
        }
        // isotropic case
        if(log_range.n_rows == 1){
          for(int j1=1; j1<bsize; j1++){
            //filling dsigma11
            for(int j2=1; j2<bsize; j2++){
              mahala_dist = sqrt( 
                ( (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) + 
                  (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2)) )
              /(rangesub(0, j1)*.5*(1 + diff_term) + rangesub(0, j2)*.5));
              dsigma11 (j1-1, j2-1) = 
              sqrt(sqrt(rangesub(0, j1) * (1 + diff_term) * rangesub(0, j2)) / 
              (rangesub(0, j1) * (1 + diff_term) *.5 + rangesub(0, j2)*.5)) * 
              maternThingy_(mahala_dist, smoothness15);
            }
            dsigma11 (j1-1, j1-1) = diag_term;
            //filling dpsigma12
            mahala_dist = sqrt( 
              ( (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)))
              /(rangesub(0, j1) * (1 + diff_term) *.5 + rangesub(0, 0)*.5));
            dpsigma12 (j1-1) = 
            sqrt(sqrt(rangesub(0, j1)  * (1 + diff_term)  *    rangesub(0, 0)) / 
            (rangesub(0, j1)  * (1 + diff_term) *.5 + rangesub(0, 0)*.5)) * 
            maternThingy_(mahala_dist, smoothness15);
            
            //filling dcsigma12
            mahala_dist = sqrt( 
              ((locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)))
              /(rangesub(0, j1) *.5 + rangesub(0, 0)*.5 * (1 + diff_term)));
            dcsigma12 (j1-1) = 
            sqrt(sqrt(rangesub(0, j1)  *    rangesub(0, 0)  * (1 + diff_term)) / 
            (rangesub(0, j1) *.5 + rangesub(0, 0)  * (1 + diff_term) *.5)) * 
            maternThingy_(mahala_dist, smoothness15);
            
          }
        }
        
        dsigma11 = (dsigma11 - sigma11)/diff_term;
        dpsigma12 = (dpsigma12 - sigma12)/diff_term;
        dcsigma12 = -(dcsigma12 - sigma12)/diff_term;
        // computing derivative of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        out(0, 1 + m*(d_idx)) = 
          //a = 0, see article's annex
          - arma::dot(salt, dcsigma12) * pow_invcondsd_3; // b
          //c = 0
          //rest of the coefficients
          the_rest = 
          agmis11 *  dcsigma12 * inverse_cond_sd // d
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
            //                d(sigma11)/d(range(j1)) 
            //                     |     .    | 
            //                     |     .    | 
            //                     |     .    | 
            //                     | .........|     
            //         (--------)  |     .    |
            //           salt
            salt_dsigma11.fill(0);
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)  += salt(j1-1)*dsigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1)+= salt(j2  )*dsigma11(j1-1, j2); //  vertical bar of the cross
            }
            
            //diagonal coefficient
            out(0, 1 + j1+ m*d_idx) = 
              // a = 0
              dpsigma12(j1-1) * salt(j1-1)  * pow_invcondsd_3 // b 
            - arma::dot(salt_dsigma11, salt) * pow_invcondsd_3/2 ; // c
            
            salt_dsigma11_agmis11 =   agmis11 * salt_dsigma11 ; 
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
    //# pragma omp critical
    result.col(i) = out;
  }
}
  //return(result);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::mat derivativeSandwiches_
(
    arma::cube vecchia, 
    arma::vec left_vector, 
    arma::vec right_vector, 
    arma::mat NNarray, 
    bool sauce_determinant_chef,
    int num_threads
)
{
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  int  n = NNarray.n_cols;
  int  m = NNarray.n_rows;
  int d = (vecchia.n_slices - 1)/vecchia.n_rows;
  arma::mat res(d,  n);
  
#pragma omp parallel for
  for(int row_idx = 0; row_idx<n; row_idx++){
    arma::vec NNarray_col = NNarray.col(row_idx);
    int bsize = std::min(row_idx+1,m);
    // anisotropy index, from 0 to 2 if aniso, from 0 to 0 if iso
    // looping over the index of the covariance parameter wrt which the derivative has been computed  (retrieved through NNarray)
    arma::rowvec right_vector_sub(m);
    for(int col_idx = 0; col_idx < bsize; col_idx++){
      right_vector_sub(col_idx) = right_vector(NNarray_col(col_idx)-1)
      ; 
    }
    arma::mat vecchia_col= vecchia.col(row_idx);
    arma::rowvec mini_sandwich = right_vector_sub * vecchia_col * left_vector(row_idx)  ;
     if(!sauce_determinant_chef){
 # pragma omp critical 
       for(int aniso_idx = 0; aniso_idx <d; aniso_idx++){
         for(int col_idx = 0; col_idx < bsize; col_idx++){
           res(aniso_idx, NNarray_col(col_idx)-1) += mini_sandwich(1 + col_idx + aniso_idx*m);
         }
       }
     }else{
 # pragma omp critical 
       for(int aniso_idx = 0; aniso_idx <d; aniso_idx++){
         for(int col_idx = 0; col_idx < bsize; col_idx++){
           res(aniso_idx, NNarray_col(col_idx)-1) += 
             mini_sandwich(   1 + col_idx + aniso_idx*m) - 
             vecchia_col(0, 1 + col_idx + aniso_idx*m)/vecchia_col(0, 0);
         }
       }
     }
   }
  return(res);
}








using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
void vecchia_stoopid(
    arma::mat log_range,
    arma::mat locs,
    arma::mat NNarray, 
    int num_threads,
    double smoothness,
    bool compute_derivative, 
    arma::cube &result){
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  bool smoothness15 = false; 
  if (smoothness == 1.5)smoothness15 = true;
  double diag_term = 1.0001;
  double diff_term = 0.0000001;
  double sqrt2 = sqrt(2); 
  // initializing stuff
  int n = locs.n_cols;
  int m = NNarray.n_rows;
  //arma::cube result(m, 1 + log_range.n_rows * m * compute_derivative, n);
  result(0,0,0) = 1;
  
  // pre-compute expensive matrix exponentials in case of aniso 
  int exprange_nrow = 1;
  int exprange_ncol = n;
  if(log_range.n_rows == 3){
    exprange_nrow = 3 + 9*compute_derivative;
  }
  arma::mat exprange(exprange_nrow, exprange_ncol);
  if(log_range.n_rows == 1){
    // loop over every observation  
    exprange = exp(log_range);
  }
  if(log_range.n_rows == 3){
    // loop over every observation  
    
#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; i++){
      arma::mat lograngemat(2,2);
      arma::mat rangemat(2,2);
      arma::vec log_range_col = log_range.col(i); 
      
      lograngemat(0,0) = log_range_col(0);
      lograngemat(1,1) = log_range_col(1);
      lograngemat(1,0) = log_range_col(2)/sqrt2;
      lograngemat(0,1) = lograngemat(1,0);
      rangemat = expmat_sym(lograngemat);
      
      exprange(0, i) = rangemat(0,0); 
      exprange(1, i) = rangemat(1,1); 
      exprange(2, i) = rangemat(0,1); 
      if(compute_derivative){
        lograngemat(0,0) = lograngemat(0,0) + diff_term;
        rangemat = expmat_sym(lograngemat);
        exprange(3, i) = rangemat(0,0); 
        exprange(4, i) = rangemat(1,1); 
        exprange(5, i) = rangemat(0,1); 
        lograngemat(0,0) = lograngemat(0,0) - diff_term;
        lograngemat(1,1) = lograngemat(1,1) + diff_term;
        rangemat = expmat_sym(lograngemat);
        exprange(6, i) = rangemat(0,0); 
        exprange(7, i) = rangemat(1,1); 
        exprange(8, i) = rangemat(0,1); 
        lograngemat(1,1) = lograngemat(1,1) - diff_term;
        lograngemat(0,1) = lograngemat(0,1) + diff_term/sqrt2;
        lograngemat(1,0) = lograngemat(0,1);
        rangemat = expmat_sym(lograngemat);
        exprange(9, i)  = rangemat(0,0); 
        exprange(10, i) = rangemat(1,1); 
        exprange(11, i) = rangemat(0,1); 
      }
    }
    
  }
  
  // Computing Vecchia approx
#pragma omp parallel 
{
  arma::mat out(m, 1 + log_range.n_rows * m * compute_derivative);
#pragma omp for
  for(int i=1; i<n; i++){
    // Fill subsets of locations, range
    int bsize = std::min(i+1,m);
    arma::mat locsub(2, bsize);
    arma::mat rangesub(exprange_nrow, bsize);
    for(int j=0; j<bsize; j++){
      int idx = NNarray(j,i);
      // subsets of locations
      locsub.col(j) = locs.col(idx - 1); 
      // subsets of range
      rangesub.col(j) = exprange.col( NNarray(0,i) - 1); 
    }
    // initialize recipient for result
    // computing covariance 
    double mahala_dist;
    double det_inv_hybrid_range[4];
    arma::mat  sigma11(bsize-1, bsize-1);
    arma::vec  sigma12(bsize-1);
    // anisotropic case
    if(log_range.n_rows == 3){
      for(int j1=1; j1<bsize; j1++){
        //filling sigma11
        for(int j2=1; j2<j1; j2++){
          DetInvSym22((rangesub(0, j1) + rangesub(0, j2))*.5,  
                      (rangesub(1, j1) + rangesub(1, j2))*.5,  
                      (rangesub(2, j1) + rangesub(2, j2))*.5, 
                      det_inv_hybrid_range);
          
          mahala_dist = sqrt(
            (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) * det_inv_hybrid_range[0]
          + (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[1]
          + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[2]
          );
          sigma11 (j1-1, j2-1) = 
            sqrt(sqrt(
                (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                  (rangesub(0, j2)*rangesub(1, j2) - rangesub(2, j2)*rangesub(2, j2))) * //determinants
                  (det_inv_hybrid_range[3])) * 
                  maternThingy_(mahala_dist, smoothness15);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
          
        }
        sigma11 (j1-1, j1-1) = diag_term;
        //filling sigma12
        DetInvSym22((rangesub(0, j1) + rangesub(0, 0))*.5, 
                    (rangesub(1, j1) + rangesub(1, 0))*.5, 
                    (rangesub(2, j1) + rangesub(2, 0))*.5, det_inv_hybrid_range);
        mahala_dist = sqrt(
          (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) * det_inv_hybrid_range[0]
        +(locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[1]
        + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
        );
        sigma12 (j1-1) = 
          sqrt(sqrt(
              (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0))) * //determinants
                (det_inv_hybrid_range[3])) * 
                maternThingy_(mahala_dist, smoothness15);
      }
    }
    // isotropic case
    if(log_range.n_rows == 1){
      for(int j1=1; j1<bsize; j1++){
        //filling sigma11
        for(int j2=1; j2<j1; j2++){
          mahala_dist = sqrt( 
            ( 
                (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) + 
                  (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2))
            )
          /(rangesub(0, j1)*.5 + rangesub(0, j2)*.5));
          sigma11 (j1-1, j2-1) = 
          sqrt(sqrt(rangesub(0, j1) *    rangesub(0, j2))/ 
          (rangesub(0, j1)*.5 + rangesub(0, j2)*.5)) * 
          maternThingy_(mahala_dist, smoothness15);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
        }
        sigma11 (j1-1, j1-1) = diag_term;
        mahala_dist = sqrt( 
          ( (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) )
          /(rangesub(0, j1)*.5 + rangesub(0, 0)*.5));
        sigma12 (j1-1) = 
        sqrt(sqrt(rangesub(0, j1) *    rangesub(0, 0)) / 
        (rangesub(0, j1)*.5 + rangesub(0, 0)*.5)) * 
        maternThingy_(mahala_dist, smoothness15);
      }
    }
    
    // solving sigma11
    arma::mat agmis11  = arma::inv_sympd(sigma11) ;
    // computing a vector used everyvhere
    arma::mat salt  = agmis11 * sigma12 ;
    //computing Vecchia approx itself
    double inverse_cond_sd = 1/sqrt(diag_term- arma::dot(salt, sigma12));
    double pow_invcondsd_3  = inverse_cond_sd*inverse_cond_sd*inverse_cond_sd;
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
      //for(int d_idx=0; d_idx<log_range.n_rows; d_idx++){
      for(int d_idx=0; d_idx<log_range.n_rows; d_idx++){
        
        // computing derivative of covariance
        // column j1 of dsigma11 contains d(sigma11(,j1))/d(range(j1))
        
        // anisotropic case
        if(log_range.n_rows == 3){
          for(int j1=1; j1<bsize; j1++){
            //filling dsigma11
            for(int j2=1; j2<bsize; j2++){
              DetInvSym22(
                (rangesub((d_idx + 1) * 3 + 0, j1) + rangesub(0, j2))*.5,
                (rangesub((d_idx + 1) * 3 + 1, j1) + rangesub(1, j2))*.5,
                (rangesub((d_idx + 1) * 3 + 2, j1) + rangesub(2, j2))*.5,
                det_inv_hybrid_range);
              double mahala_dist = sqrt(
                ((locsub(0, j1)-locsub(0, j2))) * ((locsub(0, j1)-locsub(0, j2))) * det_inv_hybrid_range[0]
              +((locsub(1, j1)-locsub(1, j2))) * ((locsub(1, j1)-locsub(1, j2))) * det_inv_hybrid_range[1]
              + 2*( locsub(0, j1)-locsub(0, j2)) * (locsub(1, j1)-locsub(1, j2)) * det_inv_hybrid_range[2]
              );
              dsigma11 (j1-1, j2-1) = 
                sqrt(sqrt(
                    (rangesub((d_idx + 1) * 3 + 0, j1)*rangesub((d_idx + 1) * 3 + 1, j1) - rangesub((d_idx + 1) * 3 + 2, j1)*rangesub((d_idx + 1) * 3 + 2, j1)) * 
                      (rangesub(0, j2)*rangesub(1, j2) - rangesub(2, j2)*rangesub(2, j2))) * //determinants
                      (det_inv_hybrid_range[3])) * 
                      maternThingy_(mahala_dist, smoothness15);
            }
            dsigma11 (j1-1, j1-1) = diag_term;
            //filling dpsigma12
            DetInvSym22(
              (rangesub((d_idx + 1) * 3 + 0, j1) + rangesub(0, 0))*.5,  
              (rangesub((d_idx + 1) * 3 + 1, j1) + rangesub(1, 0))*.5, 
              (rangesub((d_idx + 1) * 3 + 2, j1) + rangesub(2, 0))*.5, 
              det_inv_hybrid_range);
            mahala_dist = sqrt(
              (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) * det_inv_hybrid_range[0]
            +(locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[1]
            + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
            );
            dpsigma12 (j1-1) = 
              sqrt(sqrt(
                  (rangesub((d_idx + 1) * 3 + 0, j1)*rangesub((d_idx + 1) * 3 + 1, j1) - rangesub((d_idx + 1) * 3 + 2, j1)*rangesub((d_idx + 1) * 3 + 2, j1)) * 
                    (rangesub(0, 0)*rangesub(1, 0) - rangesub(2, 0)*rangesub(2, 0))) * //determinants
                    (det_inv_hybrid_range[3])) * 
                    maternThingy_(mahala_dist, smoothness15);
            //filling dcsigma12
            DetInvSym22((rangesub(0, j1) + rangesub((d_idx + 1) * 3 + 0, 0))*.5, 
                        (rangesub(1, j1) + rangesub((d_idx + 1) * 3 + 1, 0))*.5, 
                        (rangesub(2, j1) + rangesub((d_idx + 1) * 3 + 2, 0))*.5,
                        det_inv_hybrid_range);
            mahala_dist = sqrt(
              ((locsub(0, j1)-locsub(0, 0))) * ((locsub(0, j1)-locsub(0, 0))) * det_inv_hybrid_range[0]
            + ((locsub(1, j1)-locsub(1, 0))) * ((locsub(1, j1)-locsub(1, 0))) * det_inv_hybrid_range[1]
            + 2*( locsub(0, j1)-locsub(0, 0)) * (locsub(1, j1)-locsub(1, 0)) * det_inv_hybrid_range[2]
            );
            dcsigma12 (j1-1) = 
              sqrt(sqrt(
                  (rangesub(0, j1)*rangesub(1, j1) - rangesub(2, j1)*rangesub(2, j1)) * 
                    (rangesub((d_idx + 1) * 3 + 0, 0 )*rangesub((d_idx + 1) * 3 + 1, 0 ) - 
                    rangesub((d_idx + 1) * 3 + 2, 0 )*rangesub((d_idx + 1) * 3 + 2, 0 ))) * //determinants
                    (det_inv_hybrid_range[3])) * 
                    maternThingy_(mahala_dist, smoothness15);
            
          }
        }
        // isotropic case
        if(log_range.n_rows == 1){
          for(int j1=1; j1<bsize; j1++){
            //filling dsigma11
            for(int j2=1; j2<bsize; j2++){
              mahala_dist = sqrt( 
                ( (locsub(0, j1)-locsub(0, j2)) * (locsub(0, j1)-locsub(0, j2)) + 
                  (locsub(1, j1)-locsub(1, j2)) * (locsub(1, j1)-locsub(1, j2)) )
              /(rangesub(0, j1)*.5*(1 + diff_term) + rangesub(0, j2)*.5));
              dsigma11 (j1-1, j2-1) = 
              sqrt(sqrt(rangesub(0, j1) * (1 + diff_term) * rangesub(0, j2)) / 
              (rangesub(0, j1) * (1 + diff_term) *.5 + rangesub(0, j2)*.5)) * 
              maternThingy_(mahala_dist, smoothness15);
            }
            dsigma11 (j1-1, j1-1) = diag_term;
            //filling dpsigma12
            mahala_dist = sqrt( 
              ( (locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)))
              /(rangesub(0, j1) * (1 + diff_term) *.5 + rangesub(0, 0)*.5));
            dpsigma12 (j1-1) = 
            sqrt(sqrt(rangesub(0, j1)  * (1 + diff_term)  *    rangesub(0, 0)) / 
            (rangesub(0, j1)  * (1 + diff_term) *.5 + rangesub(0, 0)*.5)) * 
            maternThingy_(mahala_dist, smoothness15);
            
            //filling dcsigma12
            mahala_dist = sqrt( 
              ((locsub(0, j1)-locsub(0, 0)) * (locsub(0, j1)-locsub(0, 0)) + (locsub(1, j1)-locsub(1, 0)) * (locsub(1, j1)-locsub(1, 0)))
              /(rangesub(0, j1) *.5 + rangesub(0, 0)*.5 * (1 + diff_term)));
            dcsigma12 (j1-1) = 
            sqrt(sqrt(rangesub(0, j1)  *    rangesub(0, 0)  * (1 + diff_term)) / 
            (rangesub(0, j1) *.5 + rangesub(0, 0)  * (1 + diff_term) *.5)) * 
            maternThingy_(mahala_dist, smoothness15);
            
          }
        }
        
        dsigma11 = (dsigma11 - sigma11)/diff_term;
        dpsigma12 = (dpsigma12 - sigma12)/diff_term;
        dcsigma12 = -(dcsigma12 - sigma12)/diff_term;
        // computing derivative of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        out(0, 1 + m*(d_idx)) = 
          //a = 0, see article's annex
          - arma::dot(salt, dcsigma12) * pow_invcondsd_3; // b
          //c = 0
          //rest of the coefficients
          the_rest = 
          agmis11 *  dcsigma12 * inverse_cond_sd // d
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
            //                d(sigma11)/d(range(j1)) 
            //                     |     .    | 
            //                     |     .    | 
            //                     |     .    | 
            //                     | .........|     
            //         (--------)  |     .    |
            //           salt
            salt_dsigma11.fill(0);
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)  += salt(j1-1)*dsigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1)+= salt(j2  )*dsigma11(j1-1, j2); //  vertical bar of the cross
            }
            
            //diagonal coefficient
            out(0, 1 + j1+ m*d_idx) = 
              // a = 0
              dpsigma12(j1-1) * salt(j1-1)  * pow_invcondsd_3 // b 
            - arma::dot(salt_dsigma11, salt) * pow_invcondsd_3/2 ; // c
            
            salt_dsigma11_agmis11 =   agmis11 * salt_dsigma11 ; 
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
    //# pragma omp critical
    result.col(i) = out;
  }
}
//return(result);
}



