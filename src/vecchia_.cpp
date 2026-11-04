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


//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;
using namespace Rcpp;

// ---------------------------------------------------------------------------
// expRange_ : matrix exponentials of the log range parameters
//   isotropic   (1 row)  -> 1 x n  : exp(log_range)
//   anisotropic (3 rows) -> 3 x n  (12 x n if compute_derivative)
//     rows 0-2  : (0,0), (1,1), (0,1) entries of expm(log range matrix)
//     rows 3-5  : same, diff_term added to log_range(0)
//     rows 6-8  : same, diff_term added to log_range(1)
//     rows 9-11 : same, diff_term/sqrt2 added to log_range(2)
// ---------------------------------------------------------------------------
//[[Rcpp::export]]
arma::mat expRange_(const arma::mat& log_range,
                    bool compute_derivative = false,
                    int num_threads = 1) {
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif
  const double diff_term = 0.0001;
  const double sqrt2 = std::sqrt(2.0);
  const int n = log_range.n_cols;
  
  if (log_range.n_rows != 1 && log_range.n_rows != 3)
    stop("log_range must have 1 (isotropic) or 3 (anisotropic) rows");
  
  if (log_range.n_rows == 1) return arma::exp(log_range);
  
  const int exprange_nrow = 3 + 9 * compute_derivative;
  arma::mat exprange(exprange_nrow, n);
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++) {
    arma::mat lograngemat(2, 2);
    arma::mat rangemat(2, 2);
    arma::vec log_range_col = log_range.col(i);
    
    lograngemat(0, 0) = log_range_col(0);
    lograngemat(1, 1) = log_range_col(1);
    lograngemat(1, 0) = log_range_col(2) / sqrt2;
    lograngemat(0, 1) = lograngemat(1, 0);
    rangemat = expmat_sym(lograngemat);
    
    exprange(0, i) = rangemat(0, 0);
    exprange(1, i) = rangemat(1, 1);
    exprange(2, i) = rangemat(0, 1);
    if (compute_derivative) {
      lograngemat(0, 0) += diff_term;
      rangemat = expmat_sym(lograngemat);
      exprange(3, i) = rangemat(0, 0);
      exprange(4, i) = rangemat(1, 1);
      exprange(5, i) = rangemat(0, 1);
      lograngemat(0, 0) -= diff_term;
      
      lograngemat(1, 1) += diff_term;
      rangemat = expmat_sym(lograngemat);
      exprange(6, i) = rangemat(0, 0);
      exprange(7, i) = rangemat(1, 1);
      exprange(8, i) = rangemat(0, 1);
      lograngemat(1, 1) -= diff_term;
      
      lograngemat(0, 1) += diff_term / sqrt2;
      lograngemat(1, 0) = lograngemat(0, 1);
      rangemat = expmat_sym(lograngemat);
      exprange(9, i)  = rangemat(0, 0);
      exprange(10, i) = rangemat(1, 1);
      exprange(11, i) = rangemat(0, 1);
    }
  }
  return exprange;
}

// ---------------------------------------------------------------------------
// vecchia_ : Vecchia approximation (and derivatives); computes exprange itself
// ---------------------------------------------------------------------------
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
  const bool smoothness15 = (smoothness == 1.5);
  const double diag_term = 1 + 1e-4;
  const double diff_term = 0.0001;
  const int n = locs.n_cols;
  const int m = NNarray.n_rows;
  const int n_params = log_range.n_rows;   // 1 (isotropic) or 3 (anisotropic)
  
  result(0,0,0) = 1;
  
  // matrix exponentials (done once, before the parallel region below)
  arma::mat exprange = expRange_(log_range, compute_derivative, num_threads);
  const int exprange_nrow = exprange.n_rows;
  
  // Computing Vecchia approx
#pragma omp parallel
{
  arma::mat out(m, 1 + n_params * m * compute_derivative);
#pragma omp for
  for(int i=start_idx; i<n; i++){
    // Fill subsets of locations, range
    int bsize = std::min(i+1,m);
    arma::mat locsub(2, bsize);
    arma::mat rangesub(exprange_nrow, bsize);
    for(int j=0; j<bsize; j++){
      int idx = NNarray(j,i);
      locsub.col(j) = locs.col(idx - 1);
      rangesub.col(j) = exprange.col(idx - 1);
    }
    // computing covariance
    double mahala_dist;
    double det_inv_hybrid_range[4];
    arma::mat  sigma11(bsize-1, bsize-1);
    arma::vec  sigma12(bsize-1);
    // anisotropic case
    if(n_params == 3){
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
    if(n_params == 1){
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
    // computing a vector used everywhere
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
      arma::mat dsigma11(bsize-1, bsize-1);
      arma::vec dpsigma12(bsize-1);
      arma::vec dcsigma12(bsize-1);
      arma::vec salt_dsigma11(bsize-1);
      arma::vec the_rest(bsize-1);
      arma::vec salt_dsigma11_agmis11(bsize-1);
      for(int d_idx=0; d_idx<n_params; d_idx++){
        
        // anisotropic case
        if(n_params == 3){
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
        if(n_params == 1){
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
        
        // case child range is differentiated
        // diagonal coefficient
        out(0, 1 + m*(d_idx)) =
          - arma::dot(salt, dcsigma12) * pow_invcondsd_3; // b
          // rest of the coefficients
          the_rest =
          agmis11 *  dcsigma12 * inverse_cond_sd // d
          - salt * out(0, 1 + m*d_idx); //f
          for(int j2 = 1; j2<bsize; j2++)
          {
            out(j2, 1 + m*d_idx) = the_rest(j2-1);
          }
          
          // case parents' ranges are differentiated
          for(int j1 = 1; j1<bsize;  j1++)
          {
            salt_dsigma11.fill(0);
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)  += salt(j1-1)*dsigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1)+= salt(j2  )*dsigma11(j1-1, j2); //  vertical bar of the cross
            }
            
            //diagonal coefficient
            out(0, 1 + j1+ m*d_idx) =
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
    result.col(i) = out;
  }
}
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





// PREDICT




// Closed-form exp of the symmetric 2x2 matrix [[a, b], [b, d]]
inline void expSym22(double a, double d, double b,
                     double &e11, double &e22, double &e12) {
  const double mu = 0.5 * (a + d), h = 0.5 * (a - d);
  const double s  = std::sqrt(h * h + b * b);
  const double ch = std::cosh(s);
  const double sh = (s < 1e-8) ? 1.0 : std::sinh(s) / s;
  const double em = std::exp(mu);
  e11 = em * (ch + sh * h);
  e22 = em * (ch - sh * h);
  e12 = em * sh * b;
}

// exp(range) at every location: 1 x n (isotropic) or 3 x n (aniso: 11, 22, 12)
// [[Rcpp::export]]
arma::mat expRanges(const arma::mat& log_range) {
  const int n = log_range.n_cols;
  if (log_range.n_rows == 1) return arma::exp(log_range);
  const double sqrt2 = std::sqrt(2.0);
  arma::mat e(3, n);
#pragma omp parallel for
  for (int i = 0; i < n; i++)
    expSym22(log_range(0, i), log_range(1, i), log_range(2, i) / sqrt2,
             e(0, i), e(1, i), e(2, i));
  return e;
}

// [[Rcpp::export]]
void vecchia_predict(
    const arma::mat& log_range,            // r x n_obs   (r = 1 iso, 3 aniso)
    const arma::mat& predicted_log_range,  // r x n_pred
    const arma::mat& locs,                 // 2 x n_obs
    const arma::mat& predicted_locs,       // 2 x n_pred
    const arma::mat& NNarray,              // m x n_pred, 1-based indices into locs
    int num_threads,
    arma::mat &result,
    double smoothness) {
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  const bool smoothness15 = (smoothness == 1.5);
  const double diag_term = 1.001;
  
  const int n_obs  = locs.n_cols;
  const int n_pred = predicted_locs.n_cols;
  const int m      = NNarray.n_rows;
  const int r      = log_range.n_rows;
  const bool aniso = (r == 3);
  
  if (r != 1 && r != 3) Rcpp::stop("log_range must have 1 or 3 rows");
  if (predicted_log_range.n_rows != (arma::uword) r) Rcpp::stop("predicted_log_range must have the same number of rows as log_range");
  if (log_range.n_cols != (arma::uword) n_obs) Rcpp::stop("log_range must have ncol(locs) columns");
  if (predicted_log_range.n_cols != (arma::uword) n_pred) Rcpp::stop("predicted_log_range must have ncol(predicted_locs) columns");
  if (NNarray.n_cols != (arma::uword) n_pred) Rcpp::stop("NNarray must have ncol(predicted_locs) columns");
  // validated here, serially: exceptions must not escape the OpenMP region
  for (int i = 0; i < n_pred; i++)
    for (int j = 0; j < m; j++) {
      const double v = NNarray(j, i);
      if (!(v >= 1 && v <= n_obs)) Rcpp::stop("NNarray entries must be in 1..ncol(locs)");
    }
  
  const arma::mat exprange      = expRanges(log_range);
  const arma::mat exprange_pred = expRanges(predicted_log_range);
  
  
#pragma omp parallel for
  for (int i = 0; i < n_pred; i++) {
    // index 0 = the predicted location (child), 1..m = its neighbours in locs
    const int bsize = m + 1;
    arma::mat locsub(2, bsize);
    arma::mat rs(exprange.n_rows, bsize);
    locsub.col(0) = predicted_locs.col(i);
    rs.col(0)     = exprange_pred.col(i);
    for (int j = 1; j < bsize; j++) {
      const int idx = (int) NNarray(j - 1, i) - 1;
      locsub.col(j) = locs.col(idx);
      rs.col(j)     = exprange.col(idx);
    }

    // correlation between entries a and b of the local set
    auto cov = [&](int a, int b) -> double {
      const double dx = locsub(0, a) - locsub(0, b);
      const double dy = locsub(1, a) - locsub(1, b);
      if (aniso) {
        double dih[4];
        DetInvSym22((rs(0, a) + rs(0, b)) * .5,
                    (rs(1, a) + rs(1, b)) * .5,
                    (rs(2, a) + rs(2, b)) * .5, dih);
        const double mahala = std::sqrt(dx * dx * dih[0] + dy * dy * dih[1] + 2 * dx * dy * dih[2]);
        return std::sqrt(std::sqrt(
            (rs(0, a) * rs(1, a) - rs(2, a) * rs(2, a)) *
              (rs(0, b) * rs(1, b) - rs(2, b) * rs(2, b))) * dih[3]) *
              maternThingy_(mahala, smoothness15);
      }
      const double ra = rs(0, a), rb = rs(0, b);
      const double mahala = std::sqrt((dx * dx + dy * dy) / (ra * .5 + rb * .5));
      return std::sqrt(std::sqrt(ra * rb) / (ra * .5 + rb * .5)) *
        maternThingy_(mahala, smoothness15);
    };

    arma::mat sigma11(m, m);
    arma::vec sigma12(m);
    for (int j1 = 1; j1 < bsize; j1++) {
      for (int j2 = 1; j2 < j1; j2++) {
        sigma11(j1 - 1, j2 - 1) = cov(j1, j2);
        sigma11(j2 - 1, j1 - 1) = sigma11(j1 - 1, j2 - 1);
      }
      sigma11(j1 - 1, j1 - 1) = diag_term;
      sigma12(j1 - 1) = cov(j1, 0);
    }

    // salt = sigma11^{-1} sigma12 through a Cholesky solve (no explicit inverse)
    arma::mat L;
    arma::chol(L, sigma11, "lower");
    const arma::vec salt = arma::solve(arma::trimatu(L.t()),
                                       arma::solve(arma::trimatl(L), sigma12));
    const double cond_var = diag_term - arma::dot(salt, sigma12);
    const double inverse_cond_sd = 1.0 / std::sqrt(cond_var);
    result(i, 0) = inverse_cond_sd;
    for (int j = 1; j < bsize; j++) result(i, j) = -salt(j - 1) * inverse_cond_sd;
  }
}

