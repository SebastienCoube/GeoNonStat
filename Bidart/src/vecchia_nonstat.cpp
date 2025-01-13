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
Rcpp::List nonstat_vecchia_Linv(
    arma::mat log_range,
    std::string covfun_name,
    bool sphere,
    arma::mat locs,
    arma::mat NNarray, 
    int num_threads,
    bool compute_derivative){
#if _OPENMP
  omp_set_num_threads(num_threads);
#endif
  // data dimensions
  int n = locs.n_rows;
  int m = NNarray.n_cols;
  int dim = locs.n_cols;
  // objects of interest
  Rcpp::NumericMatrix Linv(n,m);
  Linv(0, 0) = 1;
  Rcpp::List grad(log_range.n_cols);
  // locs in 3d if on sphere
  arma::mat locs_3D(n, 3);
  if (sphere)
  {
    for(int i=0; i<n; i++){
      locs_3D(i, 0) = cos(M_PI /180 *locs(i, 0)) * cos(M_PI /180 *locs(i, 1));
      locs_3D(i, 1) = sin(M_PI /180 *locs(i, 0)) * cos(M_PI /180 *locs(i, 1));
      locs_3D(i, 2) = sin(M_PI /180 *locs(i, 1));
    }
  }
  
  if(covfun_name == "nonstationary_exponential_isotropic")
  {
    arma::vec range = exp(log_range) ; 
    arma::cube grad1(n,m, m);
    grad1.fill(0);
    // loop over every observation   
    #pragma omp parallel for
    for(int i=1; i<n; i++){
      int bsize = std::min(i+1,m);
      // first, fill in ysub, locsub, and X0 NOT in reverse order
      arma::mat locsub(bsize, dim);
      if(!sphere)
      {
        for(int j=0; j<bsize; j++){
          for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
        }
      }
      if(sphere)
      {
        // getting basis of tangent plane
        arma::mat tangent_plane_basis (3, 2);
        // parallel
        tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
        tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
        tangent_plane_basis (2, 0) = 0;
        // meridian
        tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
        tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
        tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
        
        arma::mat locsub_3D (bsize, 3);
        for(int j=0; j<bsize; j++){
          for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
        }
        locsub = locsub_3D * tangent_plane_basis ; 
        
      }
      // compute euclidean distance
      arma::mat sqeucdist(bsize, bsize);
      for(int j1=0; j1<bsize; j1++){
        for(int j2=0; j2<=j1; j2++){
          sqeucdist(j1, j2) = 0  ; 
          for(int k=0;k<dim;k++){ 
            sqeucdist(j1, j2) +=  pow((locsub(j1, k)-locsub(j2, k)), 2) ; 
          }
          sqeucdist(j2, j1) = sqeucdist(j1, j2) ; 
        }
        sqeucdist(j1, j1) = 0  ; 
      }
      
      // compute squexp covariance 
      arma::mat  sigma11(bsize-1, bsize-1);
      arma::vec  sigma12(bsize-1);
      for(int j1=1; j1<bsize; j1++){
        for(int j2=1; j2<=j1; j2++){
          sigma11 (j1-1, j2-1) = 
            pow(range(NNarray(i,j1) - 1) *    range(NNarray(i,j2) -1)    , dim*.25  ) * 
            pow(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5 , dim*(-.5)) * 
            exp(- pow(sqeucdist(j1, j2)/(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5) , .5));
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
        }
        sigma12 (j1-1) =                         
          pow(range(NNarray(i,j1) - 1)    * range(NNarray(i,0) - 1)    , dim*  .25) * 
          pow(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
          exp(-pow(sqeucdist(j1, 0)/ (range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5), .5)) ;
      }
      // solving sigma11
      arma::mat agmis11  = arma::inv(sigma11) ;
      // computing a vector used everyvhere
      arma::mat salt  = agmis11 * sigma12 ;
      
      // computing Vecchia approx itself
      double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
      Linv(i, 0) = inverse_cond_sd ;
      for(int j=1; j<bsize; j++){
        Linv(i, j) = - salt(j-1) * inverse_cond_sd;
      }
      
      if(compute_derivative)
      {
        // compute squexp covariance derivatives using finite differences
        arma::mat dsigma11(bsize-1, bsize-1);
        arma::vec dsigma12_parents(bsize-1);
        arma::vec dsigma12_child(bsize-1);
        for(int j1=1; j1<bsize; j1++){
          for(int j2=1; j2<bsize; j2++){
            //the derivative of sigma11 is a symmetric cross matrix with non-null coefficients only on one row and column. Instead of computing the cross matrix, just one row is computed. The rows are stacked in a matrix. 
            dsigma11 (j1-1, j2-1) = (
              pow(1.000001 * range(NNarray(i,j1) - 1) *    range(NNarray(i,j2) -1)    , dim*.25  ) * 
                pow(1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5 , dim*(-.5)) * 
                exp(- pow(sqeucdist(j1, j2)/ (1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5), .5))
            - sigma11 (j1-1, j2-1)) * 1000000 
            ;
          }
          dsigma11(j1-1, j1-1)=0;
          //The derivative  of sigma 12 wrt to one parent range has only one non null coefficient. They are stacked in a vector. 
          dsigma12_parents (j1-1) = (
            pow(1.000001 * range(NNarray(i,j1) - 1)    * range(NNarray(i,0) - 1)    , dim*  .25) * 
              pow(1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
              exp(-pow(sqeucdist(j1, 0)/ (1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5), .5)) 
            - sigma12(j1-1))  * 1000000 ;
          //The derivative  of sigma 12 wrt to child range is a dense vector 
          dsigma12_child (j1-1) = (
            pow(range(NNarray(i,j1) - 1)    * 1.000001 * range(NNarray(i,0) - 1)    , dim*  .25) * 
              pow(range(NNarray(i,j1) - 1)*.5 + 1.000001 * range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
              exp(-pow(sqeucdist(j1, 0)/ (range(NNarray(i,j1) - 1)*.5 + 1.000001 * range(NNarray(i,0) - 1)*.5), .5)) 
            - sigma12(j1-1))  * 1000000 ;
        }
        // computing gradient of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        grad1(i, 0, 0) = 
          //. a = 0, see article's annex
          arma::sum(salt % dsigma12_child) * pow(inverse_cond_sd, 3); // b
        //c = 0
        //rest of the coefficients
        arma::vec the_rest = 
          - agmis11 *  dsigma12_child * Linv(i, 0) // d
          //e = 0
          - salt * grad1(i, 0, 0); //f
          for(int j2 = 1; j2<bsize; j2++)
          {
            grad1(i, 0, j2) = the_rest(j2-1);
          }
          
          // case parents' ranges are differentiated 
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
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
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad1(i, j1, 0) = 
              // a = 0
              dsigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)     * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad1(i, j1, j2) = - agmis11(j2-1, j1-1) *  dsigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad1(i, j1, 0) ; //f
            }
          }
      }
    }
    grad(0) = grad1 ;
  }  
  if((covfun_name == "nonstationary_matern_isotropic"))
  {
    arma::vec range = exp(log_range) ; 
    arma::cube grad1(n,m, m);
    grad1.fill(0);
    // loop over every observation   
    #pragma omp parallel for
    for(int i=1; i<n; i++){
      int bsize = std::min(i+1,m);
      // first, fill in ysub, locsub, and X0 NOT in reverse order
      arma::mat locsub(bsize, dim);
      if(!sphere)
      {
        for(int j=0; j<bsize; j++){
          for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
        }
      }
      if(sphere)
      {
        // getting basis of tangent plane
        arma::mat tangent_plane_basis (3, 2);
        // parallel
        tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
        tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
        tangent_plane_basis (2, 0) = 0;
        // meridian
        tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
        tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
        tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
        
        arma::mat locsub_3D (bsize, 3);
        for(int j=0; j<bsize; j++){
          for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
        }
        locsub = locsub_3D * tangent_plane_basis ; 
        
      }
      // compute euclidean distance
      arma::mat sqeucdist(bsize, bsize);
      for(int j1=0; j1<bsize; j1++){
        for(int j2=0; j2<=j1; j2++){
          sqeucdist(j1, j2) = 0  ; 
          for(int k=0;k<dim;k++){ 
            sqeucdist(j1, j2) +=  pow((locsub(j1, k)-locsub(j2, k)), 2) ; 
          }
          sqeucdist(j2, j1) = sqeucdist(j1, j2) ; 
        }
        sqeucdist(j1, j1) = 0  ; 
      }
      
      // compute squexp covariance 
      arma::mat  sigma11(bsize-1, bsize-1);
      arma::vec  sigma12(bsize-1);
      for(int j1=1; j1<bsize; j1++){
        for(int j2=1; j2<=j1; j2++){
          double paciorek_dist = pow(sqeucdist(j1, j2)/(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5) , .5);
          sigma11 (j1-1, j2-1) = 
            pow(range(NNarray(i,j1) - 1) *    range(NNarray(i,j2) -1)    , dim*.25  ) * 
            pow(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5 , dim*(-.5)) * 
            exp(- paciorek_dist) * (paciorek_dist + 1.0);
          sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
          if(j1 == j2){sigma11 (j1-1, j2-1) = 1.0005;}
        }
        double paciorek_dist = pow(sqeucdist(j1, 0)/ (range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5), .5);
        sigma12 (j1-1) =                         
          pow(range(NNarray(i,j1) - 1)    * range(NNarray(i,0) - 1)    , dim*  .25) * 
          pow(range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
          exp(- paciorek_dist) * (paciorek_dist + 1.0);
      }
      // solving sigma11
      arma::mat agmis11  = arma::inv(sigma11) ;
      // computing a vector used everyvhere
      arma::mat salt  = agmis11 * sigma12 ;
      
      // computing Vecchia approx itself
      double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
      Linv(i, 0) = inverse_cond_sd ;
      for(int j=1; j<bsize; j++){
        Linv(i, j) = - salt(j-1) * inverse_cond_sd;
      }
      
      if(compute_derivative)
      {
        // compute squexp covariance derivatives using finite differences
        arma::mat dsigma11(bsize-1, bsize-1);
        arma::vec dsigma12_parents(bsize-1);
        arma::vec dsigma12_child(bsize-1);
        for(int j1=1; j1<bsize; j1++){
          for(int j2=1; j2<bsize; j2++){
            //the derivative of sigma11 is a symmetric cross matrix with non-null coefficients only on one row and column. Instead of computing the cross matrix, just one row is computed. The rows are stacked in a matrix. 
            if(j1 == j2)
            {
              dsigma11 (j1-1, j2-1) = 0;
            }
            if(j1 != j2)
            {
              double paciorek_dist = pow(sqeucdist(j1, j2)/ (1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5), .5);
              dsigma11 (j1-1, j2-1) = (
                pow(1.000001 * range(NNarray(i,j1) - 1) *    range(NNarray(i,j2) -1)    , dim*.25  ) * 
                  pow(1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,j2) -1)*.5 , dim*(-.5)) * 
                  exp(- paciorek_dist) * (paciorek_dist + 1.0)
                - sigma11 (j1-1, j2-1)) * 1000000 
              ;
            }
          }
          dsigma11(j1-1, j1-1)=0;
          //The derivative  of sigma 12 wrt to one parent range has only one non null coefficient. They are stacked in a vector. 
          double paciorek_dist = pow(sqeucdist(j1, 0)/ (1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5), .5); 
          dsigma12_parents (j1-1) = (
            pow(1.000001 * range(NNarray(i,j1) - 1)    * range(NNarray(i,0) - 1)    , dim*  .25) * 
              pow(1.000001 * range(NNarray(i,j1) - 1)*.5 + range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
              exp(- paciorek_dist) * (paciorek_dist + 1.0) 
            - sigma12(j1-1))  * 1000000 ;
          //The derivative  of sigma 12 wrt to child range is a dense vector 
          paciorek_dist = pow(sqeucdist(j1, 0)/ (range(NNarray(i,j1) - 1)*.5 + 1.000001 * range(NNarray(i,0) - 1)*.5), .5);
          dsigma12_child (j1-1) = (
            pow(range(NNarray(i,j1) - 1)    * 1.000001 * range(NNarray(i,0) - 1)    , dim*  .25) * 
              pow(range(NNarray(i,j1) - 1)*.5 + 1.000001 * range(NNarray(i,0) - 1)*.5 , dim*(-.5)) * 
              exp(- paciorek_dist) * (paciorek_dist + 1.0) 
            - sigma12(j1-1))  * 1000000 ;
        }
        // computing gradient of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        grad1(i, 0, 0) = 
          //. a = 0, see article's annex
          arma::sum(salt % dsigma12_child) * pow(inverse_cond_sd, 3); // b
        //c = 0
        //rest of the coefficients
        arma::vec the_rest = 
          - agmis11 *  dsigma12_child * Linv(i, 0) // d
          //e = 0
          - salt * grad1(i, 0, 0); //f
          for(int j2 = 1; j2<bsize; j2++)
          {
            grad1(i, 0, j2) = the_rest(j2-1);
          }
          
          // case parents' ranges are differentiated 
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
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
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad1(i, j1, 0) = 
              // a = 0
              dsigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)     * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad1(i, j1, j2) = - agmis11(j2-1, j1-1) *  dsigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad1(i, j1, 0) ; //f
            }
          }
      }
    }
    grad(0) = grad1 ;
  }
  
  //case of anisotropic range parameters
  if(covfun_name == "nonstationary_exponential_anisotropic")
  {
    if(compute_derivative)
    {
      // range matrices computed from log_range
      arma::cube range_matrices(2,2,n); 
      // computed from log_range + epsilon, used in differentiation
      arma::cube d1range_matrices(2,2,n); 
      arma::cube d2range_matrices(2,2,n); 
      arma::cube d3range_matrices(2,2,n); 
      d1range_matrices.fill(0); 
      d2range_matrices.fill(0); 
      d3range_matrices.fill(0); 
      arma::mat determinants(n, 4);
      
      for(int i=0; i<n; i++){
        arma::mat logmat(2, 2) ; 
        logmat(0, 0) = log_range(i, 0);
        logmat(1, 1) = log_range(i, 1);
        logmat(1, 0) = log_range(i, 2)/pow(2, .5);
        logmat(0, 1) = log_range(i, 2)/pow(2, .5);
        range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) + 0.000001 ; 
        d1range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) - 0.000001 ; 
        logmat(1, 1) = logmat(1, 1) + 0.000001 ; 
        d2range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(1, 1) = logmat(1, 1) - 0.000001 ; 
        logmat(0, 1) = logmat(0, 1) + 0.000001 * pow(2, -.5) ; 
        logmat(1, 0) = logmat(0, 1); 
        d3range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        //determinants(i, 0) = arma::det(range_matrices.slice(i)  )   ; 
        //determinants(i, 1) = arma::det(d1range_matrices.slice(i))  ; 
        //determinants(i, 2) = arma::det(d2range_matrices.slice(i))  ; 
        //determinants(i, 3) = arma::det(d3range_matrices.slice(i))  ; 
      }
      
      // the 3 gradients
      arma::cube grad1(n,m, m);
      arma::cube grad2(n,m, m);
      arma::cube grad3(n,m, m);
      grad1.fill(0);
      grad2.fill(0);
      grad3.fill(0);
      
      // loop over every observation   
#pragma omp parallel for
      for(int i=1; i<n; i++){
        //Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);
        
        // first, fill in ysub, locsub, and X0 NOT in reverse order
        arma::mat locsub(bsize, dim);
        if(!sphere)
        {
          for(int j=0; j<bsize; j++){
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
          }
        }
        if(sphere)
        {
          // getting basis of tangent plane
          arma::mat tangent_plane_basis (3, 2);
          // parallel
          tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 0) = 0;
          // meridian
          tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
          
          arma::mat locsub_3D (bsize, 3);
          for(int j=0; j<bsize; j++){
            for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
          }
          locsub = locsub_3D * tangent_plane_basis ; 
        }
        
        
        
        // covariance 
        arma::mat  sigma11(bsize-1, bsize-1);
        arma::vec  sigma12(bsize-1);
        // covariance derivatives using finite differences
        arma::mat d1sigma11(bsize-1, bsize-1);
        arma::mat d2sigma11(bsize-1, bsize-1);
        arma::mat d3sigma11(bsize-1, bsize-1);
        arma::vec d1sigma12_parents (bsize-1);
        arma::vec d2sigma12_parents (bsize-1);
        arma::vec d3sigma12_parents (bsize-1);
        arma::vec d1sigma12_child   (bsize-1);
        arma::vec d2sigma12_child   (bsize-1);
        arma::vec d3sigma12_child   (bsize-1);
        
        
        for(int j1=1; j1<bsize; j1++){
          //derivative of sigma11
          for(int j2=1; j2<bsize; j2++){
            arma::mat locdiff = locsub.row(j2) - locsub.row(j1);
            arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            sigma11 (j1-1, j2-1) = 
              pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,j2) -1, 0), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-pow(mahala_dist(0, 0), .5));
            sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
            
            //the derivative of sigma11 is a symmetric cross matrix with non-null coefficients only on one row and column. Instead of computing the cross matrix, just one column is computed. The columns are stacked in a matrix. 
            // first derivative
            hybrid_range = d1range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            d1sigma11 (j1-1, j2-1) = 
              pow(arma::det(d1range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 1) * determinants(NNarray(i,j2) -1, 1), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-pow(mahala_dist(0, 0), .5));
            // second derivative
            hybrid_range = d2range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            d2sigma11 (j1-1, j2-1) = 
              pow(arma::det(d2range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 2) * determinants(NNarray(i,j2) -1, 2), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-pow(mahala_dist(0, 0), .5));
            // third derivative
            hybrid_range = d3range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            d3sigma11 (j1-1, j2-1) = 
              pow(arma::det(d3range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 3) * determinants(NNarray(i,j2) -1, 3), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-pow(mahala_dist(0, 0), .5));
          }
          arma::mat locdiff = locsub.row(0) - locsub.row(j1);
          
          // sigma 12
          arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          sigma12 (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,0) -1, 1), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          
          //The derivative  of sigma 12 wrt to one parent range has only one non null coefficient. They are stacked in a vector. 
          // first coordinate of log matrix range 
          hybrid_range = d1range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d1sigma12_parents (j1-1) =                         
            pow(arma::det(d1range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          // second coordinate of log matrix range 
          hybrid_range = d2range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d2sigma12_parents (j1-1) =                         
            pow(arma::det(d2range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          // third coordinate of log matrix range 
          hybrid_range = d3range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d3sigma12_parents (j1-1) =                         
            pow(arma::det(d3range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          
          //The derivative  of sigma 12 wrt to child range is dense
          // first coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d1range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d1sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d1range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          // second coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d2range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d2sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d2range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          // third coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d3range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          d3sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d3range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
          
        }
        // removing initial state and dividing by finite difference step
        d1sigma11 = (d1sigma11 - sigma11)/0.000001;
        d2sigma11 = (d2sigma11 - sigma11)/0.000001;
        d3sigma11 = (d3sigma11 - sigma11)/0.000001;
        d1sigma12_parents = (d1sigma12_parents - sigma12)/0.000001;
        d2sigma12_parents = (d2sigma12_parents - sigma12)/0.000001;
        d3sigma12_parents = (d3sigma12_parents - sigma12)/0.000001;
        d1sigma12_child   = (d1sigma12_child   - sigma12)/0.000001;
        d2sigma12_child   = (d2sigma12_child   - sigma12)/0.000001;
        d3sigma12_child   = (d3sigma12_child   - sigma12)/0.000001;
        
        // solving sigma11
        arma::mat agmis11  = arma::inv(sigma11) ;
        // computing a vector used everyvhere
        arma::mat salt  = agmis11 * sigma12 ;
        
        // computing Vecchia approx itself
        double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
        Linv(i, 0) = inverse_cond_sd ;
        for(int j=1; j<bsize; j++){
          Linv(i, j) = - salt(j-1) * inverse_cond_sd;
        }
        
        // computing grad1ient of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        grad1(i, 0, 0) = 
          //. a = 0, see article's annex
          arma::sum(salt % d1sigma12_child) * pow(inverse_cond_sd, 3); // b
        //c = 0
        grad2(i, 0, 0) =  arma::sum(salt % d2sigma12_child) * pow(inverse_cond_sd, 3);
        grad3(i, 0, 0) =  arma::sum(salt % d3sigma12_child) * pow(inverse_cond_sd, 3);
        //rest of the coefficients
        // first gradient
        arma::vec the_rest = 
          - agmis11 *  d1sigma12_child * Linv(i, 0) // d  //e = 0
          - salt * grad1(i, 0, 0); //f
          for(int j2 = 1; j2<bsize; j2++){grad1(i, 0, j2) = the_rest(j2-1);}
          // second gradient
          the_rest = - agmis11 *  d2sigma12_child * Linv(i, 0)  - salt * grad2(i, 0, 0); 
          for(int j2 = 1; j2<bsize; j2++){grad2(i, 0, j2) = the_rest(j2-1);}
          // third gradient
          the_rest = - agmis11 *  d3sigma12_child * Linv(i, 0)  - salt * grad3(i, 0, 0); 
          for(int j2 = 1; j2<bsize; j2++){grad3(i, 0, j2) = the_rest(j2-1);}
          
          // case parents' ranges are differentiated 
          // first gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              
              salt_dsigma11(j2)   += salt(j1-1)*d1sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d1sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad1(i, j1, 0) = // a = 0
              d1sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad1(i, j1, j2) = - agmis11(j2-1, j1-1) *  d1sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad1(i, j1, 0) ; //f
            }
          }
          // second gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)   += salt(j1-1)*d2sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d2sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad2(i, j1, 0) = // a = 0
              d2sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad2(i, j1, j2) = - agmis11(j2-1, j1-1) *  d2sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad2(i, j1, 0) ; //f
            }
          }
          // third gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)   += salt(j1-1)*d3sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d3sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad3(i, j1, 0) = // a = 0
              d3sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad3(i, j1, j2) = - agmis11(j2-1, j1-1) *  d3sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad3(i, j1, 0) ; //f
            }
          }
      }
      grad(0) = grad1 ;
      grad(1) = grad2 ;
      grad(2) = grad3 ;
    }
    if(!compute_derivative)
    {
      // range matrices computed from log_range
      arma::cube range_matrices(2,2,n); 
      // computed from log_range + epsilon, used in differentiation
      arma::cube d1range_matrices(2,2,n); 
      arma::cube d2range_matrices(2,2,n); 
      arma::cube d3range_matrices(2,2,n); 
      d1range_matrices.fill(0); 
      d2range_matrices.fill(0); 
      d3range_matrices.fill(0); 
      arma::mat determinants(n, 4);
      
      for(int i=0; i<n; i++){
        arma::mat logmat(2, 2) ; 
        logmat(0, 0) = log_range(i, 0);
        logmat(1, 1) = log_range(i, 1);
        logmat(1, 0) = log_range(i, 2)/pow(2, .5);
        logmat(0, 1) = log_range(i, 2)/pow(2, .5);
        range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        //determinants(i, 0) = arma::det(range_matrices.slice(i)  )   ; 
        //determinants(i, 1) = arma::det(d1range_matrices.slice(i))  ; 
        //determinants(i, 2) = arma::det(d2range_matrices.slice(i))  ; 
        //determinants(i, 3) = arma::det(d3range_matrices.slice(i))  ; 
      }
      
      // the 3 gradients
      arma::cube grad1(n,m, m);
      arma::cube grad2(n,m, m);
      arma::cube grad3(n,m, m);
      grad1.fill(0);
      grad2.fill(0);
      grad3.fill(0);
      
      // loop over every observation   
#pragma omp parallel for
      for(int i=1; i<n; i++){
        //Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);
        
        // first, fill in ysub, locsub, and X0 NOT in reverse order
        arma::mat locsub(bsize, dim);
        if(!sphere)
        {
          for(int j=0; j<bsize; j++){
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
          }
        }
        if(sphere)
        {
          // getting basis of tangent plane
          arma::mat tangent_plane_basis (3, 2);
          // parallel
          tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 0) = 0;
          // meridian
          tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
          
          arma::mat locsub_3D (bsize, 3);
          for(int j=0; j<bsize; j++){
            for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
          }
          locsub = locsub_3D * tangent_plane_basis ; 
        }
        
        
        
        // covariance 
        arma::mat  sigma11(bsize-1, bsize-1);
        arma::vec  sigma12(bsize-1);
        
        
        for(int j1=1; j1<bsize; j1++){
          //derivative of sigma11
          for(int j2=1; j2<bsize; j2++){
            arma::mat locdiff = locsub.row(j2) - locsub.row(j1);
            arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            sigma11 (j1-1, j2-1) = 
              pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,j2) -1, 0), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-pow(mahala_dist(0, 0), .5));
            sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
          }
          arma::mat locdiff = locsub.row(0) - locsub.row(j1);
          
          // sigma 12
          arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          sigma12 (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,0) -1, 1), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-pow(mahala_dist(0, 0), .5));
        }
        // solving sigma11
        arma::mat agmis11  = arma::inv(sigma11) ;
        // computing a vector used everyvhere
        arma::mat salt  = agmis11 * sigma12 ;
        
        // computing Vecchia approx itself
        double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
        Linv(i, 0) = inverse_cond_sd ;
        for(int j=1; j<bsize; j++){
          Linv(i, j) = - salt(j-1) * inverse_cond_sd;
        }
        
      }
      grad(0) = grad1 ;
      grad(1) = grad2 ;
      grad(2) = grad3 ;
    }
  }
  if((covfun_name == "nonstationary_matern_anisotropic")){
    if(compute_derivative)
    {
      
      
      // range matrices computed from log_range
      arma::cube range_matrices(2,2,n); 
      // computed from log_range + epsilon, used in differentiation
      arma::cube d1range_matrices(2,2,n); 
      arma::cube d2range_matrices(2,2,n); 
      arma::cube d3range_matrices(2,2,n); 
      d1range_matrices.fill(0); 
      d2range_matrices.fill(0); 
      d3range_matrices.fill(0); 
      arma::mat determinants(n, 4);
      
      for(int i=0; i<n; i++){
        arma::mat logmat(2, 2) ; 
        logmat(0, 0) = log_range(i, 0);
        logmat(1, 1) = log_range(i, 1);
        logmat(1, 0) = log_range(i, 2)/pow(2, .5);
        logmat(0, 1) = log_range(i, 2)/pow(2, .5);
        range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) + 0.000001 ; 
        d1range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) - 0.000001 ; 
        logmat(1, 1) = logmat(1, 1) + 0.000001 ; 
        d2range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(1, 1) = logmat(1, 1) - 0.000001 ; 
        logmat(0, 1) = logmat(0, 1) + 0.000001 * pow(2, -.5) ; 
        logmat(1, 0) = logmat(0, 1); 
        d3range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        //determinants(i, 0) = arma::det(range_matrices.slice(i)  )   ; 
        //determinants(i, 1) = arma::det(d1range_matrices.slice(i))  ; 
        //determinants(i, 2) = arma::det(d2range_matrices.slice(i))  ; 
        //determinants(i, 3) = arma::det(d3range_matrices.slice(i))  ; 
      }
      
      // the 3 gradients
      arma::cube grad1(n,m, m);
      arma::cube grad2(n,m, m);
      arma::cube grad3(n,m, m);
      grad1.fill(0);
      grad2.fill(0);
      grad3.fill(0);
      
      // loop over every observation   
#pragma omp parallel for
      for(int i=1; i<n; i++){
        //Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);
        
        // first, fill in ysub, locsub, and X0 NOT in reverse order
        arma::mat locsub(bsize, dim);
        if(!sphere)
        {
          for(int j=0; j<bsize; j++){
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
          }
        }
        if(sphere)
        {
          // getting basis of tangent plane
          arma::mat tangent_plane_basis (3, 2);
          // parallel
          tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 0) = 0;
          // meridian
          tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
          
          arma::mat locsub_3D (bsize, 3);
          for(int j=0; j<bsize; j++){
            for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
          }
          locsub = locsub_3D * tangent_plane_basis ; 
        }
        
        
        
        // covariance 
        arma::mat  sigma11(bsize-1, bsize-1);
        arma::vec  sigma12(bsize-1);
        // covariance derivatives using finite differences
        arma::mat d1sigma11(bsize-1, bsize-1);
        arma::mat d2sigma11(bsize-1, bsize-1);
        arma::mat d3sigma11(bsize-1, bsize-1);
        arma::vec d1sigma12_parents (bsize-1);
        arma::vec d2sigma12_parents (bsize-1);
        arma::vec d3sigma12_parents (bsize-1);
        arma::vec d1sigma12_child   (bsize-1);
        arma::vec d2sigma12_child   (bsize-1);
        arma::vec d3sigma12_child   (bsize-1);
        
        
        for(int j1=1; j1<bsize; j1++){
          //derivative of sigma11
          for(int j2=1; j2<bsize; j2++){
            arma::mat locdiff = locsub.row(j2) - locsub.row(j1);
            arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
            sigma11 (j1-1, j2-1) = 
              pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,j2) -1, 0), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
            sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
            
            //the derivative of sigma11 is a symmetric cross matrix with non-null coefficients only on one row and column. Instead of computing the cross matrix, just one column is computed. The columns are stacked in a matrix. 
            // first derivative
            hybrid_range = d1range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
            d1sigma11 (j1-1, j2-1) = 
              pow(arma::det(d1range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 1) * determinants(NNarray(i,j2) -1, 1), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
            // second derivative
            hybrid_range = d2range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
            d2sigma11 (j1-1, j2-1) = 
              pow(arma::det(d2range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 2) * determinants(NNarray(i,j2) -1, 2), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
            // third derivative
            hybrid_range = d3range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
            d3sigma11 (j1-1, j2-1) = 
              pow(arma::det(d3range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 3) * determinants(NNarray(i,j2) -1, 3), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          }
          arma::mat locdiff = locsub.row(0) - locsub.row(j1);
          
          // sigma 12
          arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          sigma12 (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,0) -1, 1), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          
          //The derivative  of sigma 12 wrt to one parent range has only one non null coefficient. They are stacked in a vector. 
          // first coordinate of log matrix range 
          hybrid_range = d1range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d1sigma12_parents (j1-1) =                         
            pow(arma::det(d1range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          // second coordinate of log matrix range 
          hybrid_range = d2range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d2sigma12_parents (j1-1) =                         
            pow(arma::det(d2range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          // third coordinate of log matrix range 
          hybrid_range = d3range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d3sigma12_parents (j1-1) =                         
            pow(arma::det(d3range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          
          //The derivative  of sigma 12 wrt to child range is dense
          // first coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d1range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d1sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d1range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          // second coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d2range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d2sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d2range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          // third coordinate of log matrix range 
          hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + d3range_matrices.slice(NNarray(i,0) -1)*.5 ;
          mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          d3sigma12_child (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(d3range_matrices.slice(NNarray(i,0) -1)), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          
        }
        // removing initial state and dividing by finite difference step
        d1sigma11 = (d1sigma11 - sigma11)/0.000001;
        d2sigma11 = (d2sigma11 - sigma11)/0.000001;
        d3sigma11 = (d3sigma11 - sigma11)/0.000001;
        d1sigma12_parents = (d1sigma12_parents - sigma12)/0.000001;
        d2sigma12_parents = (d2sigma12_parents - sigma12)/0.000001;
        d3sigma12_parents = (d3sigma12_parents - sigma12)/0.000001;
        d1sigma12_child   = (d1sigma12_child   - sigma12)/0.000001;
        d2sigma12_child   = (d2sigma12_child   - sigma12)/0.000001;
        d3sigma12_child   = (d3sigma12_child   - sigma12)/0.000001;
        
        // solving sigma11
        arma::mat agmis11  = arma::inv(sigma11) ;
        // computing a vector used everyvhere
        arma::mat salt  = agmis11 * sigma12 ;
        
        // computing Vecchia approx itself
        double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
        Linv(i, 0) = inverse_cond_sd ;
        for(int j=1; j<bsize; j++){
          Linv(i, j) = - salt(j-1) * inverse_cond_sd;
        }
        
        // computing grad1ient of Vecchia approx 
        
        // case child range is differentiated 
        //diagonal coefficient
        grad1(i, 0, 0) = 
          //. a = 0, see article's annex
          arma::sum(salt % d1sigma12_child) * pow(inverse_cond_sd, 3); // b
        //c = 0
        grad2(i, 0, 0) =  arma::sum(salt % d2sigma12_child) * pow(inverse_cond_sd, 3);
        grad3(i, 0, 0) =  arma::sum(salt % d3sigma12_child) * pow(inverse_cond_sd, 3);
        //rest of the coefficients
        // first gradient
        arma::vec the_rest = 
          - agmis11 *  d1sigma12_child * Linv(i, 0) // d  //e = 0
          - salt * grad1(i, 0, 0); //f
          for(int j2 = 1; j2<bsize; j2++){grad1(i, 0, j2) = the_rest(j2-1);}
          // second gradient
          the_rest = - agmis11 *  d2sigma12_child * Linv(i, 0)  - salt * grad2(i, 0, 0); 
          for(int j2 = 1; j2<bsize; j2++){grad2(i, 0, j2) = the_rest(j2-1);}
          // third gradient
          the_rest = - agmis11 *  d3sigma12_child * Linv(i, 0)  - salt * grad3(i, 0, 0); 
          for(int j2 = 1; j2<bsize; j2++){grad3(i, 0, j2) = the_rest(j2-1);}
          
          // case parents' ranges are differentiated 
          // first gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              
              salt_dsigma11(j2)   += salt(j1-1)*d1sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d1sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad1(i, j1, 0) = // a = 0
              d1sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad1(i, j1, j2) = - agmis11(j2-1, j1-1) *  d1sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad1(i, j1, 0) ; //f
            }
          }
          // second gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)   += salt(j1-1)*d2sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d2sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad2(i, j1, 0) = // a = 0
              d2sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad2(i, j1, j2) = - agmis11(j2-1, j1-1) *  d2sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad2(i, j1, 0) ; //f
            }
          }
          // third gradient
          for(int j1 = 1; j1<bsize;  j1++)
          {
            arma::vec salt_dsigma11(bsize-1);
            salt_dsigma11.fill(0);
            // computing stuff 
            for(int j2 = 0; j2<bsize-1; j2++)
            {
              salt_dsigma11(j2)   += salt(j1-1)*d3sigma11(j1-1, j2); //  horizontal bar of the cross
              salt_dsigma11(j1-1) += salt(j2  )*d3sigma11(j1-1, j2); //  vertical bar of the cross
            }
            arma::vec salt_dsigma11_agmis11 = agmis11 * salt_dsigma11 ; 
            //diagonal coefficient
            grad3(i, j1, 0) = // a = 0
              d3sigma12_parents(j1-1) * salt(j1-1)  * pow(inverse_cond_sd, 3) // b 
              -arma::sum(salt_dsigma11 % salt)      * pow(inverse_cond_sd, 3)/2 ; // c
            
            //rest of the coefficients
            for(int j2 = 1; j2<bsize; j2++)
            {
              grad3(i, j1, j2) = - agmis11(j2-1, j1-1) *  d3sigma12_parents(j1-1) * Linv(i, 0) // d
              + salt_dsigma11_agmis11(j2-1)                   * Linv(i, 0) // e
              - salt(j2-1)                                    * grad3(i, j1, 0) ; //f
            }
          }
      }
      grad(0) = grad1 ;
      grad(1) = grad2 ;
      grad(2) = grad3 ;
      
      
    }
    if(!compute_derivative)
    {
      // range matrices computed from log_range
      arma::cube range_matrices(2,2,n); 
      // computed from log_range + epsilon, used in differentiation
      arma::cube d1range_matrices(2,2,n); 
      arma::cube d2range_matrices(2,2,n); 
      arma::cube d3range_matrices(2,2,n); 
      d1range_matrices.fill(0); 
      d2range_matrices.fill(0); 
      d3range_matrices.fill(0); 
      arma::mat determinants(n, 4);
      
      for(int i=0; i<n; i++){
        arma::mat logmat(2, 2) ; 
        logmat(0, 0) = log_range(i, 0);
        logmat(1, 1) = log_range(i, 1);
        logmat(1, 0) = log_range(i, 2)/pow(2, .5);
        logmat(0, 1) = log_range(i, 2)/pow(2, .5);
        range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) + 0.000001 ; 
        d1range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(0, 0) = logmat(0, 0) - 0.000001 ; 
        logmat(1, 1) = logmat(1, 1) + 0.000001 ; 
        d2range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        logmat(1, 1) = logmat(1, 1) - 0.000001 ; 
        logmat(0, 1) = logmat(0, 1) + 0.000001 * pow(2, -.5) ; 
        logmat(1, 0) = logmat(0, 1); 
        d3range_matrices.slice(i) = expmat_sym(logmat) ; 
        
        //determinants(i, 0) = arma::det(range_matrices.slice(i)  )   ; 
        //determinants(i, 1) = arma::det(d1range_matrices.slice(i))  ; 
        //determinants(i, 2) = arma::det(d2range_matrices.slice(i))  ; 
        //determinants(i, 3) = arma::det(d3range_matrices.slice(i))  ; 
      }
      
      // the 3 gradients
      arma::cube grad1(n,m, m);
      arma::cube grad2(n,m, m);
      arma::cube grad3(n,m, m);
      grad1.fill(0);
      grad2.fill(0);
      grad3.fill(0);
      
      // loop over every observation   
#pragma omp parallel for
      for(int i=1; i<n; i++){
        //Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);
        
        // first, fill in ysub, locsub, and X0 NOT in reverse order
        arma::mat locsub(bsize, dim);
        if(!sphere)
        {
          for(int j=0; j<bsize; j++){
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,j) - 1, k); }
          }
        }
        if(sphere)
        {
          // getting basis of tangent plane
          arma::mat tangent_plane_basis (3, 2);
          // parallel
          tangent_plane_basis (0, 0) =  sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 0) = -cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 0) = 0;
          // meridian
          tangent_plane_basis (0, 1) =  - sin(M_PI /180 *locs(i, 1)) * cos(M_PI /180 *locs(i, 0));
          tangent_plane_basis (1, 1) =  - sin(M_PI /180 *locs(i, 1)) * sin(M_PI /180 *locs(i, 0));
          tangent_plane_basis (2, 1) =    cos(M_PI /180 *locs(i, 1)) ;
          
          arma::mat locsub_3D (bsize, 3);
          for(int j=0; j<bsize; j++){
            for(int k=0;k<3;k++){ locsub_3D(j,k) = locs_3D( NNarray(i,j) - 1, k); }
          }
          locsub = locsub_3D * tangent_plane_basis ; 
        }
        
        
        
        // covariance 
        arma::mat  sigma11(bsize-1, bsize-1);
        arma::vec  sigma12(bsize-1);
        
        for(int j1=1; j1<bsize; j1++){
          //derivative of sigma11
          for(int j2=1; j2<bsize; j2++){
            arma::mat locdiff = locsub.row(j2) - locsub.row(j1);
            arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,j2) -1)*.5 ;
            
            arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
            mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
            sigma11 (j1-1, j2-1) = 
              pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,j2) -1)), .25) * 
              //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,j2) -1, 0), .25) * 
              pow(arma::det(hybrid_range), -.5) * 
              exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
            sigma11 (j2-1, j1-1) = sigma11 (j1-1, j2-1) ; 
            
          }
          arma::mat locdiff = locsub.row(0) - locsub.row(j1);
          
          // sigma 12
          arma::mat hybrid_range = range_matrices.slice(NNarray(i,j1) - 1)*.5 + range_matrices.slice(NNarray(i,0) -1)*.5 ;
          arma::mat mahala_dist = locdiff * arma::inv(hybrid_range) * locdiff.t();
          mahala_dist(0, 0) = pow(mahala_dist(0, 0), .5);
          sigma12 (j1-1) =                         
            pow(arma::det(range_matrices.slice(NNarray(i,j1) - 1)) * arma::det(range_matrices.slice(NNarray(i,0) -1)), .25) * 
            //pow(determinants(NNarray(i,j1) - 1, 0) * determinants(NNarray(i,0) -1, 1), .25) * 
            pow(arma::det(hybrid_range), -.5) * 
            exp(-mahala_dist(0, 0)) * (1+ mahala_dist(0,0));
          
        }
        // solving sigma11
        arma::mat agmis11  = arma::inv(sigma11) ;
        // computing a vector used everyvhere
        arma::mat salt  = agmis11 * sigma12 ;
        
        // computing Vecchia approx itself
        double inverse_cond_sd = pow(1- sum(salt % sigma12), -.5);
        Linv(i, 0) = inverse_cond_sd ;
        for(int j=1; j<bsize; j++){
          Linv(i, j) = - salt(j-1) * inverse_cond_sd;
        }
      }
      grad(0) = grad1 ;
      grad(1) = grad2 ;
      grad(2) = grad3 ;
    }
  }
  
  Rcpp::List out(2);    
  out(0)= Linv;
  out(1)= grad;
  return(out);
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



/// #include <RcppArmadillo.h>
/// //[[Rcpp::depends(RcppArmadillo)]]
/// //[[Rcpp::export]]
/// arma::mat derivative_sandwiches_cpp 
/// (
///     List derivatives, 
///     arma::vec left_vector, 
///     arma::vec right_vector, 
///     Rcpp::IntegerMatrix NNarray 
/// )
/// {
///   int  n = NNarray.nrow() ;
///   int  m = NNarray.ncol() ;
///   int  d = derivatives.length() ;
///   arma::mat res(n, d);
///   res.fill(0);
/// # pragma omp parallel for
///   for(int d_idx = 0; d_idx<d; d_idx++){
///     arma::cube derivative = derivatives(d_idx);
///     // looping over the row idx of the derivative of tilde R
///     // # pragma omp parallel
///     for(int row_idx = 0; row_idx<n; row_idx++){
///       int bsize = std::min(row_idx+1,m);
///       // looping over the index of the covariance parameter wrt which the derivative has been computed  (retrieved through NNarray)
///       for(int covparm_idx = 0; covparm_idx < bsize; covparm_idx++)
///       {
///         // looping over the column idx of the derivative of tilde R (column idx retrieved through NNarray)
///         for(int col_idx = 0; col_idx < bsize; col_idx++)
///         {
///            derivative(row_idx, covparm_idx, col_idx)
///           * left_vector(row_idx) 
///           * right_vector(NNarray(row_idx, col_idx)-1); 
///         }
///         res(NNarray(row_idx, covparm_idx)-1,d_idx)+= increment
///       }
///       
///     }
///   }
/// 
///   return(res);
/// }



