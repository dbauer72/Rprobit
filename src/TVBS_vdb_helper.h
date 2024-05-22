#include <Rcpp.h>
#include <RcppEigen.h>



using namespace Rcpp;
//TVBS Approximation
//INPUTS
//- x: a vector containing the upper integration limits
//- mu: a vector containing the means of the distribution
//- Sigma: a covariance matrix
//- ll: 0 if Probability should be the output, 1 if the logarithm should be returned
//OUTPUT
//- TVBS: (log)-probability


//[[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies

//#ifndef __TOMS462_INCLUDED__
//#define __TOMS462_INCLUDED__
#include "toms462.h"                                          // allows bivariate normal cdf calculations
#include "distrib.h"
//#endif


// [[Rcpp::plugins(cpp11)]]


double cal_delta(double w0, double w1, double rho);
Eigen::MatrixXd Hessian_delta(double w0, double w1, double rho);
double TVBS_std_normal_cdf_inv(Eigen::VectorXd w, Eigen::VectorXd b, Eigen::MatrixXd L);
struct TVBS_Vec_Mat TVBS_reorder_chol_cpp(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
// reorder variables for better numerical results
struct TVBS_Vec_Mat_Vec TVBS_reorder_sig_cpp(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::VectorXd TVBS_biv_norm_trunc(Eigen::VectorXd w, double rho);
Eigen::VectorXd TVBS_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);
double tildePhithree(Eigen::VectorXd x, Eigen::MatrixXd Sigma);
double tildePhifour(Eigen::VectorXd x, Eigen::MatrixXd Sigma);
double TVBS_pmvnorm_cpp(Eigen::VectorXd upper, Eigen::VectorXd mu, Eigen::MatrixXd Sigma);

struct TVBS_Matrix_two TVBS_LDLT_decomp_cpp(Eigen::MatrixXd Sigma);                        // create a function to calculate two-block LDL decomposition
struct TVBS_Matrix_two TVBS_grad_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind);                        // create a function to calculate two-block LDL decomposition
struct TVBS_Matrix_two TVBS_Hess_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind1, Eigen::VectorXd ind2);                        // create a function to calculate two-block LDL decomposition  

struct TVBS_Matrix_two TVBS_LDLT_update_cpp(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, int m);
struct TVBS_Matrix_two TVBS_LDLT_update_grad(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd dL_k, Eigen::MatrixXd dD_k, Eigen::MatrixXd dOmega, int m);   
struct TVBS_Matrix_two TVBS_LDLT_update_Hess(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd d1L_k, Eigen::MatrixXd d1D_k, Eigen::MatrixXd d1Omega, Eigen::MatrixXd d2L_k, Eigen::MatrixXd d2D_k, Eigen::MatrixXd d2Omega, Eigen::MatrixXd hL_k, Eigen::MatrixXd hD_k, Eigen::MatrixXd hOmega, int m);  



