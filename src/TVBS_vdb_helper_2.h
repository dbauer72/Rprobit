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

struct TVBS_Matrix_two TVBS_truncate_bi_normal_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point, double tol);

// gradient of truncation function
struct TVBS_Matrix_two TVBS_truncate_bi_normal_gradient_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point);

// function to truncate mu for quadrovariate normal distribution
struct TVBS_Matrix_two TVBS_mu_l_trunc_bivariate_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd tr_point);
Eigen::MatrixXd TVBS_Jacobian_biv_norm_trunc(Eigen::VectorXd w, double rho);

Eigen::VectorXi TVBS_logic2position_cpp(Eigen::VectorXd vec_logical); 
Eigen::MatrixXd TVBS_vec_2_upper_diag_cpp(Eigen::VectorXd vec_in, int diag_inc);

Eigen::VectorXd TVBS_lower_tri_entries_cpp(Eigen::MatrixXd x, int diag_inc);
Eigen::VectorXd TVBS_sincs_grad_cpp(double x);

Eigen::VectorXd TVBS_pntgnd_grad_cpp(double ba, double bb, double bc, double ra, double rb, double r, double rr); 

Eigen::MatrixXd TVBS_vec_symmetry_cpp(int r); 



Eigen::MatrixXd TVBS_g_asym_to_sym_cpp(Eigen::MatrixXd x); 

Eigen::MatrixXd TVBS_g_inverse_cpp(Eigen::MatrixXd x);
