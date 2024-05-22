#include <Rcpp.h>
#include <RcppEigen.h>



using namespace Rcpp;


//[[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies

//#ifndef __TOMS462_INCLUDED__
//#define __TOMS462_INCLUDED__
#include "toms462.h"                                          // allows bivariate normal cdf calculations
#include "distrib.h"
//#endif



Eigen::MatrixXd TVBS_g_a_omega_b_cpp(Eigen::MatrixXd x1, Eigen::MatrixXd x2); //, int omsymmetric, int omdiagonal); 

Eigen::MatrixXd TVBS_g_cholesky_cov_cpp(Eigen::MatrixXd Sigma);

Eigen::MatrixXd TVBS_g_cholesky_cor_cpp(Eigen::MatrixXd Sigma);

// input: x_norm, Sigma_diag_sqrt, Cor_mat
struct TVBS_Matrix_two TVBS_grad_cor_cov_cpp(Eigen::VectorXd x_norm, Eigen::VectorXd Sigma_diag_sqrt, Eigen::MatrixXd Cor_mat);

// issues arrise for some cases due to divergence to behaviour compared to R code (when tempselmat length is not equal to number of rows/columns to select from)
struct TVBS_Matrix_two TVBS_g_both_x_omega_x_cpp(Eigen::MatrixXd Mat_1, Eigen::MatrixXd Mat_2);

struct TVBS_Matrix_two TVBS_g_cond_cov_cpp(Eigen::MatrixXd Id_mat, Eigen::MatrixXd Sigma);

struct TVBS_Matrix_three TVBS_g_cond_cov_trunc_cpp(Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x);

struct TVBS_Matrix_four TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd Id_mat, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x);



struct TVBS_double_three TVBS_grad_cdf_bvn_cpp(Eigen::VectorXd x_temp, Eigen::MatrixXd Rho);
struct TVBS_double_three TVBS_grad_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd w, Eigen::MatrixXd Rho);
struct TVBS_Matrix_three TVBS_grad_non_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);

struct TVBS_Matrix_two TVBS_grad_cdf_tvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr);

// [[Rcpp::plugins(cpp11)]]

struct Hessian_trunc_struct TVBS_Hessian_biv_norm_trunc(Eigen::VectorXd w, double rho);
Eigen::MatrixXd TVBS_Hessian_biv_norm_trunc_R(Eigen::VectorXd w, double rho, int out);
Eigen::MatrixXd TVBS_Jacobian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);


struct Hessian_trunc_struct TVBS_Hessian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);
