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
double TVBS_pmvnorm_cpp(Eigen::VectorXd upper, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, int N, int h, int s);

struct TVBS_Matrix_two TVBS_LDLT_decomp_cpp(Eigen::MatrixXd Sigma);                        // create a function to calculate two-block LDL decomposition
struct TVBS_Matrix_two TVBS_grad_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind);                        // create a function to calculate two-block LDL decomposition
struct TVBS_Matrix_two TVBS_Hess_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind1, Eigen::VectorXd ind2);                        // create a function to calculate two-block LDL decomposition  

struct TVBS_Matrix_two TVBS_LDLT_update_cpp(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, int m);
struct TVBS_Matrix_two TVBS_LDLT_update_grad(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd dL_k, Eigen::MatrixXd dD_k, Eigen::MatrixXd dOmega, int m);   
struct TVBS_Matrix_two TVBS_LDLT_update_Hess(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd d1L_k, Eigen::MatrixXd d1D_k, Eigen::MatrixXd d1Omega, Eigen::MatrixXd d2L_k, Eigen::MatrixXd d2D_k, Eigen::MatrixXd d2Omega, Eigen::MatrixXd hL_k, Eigen::MatrixXd hD_k, Eigen::MatrixXd hOmega, int m);  
  
struct TVBS_Matrix_two TVBS_truncate_bi_normal_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point, double tol);

// gradient of truncation function
struct TVBS_Matrix_two TVBS_truncate_bi_normal_gradient_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point);

// function to truncate mu for quadrovariate normal distribution
struct TVBS_Matrix_two TVBS_mu_l_trunc_bivariate_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd tr_point);
Eigen::MatrixXd TVBS_Jacobian_biv_norm_trunc(Eigen::VectorXd w, double rho);

struct Hessian_trunc_struct TVBS_Hessian_biv_norm_trunc(Eigen::VectorXd w, double rho);
Eigen::MatrixXd TVBS_Hessian_biv_norm_trunc_R(Eigen::VectorXd w, double rho, int out);
Eigen::MatrixXd TVBS_Jacobian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);
struct Hessian_trunc_struct TVBS_Hessian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);
Eigen::MatrixXd TVBS_Hessian_trunc_gen_R(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma, int out);

Eigen::VectorXi TVBS_logic2position_cpp(Eigen::VectorXd vec_logical); 
Eigen::MatrixXd TVBS_vec_2_upper_diag_cpp(Eigen::VectorXd vec_in, int diag_inc);
Eigen::VectorXd TVBS_lower_tri_entries_cpp(Eigen::MatrixXd x, int diag_inc);

Eigen::VectorXd TVBS_sincs_grad_cpp(double x);

Eigen::VectorXd TVBS_pntgnd_grad_cpp(double ba, double bb, double bc, double ra, double rb, double r, double rr); 

Eigen::MatrixXd TVBS_vec_symmetry_cpp(int r); 

Eigen::MatrixXd TVBS_g_asym_to_sym_cpp(Eigen::MatrixXd x); 

Eigen::MatrixXd TVBS_g_inverse_cpp(Eigen::MatrixXd x);


Eigen::MatrixXd TVBS_g_a_omega_b_cpp(Eigen::MatrixXd x1, Eigen::MatrixXd x2, int omsymmetric, int omdiagonal); 


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
struct TVBS_Matrix_two TVBS_grad_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);
struct TVBS_Matrix_three TVBS_grad_non_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);
struct TVBS_Matrix_two TVBS_grad_cdf_qvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr);
struct TVBS_Matrix_two TVBS_grad_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);
struct TVBS_Matrix_three TVBS_grad_non_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);
Eigen::VectorXd cal_grad_rhokl(Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::VectorXd He_onetofour);

Eigen::VectorXd TVBS_grad_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::MatrixXd Hessian_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);

struct TVBS_Vector_three TVBS_pdf_mvna_tvbs_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
struct TVBS_Vector_four TVBS_pdf_mvn_analytic_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x, int log_out);

Eigen::VectorXd TVBS_v2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
double TVBS_pv2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);

Eigen::VectorXd TVBS_grad(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
Eigen::MatrixXd TVBS_hess_new(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat);



