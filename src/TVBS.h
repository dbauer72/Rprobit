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
#include "toms462.h" // allows bivariate normal cdf calculations
#include "TVBS_vdb.h"
//#endif


// [[Rcpp::plugins(cpp11)]]



double TVBS_std_normal_cdf_cpp(double x);                                 // create a function for the standard normal cumulative distribution function
double TVBS_std_normal_pdf_cpp(double x);                                 // create a function for the standard normal probability distribution function
double TVBS_std_normal_cdf_inv_1(double y);
Eigen::VectorXd TVBS_std_normal_cdf_inv(Eigen::VectorXd y);

double TVBS_std_normal_cdf_inv(Eigen::VectorXd w, Eigen::VectorXd b, Eigen::MatrixXd L);
Eigen::VectorXi TVBS_v_Korobov_cpp(int h, int K, int N);
Eigen::MatrixXd TVBS_L_N_shift_cpp(int h, int K, int N, int s);

double TVBS_pmvnorm_old_cpp(Eigen::VectorXd upper, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, int N, int h, int s);


struct TVBS_Matrix_two TVBSo_grad_cdf_tvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr);
struct TVBS_Matrix_two TVBSo_grad_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);
struct TVBS_Matrix_three TVBSo_grad_non_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);
struct TVBS_Matrix_two TVBSo_grad_cdf_qvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr);
struct TVBS_Matrix_two TVBSo_grad_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);
struct TVBS_Matrix_three TVBSo_grad_non_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);
struct TVBS_Vector_three TVBSo_pdf_mvna_tvbs_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
struct TVBS_Vector_four TVBSo_pdf_mvn_analytic_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x, int log_out);
Eigen::VectorXd TVBS(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
double TVBS_p(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);




