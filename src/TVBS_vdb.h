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

//struct Hessian_trunc_struct TVBS_Hessian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma);
// Eigen::MatrixXd TVBS_Hessian_trunc_gen_R(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma, int out);




struct TVBS_Matrix_two TVBS_grad_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);
struct TVBS_Matrix_three TVBS_grad_non_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);
struct TVBS_Matrix_two TVBS_grad_cdf_qvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr);
struct TVBS_Matrix_two TVBS_grad_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor);

struct TVBS_Matrix_three TVBS_grad_non_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x);

Eigen::VectorXd TVBS_grad_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);


Eigen::VectorXd cal_grad_rhokl(Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::VectorXd He_onetofour);

Eigen::MatrixXd Hessian_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);

struct TVBS_Vector_three TVBS_pdf_mvna_tvbs_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
struct TVBS_Vector_four TVBS_pdf_mvn_analytic_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x, int log_out);

Eigen::VectorXd TVBS_v2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
double TVBS_pv2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);

Eigen::VectorXd TVBS_grad(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out);
Eigen::MatrixXd TVBS_hess_new(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat);






