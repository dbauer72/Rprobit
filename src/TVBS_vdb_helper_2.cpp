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

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies


#include "toms462.h"  // allows bivariate normal cdf calculations
#include "distrib.h"
#include "lin_alg.h"
#include "TVBS_vdb_Hessian.h"
#include "TVBS_vdb_helper.h"

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <random>       // for random numbers


extern double tol;

// // // // // // // // // // // // // // // // //
//                                              //
//    DECLARATION OF GLOBAL VARIABLES           //
//                                              //
// // // // // // // // // // // // // // // // //

extern int x1symmetric_global;
extern int x2symmetric_global;
extern int x1diagonal_global;
extern int x2diagonal_global;    // = 0 always
extern int x2correlation_global;
extern int xinvsymmetric_global;
extern int xinvcorrelation_global;
extern int xinvdiagonal_global;    // = 0 always
extern int omsymmetric_global;
extern int omdiagonal_global;    // = 0 always
extern int condcov_global;
extern int cholesky_global;    // = 0 always
extern int condcovmeantrunc_global;
extern int condcovsigtrunc_global;
extern int optimal_global;
extern int covarr_global;    // = 1 always

extern int counter_check;

extern double pi_global;
extern double tol;
// 

// // // // // // // // // // // // // // // // //
//                                              //
//    STRUCTURES                                //
//                                              //
// // // // // // // // // // // // // // // // //

struct TVBS_Matrix_two
{
  Eigen::MatrixXd mat1, mat2;
};

struct TVBS_Matrix_three
{
  Eigen::MatrixXd mat1, mat2, mat3;
};

struct TVBS_Matrix_four
{
  Eigen::MatrixXd mat1, mat2, mat3, mat4;
};

struct TVBS_double_three
{
  double d1, d2, d3;
};

struct TVBS_double_Vector_two
{
  // not used
  double d1;
  Eigen::VectorXd v1, v2;
};

struct TVBS_Vector_three
{
  Eigen::VectorXd v1, v2, v3;
};

struct TVBS_Vector_four
{
  Eigen::VectorXd v1, v2, v3, v4;
};

struct TVBS_Vec_Mat
{
  Eigen::VectorXd vec1;
  Eigen::MatrixXd mat1;
};

struct TVBS_Vec_Mat_Vec
{
  Eigen::VectorXd vec1;
  Eigen::MatrixXd mat1;
  Eigen::VectorXd vec2;
};


struct Hessian_trunc_struct
{
  Eigen::MatrixXd mat1, mat2, mat3, mat4, mat5;
};



/////////////////////////////////////////////////
/// functions //
////////////////////////////////////////////////


// // // // // // // // // // // // // // // // //
//                                              //
//    TRUNCATIONS                               //
//                                              //
// // // // // // // // // // // // // // // // //

struct TVBS_Matrix_two TVBS_truncate_bi_normal_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point, double tol = 0.000001)
{
  // extract relevant variables
  double mu_2km1 = mu_untruncated(0);
  double mu_2k = mu_untruncated(1);
  double sigma_1s = Covariance_untruncated(0,0);
  double sigma_2s = Covariance_untruncated(1,1);
  double sigma_12 = Covariance_untruncated(0,1);
  
  // calculate standardized w and correlation
  double rho = sigma_12/sqrt(sigma_1s * sigma_2s);
  Eigen::MatrixXd Rho(2,2);
  Rho.setOnes();
  Rho(0,1) = rho;
  Rho(1,0) = rho;
  double w_2km1 = (truncation_point(0)-mu_2km1)/sqrt(sigma_1s);
  double w_2k = (truncation_point(1)-mu_2k)/sqrt(sigma_2s);
  Eigen::VectorXd w_temp(2);
  w_temp(0) = w_2km1;
  w_temp(1) = w_2k;
  
  
  // calculate mu_tilde_k
  double delta_2km1 = std_normal_pdf(w_2km1)*std_normal_cdf( (w_2k - rho * w_2km1)/sqrt(1.0-pow(rho,2)) );
  delta_2km1 = std::max(delta_2km1, pow(tol,2));
  
  
  double delta_2k = std_normal_pdf( w_2k ) * std_normal_cdf( (w_2km1 - rho * w_2k)/sqrt(1.0-pow(rho,2)) );
  delta_2k = std::max(delta_2k, pow(tol,2));
  
  
  Eigen::VectorXd mu_zero_2(2);
  mu_zero_2.setZero();
  double Phi_2_w = TVBS_pmvnorm_cpp(w_temp, mu_zero_2, Rho);
  Phi_2_w = std::max(Phi_2_w, pow(tol,2));
  
  double lambda_2km1 = -1.0*(delta_2km1 + rho*delta_2k)/Phi_2_w;
  double lambda_2k = -1.0*(delta_2k + rho*delta_2km1)/Phi_2_w;
  
  Eigen::VectorXd mu_tilde_k(2);
  mu_tilde_k(0) = mu_2km1+sqrt(sigma_1s)*lambda_2km1;
  mu_tilde_k(1) = mu_2k+sqrt(sigma_2s)*lambda_2k;
  
  
  // calculate Omega_k
  double phi_2_w = TVBS_biv_std_norm_pdf(w_temp, rho);
  phi_2_w = std::max(phi_2_w, pow(tol,2));
  
  double sigma_tilde_1s   = 1.0 - ( w_2km1 * delta_2km1 + pow(rho,2) * w_2k   * delta_2k   - (1.0 - pow(rho,2))*rho * phi_2_w)/Phi_2_w - pow(lambda_2km1,2);
  double sigma_tilde_2s   = 1.0 - ( w_2k   * delta_2k   + pow(rho,2) * w_2km1 * delta_2km1 - (1.0 - pow(rho,2))*rho * phi_2_w)/Phi_2_w - pow(lambda_2k,2);
  double sigma_tilde_12   = rho - ( rho * w_2km1 * delta_2km1 + rho * w_2k * delta_2k - (1.0-pow(rho,2)) * phi_2_w )/Phi_2_w - lambda_2km1 * lambda_2k;
  
  
  
  // make sure that the result is a Covariance matrix!
  if(sigma_tilde_1s<tol){
    sigma_tilde_1s = tol;
  }
  if(sigma_tilde_2s<tol){
    sigma_tilde_2s = tol;
  }
  double cor_sigma = sigma_tilde_12/(sqrt(sigma_tilde_1s*sigma_tilde_2s));
  if(abs(cor_sigma)>1){
    sigma_tilde_12 = sigma_tilde_12/abs(sigma_tilde_12)*(1-tol)*(sqrt(sigma_tilde_1s*sigma_tilde_2s));
    cor_sigma = sigma_tilde_12/(sqrt(sigma_tilde_1s*sigma_tilde_2s));
  }
  
  
  
  Eigen::MatrixXd Cov_Z(2,2);
  Cov_Z(0,0) = sigma_tilde_1s;
  Cov_Z(0,1) = sigma_tilde_12;
  Cov_Z(1,0) = sigma_tilde_12;
  Cov_Z(1,1) = sigma_tilde_2s;
  
  Eigen::MatrixXd Gamma_Sigma(2,2);
  Gamma_Sigma.setZero();
  Gamma_Sigma(0,0) = sqrt(sigma_1s);
  Gamma_Sigma(1,1) = sqrt(sigma_2s);
  
  Eigen::MatrixXd Omega_k = Gamma_Sigma*Cov_Z*Gamma_Sigma;
  
  struct TVBS_Matrix_two output = {mu_tilde_k, Omega_k};
  return output;
  
}



// gradient of truncation function
struct TVBS_Matrix_two TVBS_truncate_bi_normal_gradient_cpp(Eigen::VectorXd mu_untruncated, Eigen::MatrixXd Covariance_untruncated, Eigen::VectorXd truncation_point)
{
  Eigen::VectorXd mu_zero_2(2), mu_zero_3(3), mu_zero_4(4);
  mu_zero_2.setZero();
  mu_zero_3.setZero();
  mu_zero_4.setZero();
  
  int m = truncation_point.size();
  Eigen::VectorXd sigma_1_2_s = Covariance_untruncated.diagonal();
  Eigen::VectorXd sigma_1_2 = sigma_1_2_s.array().sqrt();
  Eigen::VectorXd sigma_1_2_inv = 1.0/sigma_1_2.array();
  Eigen::VectorXd sigma_1_2_inv_rev = sigma_1_2_inv.reverse();
  
  // Calculate normalized truncation point
  Eigen::VectorXd truncation_point_normalized = (truncation_point.array() - mu_untruncated.array())/sigma_1_2.array();
  
  Eigen::VectorXd truncation_point_normalized_rev = truncation_point_normalized.reverse();
  
  // calculate correlation matrix from covariance matrix
  Eigen::MatrixXd Sigma_1_2_inv(sigma_1_2_inv.size(), sigma_1_2_inv.size());
  Sigma_1_2_inv.setZero();
  Sigma_1_2_inv.diagonal() = sigma_1_2_inv;
  Eigen::MatrixXd Rho = Sigma_1_2_inv*Covariance_untruncated*Sigma_1_2_inv;
  double rho = Rho(0,1);
  
  // calculate the "probability of the truncation points"
  double p = TVBS_pmvnorm_cpp(truncation_point_normalized, mu_zero_2, Rho);
  
  // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  // 
  // THIS IS TO MAKE SURE WE DO NOT DIVIDE BY ZERO
  //
  p = std::max(p,pow(tol,2));
  //
  //
  // // // // // // // // // // // // // // // // // // // // // // // // // // // 
  
  
  
  // auxiliary variables
  double rhotilde2 = 1.0-pow(rho,2);
  double rhotilde = sqrt(rhotilde2);
  
  
  Eigen::VectorXd tr1 = (truncation_point_normalized.array() - rho*truncation_point_normalized_rev.array())/rhotilde;
  Eigen::VectorXd tr2 = tr1.reverse();
  
  Eigen::VectorXd trcomp = rho*truncation_point_normalized.array()-truncation_point_normalized_rev.array();
  Eigen::VectorXd trcomp_rev = trcomp.reverse();
  
  Eigen::VectorXd pd1 = truncation_point_normalized.unaryExpr(&std_normal_pdf);
  Eigen::VectorXd pd2 = pd1.reverse();
  
  Eigen::VectorXd cd1 = tr1.unaryExpr(&std_normal_cdf);
  Eigen::VectorXd cd2 = cd1.reverse();
  
  Eigen::VectorXd del1 = pd1.array()*cd2.array();
  Eigen::VectorXd del2 = del1.reverse();
  
  double pdf2 = 1.0/rhotilde*pd1(0)*std_normal_pdf(tr1(1));
  
  Eigen::VectorXd mu = (-1.0*rho*del2.array()-del1.array())/p;
  
  
  // calculate derivatives of the mean mu using chain rule
  Eigen::VectorXd dmutilde1_dw1 = (1.0/p)*sigma_1_2.array()*(truncation_point_normalized.array()*del1.array() + del1.array()*(-1.0*mu.array()));
  Eigen::VectorXd dmutilde1_dw2 = (1.0/p) * (-1.0*sigma_1_2.array()) * (rhotilde2*pdf2 + (mu.array() - rho*truncation_point_normalized_rev.array())*del2.array() );
  Eigen::VectorXd dmutilde_drho = (1.0/p) * (-sigma_1_2.array())*( del2.array() + pdf2*(mu.array()-truncation_point_normalized.array()) );
  Eigen::VectorXd dmutilde1_da1 = 1.0- sigma_1_2_inv.array()*dmutilde1_dw1.array();
  Eigen::VectorXd dmutilde1_da2 = -1.0*sigma_1_2_inv_rev.array() * dmutilde1_dw2.array();
  Eigen::VectorXd dmutilde1_dsig1 = -1.0*(sigma_1_2_inv.array()*truncation_point_normalized.array()*dmutilde1_dw1.array() - mu.array() + rho*sigma_1_2_inv.array()*dmutilde_drho.array());
  Eigen::VectorXd dmutilde1_dsig2 = -1.0*(sigma_1_2_inv_rev.array()*truncation_point_normalized_rev.array()*dmutilde1_dw2.array() + rho*sigma_1_2_inv_rev.array()*dmutilde_drho.array());
  Eigen::VectorXd dmutilde1_dsig1sq = 0.5*sigma_1_2_inv.array()*dmutilde1_dsig1.array();
  Eigen::VectorXd dmutilde1_dsig2sq = 0.5*sigma_1_2_inv_rev.array()*dmutilde1_dsig2.array();
  Eigen::VectorXd dmutilde_dsig12 = sigma_1_2_inv.prod()*dmutilde_drho.array();
  Eigen::VectorXd dmutilde1_dtr1 = sigma_1_2_inv.array()*dmutilde1_dw1.array();
  Eigen::VectorXd dmutilde1_dtr2 = sigma_1_2_inv_rev.array()*dmutilde1_dw2.array();
  
  Eigen::VectorXd dmu_deriv1(7), dmu_deriv2(7);
  dmu_deriv1.setZero();
  dmu_deriv2.setZero();
  dmu_deriv1 << dmutilde1_da1(0), dmutilde1_da2(0), dmutilde1_dsig1sq(0), dmutilde_dsig12(0), dmutilde1_dsig2sq(0), dmutilde1_dtr1(0), dmutilde1_dtr2(0);
  dmu_deriv2 << dmutilde1_da2(1), dmutilde1_da1(1), dmutilde1_dsig2sq(1), dmutilde_dsig12(1), dmutilde1_dsig1sq(1), dmutilde1_dtr2(1), dmutilde1_dtr1(1);
  Eigen::MatrixXd dmu_deriv(7, 2);
  dmu_deriv.setZero();
  dmu_deriv.col(0) = dmu_deriv1;
  dmu_deriv.col(1) = dmu_deriv2;
  
  
  // calculate derivatives of the covariance matrix entries using chain rule
  Eigen::VectorXd sig = ((p - (truncation_point_normalized.array()*pd1.array()*cd2.array()) - (pow(rho,2)*truncation_point_normalized_rev.array()*pd2.array()*cd1.array()) + (rhotilde*rho*pd2.array()*(tr1.unaryExpr(&std_normal_pdf)).array() ) ) / p) - (mu.array()*mu.array());
  double sig_12 = ((rho*p - rho*truncation_point_normalized(0)*pd1(0)*cd2(0) + rhotilde*pd1(0)*std_normal_pdf(tr2(0)) - rho*truncation_point_normalized_rev(0)*pd2(0)*cd1(0))/p) - mu.prod();
  
  Eigen::VectorXd dsigtilde1_dw1 = -1.0*(del1.array()/p) * ((-1.0*truncation_point_normalized.array()*truncation_point_normalized.array()) + sig.array() + mu.array()*mu.array()) - 2* mu.array() * (dmutilde1_dw1.array() / sigma_1_2.array());
  Eigen::VectorXd dsigtilde1_dw2 = -1.0*(1.0/p) * ((rhotilde2 * ( rho*pdf2*truncation_point_normalized_rev.array() + truncation_point_normalized.array()*pdf2 - del2.array())) - pow(rho,2)*(truncation_point_normalized_rev.array()*truncation_point_normalized_rev.array())*del2.array() + (sig.array()+mu.array()*mu.array())*del2.array()) - 2*mu.array()*(dmutilde1_dw2.array()/sigma_1_2.array());
  Eigen::VectorXd dsigtilde_drho = -(1.0/p) * ( pdf2*(-1.0*truncation_point_normalized.array()*truncation_point_normalized.array()) + 2*rho*truncation_point_normalized_rev.array()*del2.array() - (2*pdf2*rhotilde2) + ((sig.array()+mu.array()*mu.array())*pdf2)) - 2*mu.array()*(dmutilde_drho.array()/sigma_1_2.array());
  Eigen::VectorXd dsigtilde = (1.0/p)*( (rho*truncation_point_normalized.array())*((rho*pdf2) + truncation_point_normalized.array()*del1.array()) - pdf2*truncation_point_normalized.array() - del1.array()*(sig_12+mu.prod()));
  
  Eigen::VectorXd foo_vector_1(2), foo_vector_2(2), foo_vector_3(2);
  foo_vector_1 << dmutilde1_dw2(1), dmutilde1_dw1(0);
  foo_vector_2 << dmutilde1_dw1(1), dmutilde1_dw2(0);
  
  double dsigtilde12_dw1 = dsigtilde(0) - ( mu.array()*( foo_vector_1.array() / (sigma_1_2.reverse()).array() ) ).sum();
  double dsigtilde12_dw2 = dsigtilde(1) - ( mu.array()*( foo_vector_2.array() / (sigma_1_2.reverse()).array() ) ).sum();
  double dsigtilde12_drho = 1.0 - (1.0/p) * ( -1.0*(truncation_point_normalized.prod())*pdf2 + (truncation_point_normalized.array()*del1.array()).sum() + ((mu.prod())+sig_12)*pdf2);
  dsigtilde12_drho = dsigtilde12_drho - (mu.array()*(dmutilde_drho.reverse()).array()*(sigma_1_2_inv_rev.array())).sum() ;
  
  foo_vector_3 << dsigtilde12_dw1, dsigtilde12_dw2;
  Eigen::VectorXd domg1_dtr1 = (sigma_1_2.array())*dsigtilde1_dw1.array();
  Eigen::VectorXd domg12_dtr = (sigma_1_2.reverse()).array()*foo_vector_3.array();
  Eigen::VectorXd domg1_dtr2 = (sigma_1_2.array()*sigma_1_2.array())*sigma_1_2_inv_rev.array()*dsigtilde1_dw2.array();
  Eigen::VectorXd domg1_da1 = -1.0*domg1_dtr1;
  Eigen::VectorXd domg12_da = -1.0*domg12_dtr;
  Eigen::VectorXd domg1_da2 = -1.0*domg1_dtr2;
  Eigen::VectorXd domg1_s1 = sig.array()+dsigtilde1_dw1.array()*truncation_point_normalized.array()*(-0.5)+dsigtilde_drho.array()*rho*(-0.5);
  Eigen::VectorXd domg1_s12 = (sigma_1_2.array()/(sigma_1_2.reverse()).array())*dsigtilde_drho.array();
  Eigen::VectorXd domg1_s2 = (-0.5)*(sigma_1_2_s.array()/(sigma_1_2_s.reverse()).array())*(dsigtilde1_dw2.array()*truncation_point_normalized_rev.array()+dsigtilde_drho.array()*rho);
  Eigen::VectorXd domg12_s1 = 0.5*((sigma_1_2.reverse()).array()/sigma_1_2.array())*(sig_12-truncation_point_normalized.array()*foo_vector_3.array()-rho*dsigtilde12_drho);
  
  double domg12_s12 = dsigtilde12_drho;
  
  Eigen::VectorXd domg_deriv_11(7), domg_deriv_12(7), domg_deriv_22(7);
  domg_deriv_11 << domg1_da1(0), domg1_da2(0), domg1_s1(0), domg1_s12(0), domg1_s2(0), domg1_dtr1(0), domg1_dtr2(0);
  domg_deriv_12 << domg12_da(0), domg12_da(1), domg12_s1(0), domg12_s12, domg12_s1(1), domg12_dtr(0), domg12_dtr(1);
  domg_deriv_22 << domg1_da2(1), domg1_da1(1), domg1_s2(1), domg1_s12(1), domg1_s1(1), domg1_dtr2(1), domg1_dtr1(1);
  
  Eigen::MatrixXd domg_deriv(7,3);
  domg_deriv.col(0) = domg_deriv_11;
  domg_deriv.col(1) = domg_deriv_12;
  domg_deriv.col(2) = domg_deriv_22;
  
  struct TVBS_Matrix_two output_final = {dmu_deriv, domg_deriv};
  return output_final;
}



// function to truncate mu for quadrovariate normal distribution
struct TVBS_Matrix_two TVBS_mu_l_trunc_bivariate_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd tr_point)
{
  int m = mu.size();
  Eigen::VectorXd mu_untrunc = mu.head(2);
  Eigen::MatrixXd v11 = cov.block(0,0,2,2);
  Eigen::MatrixXd v11inv = v11.inverse();
  Eigen::MatrixXd v12 = cov.block(0,2,2,m-2);
  Eigen::MatrixXd v21 = v12.transpose();
  Eigen::MatrixXd v22 = cov.block(2,2,m-2,m-2);
  
  struct TVBS_Matrix_two output = TVBS_truncate_bi_normal_cpp(mu_untrunc, v11, tr_point);
  Eigen::VectorXd mu1 = output.mat1;
  Eigen::MatrixXd omega11 = output.mat2;
  
  Eigen::VectorXd mu2(4);
  mu2 << mu1, mu.tail(m-2) + (v21*v11inv)*(mu1-mu_untrunc);
  Eigen::MatrixXd omega12 = omega11*v11inv*v12;
  Eigen::MatrixXd omega22 = v22-v21*(v11inv-v11inv*omega11*v11inv)*v12;
  
  Eigen::MatrixXd omg(omega11.rows()+omega22.rows(), omega11.cols()+omega12.cols());
  omg.block(0,0,omega11.rows(),omega11.cols()) = omega11;
  omg.block(0,omega11.cols(),omega12.rows(),omega12.cols()) = omega12;
  omg.block(omega11.rows(),0,omega12.cols(),omega12.rows()) = omega12.transpose();
  omg.block(omega11.rows(),omega11.cols(),omega22.rows(),omega22.cols()) = omega22;
  
  struct TVBS_Matrix_two output_final = {mu2, omg};
  return output_final;
}

Eigen::MatrixXd TVBS_Jacobian_biv_norm_trunc(Eigen::VectorXd w, double rho)
{
  Eigen::MatrixXd Jacob(5,3); // Jacobian of lambda_{0,1}, Omega_{0,0},Omega_{0,1},Omega_{1,1} as a function of w_0, w_1, rho. (indexation according to CPP to match code) 
  Jacob.setZero();
  
  // calculate deltas 
  double cdf_0 = std_normal_cdf((w(1) - rho*w(0))/sqrt(1-rho*rho));
  double delta_0 = std_normal_pdf(w(0))*cdf_0;
  double cdf_1 = std_normal_cdf((w(0) - rho*w(1))/sqrt(1-rho*rho));
  double delta_1 = std_normal_pdf(w(1))*cdf_1;
  
  // now for derivatives 
  // gradients of delta_0, 
  Eigen::MatrixXd grad_delta_0(1,3);
  grad_delta_0.setZero();
  
  grad_delta_0(0,0) = -std_normal_pdf(w(0))*cdf_0*w(0) - std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_0(0,1) = std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_0(0,2) = std_normal_pdf(w(0))*normal_pdf((w(1) - rho*w(0)),(1-rho*rho))*(-w(0)+(w(1)-rho*w(0))*rho/(1-rho*rho));
  
  // gradients of delta_1. 
  Eigen::MatrixXd grad_delta_1(1,3);
  grad_delta_1.setZero();
  grad_delta_1(0,1) = -std_normal_pdf(w(1))*cdf_1*w(1) - std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_1(0,0) = std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_1(0,2) = std_normal_pdf(w(1))*normal_pdf((w(0) - rho*w(1)),(1-rho*rho))*(-w(1)+(w(0)-rho*w(1))*rho/(1-rho*rho));  
  // now for the lambdas. 
  // calculate lambdas
  double Phi2 = biv_normal_cdf(w(0),w(1),rho);
  double lambda_0 = -(delta_0+rho*delta_1)/Phi2;
  double lambda_1 = -(delta_1+rho*delta_0)/Phi2;
  
  Eigen::VectorXd dPhi2 = grad_cdf(w(0),w(1),rho);
  
  // derivatives of lambda_0 and lambda_1. 
  Jacob(0,0) = -(grad_delta_0(0,0)+rho*grad_delta_1(0,0))/Phi2 + (delta_0+rho*delta_1)/pow(Phi2,2)*dPhi2(0);
  Jacob(0,1) = -(grad_delta_0(0,1)+rho*grad_delta_1(0,1))/Phi2 + (delta_0+rho*delta_1)/pow(Phi2,2)*dPhi2(1);
  Jacob(0,2) = -(grad_delta_0(0,2)+rho*grad_delta_1(0,2))/Phi2 + (delta_0+rho*delta_1)/pow(Phi2,2)*dPhi2(2);
  Jacob(0,2) += -delta_1/Phi2; 
  
  Jacob(1,0) = -(grad_delta_1(0,0)+rho*grad_delta_0(0,0))/Phi2 + (delta_1+rho*delta_0)/pow(Phi2,2)*dPhi2(0);
  Jacob(1,1) = -(grad_delta_1(0,1)+rho*grad_delta_0(0,1))/Phi2 + (delta_1+rho*delta_0)/pow(Phi2,2)*dPhi2(1);
  Jacob(1,2) = -(grad_delta_1(0,2)+rho*grad_delta_0(0,2))/Phi2 + (delta_1+rho*delta_0)/pow(Phi2,2)*dPhi2(2);
  Jacob(1,2) += -delta_0/Phi2; 
  
  
  // turn to the Omegas. 
  double phi2 = biv_normal_pdf(w(0),w(1),rho);
  
  Eigen::VectorXd dphi2(3);
  dphi2.setZero();
  
  double detSigma = 1 - rho*rho;
  double bexp = w(0)*w(0) + w(1)*w(1) - 2*w(0)*w(1)*rho;
  
  Eigen::VectorXd ddetSigma(3);
  ddetSigma.setZero();
  ddetSigma(2)= -2*rho;
  
  dphi2(0) = (-phi2)*(w(0)-rho*w(1))/(detSigma);
  dphi2(1) = (-phi2)*(w(1)-rho*w(0))/(detSigma);
  dphi2(2) = phi2*(rho/(detSigma)+ (w(0)*w(1)*detSigma-(bexp*rho))/pow(detSigma,2) );
  
  
  
  // Omega_{00}
  double exp_0 = w(0)*delta_0 + rho*rho*w(1)*delta_1 - rho*(1-rho*rho)*phi2;
  Eigen::VectorXd grad_exp_0(3);
  grad_exp_0.setZero();
  for (int k=0;k<3;k++){
    grad_exp_0(k) = w(0)*grad_delta_0(k) + rho*rho*w(1)*grad_delta_1(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_0(0) += delta_0; 
  grad_exp_0(1) += rho*rho*delta_1; 
  grad_exp_0(2) += 2*rho*w(1)*delta_1 - phi2*(1-3*rho*rho);
  
  
  //double omega_00 = 1 - exp_0/Phi2 - lambda_0*lambda_0; 
  for (int k=0;k<3;k++){
    Jacob(2,k) =  -grad_exp_0(k)/Phi2 + exp_0*dPhi2(k)/pow(Phi2,2) - 2*lambda_0 * Jacob(0,k);
  }
  
  // Omega_{11}
  double exp_1 = w(1)*delta_1 + rho*rho*w(0)*delta_0 - rho*(1-rho*rho)*phi2;
  Eigen::VectorXd grad_exp_1(3);
  grad_exp_1.setZero();
  for (int k=0;k<3;k++){
    grad_exp_1(k) = w(1)*grad_delta_1(k) + rho*rho*w(0)*grad_delta_0(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_1(0) += rho*rho*delta_0; 
  grad_exp_1(1) += delta_1; 
  grad_exp_1(2) += 2*rho*w(0)*delta_0 - phi2*(1-3*rho*rho);
  
  //double omega_11 = 1 - exp_1/Phi2 - lambda_1*lambda_1; 
  for (int k=0;k<3;k++){
    Jacob(4,k) =  -grad_exp_1(k)/Phi2 + exp_1*dPhi2(k)/pow(Phi2,2) - 2*lambda_1 * Jacob(1,k);  
  }
  
  // Omega_{01}
  double exp_01 = rho*w(1)*delta_1 + rho*w(0)*delta_0  - (1-rho*rho)*phi2;
  Eigen::VectorXd grad_exp_01(3);
  grad_exp_01.setZero();
  for (int k=0;k<3;k++){
    grad_exp_01(k) = rho*w(1)*grad_delta_1(k) + rho*w(0)*grad_delta_0(k) - (1-rho*rho)*dphi2(k);
  }
  grad_exp_01(0) += rho*delta_0; 
  grad_exp_01(1) += rho*delta_1; 
  grad_exp_01(2) += w(0)*delta_0 + w(1)*delta_1 + (2*rho* phi2);
  
  //double omega_01 = rho - exp_01/Phi2 - lambda_0*lambda_1; 
  for (int k=0;k<3;k++){
    Jacob(3,k) =  -grad_exp_01(k)/Phi2 + exp_01*dPhi2(k)/pow(Phi2,2) - lambda_1 * Jacob(0,k)- lambda_0 * Jacob(1,k);
  }
  Jacob(3,2) += 1;
  
  // return matrix. 
  return Jacob;
}



// // // // // // // // // // // // // // // // //
//                                              //
//    NEW AUXILLARY FUNCTIONS                   //
//                                              //
// // // // // // // // // // // // // // // // //

// own functions

Eigen::VectorXi TVBS_logic2position_cpp(Eigen::VectorXd vec_logical) 
{
  int n = vec_logical.rows();
  int n1 = vec_logical.sum();
  
  Eigen::VectorXi vec_position(n1);
  
  int current_position = 0;
  for(int j = 0; j < n; j++ )
  {
    if(vec_logical(j)>0)
    {
      vec_position(current_position) = j;
      current_position++;
    }
  }
  return vec_position;
}


Eigen::MatrixXd TVBS_vec_2_upper_diag_cpp(Eigen::VectorXd vec_in, int diag_inc = 1)
{
  int dim = (-1+sqrt(1+8*vec_in.size()))/2+(1-diag_inc);
  Eigen::MatrixXd mat_out(dim,dim);
  mat_out.setZero();
  
  int i_entry = 0;
  for(int i_row = 0; i_row<dim; i_row++)
  {
    for(int i_col = i_row+(1-diag_inc); i_col<dim; i_col++)
    {
      mat_out(i_row, i_col) = vec_in(i_entry);
      i_entry++;
    }
  }
  return mat_out;
}



Eigen::VectorXd TVBS_lower_tri_entries_cpp(Eigen::MatrixXd x, int diag_inc = 1)
{
  int dimdiff = x.rows()-x.cols();
  int vec_size = (x.rows()-1+diag_inc)*(x.rows()+diag_inc)/2 + std::max(0, dimdiff)*x.cols();
  
  Eigen::VectorXd output(vec_size);
  int current_position = 0;
  for(int j_col = 0; j_col < x.cols(); j_col++)
  {
    for(int j_row = (j_col+1-diag_inc); j_row < x.rows(); j_row++)
    {
      output(current_position) = x(j_row,j_col);
      current_position++;
    }
  }
  return output;
}



//Eigen::MatrixXd TVBS_Hessian_trunc_gen_R(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma, int out)
//{
//  struct Hessian_trunc_struct  output = TVBS_Hessian_biv_gen_trunc(x, alpha, Sigma);
//  
//  Eigen::MatrixXd outp = output.mat1;
//  if (out==2){ outp = output.mat2;}
//  if (out==3){ outp = output.mat3;}
//  if (out==4){ outp = output.mat4;}
//  if (out==5){ outp = output.mat5;}
//  
//  return outp;
//}


// "translated functions"
Eigen::VectorXd TVBS_sincs_grad_cpp(double x) 
{
  double pi     = 3.14159265358979323846;
  double ee, d_ee, cs, d_cs, sign_x, sx, d_sx;
  
  ee      = pow(pi/2 - abs(x),2);
  if(x<0)
  {
    d_ee    = 2*( pi/2 + x );
  }
  else
  {
    d_ee    = 2*( x - pi/2  );
  }
  
  if(ee<0.00005)
  {
    cs      = ee*(1.0-ee*(1.0-2*ee/15)/3);
    d_cs    = (1.0 - ee*(1.0-2*ee/15 )/3) + ee*( -(( 1.0 - 2*ee/15 )/3 - ee * (2/45) ) );
    d_cs    = d_ee*d_cs;
    sign_x  = x/abs(x);
    sx      = (1.0-ee*(1.0-ee/12)/2)*sign_x;
    
    if(x>0)
    {
      d_sx    = -((1.0 - ee/12)/2 + ee*(-1.0/12)/2);
    }
    else
    {
      d_sx    = ( 1.0 - ee/12 )/2 + ee*( -1.0/12 )/2;
    }
    d_sx    = d_sx*d_ee;
  }
  else
  {
    sx      = sin(x);
    d_sx    = cos(x);
    cs      = 1.0 - pow(sx,2);
    d_cs    = -2*sx*d_sx;
  }
  
  Eigen::VectorXd output_vec(4);
  output_vec << d_sx, d_cs, sx, cs;
  return output_vec;
}



Eigen::VectorXd TVBS_pntgnd_grad_cpp(double ba, double bb, double bc, double ra, double rb, double r, double rr) 
{
  double dt, d_dt_ra, d_dt_rb, d_dt_r, d_dt_rr, bt, ft, d_bt_dt, d_bt_ba, d_ft_ba, d_bt_bb, d_ft_bb, d_bt_bc, d_bt_ra, d_bt_rb, d_bt_r, d_ft_r, d_bt_rr, d_ft_rr;
  
  double f        = 0;
  double d_ba     = 0;
  double d_bb     = 0;
  double d_bc     = 0;
  double d_ra     = 0;
  double d_rb     = 0;
  double d_r      = 0;
  double d_rr     = 0;
  
  dt              = rr*(rr - pow(ra - rb,2) - 2*ra*rb*(1.0 - r));
  
  d_dt_ra         = rr*( -2*(ra - rb) - 2*rb*(1.0 - r));
  d_dt_rb         = rr*(2*(ra - rb) - 2*ra*(1.0 - r));
  d_dt_r          = rr*2*ra*rb;
  d_dt_rr         = 2*rr - pow(ra - rb, 2) - 2*ra*rb*(1.0 - r);
  
  
  if(dt>0)
  {
    bt              = (bc*rr + ba*(r*rb - ra) + bb*(r*ra - rb))/sqrt(dt);
    ft              = pow(ba - r*bb,2)/rr + bb*bb;
    
    d_bt_dt         = -0.5*bt/dt;
    
    d_bt_ba         = (r*rb - ra)/sqrt(dt);
    d_ft_ba         = 2*(ba - r*bb)/rr;
    d_bt_bb         = (r*ra - rb)/sqrt(dt);
    d_ft_bb         = 2*(ba - r*bb)*(-r)/rr + 2*bb;
    d_bt_bc         = rr/sqrt(dt);
    
    d_bt_ra         = (ba*(-1.0) + bb*(r) )/sqrt(dt)  +  d_bt_dt*d_dt_ra ;
    d_bt_rb         = (ba*(r)  + bb*(-1.0))/sqrt(dt)  +  d_bt_dt*d_dt_rb;
    d_bt_r          = (ba*(rb) + bb*(ra))/sqrt(dt)  +  d_bt_dt*d_dt_r;
    d_ft_r          = 2*( ba - r*bb )*(-bb)/rr;
    d_bt_rr         = bc/sqrt(dt)                   +  d_bt_dt*d_dt_rr;
    d_ft_rr         = -pow(ba - r*bb, 2)/pow(rr, 2);
    
    
    if( (bt>-10) && (ft<100) )
    {
      f               = std::exp( -ft/2 );
      d_ba            = -0.5*f*d_ft_ba;
      d_bb            = -0.5*f*d_ft_bb;
      d_r             = -0.5*f*d_ft_r;
      d_rr            = -0.5*f*d_ft_rr;
      
      
      if(bt<10)
      {
        double d_norm_bt, p_norm_bt;
        d_norm_bt       = std_normal_pdf(bt);
        p_norm_bt       = std_normal_cdf(bt);
        
        
        d_ba            = d_ba*p_norm_bt + f*d_norm_bt*d_bt_ba;
        d_bb            = d_bb*p_norm_bt + f*d_norm_bt*d_bt_bb;
        d_bc            =                  f*d_norm_bt*d_bt_bc;
        d_ra            =                  f*d_norm_bt*d_bt_ra;
        d_rb            =                  f*d_norm_bt*d_bt_rb;
        d_r             = d_r*p_norm_bt  + f*d_norm_bt*d_bt_r;
        d_rr            = d_rr*p_norm_bt + f*d_norm_bt*d_bt_rr;
        f               = f*p_norm_bt;
      }
    }
  }
  
  Eigen::VectorXd output_vec(8);
  output_vec.setZero();
  output_vec << d_ba, d_bb, d_bc, d_ra, d_rb, d_r, d_rr, f;
  return output_vec;
}



Eigen::MatrixXd TVBS_vec_symmetry_cpp(int r) 
{
  Eigen::MatrixXd temp1((r*(r+1)/2),(r*r));
  
  int current_row = 0;
  for(int j1 = 0; j1 < r; j1++)
  {
    for(int j2 = j1; j2 < r; j2++)
    {
      Eigen::MatrixXd temp(r,r);
      temp.setZero();
      temp(j1,j2) = 1.0;
      temp(j2,j1) = 1.0;
      temp.resize(1, r*r);
      
      temp1.row(current_row) = temp;
      current_row++;
    }
  }
  return temp1;
}


Eigen::MatrixXd TVBS_g_asym_to_sym_cpp(Eigen::MatrixXd x) 
{
  
  int k = sqrt(x.rows());
  int n = sqrt(x.cols());
  
  Eigen::MatrixXd temp_1 = TVBS_vec_symmetry_cpp(k);
  
  Eigen::MatrixXd temp_sel_matrix(n,n), temp_sel_matrix_ones(n,n);
  temp_sel_matrix.setZero();
  temp_sel_matrix_ones.setOnes();
  temp_sel_matrix.triangularView<Eigen::Lower>() = temp_sel_matrix_ones.triangularView<Eigen::Lower>();
  temp_sel_matrix.resize(n*n,1);
  
  
  Eigen::MatrixXd newx(k*k,(n*(n+1))/2);
  newx.setZero();
  
  Eigen::VectorXi col_select = TVBS_logic2position_cpp(temp_sel_matrix);
  
  // THERE MUST BE A BETTER WAY SELECTING THESE for example with row_select
  for(int j = 0; j<(n*(n+1))/2; j++)
  {
    newx.col(j) = x.col(col_select(j));
  }
  
  Eigen::MatrixXd gsym = temp_1*newx;
  return gsym;
}




Eigen::MatrixXd TVBS_g_inverse_cpp(Eigen::MatrixXd x) 
{
  int k = x.rows();
  int nrows, ncols;
  
  Eigen::MatrixXd g1(k*k,k*k), g1_Lhs(k,k), g1_Rhs(k,k);
  g1_Lhs = (x.transpose()).inverse();
  g1_Rhs = x.inverse();
  g1 = -1.0*(kroneckerProduct(g1_Lhs, g1_Rhs));
  
  nrows = g1.rows();
  ncols = g1.cols();
  
  if(xinvsymmetric_global==1)
  {
    int diag_inc = 1;
    Eigen::VectorXd temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k) ,diag_inc);
    g1 = TVBS_g_asym_to_sym_cpp(g1);
    nrows = g1.rows();
    ncols = g1.cols();
    
    if(xinvdiagonal_global==1)
    {
      Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix);
      nrows = row_select.size();
      ncols = row_select.size();
      for(int j = 0; j<nrows; j++)
      {
        g1.row(j) = g1.row(row_select(j));
      }
      for(int j = 0; j<ncols; j++)
      {
        g1.col(j) = g1.col(row_select(j));
      }
    }
    if( (xinvcorrelation_global==1) && (xinvdiagonal_global==0) )
    {
      Eigen::VectorXd temp_sel_matrix_2 = (1.0 - temp_sel_matrix.array()).matrix();
      Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix_2);
      nrows = row_select.size();
      for(int j = 0; j<nrows; j++)
      {
        g1.row(j) = g1.row(row_select(j));
      }
    }
    if((xinvcorrelation_global==1) && (xinvdiagonal_global==1))
    {
      g1.setZero();
      nrows = 1;
      ncols = 1;
    }
  }
  
  Eigen::MatrixXd g1_out = g1.block(0,0,nrows,ncols);
  return g1_out;
}

