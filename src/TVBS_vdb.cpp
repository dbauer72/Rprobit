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
#include "TVBS_vdb_helper.h"
#include "TVBS_vdb_helper_2.h"
#include "TVBS_vdb_Hessian.h"

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <random>       // for random numbers



// // // // // // // // // // // // // // // // //
//                                              //
//    DECLARATION OF GLOBAL VARIABLES           //
//                                              //
// // // // // // // // // // // // // // // // //

int x1symmetric_global        = 0;
int x2symmetric_global        = 0;
int x1diagonal_global         = 0;
int x2diagonal_global         = 0;    // = 0 always
int x2correlation_global      = 0;
int xinvsymmetric_global      = 0;
int xinvcorrelation_global    = 0;
int xinvdiagonal_global       = 0;    // = 0 always
int omsymmetric_global        = 0;
int omdiagonal_global         = 0;    // = 0 always
int condcov_global            = 0;
int cholesky_global           = 0;    // = 0 always
int condcovmeantrunc_global   = 0;
int condcovsigtrunc_global    = 0;
int optimal_global            = 3;
int covarr_global             = 1;    // = 1 always

//int counter_check = 0;

double pi_global              = 3.14159265358979323846;
extern double tol;


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


///////////////////////////////////////////////
///
//////////////////////////////////////////////


// // // // // // // // //
//  trivariate          //
// // // // // // // // //










// some differences to R due to trivariate CDF approximation
struct TVBS_Matrix_two TVBS_grad_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor)
{
  Eigen::VectorXd mu_3_0(3), mu_2_0(2);
  mu_3_0.setZero();
  mu_2_0.setZero();
  double tri_var_cdf = TVBS_pmvnorm_cpp(w, mu_3_0, cor);
  double bi_var_cdf = TVBS_pmvnorm_cpp(w.head(2), mu_2_0, cor.block(0,0,2,2));
  
  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(w.head(2), cor.block(0,0,2,2));
  double g_w_1 = output_3.d1;
  double g_w_2 = output_3.d2;
  double g_rho = output_3.d3;
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_tvn_cpp(w, cor);
  Eigen::VectorXd g_w_tri = output_2.mat1;
  Eigen::VectorXd g_rho_tri = output_2.mat2;
  
  Eigen::VectorXd g_w_12(2);
  g_w_12 << g_w_1, g_w_2;
  Eigen::VectorXd g_w1_w2_new = (bi_var_cdf*g_w_tri.head(2).array() - tri_var_cdf*g_w_12.array())/pow(bi_var_cdf,2);
  
  
  double g_rho12_new = (bi_var_cdf*g_rho_tri(0) - tri_var_cdf*g_rho)/pow(bi_var_cdf,2);
  
  Eigen::VectorXd output_final_1(3), output_final_2(3);
  output_final_1 << g_w1_w2_new, g_w_tri(2)/bi_var_cdf;
  output_final_2 << g_rho12_new, g_rho_tri.segment(1,2).array()/bi_var_cdf;
  
  struct TVBS_Matrix_two output_final = {output_final_1, output_final_2};
  
  return output_final;
}




// numerical issues might arise
struct TVBS_Matrix_three TVBS_grad_non_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  // cap the x values over 6 and under -6
  
  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  Eigen::MatrixXd Sigma_diag_sqrt_inv(3,3);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_tvn_by_cdf_bvn_cpp(x_norm.head(3), Cor_mat.block(0,0,3,3));
  Eigen::VectorXd g_w = output_2.mat1;
  Eigen::VectorXd g_rho = output_2.mat2;
  
  struct TVBS_Matrix_two output_3 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  Eigen::MatrixXd g_b_cor_cov = output_3.mat1;
  Eigen::MatrixXd g_omega_cor_cov = output_3.mat2;
  
  
  Eigen::VectorXd g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  Eigen::VectorXd g_cov = g_b_cor_cov*g_w;
  g_cov += g_omega_cor_cov*g_rho;
  Eigen::VectorXd g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  
  return output_final;
}





// issues arrise with approximated trivariate cdf
struct TVBS_Matrix_two TVBS_grad_cdf_qvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr)
{
  // initiate probability vector
  Eigen::VectorXd p(2);
  p.setZero();
  
  // include reordering here if wanted
  Eigen::VectorXd x_temp = x_norm;
  Eigen::MatrixXd Corr_temp = Corr;
  Eigen::VectorXd mu_temp(x_norm.size());
  mu_temp.setZero();
  
  // calculate first part of probability using trivariate normal distribution
  p(0) = TVBS_pmvnorm_cpp(x_temp.head(3), mu_temp.head(3), Corr_temp.block(0,0,3,3));
  
  // truncate on first two variables
  struct TVBS_Matrix_two output_2 = TVBS_mu_l_trunc_bivariate_cpp(mu_temp, Corr_temp, x_temp.head(2));
  Eigen::VectorXd mu_temp_1 = output_2.mat1;
  Eigen::MatrixXd sigma_temp_1 = output_2.mat2;
  
  // calculate second part by calculating P(4)=P(3,4)/P(3)
  double p_2_1 = TVBS_pmvnorm_cpp(x_temp.segment(2,2), mu_temp_1.segment(2,2), sigma_temp_1.block(2,2,2,2));
  double p_2_2 = std_normal_cdf((x_temp(2)-mu_temp_1(2))/sqrt(sigma_temp_1(2,2)));
  if (p_2_2<tol) { p_2_2 = tol; }
  p(1) = p_2_1/p_2_2;
  
  // update global condition numbers
  condcovsigtrunc_global = 0;
  condcovmeantrunc_global = 0;
  
  struct TVBS_Matrix_four output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(2,2), mu_temp, Corr_temp, x_temp.head(2));
  Eigen::MatrixXd g_y = output_4.mat1;
  Eigen::MatrixXd g_mu_mean = output_4.mat2;
  Eigen::MatrixXd g_x_mean = output_4.mat3;
  Eigen::MatrixXd g_c_mean = output_4.mat4;
  
  struct TVBS_Matrix_three output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Corr_temp, x_temp.head(2));
  Eigen::MatrixXd g_mu_cov = output_3.mat1;
  Eigen::MatrixXd g_x_cov = output_3.mat2;
  Eigen::MatrixXd g_c_cov = output_3.mat3;
  
  
  Eigen::MatrixXd g_cu_mu_l_mu_sig(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
  g_cu_mu_l_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
  g_cu_mu_l_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
  
  Eigen::MatrixXd g_cu_mu_l_c(g_c_mean.rows()+2,5);
  g_cu_mu_l_c.setZero();
  g_cu_mu_l_c.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
  g_cu_mu_l_c.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
  
  
  Eigen::VectorXd g_c(4);
  g_c.setZero();
  
  output_3 = TVBS_grad_non_cdf_bvn_by_cdfn_cpp(mu_temp_1.segment(2,2), sigma_temp_1.block(2,2,2,2), x_temp.segment(2,2));
  Eigen::VectorXd g_from_p_mu = output_3.mat1;
  Eigen::VectorXd g_from_p_cov = output_3.mat2;
  Eigen::VectorXd g_from_p_c = output_3.mat3;
  
  
  Eigen::VectorXd g_rho_rho = p(0)*(g_cu_mu_l_mu_sig.block(0,0,g_cu_mu_l_mu_sig.rows(),2)*g_from_p_mu + g_cu_mu_l_mu_sig.block(0,2,g_cu_mu_l_mu_sig.rows(),3)*g_from_p_cov);
  
  g_c = p(0)*(g_cu_mu_l_c.block(0,0,g_cu_mu_l_c.rows(),2)*g_from_p_mu + g_cu_mu_l_c.block(0,2,g_cu_mu_l_c.rows(),3)*g_from_p_cov);
  
  // getting gradient contribution directly from noncdfbvn for absiccae except the first two
  g_c.segment(2,2) = p(0)*g_from_p_c;
  
  output_2 = TVBS_grad_cdf_tvn_cpp(x_temp.segment(0,3), Corr_temp.block(0,0,3,3));
  Eigen::VectorXd g_c_1 = output_2.mat1;
  Eigen::VectorXd g_rho_1 = output_2.mat2;
  
  Eigen::VectorXd g_rho_rho_add(6);
  g_rho_rho_add << g_rho_1.head(2), 0, g_rho_1(2), 0, 0;
  g_rho_rho += p(1)*g_rho_rho_add;
  
  g_c.head(g_c_1.size()) += p(1)*g_c_1;
  // restore order, if reordered earlier, according to temp1
  
  Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
  g_rho_rho_mat += g_rho_rho_mat.transpose();
  
  // reorder the g_rho_rho_mat according to temp1
  g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat, 0);
  
  struct TVBS_Matrix_two output_final = {g_c, g_rho_rho};
  return output_final;
}


struct TVBS_Matrix_two TVBS_grad_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor)
{
  Eigen::VectorXd mu_temp_4(4);
  mu_temp_4.setZero();
  
  double quad_var_cdf = TVBS_pmvnorm_cpp(w.head(4), mu_temp_4, cor.block(0,0,4,4));
  double bi_var_cdf = TVBS_pmvnorm_cpp(w.head(2), mu_temp_4.head(2), cor.block(0,0,2,2));
  
  
  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(w.head(2),cor.block(0,0,2,2));
  double g_w_1 = output_3.d1;
  double g_w_2 = output_3.d2;
  double g_rho = output_3.d3;
  Eigen::VectorXd g_w(2);
  g_w << g_w_1, g_w_2;
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_qvn_cpp(w, cor);
  Eigen::VectorXd g_w_quad = output_2.mat1;
  Eigen::VectorXd g_rho_quad = output_2.mat2;
  
  
  Eigen::VectorXd g_w1_w2_new = (bi_var_cdf*g_w_quad.head(2) - quad_var_cdf*g_w)/pow(bi_var_cdf,2);
  double g_rho12_new = (bi_var_cdf*g_rho_quad(0) - quad_var_cdf*g_rho)/pow(bi_var_cdf,2);
  
  Eigen::VectorXd output_final_1(4), output_final_2(6);
  output_final_1 << g_w1_w2_new, g_w_quad.segment(2,2)/bi_var_cdf;
  output_final_2 << g_rho12_new, g_rho_quad.segment(1,5)/bi_var_cdf;
  
  struct TVBS_Matrix_two output_final = {output_final_1, output_final_2};
  return output_final;
}











// strange behaviour with cholesky_global==1
struct TVBS_Matrix_three TVBS_grad_non_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  // normalize x and calculate correlation matrix from covariance matrix
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Sigma_diag_sqrt_inv(4,4);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  
  // cap the x values over 6 and under -6
  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_qvn_by_cdf_bvn_cpp(x_norm.head(4), Cor_mat.block(0,0,4,4));
  Eigen::VectorXd g_w = output_2.mat1;
  Eigen::VectorXd g_rho = output_2.mat2;
  
  output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  Eigen::MatrixXd g_b_cor_cov = output_2.mat1;
  Eigen::MatrixXd g_omega_cor_cov = output_2.mat2;
  
  Eigen::VectorXd g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  Eigen::VectorXd g_cov = g_b_cor_cov*g_w + g_omega_cor_cov*g_rho;
  Eigen::VectorXd g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  
  return output_final;
  
}

Eigen::VectorXd TVBS_grad_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd grad(14);
  grad.setZero();
  
  // derivative with respect to b_i. 
  
  Eigen::VectorXd bh=b;
  Eigen::MatrixXd Sigmah = Sigma;
  
  for (int i=0;i<4;i++){
    double f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
    Eigen::MatrixXd R21 = Sigmah.block(1,0,3,1) / Sigmah(0,0); 
    Eigen::VectorXd mub = bh.block(1,0,3,1) - R21 * bh(0);
    Eigen::MatrixXd Sigb = Sigmah.block(1,1,3,3) - R21 * Sigmah(0,0) * R21.transpose();  
    
    double f2 = TVBS_pmvnorm_cpp(mub,mub*0,Sigb);
    
    grad(i) = f1* f2 /sqrt(Sigmah(0,0));
    
    // rotate 
    Eigen::MatrixXd ROT(4,4);
    ROT.setZero();
    
    ROT(0,1) =1;
    ROT(1,2) =1;
    ROT(2,3)= 1;
    ROT(3,0)=1;
    
    bh = ROT * bh;
    Sigmah = ROT*Sigmah * ROT.transpose();
    
  }
  
  // derivative with respect to rho_{kl}
  Eigen::MatrixXd Gr(4,4);
  Gr.setZero();
  
  Eigen::MatrixXd D(4,4);
  D.setZero();
  for (int d=0;d<4;d++){ 
    D(d,d)= sqrt(Sigma(d,d));
  }
  bh=D.inverse()* b;
  Sigmah = D.inverse()* Sigma * D.inverse();
  
  Eigen::VectorXd ord(4);
  ord(0)=1;
  ord(1)=2;
  ord(2)=3;
  ord(3)=4;
  
  for (int a=0;a<4;a++){
    for (int b=a+1;b<4;b++){
      
      
      Eigen::MatrixXd Sigr(2,2);
      Sigr = Sigmah.block(2,2,2,2);
      
      Eigen::MatrixXd Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,2) * Sigr.inverse() * Sigmah.block(2,0,2,2);
      Eigen::VectorXd bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,2) * Sigr.inverse() *bh.block(2,0,2,1);
      double f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
      
      Eigen::MatrixXd Sh = Sigmah.block(0,0,2,2);
      Eigen::MatrixXd R12 = Sigmah.block(2,0,2,2) * Sh.inverse(); 
      Eigen::VectorXd mub = bh.block(2,0,2,1) - (R12.block(0,0,2,1) * bh(0)) - (R12.block(0,1,2,1)* bh(1));
      Eigen::MatrixXd Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
      
      double f2 = biv_gen_cdf(mub, Sigr - Sigb);
      
      
      Gr(a,b) = f1* f2;
      
      
      
      // rotate second entry to the back. 
      Eigen::MatrixXd ROT(4,4);
      ROT.setZero();
      
      ROT(0,0) =1;
      ROT(1,2) =1;
      ROT(2,3)=1;
      ROT(3,1) =1; 
      
      bh = ROT * bh;
      Sigmah = ROT*Sigmah * ROT.transpose();
      
      ord = ROT * ord;
    }
    
    // restore original order
    bh=D.inverse()* b;
    Sigmah = D.inverse()* Sigma * D.inverse();
    ord(0)=1;
    ord(1)=2;
    ord(2)=3;
    ord(3)=4;
    
    // rotate first entry to last entry. 
    Eigen::MatrixXd ROT(4,4);
    ROT.setZero();
    
    ROT(0,1) =1;
    ROT(1,2) =1;
    ROT(2,3)=1;
    ROT(3,0)=1; 
    // now rotate a+1 times 
    for (int c=0;c<a+1;c++){
      bh = ROT * bh;
      Sigmah = ROT*Sigmah * ROT.transpose();
      ord = ROT * ord;
      
    }
  }
  
  // fill in results into grad    
  grad(5) = Gr(0,1) / (D(0,0)* D(1,1)); // entry sigma(1,2)
  grad(6) = Gr(0,2) / (D(0,0)* D(2,2)); // entry sigma(1,3)
  grad(7) = Gr(0,3) / (D(0,0)* D(3,3)); // entry sigma(1,4)
  grad(9) = Gr(1,2) / (D(1,1)* D(2,2)); // entry sigma(2,3)
  grad(10) = Gr(1,3) / (D(1,1)* D(3,3)); // entry sigma(2,4)
  grad(12) = Gr(2,3) / (D(2,2)* D(3,3)); // entry sigma(3,4)
  
  // finally diagonal entries
  grad(4) =  (-grad(0)*b(0)-grad(5)*Sigma(0,1) - grad(6)*Sigma(0,2) - grad(7)*Sigma(0,3))/(2*Sigma(0,0)); 
  grad(8) =  (-grad(1)*b(1)-grad(5)*Sigma(0,1) - grad(9)*Sigma(1,2) - grad(10)*Sigma(1,3))/(2*Sigma(1,1)); 
  grad(11) = (-grad(2)*b(2)-grad(6)*Sigma(0,2) - grad(9)*Sigma(1,2) - grad(12)*Sigma(2,3))/(2*Sigma(2,2)); 
  grad(13) = (-grad(3)*b(3)-grad(7)*Sigma(0,3) - grad(10)*Sigma(1,3)- grad(12)*Sigma(2,3))/(2*Sigma(3,3));
  
  return grad;
}






// calculates the second derivative with respect to twice elements sigma_{kl}. 
// is called several times in the Hessian calculation.
Eigen::VectorXd cal_grad_rhokl(Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::VectorXd He_onetofour)
{
  
  double Gr=0;
  
  Eigen::MatrixXd D(4,4);
  D.setZero();
  for (int d=0;d<4;d++){ 
    D(d,d)= sqrt(Sigma(d,d));
    if (D(d,d)<0.00001)
    {   
      D(d,d)= 0.00001;
    }
  }
  
  Eigen::VectorXd bh(4);
  bh.setZero();
  Eigen::MatrixXd Sigmah(4,4);
  Sigmah.setZero();
  
  bh=D.inverse()* b;
  Sigmah = D.inverse()* Sigma * D.inverse();
  
  
  // order (0,1,2,3)
  Eigen::MatrixXd Sigr(2,2);
  Sigr = Sigmah.block(2,2,2,2);
  
  Eigen::MatrixXd Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,2) * Sigr.inverse() * Sigmah.block(2,0,2,2);
  Eigen::VectorXd bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,2) * Sigr.inverse() *bh.block(2,0,2,1);
  double f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
  
  Eigen::MatrixXd Sh = Sigmah.block(0,0,2,2);
  Eigen::MatrixXd R12 = Sigmah.block(2,0,2,2) * Sh.inverse(); 
  Eigen::VectorXd mub = bh.block(2,0,2,1) - (R12.block(0,0,2,1) * bh(0)) - (R12.block(0,1,2,1)* bh(1));
  Eigen::MatrixXd Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
  double f2 = biv_gen_cdf(mub, Sigr - Sigb);
  
  Gr = f1* f2;
  
  // now start Hessian calculations. 
  // f1 derived with respect to (s_00,s_01,s_02,s_03,s_11,s_12,s_13,s_22,s_23,s_33)
  Eigen::MatrixXd df1(1,10);
  df1.setZero();
  
  double detSigma = Sigmah(0,0)*Sigmah(1,1) - Sigmah(0,1)*Sigmah(1,0);
  if (detSigma<0.00001){ detSigma = 0.00001;}
  double bexp = bh(0)*bh(0)*Sigmah(1,1) + bh(1)*bh(1)*Sigmah(0,0) - 2*bh(0)*bh(1)*Sigmah(0,1);
  
  Eigen::MatrixXd ddetSigma(1,10);
  ddetSigma.setZero();
  
  ddetSigma(0,0)= Sigmah(1,1);
  ddetSigma(0,1)= -2*Sigmah(0,1);
  ddetSigma(0,4) = Sigmah(0,0);
  
  df1 = f1*(-1/(2*detSigma)+bexp/(2*detSigma*detSigma)) * ddetSigma;
  df1(0,0) += (-f1)*(bh(1)*bh(1)/(2*(detSigma)));
  df1(0,1) += (f1)*(bh(0)*bh(1)/((detSigma)));
  df1(0,4) += (-f1)*(bh(0)*bh(0)/(2*(detSigma)));
  
  
  // f2 derived with respect to (s_00,s_01,s_02,s_03,s_11,s_12,s_13,s_22,s_23,s_33)
  Eigen::MatrixXd sig_bar = Sigr - Sigb;
  Eigen::MatrixXd df2(1,10);
  df2.setZero();
  
  Eigen::MatrixXd db_bar(2,10);
  db_bar.setZero();
  
  Eigen::MatrixXd dsig_bar(3,10);
  dsig_bar.setZero();
  
  // f2 derived with respect to (s_00,s_01,s_02,s_03,s_11,s_12,s_13,s_22,s_23,s_33) 
  // ordering: R_12(0,0),R_12(0,1);R_12(1,0),R_12(1,1)
  Eigen::MatrixXd dR_12(4,10);
  dR_12.setZero();
  
  dR_12(0,1) = 2*(Sigmah(0,2)-Sigmah(1,2)*Sigmah(0,1))/(detSigma*detSigma)*Sigmah(0,1);
  dR_12(0,1) += (-1/detSigma)*Sigmah(1,2);
  dR_12(0,2) += 1/detSigma;
  dR_12(0,5) = (-1/detSigma)*Sigmah(0,1);
  
  dR_12(1,1) = 2*(-Sigmah(0,2)*Sigmah(0,1)+Sigmah(1,2))/(detSigma*detSigma)*Sigmah(0,1);
  dR_12(1,1) += (-1/detSigma)*Sigmah(0,2);
  dR_12(1,2) = -Sigmah(0,1)/detSigma;
  dR_12(1,5) += 1/detSigma;
  
  // R_12(1,0) = (s_03s_11 - s_31 s_01)/detSigma 
  dR_12(2,1)= 2*(Sigmah(0,3)-Sigmah(1,3)*Sigmah(0,1))/(detSigma*detSigma)*Sigmah(0,1);
  dR_12(2,1) += -Sigmah(1,3)/detSigma;
  dR_12(2,3) = Sigmah(1,1)/detSigma;
  dR_12(2,6) += -Sigmah(0,1)/detSigma;
  
  // R_12(1,1) = (-s_03s_01 + s_31 s_00)/detSigma 
  dR_12(3,1) += -Sigmah(0,3)/detSigma;
  dR_12(3,1) += 2*(-Sigmah(0,3)*Sigmah(0,1)+Sigmah(1,3)*Sigmah(0,0))/(detSigma*detSigma)*Sigmah(0,1);
  dR_12(3,3) = -Sigmah(0,1)/detSigma;
  dR_12(3,6) += 1/detSigma;
  
  
  db_bar.block(0,0,1,10) = -bh(0)*dR_12.block(0,0,1,10) - bh(1)* dR_12.block(1,0,1,10);
  db_bar.block(1,0,1,10) = -bh(0)*dR_12.block(2,0,1,10) - bh(1)* dR_12.block(3,0,1,10);
  
  // element (0,0)
  dsig_bar.block(0,0,1,10) = -Sigmah(0,2)*dR_12.block(0,0,1,10) - Sigmah(1,2)*dR_12.block(1,0,1,10);
  dsig_bar(0,2) += -R12(0,0);
  dsig_bar(0,5) += -R12(0,1);
  dsig_bar(0,7) = 1; 
  
  // element (0,1)
  dsig_bar.block(1,0,1,10) = -Sigmah(0,3)*dR_12.block(0,0,1,10) - Sigmah(1,3)*dR_12.block(1,0,1,10);
  dsig_bar(1,3) += -R12(0,0);
  dsig_bar(1,6) += -R12(0,1);
  dsig_bar(1,8) = 1; 
  
  // element (1,1)
  dsig_bar.block(2,0,1,10) = - Sigmah(0,3)*dR_12.block(2,0,1,10) - Sigmah(1,3)*dR_12.block(3,0,1,10);
  dsig_bar(2,3) += -R12(1,0);
  dsig_bar(2,6) += -R12(1,1);
  dsig_bar(2,9) = 1; 
  
  
  Eigen::VectorXd gen_grad = grad_gen_cdf(mub, Sigr - Sigb); //normal_pdf(mu_b,sig_bar)*(db_bar - mu_b*dsig_bar/(2*sig_bar)); 
  df2 = gen_grad(0) * db_bar.block(0,0,1,10) + gen_grad(1)* db_bar.block(1,0,1,10);
  df2 += gen_grad(2)* dsig_bar.block(0,0,1,10) +gen_grad(4)* dsig_bar.block(1,0,1,10) + gen_grad(3)* dsig_bar.block(2,0,1,10);  
  
  
  // combine the pieces 
  // grad = f1*f2 / sqrt(sigma_11 *sigma_22)
  Eigen::VectorXd He_line(14);
  He_line.setZero();
  
  He_line.head(4) = He_onetofour;
  
  He_line.tail(10) = (df1*f2 + f1*df2)/(D(0,0)*D(1,1));
  
  He_line(5) = He_line(5)/(D(0,0)*D(1,1));
  He_line(6) = He_line(6)/(D(0,0)*D(2,2));
  He_line(7) = He_line(7)/(D(0,0)*D(3,3));
  He_line(8) = He_line(8)/(D(1,1)*D(1,1));
  He_line(9) = He_line(9)/(D(1,1)*D(2,2));
  He_line(10) = He_line(10)/(D(1,1)*D(3,3));
  He_line(11) = He_line(11)/(D(2,2)*D(2,2));
  He_line(12) = He_line(12)/(D(2,2)*D(3,3));
  He_line(13) = He_line(13)/(D(3,3)*D(3,3));
  
  
  // derivative w.r.t. sigma_11 reconstructed from others. 
  He_line(4) = -(He_line(0)*b(0) + He_line(5)*Sigma(0,1)+ He_line(6)*Sigma(0,2)+ He_line(7)*Sigma(0,3)+ (Gr / (D(0,0)* D(1,1))))/(2*Sigma(0,0));
  He_line(8) = -(He_line(1)*b(1) + He_line(5)*Sigma(0,1)+ He_line(9)*Sigma(1,2)+ He_line(10)*Sigma(1,3) +(Gr / (D(0,0)* D(1,1))))/(2*Sigma(1,1));
  He_line(11) = -(He_line(2)*b(2) + He_line(6)*Sigma(0,2)+ He_line(9)*Sigma(1,2)+ He_line(12)*Sigma(2,3))/(2*Sigma(2,2));
  He_line(13) = -(He_line(3)*b(3) + He_line(7)*Sigma(0,3)+ He_line(10)*Sigma(1,3)+ He_line(12)*Sigma(2,3))/(2*Sigma(3,3));  
  
  return He_line;
}


Eigen::MatrixXd Hessian_gen_qvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  // initialize 
  Eigen::MatrixXd Hess(14,14);
  Hess.setZero();
  Eigen::VectorXd bh = b;
  Eigen::MatrixXd Sigmah = Sigma;
  
  // cycle over rows. 
  // first with respect to b1
  
  double f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  double f1sqrt = f1/sqrt(Sigmah(0,0));
  
  Eigen::MatrixXd R21 = Sigmah.block(1,0,3,1) / Sigmah(0,0); 
  Eigen::VectorXd mub = bh.block(1,0,3,1) - R21 * bh(0);
  Eigen::MatrixXd Sigb = Sigmah.block(1,1,3,3) - R21 * Sigmah(0,0) * R21.transpose();  
  double f2 = TVBS_pmvnorm_cpp(mub,mub*0,Sigb);
  Eigen::VectorXd gr = grad_gen_tvn_cdf(mub, Sigb);
  
  // second deriv w.r.t. bi 
  Hess(0,0) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2)+gr(2)*Sigmah(0,3))/Sigmah(0,0);
  
  Hess(0,1) = f1sqrt*(gr(0));
  Hess(0,2) = f1sqrt*(gr(1));
  Hess(0,3) = f1sqrt*(gr(2));
  
  
  // derivative with respect to sigma_{11}. 
  Eigen::VectorXd dSigb_dSig(9);
  dSigb_dSig.setZero();
  
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2)= bh(0)*Sigmah(0,3);
  dSigb_dSig(3) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(5) = Sigmah(0,3)*Sigmah(0,1);
  dSigb_dSig(6) = Sigmah(0,2)*Sigmah(0,2);
  dSigb_dSig(7) = Sigmah(0,2)*Sigmah(0,3);
  dSigb_dSig(8) = Sigmah(0,3)*Sigmah(0,3);
  
  
  double pr = (gr.transpose() * dSigb_dSig);
  Hess(0,4) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  
  // derivative w.r.t sigma_{1,2} 
  Hess(0,5) = -f1sqrt*( gr(0)*bh(0)+ gr(3)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2)+gr(5)*Sigmah(0,3))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,3}
  Hess(0,6) =  -f1sqrt*( gr(1)*bh(0) + gr(6)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1)+gr(7)*Sigmah(0,3))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,4}
  Hess(0,7) =  -f1sqrt*( gr(2)*bh(0) + gr(8)*2*Sigmah(0,3)+gr(5)*Sigmah(0,1)+gr(7)*Sigmah(0,2))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{2,2}
  Hess(0,8)= f1sqrt*gr(3);
  
  // derivative w.r.t sigma_{2,3} 
  Hess(0,9) = f1sqrt*gr(4);
  
  // derivative w.r.t sigma_{2,4} 
  Hess(0,10) = f1sqrt*gr(5);
  
  
  // derivative w.r.t sigma_{3,3}
  Hess(0,11) = f1sqrt*gr(6);
  
  // derivative w.r.t sigma_{3,4}
  Hess(0,12) = f1sqrt*gr(7);
  
  // derivative w.r.t sigma_{4,4}
  Hess(0,13) = f1sqrt*gr(8);
  
  
  
  // next up: derivativves with respect to b2. 
  // rotate 
  Eigen::MatrixXd ROT(4,4);
  ROT.setZero();
  
  ROT(0,1) =1;
  ROT(1,2) =1;
  ROT(2,3)= 1;
  ROT(3,0)=1;
  
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new ordering: (1,2,3,0).
  f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  f1sqrt = f1/sqrt(Sigmah(0,0));
  R21 = Sigmah.block(1,0,3,1) / Sigmah(0,0); 
  mub = bh.block(1,0,3,1) - R21 * bh(0);
  Sigb = Sigmah.block(1,1,3,3) - R21 * Sigmah(0,0) * R21.transpose();  
  f2 = TVBS_pmvnorm_cpp(mub,mub*0,Sigb);
  gr = grad_gen_tvn_cdf(mub, Sigb);
  
  // second deriv w.r.t. b1 
  Hess(1,1) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2)+gr(2)*Sigmah(0,3))/Sigmah(0,0);
  
  Hess(1,2) = f1sqrt*(gr(0));
  Hess(1,3) = f1sqrt*(gr(1));
  Hess(1,0) = f1sqrt*(gr(2));
  
  
  // derivative with respect to sigma_{22}. 
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2)= bh(0)*Sigmah(0,3);
  dSigb_dSig(3) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(5) = Sigmah(0,3)*Sigmah(0,1);
  dSigb_dSig(6) = Sigmah(0,2)*Sigmah(0,2);
  dSigb_dSig(7) = Sigmah(0,2)*Sigmah(0,3);
  dSigb_dSig(8) = Sigmah(0,3)*Sigmah(0,3);
  
  
  pr = (gr.transpose() * dSigb_dSig);
  // derivative w.r.t. sigma_{1,1}
  Hess(1,8) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  
  // derivative w.r.t sigma_{1,2} 
  Hess(1,9) = -f1sqrt*( gr(0)*bh(0)+ gr(3)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2)+gr(5)*Sigmah(0,3))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,3}
  Hess(1,10) =  -f1sqrt*( gr(1)*bh(0) + gr(6)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1)+gr(7)*Sigmah(0,3))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,0}
  Hess(1,5) =  -f1sqrt*( gr(2)*bh(0) + gr(8)*2*Sigmah(0,3)+gr(5)*Sigmah(0,1)+gr(7)*Sigmah(0,2))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{2,2}
  Hess(1,11)= f1sqrt*gr(3);
  
  // derivative w.r.t sigma_{2,3} 
  Hess(1,12) = f1sqrt*gr(4);
  
  // derivative w.r.t sigma_{2,0} 
  Hess(1,6) = f1sqrt*gr(5);
  
  
  // derivative w.r.t sigma_{3,3}
  Hess(1,13) = f1sqrt*gr(6);
  
  // derivative w.r.t sigma_{3,0}
  Hess(1,7) = f1sqrt*gr(7);
  
  // derivative w.r.t sigma_{0,0}
  Hess(1,4) = f1sqrt*gr(8);
  
  
  // derivative with respect to b3. 
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new ordering: (2,3,0,1).
  
  f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  f1sqrt = f1/sqrt(Sigmah(0,0));
  R21 = Sigmah.block(1,0,3,1) / Sigmah(0,0); 
  mub = bh.block(1,0,3,1) - R21 * bh(0);
  Sigb = Sigmah.block(1,1,3,3) - R21 * Sigmah(0,0) * R21.transpose();  
  f2 = TVBS_pmvnorm_cpp(mub,mub*0,Sigb);
  gr = grad_gen_tvn_cdf(mub, Sigb);
  
  // second deriv w.r.t. b1 
  Hess(2,2) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2)+gr(2)*Sigmah(0,3))/Sigmah(0,0);
  
  Hess(2,3) = f1sqrt*(gr(0));
  Hess(2,0) = f1sqrt*(gr(1));
  Hess(2,1) = f1sqrt*(gr(2));
  
  // derivative with respect to sigma_{22}. 
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2)= bh(0)*Sigmah(0,3);
  dSigb_dSig(3) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(5) = Sigmah(0,3)*Sigmah(0,1);
  dSigb_dSig(6) = Sigmah(0,2)*Sigmah(0,2);
  dSigb_dSig(7) = Sigmah(0,2)*Sigmah(0,3);
  dSigb_dSig(8) = Sigmah(0,3)*Sigmah(0,3);
  
  pr = (gr.transpose() * dSigb_dSig);
  // derivative w.r.t. sigma_{2,2}
  Hess(2,11) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  // de2ivative w.r.t sigma_{2,3} 
  Hess(2,12) = -f1sqrt*( gr(0)*bh(0)+ gr(3)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2)+gr(5)*Sigmah(0,3))/Sigmah(0,0);
  // derivative w.r.t. sigma_{2,0}
  Hess(2,6) =  -f1sqrt*( gr(1)*bh(0) + gr(6)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1)+gr(7)*Sigmah(0,3))/Sigmah(0,0);
  // derivative w.r.t. sigma_{2,1}
  Hess(2,9) =  -f1sqrt*( gr(2)*bh(0) + gr(8)*2*Sigmah(0,3)+gr(5)*Sigmah(0,1)+gr(7)*Sigmah(0,2))/Sigmah(0,0);
  // derivative w.r.t. sigma_{3,3}
  Hess(2,13)= f1sqrt*gr(3);
  // derivative w.r.t sigma_{3,0} 
  Hess(2,7) = f1sqrt*gr(4);
  // derivative w.r.t sigma_{3,1} 
  Hess(2,10) = f1sqrt*gr(5);
  // derivative w.r.t sigma_{0,0}
  Hess(2,4) = f1sqrt*gr(6);
  // derivative w.r.t sigma_{0,1}
  Hess(2,5) = f1sqrt*gr(7);
  // derivative w.r.t sigma_{1,1}
  Hess(2,8) = f1sqrt*gr(8);
  
  
  // derivative with respect to b4. 
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new ordering: (3,0,1,2).
  
  f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  f1sqrt = f1/sqrt(Sigmah(0,0));
  R21 = Sigmah.block(1,0,3,1) / Sigmah(0,0); 
  mub = bh.block(1,0,3,1) - R21 * bh(0);
  Sigb = Sigmah.block(1,1,3,3) - R21 * Sigmah(0,0) * R21.transpose();  
  f2 = TVBS_pmvnorm_cpp(mub,mub*0,Sigb);
  gr = grad_gen_tvn_cdf(mub, Sigb);
  
  // second deriv w.r.t. b1 
  Hess(3,3) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2)+gr(2)*Sigmah(0,3))/Sigmah(0,0);
  
  Hess(3,0) = f1sqrt*(gr(0));
  Hess(3,1) = f1sqrt*(gr(1));
  Hess(3,2) = f1sqrt*(gr(2));
  
  // derivative with respect to sigma_{33}. 
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2)= bh(0)*Sigmah(0,3);
  dSigb_dSig(3) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(5) = Sigmah(0,3)*Sigmah(0,1);
  dSigb_dSig(6) = Sigmah(0,2)*Sigmah(0,2);
  dSigb_dSig(7) = Sigmah(0,2)*Sigmah(0,3);
  dSigb_dSig(8) = Sigmah(0,3)*Sigmah(0,3);
  
  pr = (gr.transpose() * dSigb_dSig);
  // derivative w.r.t. sigma_{3,3}
  Hess(3,13) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  // de2ivative w.r.t sigma_{3,0} 
  Hess(3,7) = -f1sqrt*( gr(0)*bh(0)+ gr(3)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2)+gr(5)*Sigmah(0,3))/Sigmah(0,0);
  // derivative w.r.t. sigma_{3,1}
  Hess(3,10) =  -f1sqrt*( gr(1)*bh(0) + gr(6)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1)+gr(7)*Sigmah(0,3))/Sigmah(0,0);
  // derivative w.r.t. sigma_{3,2}
  Hess(3,12) =  -f1sqrt*( gr(2)*bh(0) + gr(8)*2*Sigmah(0,3)+gr(5)*Sigmah(0,1)+gr(7)*Sigmah(0,2))/Sigmah(0,0);
  // derivative w.r.t. sigma_{0,0}
  Hess(3,4)= f1sqrt*gr(3);
  // derivative w.r.t sigma_{0,1} 
  Hess(3,5) = f1sqrt*gr(4);
  // derivative w.r.t sigma_{0,2} 
  Hess(3,6) = f1sqrt*gr(5);
  // derivative w.r.t sigma_{1,1}
  Hess(3,8) = f1sqrt*gr(6);
  // derivative w.r.t sigma_{1,2}
  Hess(3,9) = f1sqrt*gr(7);
  // derivative w.r.t sigma_{2,2}
  Hess(3,11) = f1sqrt*gr(8);
  
  
  // use symmetry of Hessian
  
  Eigen::MatrixXd He = Hess + Hess.transpose();
  
  He.block(0,0,4,4) = He.block(0,0,4,4)/2; 
  
  ////////////////////////////////////////
  /// done with derivatives involving b //
  ////////////////////////////////////////
  
  
  ////////////////////////////
  /// sigma_01/////
  ////////////////////////////
  Eigen::VectorXd He_line(14);
  He_line.setZero();
  Eigen::VectorXd he_start(4);
  for (int j=0;j<4;j++){ he_start(j) = He(5,j);}
  
  
  He_line = cal_grad_rhokl(b, Sigma, he_start);
  
  for (int j=4;j<14;j++){
    He(5,j) = He_line(j); 
  }
  
  
  ///////////////////////////////////////////////
  //// derivative w.r.t. sigma_02           ///
  ///////////////////////////////////////////////
  
  // rotate second entry to the back. 
  ROT.setZero();
  
  ROT(0,0) =1;
  ROT(1,2) =1;
  ROT(2,3)=1;
  ROT(3,1)=1; 
  
  // new order: (1,3,4,2) 
  
  
  Eigen::VectorXd btil = ROT * b;
  Eigen::MatrixXd Sigmatil = ROT * Sigma * ROT.transpose(); 
  
  for (int j=0;j<4;j++){ he_start(j) = He(6,j);}
  
  
  He_line = cal_grad_rhokl(btil, Sigmatil, ROT* he_start);
  
  // reorder
  He(6,4) = He_line(4);
  He(6,5) = He_line(7);
  He(6,6) = He_line(5);
  He(6,7) = He_line(6);
  He(6,8) = He_line(13);
  He(6,9) = He_line(10);
  He(6,10) = He_line(12);
  He(6,11) = He_line(8);
  He(6,12) = He_line(9);
  He(6,13) = He_line(11);
  
  ///////////////////////////////////////////////
  //// derivative w.r.t. sigma_03            ///
  ///////////////////////////////////////////////
  
  // exchange second and fourth entry. 
  ROT.setZero();
  
  ROT(0,0) =1;
  ROT(1,3) =1;
  ROT(2,2)=1;
  ROT(3,1)=1; 
  
  // new order: (1,3,4,2) 
  btil = ROT * b;
  Sigmatil = ROT * Sigma * ROT.transpose(); 
  
  for (int j=0;j<4;j++){ he_start(j) = He(7,j);}
  
  
  He_line = cal_grad_rhokl(btil, Sigmatil, ROT* he_start);
  
  // reorder
  He(7,4) = He_line(4);
  He(7,5) = He_line(7);
  He(7,6) = He_line(6);
  He(7,7) = He_line(5);
  He(7,8) = He_line(13);
  He(7,9) = He_line(12);
  He(7,10) = He_line(10);
  He(7,11) = He_line(11);
  He(7,12) = He_line(9);
  He(7,13) = He_line(8);
  
  
  ///////////////////////////////////////////////
  //// derivative w.r.t. sigma_12            ///
  ///////////////////////////////////////////////
  
  // rotate first entry to last position.  
  ROT.setZero();
  
  ROT(0,1) =1;
  ROT(1,2) =1;
  ROT(2,3)=1;
  ROT(3,0)=1; 
  
  // new order: (2,3,4,1) 
  btil = ROT * b;
  Sigmatil = ROT * Sigma * ROT.transpose(); 
  
  for (int j=0;j<4;j++){ he_start(j) = He(9,j);}
  
  
  
  He_line = cal_grad_rhokl(btil, Sigmatil, ROT* he_start);
  
  // reorder
  He(9,4) = He_line(13);
  He(9,5) = He_line(7);
  He(9,6) = He_line(10);
  He(9,7) = He_line(12);
  He(9,8) = He_line(4);
  He(9,9) = He_line(5);
  He(9,10) = He_line(6);
  He(9,11) = He_line(8);
  He(9,12) = He_line(9);
  He(9,13) = He_line(11);
  
  ///////////////////////////////////////////////
  //// derivative w.r.t. sigma_13            ///
  ///////////////////////////////////////////////
  
  // exchange first and last entry.  
  ROT.setZero();
  
  ROT(0,3) =1;
  ROT(1,1) =1;
  ROT(2,2)=1;
  ROT(3,0)=1; 
  
  // new order: (4,2,3,1) 
  btil = ROT * b;
  Sigmatil = ROT * Sigma * ROT.transpose(); 
  
  for (int j=0;j<4;j++){ he_start(j) = He(10,j);}
  
  He_line = cal_grad_rhokl(btil, Sigmatil, ROT* he_start);
  
  // reorder
  He(10,4) = He_line(13);
  He(10,5) = He_line(10);
  He(10,6) = He_line(12);
  He(10,7) = He_line(7);
  He(10,8) = He_line(8);
  He(10,9) = He_line(9);
  He(10,10) = He_line(5);
  He(10,11) = He_line(11);
  He(10,12) = He_line(6);
  He(10,13) = He_line(4);
  
  ///////////////////////////////////////////////
  //// derivative w.r.t. sigma_23            ///
  ///////////////////////////////////////////////
  
  // exchange first two and last two entries.  
  ROT.setZero();
  
  ROT(0,2) =1;
  ROT(1,3) =1;
  ROT(2,0)=1;
  ROT(3,1)=1; 
  
  // new order: (3,4,1,2) 
  btil = ROT * b;
  Sigmatil = ROT * Sigma * ROT.transpose(); 
  
  for (int j=0;j<4;j++){ he_start(j) = He(12,j);}
  
  He_line = cal_grad_rhokl(btil, Sigmatil, ROT*  he_start);
  
  // reorder
  He(12,4) = He_line(11);
  He(12,5) = He_line(12);
  He(12,6) = He_line(6);
  He(12,7) = He_line(9);
  He(12,8) = He_line(13);
  He(12,9) = He_line(7);
  He(12,10) = He_line(10);
  He(12,11) = He_line(4);
  He(12,12) = He_line(5);
  He(12,13) = He_line(8);
  
  // now fill in rows for columns for non-diagonal entries
  // 5th col.
  for (int k=4;k<14;k++){
    He(k,6)=He(6,k);
    He(k,7)=He(7,k);
    He(k,5)=He(5,k);
    He(k,9)=He(9,k);
    He(k,10)=He(10,k);
    He(k,12)=He(12,k);
  }
  
  ///////////////////////////////////////////
  ////// entries w.r.t. diagonal elements  //
  ///////////////////////////////////////////
  
  // first get gradient
  Eigen::VectorXd grad = TVBS_grad_gen_qvn_cdf(b, Sigma); 
  
  //////////////
  // sigma_00 
  /////////////
  He(4,4) = (-He(0,4)*b(0)   - He(5,4)*Sigma(0,1)  - He(6,4)*Sigma(0,2) - He(7,4)*Sigma(0,3) - 2*grad(4))/(2*Sigma(0,0));
  He(4,8) = (-He(0,8)*b(0)   - He(5,8)*Sigma(0,1)  - He(6,8)*Sigma(0,2) - He(7,8)*Sigma(0,3))/(2*Sigma(0,0));
  He(4,11) = (-He(0,11)*b(0) - He(5,11)*Sigma(0,1) - He(6,11)*Sigma(0,2)- He(7,11)*Sigma(0,3))/(2*Sigma(0,0));
  He(4,13) = (-He(0,13)*b(0) - He(5,13)*Sigma(0,1) - He(6,13)*Sigma(0,2)- He(7,13)*Sigma(0,3))/(2*Sigma(0,0));
  
  He(8,4) = He(4,8);
  He(11,4) = He(4,11); 
  He(13,4) = He(4,13); 
  
  //////////////
  // sigma_11 
  /////////////
  He(8,8) = (-He(1,8)*b(1)   - He(5,8)*Sigma(0,1)  - He(9,8)*Sigma(1,2)   - He(10,8)*Sigma(1,3) - 2*grad(8))/(2*Sigma(1,1));
  He(11,8) = (-He(1,11)*b(1) - He(5,11)*Sigma(0,1) - He(9,11)*Sigma(1,2)- He(10,11)*Sigma(1,3))/(2*Sigma(1,1));
  He(13,8) = (-He(1,13)*b(1) - He(5,13)*Sigma(0,1) - He(9,13)*Sigma(1,2)- He(10,13)*Sigma(1,3))/(2*Sigma(1,1));
  
  He(8,11)= He(11,8); 
  He(8,13)= He(13,8); 
  
  //////////////
  // sigma_22 
  /////////////
  
  He(11,11) = (-He(2,11)*b(2) - He(6,11)*Sigma(0,2) - He(9,11)*Sigma(1,2) - He(12,11)*Sigma(2,3) - 2*grad(11))/(2*Sigma(2,2));
  He(11,13) = (-He(2,13)*b(2) - He(6,13)*Sigma(0,2) - He(9,13)*Sigma(1,2) - He(12,13)*Sigma(2,3) )/(2*Sigma(2,2));
  He(13,11)= He(11,13); 
  
  
  //////////////
  // sigma_33 
  /////////////
  
  He(13,13) = (-He(3,13)*b(3)-He(7,13)*Sigma(0,3) - He(10,8)*Sigma(1,3) - He(12,13)*Sigma(2,3) - 2*grad(13))/(2*Sigma(3,3)); 
  
  
  // STOPP!
  // return statement
  return He;
  
}


struct TVBS_Vector_three TVBS_pdf_mvna_tvbs_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{
  // determine dimension of the problem
  int m = x_norm.size();
  
  // calculate number of iterations
  int k_tilde = std::max((m-1)/2,1);
  
  // initialize variables for the output
  Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
  p.setZero();
  p_out.setZero();
  g_c.setZero();
  g_rho_rho.setZero();
  p_out.setZero();
  
  // initialize mu_zero variable
  Eigen::VectorXd mu_zero(m);
  mu_zero.setZero();
  
  
  // // // // // // // // // // //
  // reorder the dimensions     //
  // // // // // // // // // // //
  
  struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
  
  Eigen::MatrixXd Cor_mat_temp = output_reorder.mat1;        // we'll work with this one
  Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
  
  Eigen::VectorXd x_temp = output_reorder.vec1;
  
  Eigen::VectorXd temp1 = output_reorder.vec2;
  
  Eigen::VectorXd mu_temp(m);
  mu_temp.setZero();
  
  // // // // // // // // // // //
  // end of reordering          //
  // // // // // // // // // // //
  
  
  
  // calculate initial LDLT decomposition
  struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
  Eigen::MatrixXd L = output_2.mat1;
  Eigen::MatrixXd D = output_2.mat2;
  
  
  // deal with low dimensional cases directly
  if(m==2)
  {
    struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(x_temp, Cor_mat_temp.block(0,0,2,2));
    g_c << output_3.d1, output_3.d2;
    g_rho_rho(0) = output_3.d3;
    p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==3)
  {
    struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_tvn_cpp(x_temp, Cor_mat_temp.block(0,0,3,3));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==4)
  {
    struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_qvn_cpp(x_temp, Cor_mat_temp.block(0,0,4,4));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else 
  {
    // initiate some variables
    struct TVBS_Matrix_two output_2;
    struct TVBS_Matrix_three output_3;
    struct TVBS_Matrix_four output_4;
    Eigen::VectorXd mu_tilde;
    Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
    double p_k1_nomi, p_k1_denomi;
    Eigen::VectorXd g_from_p_mu, g_from_p_cov, g_from_p_c;
    
    
    // calculate P1 using the first four variables
    p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    Eigen::VectorXd g_cu_mul_from_pc(m);
    g_cu_mul_from_pc.setZero();
    
    
    // start here from k1=1 to make indexing easier
    for(int k1=1; k1<k_tilde; k1++)
    {
      // claculate the truncated mean and variance
      double tol = 0.000001;
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
      mu_tilde = output_2.mat1;
      omega = output_2.mat2;
      
      
      // set global variables according to the case (if they have already been truncated or not)
      if(k1==1)
      {
        condcovsigtrunc_global = 0;
        condcovmeantrunc_global = 0;
      }
      else
      {
        condcovsigtrunc_global = 1;
        condcovmeantrunc_global = 1;
      }
      
      
      output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(m-2*k1, m-2*k1), mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_y = output_4.mat1;
      g_mu_mean = output_4.mat2;
      g_x_mean = output_4.mat3;
      g_c_mean = output_4.mat4;
      
      
      output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_mu_cov = output_3.mat1;
      g_x_cov = output_3.mat2;
      g_c_cov = output_3.mat3;
      
      
      // update g_cu_mul_mu_sig and g_cu_mu_lc
      if(k1==1)
      {
        g_cu_mul_mu_sig.resize(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
        g_cu_mul_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        
        g_cu_mu_lc.resize(g_c_mean.rows()+m-2,g_c_mean.cols()+g_c_cov.cols());
        g_cu_mu_lc.setZero();
        g_cu_mu_lc.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
        
      }
      else if(k1>1)
      {
        g_cu_mul_mu_sig_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(), g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mul_mu_sig_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mul_mu_sig_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mul_mu_sig = g_cu_mul_mu_sig*g_cu_mul_mu_sig_rhs;
        
        // update g_cu_mu_lc
        g_cu_mu_lc_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(),g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mu_lc_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mu_lc_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mu_lc = g_cu_mu_lc*g_cu_mu_lc_rhs;
        
        g_cu_mu_lc.block((2*(k1-1)+1)-1,0,2,g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block((2*(k1-1)+1)-1,g_c_mean.cols(),2,g_c_cov.cols()) = g_c_cov;
      }
      
      
      // update pi_kp1
      mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
      
      
      // update the covariance matrix LDLT decomposition
      output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
      L = output_2.mat1;
      D = output_2.mat2;
      Cor_mat_temp = L*D*(L.transpose());
      
      
      if(m>=2*k1+4)
      {
        // regular case
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
        p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        
        // update gradients
        output_3 = TVBS_grad_non_cdf_qvn_by_cdf_bvn_cpp(mu_temp.head(4), Cor_mat_temp.block(0,0,4,4), x_temp.segment((2*k1+1)-1,4));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 4+3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),4) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),4);
        g_rho_rho_add.block(0,4,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,4+3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(3*m-6*k1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,4+3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,(4*m-8*k1-2)-1,g_cu_mul_mu_sig.rows(),1);
        
        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),4)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,4) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),4+3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),4) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),4);
        g_c_add.block(0,4,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,4+3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(3*m-6*k1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,4+3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,(4*m-8*k1-2)-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
      }
      else
      {
        // last iteration, if uneven dimension
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
        p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        
        
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        
        // update gradients
        
        // // // 
        // Errors of e-6 magnitude occur here
        // // //
        output_3 = TVBS_grad_non_cdf_tvn_by_cdf_bvn_cpp(mu_temp.head(3), Cor_mat_temp.block(0,0,3,3), x_temp.segment((2*k1+1)-1,3));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,3*m-6*k1-1,g_cu_mul_mu_sig.rows(),1);
        
        // // // 
        // Errors of e-4 magnitude occur here
        // // //
        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),3)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,3) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,3*m-6*k1-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
        
      }
      
    }
    
    
    
    output_2 = TVBS_grad_cdf_qvn_cpp(x_temp.head(4), Cor_mat_ordered.block(0,0,4,4));
    Eigen::VectorXd g_w = output_2.mat1;
    Eigen::VectorXd g_rho = output_2.mat2;
    
    
    // adding contribution of rho12,rho13, and rho14 from initial cdfqvn function
    g_rho_rho.head(3) += g_rho.head(3)/p(0);
    
    
    // adding contribution of rho23 & rho24 from initial cdfqvn function
    g_rho_rho.segment(m-1,2) += g_rho.segment(3,2)/p(0);
    
    
    // adding contribution of rho34 from initial cdfqvn function
    g_rho_rho(2*m-2-1) += g_rho(5)/p(0);
    
    
    // adding contribution of all abscissa originating from the probability function
    g_c += g_cu_mul_from_pc;
    
    
    // inserting gradient contribution of first four abscissae directly from the cdfqvn function
    g_c.head(4) += g_w/p(0);
    
    
    if(log_out == 0)
    {
      p_out(0) = p.prod();
      g_c = (p_out(0)*g_c);
      g_rho_rho = (p_out(0)*g_rho_rho);
    }
    else
    {
      p_out(0) = p.array().log().sum();
    }
  }
  
  
  
  
  // // // // // // // // // // // // // // //
  // invert the reordering                  //
  // // // // // // // // // // // // // // //
  
  // g_rho_rho has to be transformed into matrix form, reordered and then back into vectorization
  
  Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
  Eigen::MatrixXd g_rho_rho_mat_new = g_rho_rho_mat + g_rho_rho_mat.transpose();
  
  Eigen::VectorXd g_c_new(m);
  g_c_new.setZero();
  
  for(int i_row=0; i_row<m; i_row++)
  {
    g_c_new[temp1[i_row]] = g_c(i_row);
    for(int i_col=0; i_col<m; i_col++)
    {
      g_rho_rho_mat((int)temp1[i_row],(int)temp1[i_col]) = g_rho_rho_mat_new(i_row,i_col);
    }
  }
  g_c = g_c_new;
  g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat,0);
  
  
  // // // // // // // // // // // // // // //
  // end of reorder invertion               //
  // // // // // // // // // // // // // // //
  
  struct TVBS_Vector_three output_final = {p_out, g_c, g_rho_rho};
  return output_final;
}


struct TVBS_Vector_four TVBS_pdf_mvn_analytic_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x, int log_out = 0)
{
  // initialize output variables
  Eigen::VectorXd p, g_mu, g_cov, g_x; 
  
  // calculate normalized x vector
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  
  
  // if x values are further away from mu than 6 sd => P is in {0,1}
  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  
  // calculate correlation matrix from covariance matrix
  Eigen::MatrixXd Sigma_diag_sqrt_inv(Sigma_diag_sqrt.size(),Sigma_diag_sqrt.size());
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  
  // adjust global variable
  if(optimal_global==1)
  {
    optimal_global = 0;
  }
  
  
  // get probability and gradient of the normalized problem
  struct TVBS_Vector_three output_3 = TVBS_pdf_mvna_tvbs_cpp(x_norm, Cor_mat, log_out);
  p = output_3.v1;
  Eigen::VectorXd g_w = output_3.v2;
  Eigen::VectorXd g_rho = output_3.v3;
  
  
  
  // calculate the derivatives of the COVARIANCE matrix from the CORRELATION matrix
  if(covarr_global==1)
  {
    struct TVBS_Matrix_two output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
    Eigen::MatrixXd g_b_cor_cov = output_2.mat1;
    Eigen::MatrixXd g_omega_cor_cov = output_2.mat2;
    
    g_cov = g_b_cor_cov*g_w + g_omega_cor_cov*g_rho;
  }
  else
  {
    g_cov = g_rho;
  }
  
  
  // resize the gradient of the mean and x vector
  g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  g_x = -1.0*g_mu;
  
  
  struct TVBS_Vector_four output_final = {p, g_mu, g_cov, g_x};
  return output_final;
}


// // // // // // // // //
//  TVBS + grad norm    //
// // // // // // // // //

//' TVBS approximation to multivariate Gaussian CDF. 
//' @description
//' The function computes the  CDF for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS). 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @param log_out
//' integer; probability or log of probability? 
//' @return 
//' double; (log of) probability. 
//' @keywords internal
//'
//[[Rcpp::export]]
Eigen::VectorXd TVBS_v2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{
   // determine dimension of the problem
   int m = x_norm.size();
   
   // calculate number of iterations
   int k_tilde = std::max((m-1)/2,1);
   
   // initialize variables for the output
   Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
   p.setZero();
   p_out.setZero();
   g_c.setZero();
   g_rho_rho.setZero();
   p_out.setZero();
   
   // initialize mu_zero variable
   Eigen::VectorXd mu_zero(m);
   mu_zero.setZero();
   
   // // // // // // // // // // //
   // reorder the dimensions     //
   // // // // // // // // // // //
   
   // take out reordering. 
   struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
   
   Eigen::MatrixXd Cor_mat_temp = Cor_mat; // output_reorder.mat1;        // we'll work with this one
   Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
   
   Eigen::VectorXd x_temp = x_norm; // output_reorder.vec1;
   
   int kb = x_norm.size();
   Eigen::VectorXd temp1 = Eigen::VectorXd::LinSpaced(kb,0,kb-1); // output_reorder.vec2;
   
   Eigen::VectorXd mu_temp(m);
   mu_temp.setZero();
   
   // // // // // // // // // // //
   // end of reordering          //
   // // // // // // // // // // //
   
   
   // calculate initial LDLT decomposition
   struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
   Eigen::MatrixXd L = output_2.mat1;
   Eigen::MatrixXd D = output_2.mat2;
   
   
   // deal with low dimensional cases directly
   if(m==2){
     struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(x_temp, Cor_mat_temp.block(0,0,2,2));
     g_c << output_3.d1, output_3.d2;
     g_rho_rho(0) = output_3.d3;
     
     p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
     
     // if the output should be for log
     if(log_out==1)
     {
       g_c /= p_out(0);
       g_rho_rho /= p_out(0);
       double p = log(p_out(0));
       p_out(0) = p;
     }
   }
   if(m==3){
     struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_tvn_cpp(x_temp, Cor_mat_temp.block(0,0,3,3));
     g_c = output_2.mat1;
     g_rho_rho = output_2.mat2;
     p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
     
     // if the output should be for log
     if(log_out==1)
     {
       g_c /= p_out(0);
       g_rho_rho /= p_out(0);
       double p = log(p_out(0));
       p_out(0) = p;
     }
   }
   if(m==4){
     
     struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_qvn_cpp(x_temp, Cor_mat_temp.block(0,0,4,4));
     g_c = output_2.mat1;
     g_rho_rho = output_2.mat2;
     
     
     p_out(0) = TVBS_pmvnorm_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     // if the output should be for log
     if(log_out==1)
     {
       g_c /= p_out(0);
       g_rho_rho /= p_out(0);
       double p = log(p_out(0));
       p_out(0) = p;
     }
   }
   if (m>4){
     
     // initiate some variables
     struct TVBS_Matrix_two output_2;
     struct TVBS_Matrix_three output_3;
     struct TVBS_Matrix_four output_4;
     Eigen::VectorXd mu_tilde;
     Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
     double p_k1_nomi, p_k1_denomi;
     Eigen::VectorXd g_from_p_mu, g_from_p_cov, g_from_p_c;
     
     
     // calculate P1 using the first four variables
     p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     
     Eigen::VectorXd g_cu_mul_from_pc(m);
     g_cu_mul_from_pc.setZero();
     
     
     // start here from k1=1 to make indexing easier
     for(int k1=1; k1<k_tilde; k1++)
     {
       // claculate the truncated mean and variance// [[Rcpp::export]]
       double tol = 0.000001;
       output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
       mu_tilde = output_2.mat1;
       omega = output_2.mat2;
       
       // set global variables according to the case (if they have already been truncated or not)
       if(k1==1)
       {
         condcovsigtrunc_global = 0;
         condcovmeantrunc_global = 0;
       }
       else
       {
         condcovsigtrunc_global = 1;
         condcovmeantrunc_global = 1;
       }
       
       
       output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(m-2*k1, m-2*k1), mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
       g_y = output_4.mat1;
       g_mu_mean = output_4.mat2;
       g_x_mean = output_4.mat3;
       g_c_mean = output_4.mat4;
       
       
       output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
       g_mu_cov = output_3.mat1;
       g_x_cov = output_3.mat2;
       g_c_cov = output_3.mat3;
       
       
       // update g_cu_mul_mu_sig and g_cu_mu_lc
       if(k1==1)
       {
         g_cu_mul_mu_sig.resize(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
         g_cu_mul_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
         g_cu_mul_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
         
         g_cu_mu_lc.resize(g_c_mean.rows()+m-2,g_c_mean.cols()+g_c_cov.cols());
         g_cu_mu_lc.setZero();
         g_cu_mu_lc.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
         g_cu_mu_lc.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
         
       }
       else if(k1>1)
       {
         g_cu_mul_mu_sig_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(), g_mu_mean.cols()+g_mu_cov.cols());
         g_cu_mul_mu_sig_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
         g_cu_mul_mu_sig_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
         g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
         g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
         g_cu_mul_mu_sig = g_cu_mul_mu_sig*g_cu_mul_mu_sig_rhs;
         
         // update g_cu_mu_lc
         g_cu_mu_lc_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(),g_mu_mean.cols()+g_mu_cov.cols());
         g_cu_mu_lc_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
         g_cu_mu_lc_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
         g_cu_mu_lc_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
         g_cu_mu_lc_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
         g_cu_mu_lc = g_cu_mu_lc*g_cu_mu_lc_rhs;
         
         g_cu_mu_lc.block((2*(k1-1)+1)-1,0,2,g_c_mean.cols()) = g_c_mean;
         g_cu_mu_lc.block((2*(k1-1)+1)-1,g_c_mean.cols(),2,g_c_cov.cols()) = g_c_cov;
       }
       
       
       // update pi_kp1
       mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
       
       
       // update the covariance matrix LDLT decomposition
       output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
       L = output_2.mat1;
       D = output_2.mat2;
       Cor_mat_temp = L*D*(L.transpose());
       
       if(m>=2*k1+4)
       {
         // regular case
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
         p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
         p(k1) = p_k1_nomi/p_k1_denomi;
         
         
         // update gradients
         output_3 = TVBS_grad_non_cdf_qvn_by_cdf_bvn_cpp(mu_temp.head(4), Cor_mat_temp.block(0,0,4,4), x_temp.segment((2*k1+1)-1,4));
         g_from_p_mu = output_3.mat1;
         g_from_p_cov = output_3.mat2;
         g_from_p_c = output_3.mat3;
         
         
         // update g_rho_rho
         g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 4+3+2+1);
         g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),4) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),4);
         g_rho_rho_add.block(0,4,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
         g_rho_rho_add.block(0,4+3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(3*m-6*k1)-1,g_cu_mul_mu_sig.rows(),2);
         g_rho_rho_add.block(0,4+3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,(4*m-8*k1-2)-1,g_cu_mul_mu_sig.rows(),1);
         
         g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),4)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
         
         
         // update g_cu_mul_from_pc
         g_cu_mul_from_pc.segment((2*k1+1)-1,4) += g_from_p_c/p(k1);
         
         
         // update g_c
         g_c_add.resize(g_cu_mu_lc.rows(),4+3+2+1);
         g_c_add.block(0,0,g_cu_mu_lc.rows(),4) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),4);
         g_c_add.block(0,4,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),3);
         g_c_add.block(0,4+3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(3*m-6*k1)-1,g_cu_mu_lc.rows(),2);
         g_c_add.block(0,4+3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,(4*m-8*k1-2)-1,g_cu_mu_lc.rows(),1);
         
         g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
         
       }
       else
       {
         // last iteration, if uneven dimension
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
         p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
         
         p(k1) = p_k1_nomi/p_k1_denomi;
         
         // update gradients
         
         // // // 
         // Errors of e-6 magnitude occur here
         // // //
         
         output_3 = TVBS_grad_non_cdf_tvn_by_cdf_bvn_cpp(mu_temp.head(3), Cor_mat_temp.block(0,0,3,3), x_temp.segment((2*k1+1)-1,3));
         g_from_p_mu = output_3.mat1;
         g_from_p_cov = output_3.mat2;
         g_from_p_c = output_3.mat3;
         
         // update g_rho_rho
         g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 3+2+1);
         g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
         g_rho_rho_add.block(0,3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),2);
         g_rho_rho_add.block(0,3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,3*m-6*k1-1,g_cu_mul_mu_sig.rows(),1);
         
         // // // 
         // Errors of e-4 magnitude occur here
         // // //
         g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),3)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
         
         // update g_cu_mul_from_pc
         g_cu_mul_from_pc.segment((2*k1+1)-1,3) += g_from_p_c/p(k1);
         
         
         // update g_c
         g_c_add.resize(g_cu_mu_lc.rows(),3+2+1);
         g_c_add.block(0,0,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),3);
         g_c_add.block(0,3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),2);
         g_c_add.block(0,3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,3*m-6*k1-1,g_cu_mu_lc.rows(),1);
         
         g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
       }
       
     }
     
     
     
     output_2 = TVBS_grad_cdf_qvn_cpp(x_temp.head(4), Cor_mat_ordered.block(0,0,4,4));
     Eigen::VectorXd g_w = output_2.mat1;
     Eigen::VectorXd g_rho = output_2.mat2;
     
     
     // adding contribution of rho12,rho13, and rho14 from initial cdfqvn function
     g_rho_rho.head(3) += g_rho.head(3)/p(0);
     
     
     // adding contribution of rho23 & rho24 from initial cdfqvn function
     g_rho_rho.segment(m-1,2) += g_rho.segment(3,2)/p(0);
     
     
     // adding contribution of rho34 from initial cdfqvn function
     g_rho_rho(2*m-2-1) += g_rho(5)/p(0);
     
     
     // adding contribution of all abscissa originating from the probability function
     g_c += g_cu_mul_from_pc;
     
     
     // inserting gradient contribution of first four abscissae directly from the cdfqvn function
     g_c.head(4) += g_w/p(0);
     
     
     if(log_out == 0)
     {
       p_out(0) = p.prod();
       g_c = (p_out(0)*g_c);
       g_rho_rho = (p_out(0)*g_rho_rho);
     }
     else
     {
       p_out(0) = p.array().log().sum();
     }
     
   }
   
   
   
   
   // // // // // // // // // // // // // // //
   // invert the reordering                  //
   // // // // // // // // // // // // // // //
   
   // g_rho_rho has to be transformed into matrix form, reordered and then back into vectorization
   
   Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
   Eigen::MatrixXd g_rho_rho_mat_new = g_rho_rho_mat + g_rho_rho_mat.transpose();
   
   Eigen::VectorXd g_c_new(m);
   g_c_new.setZero();
   
   for(int i_row=0; i_row<m; i_row++)
   {
     g_c_new[temp1[i_row]] = g_c(i_row);
     for(int i_col=0; i_col<m; i_col++)
     {
       g_rho_rho_mat((int)temp1[i_row],(int)temp1[i_col]) = g_rho_rho_mat_new(i_row,i_col);
     }
   }
   g_c = g_c_new;
   g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat,0);
   
   
   // // // // // // // // // // // // // // //
   // end of reorder invertion               //
   // // // // // // // // // // // // // // //
   
   
   
   Eigen::VectorXd output_final(1+g_c.size()+g_rho_rho.size()); 
   output_final << g_c, g_rho_rho, p_out(0);
   
   
   return output_final;
 }



// // // // // // // // //
//  TVBS - prob only    //
// // // // // // // // //

//' TVBS approximation to multivariate Gaussian CDF. 
//' @description
//' The function computes the  CDF for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS). 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @param log_out
//' integer; probability or log of probability? 
//' @return 
//' double; (log of) probability. 
//' @keywords internal
//'
// [[Rcpp::export]]
double TVBS_pv2(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{

   // determine dimension of the problem
   int m = x_norm.size();
   
   // calculate number of iterations
   int k_tilde = std::max((m-1)/2,1);
   
   
   // initialize variables for the output
   Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
   p.setZero();
   p_out.setZero();
   g_c.setZero();
   g_rho_rho.setZero();
   p_out.setZero();
   
   // initialize mu_zero variable
   Eigen::VectorXd mu_zero(m);
   mu_zero.setZero();
   
   
   
   
   // // // // // // // // // // //
   // reorder the dimensions     //
   // // // // // // // // // // //
   
   struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
   
   Eigen::MatrixXd Cor_mat_temp = Cor_mat;        // we'll work with this one
   Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
   
   Eigen::VectorXd x_temp = x_norm;
   
   int kb = x_norm.size();
   Eigen::VectorXd temp1 = Eigen::VectorXd::LinSpaced(kb,0,kb-1); // output_reorder.vec2;
   
   Eigen::VectorXd mu_temp(m);
   mu_temp.setZero();
   
   // // // // // // // // // // //
   // end of reordering          //
   // // // // // // // // // // //
   
   
   
   
   // calculate initial LDLT decomposition
   struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
   Eigen::MatrixXd L = output_2.mat1;
   Eigen::MatrixXd D = output_2.mat2;
   
   
   
   // deal with low dimensional cases directly
   if(m==2){
     p_out(0) = TVBS_pmvnorm_cpp(x_temp, x_temp*0, Cor_mat_temp.block(0,0,2,2));
     
     // if the output should be for log
     if(log_out==1)
     {
       double p_temp = log(p_out(0));
       p_out(0) = p_temp;
     }
   }
   if(m==3){
     p_out(0) = TVBS_pmvnorm_cpp(x_temp, x_temp*0, Cor_mat_temp.block(0,0,3,3));
     
     // if the output should be for log
     if(log_out==1)
     {
       double p_temp = log(p_out(0));
       p_out(0) = p_temp;
     }
   }
   if(m==4){
     
     p_out(0) = TVBS_pmvnorm_cpp(x_temp,  x_temp*0, Cor_mat_temp.block(0,0,4,4));
     
     // if the output should be for log
     if(log_out==1)
     {
       double p_temp = log(p_out(0));
       p_out(0) = p_temp;
     }
   }
   if (m>4) {
     // initiate some variables
     struct TVBS_Matrix_two output_2;
     struct TVBS_Matrix_three output_3;
     struct TVBS_Matrix_four output_4;
     Eigen::VectorXd mu_tilde;
     Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
     double p_k1_nomi, p_k1_denomi;
     
     
     // calculate P1 using the first four variables
     p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     
     // start here from k1=1 to make indexing easier
     for(int k1=1; k1<k_tilde; k1++)
     {
       // calculate the truncated mean and variance
       double tol = 0.000001;
       output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
       mu_tilde = output_2.mat1;
       omega = output_2.mat2;
       
       // update pi_kp1
       mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
       
       //Rcout << mu_temp << std::endl; 
       
       // update the covariance matrix LDLT decomposition
       output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
       L = output_2.mat1;
       D = output_2.mat2;
       Cor_mat_temp = L*D*(L.transpose());
       
       
       if(m>=2*k1+4)
       {
         // regular case
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
         //p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         //p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
         p(k1) = p_k1_nomi/p_k1_denomi;
         
       }
       else
       {
         // last iteration, if uneven dimension
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         p(k1) = p_k1_nomi/p_k1_denomi;
         
       }
       
     }
     
     
     if(log_out == 0)
     {
       p_out(0) = p.prod();
     }
     else
     {
       p_out(0) = p.array().log().sum();
     }
     
     
   }
   
   double p_final = p_out(0);
   return p_final;
 }

// // // // // // // // // // //   
// TVBS gradient calculations //
// // // // // // // // // // // 

//' TVBS approximation to multivariate Gaussian CDF: calculates CDF and gradient. 
//' @description
//' The function computes the  CDF and the corresponding gradient for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS). 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @param log_out
//' integer; probability or log of probability? 
//' @return 
//' vector; gradient of log probability, double; (log of) probability.
//' @keywords internal 
//'  
// [[Rcpp::export]]
Eigen::VectorXd TVBS_grad(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{

   // determine dimension of the problem
   int m = x_norm.size();
   
   // calculate number of iterations
   int k_tilde = std::max((m-1)/2,1);
   
   int npar = m + m*(m-1)/2; // m entries in b, the rest in upper triangluar of Sigma. 
   Eigen::MatrixXd dp(k_tilde,npar),dpd(k_tilde,npar),dpn(k_tilde,npar);
   dp.setZero();
   
   // initialize variables for the output
   Eigen::VectorXd p(k_tilde), p_out(1), grad(npar);
   p.setZero();
   p_out.setZero();
   grad.setZero();
   
   Eigen::MatrixXd Cor_mat_temp = Cor_mat;        // we'll work with this one
   Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
   
   Eigen::VectorXd x_temp = x_norm;
   
   Eigen::VectorXd temp1 = x_norm;
   
   // initialize mu_zero variable
   Eigen::VectorXd mu_zero(m); 
   Eigen::MatrixXd mu_temp(m,k_tilde);
   mu_zero.setZero();
   mu_temp.setZero();
   
   Eigen::MatrixXd dpmu(3,npar);
   
   // calculate initial LDLT decomposition
   
   struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
   Eigen::MatrixXd L(m,m*k_tilde);
   L.setZero();
   L.block(0,0,m,m)= output_2.mat1;
   Eigen::MatrixXd D(m,m*k_tilde);
   D.setZero();
   D.block(0,0,m,m)= output_2.mat2;
   
   
   Eigen::VectorXd numrows(k_tilde);
   numrows.setZero();
   numrows(0)=m; 
   
   // deal with low dimensional cases directly
   if(m==2){
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2)));
     Eigen::VectorXd gr = grad_gen_cdf(x_temp - mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
     grad(0)=gr(0);
     grad(1)=gr(1);
     grad(2)=gr(4);
     grad = grad/exp(p_out(0));
     dp = grad;
   }
   if(m==3){
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3)));
     Eigen::VectorXd gr = grad_gen_tvn_cdf(x_temp - mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
     
     grad.head(3) = gr.head(3); 
     grad(3) = gr(4);
     grad(4)= gr(5);
     grad(5) = gr(7); 
     grad = grad/exp(p_out(0));
     dp = grad;
   }
   if(m==4){
     
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, x_temp*0, Cor_mat_temp.block(0,0,4,4)));
     Eigen::VectorXd gr = TVBS_grad_gen_qvn_cdf(x_temp, Cor_mat_temp.block(0,0,4,4));
     
     grad.head(4) = gr.head(4); 
     grad(4) = gr(5);
     grad(5)= gr(6);
     grad(6) = gr(7);
     grad(7) = gr(9);
     grad(8) = gr(10);
     grad(9) = gr(12);
     grad = grad/exp(p_out(0));
     dp = grad;
   }
   if (m>4){
     // initiate some variables
     struct TVBS_Matrix_two output_2;
     struct TVBS_Matrix_three output_3;
     struct TVBS_Matrix_four output_4;
     Eigen::MatrixXd mu_tilde(m,k_tilde);
     mu_tilde.setZero();
     Eigen::VectorXd out_v(5);
     
     Eigen::MatrixXd Cor_mat_temp(m,m*k_tilde);
     Cor_mat_temp.setZero();
     
     Eigen::MatrixXd Omega(2,2*k_tilde), dOmega(2,2), dL(m,m), dD(m,m), dCor_mat(m,m);
     Eigen::VectorXd dmu_temp(m), dmu_tilde(m), dx_temp(m); 
     
     Eigen::MatrixXd doutput_2, doutput_2_all(5*k_tilde,7);
     doutput_2_all.setZero(); 
     doutput_2.setZero();
     
     Eigen::VectorXd gr(14);
     gr.setZero(); 
     
     Eigen::MatrixXd gr_all(k_tilde,14);
     gr_all.setZero();
     
     Eigen::VectorXd grd(5);
     grd.setZero(); 
     
     Eigen::MatrixXd grd_all(k_tilde,5);
     grd_all.setZero();
     
     
     struct TVBS_Matrix_three doutput_3;
     struct TVBS_Matrix_four doutput_4;
     
     Omega.setZero(); 
     
     
     Eigen:VectorXd p_k1_nomi(k_tilde), p_k1_denomi(k_tilde);
     
     //////////////////////////////////////////////////////////////////////////////////
     //  gradient calculations: iterate over the elements                          ///
     //  this allows to reuse matrices of the correct size and makes coding easier. //
     /////////////////////////////////////////////////////////////////////////////////
     
     
     
     dp.setZero();
     dpn.setZero();
     dpd.setZero();
     
     
     // find the indices. 
     // ind_mat has three columns: first column: index in b, if applicable, -1 else; second column: row index; third column col index.
     Eigen::MatrixXd ind_mat (npar,3);
     ind_mat.setZero();
     
     int cur = -1;
     for (int nm=0; nm< m; nm++){
       cur+=1;
       ind_mat(nm,0)=nm;
     }
     for (int nm=0; nm< m; nm++){
       for (int nk=0;nk<m-nm-1;nk++){
         cur+=1;
         ind_mat(cur,0)=-1;
         ind_mat(cur,1)=nm;
         ind_mat(cur,2)=nm+nk+1;
       }
     }
     
     /////////////////////////////////////////////////////////////////
     // first iterate over k1 to perform all calculations that need to be done only once. 
     /////////////////////////////////////////////////////////////////
     // calculate P1 using the first four variables
     
     Cor_mat_temp.block(0,0,m,m) = Cor_mat;     
     
     p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     gr_all.row(0) = TVBS_grad_gen_qvn_cdf(x_temp.head(4) - mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     // start here from k1=1 to make indexing easier
     for(int k1=1; k1<k_tilde; k1++)
     {
       // calculate the truncated mean and variance
       double tol = 0.000001;
       struct TVBS_Matrix_two output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.block(0,k1-1,2,1), D.block(0,(k1-1)*m,2,2), x_temp.segment((2*k1-1)-1,2),tol);
       mu_tilde.block(0,k1,2,1) = output_2.mat1;
       Omega.block(0,k1*2,2,2) = output_2.mat2;
       
       doutput_2_all.block(5*(k1-1),0,5,7) = TVBS_Jacobian_biv_gen_trunc(x_temp.segment((2*k1-1)-1,2), mu_temp.block(0,k1-1,2,1), D.block(0,(k1-1)*m,2,2));
       
       
       // update pi_kp1
       mu_temp.block(0,k1,m-2*k1,1) = mu_temp.block(2,k1-1,numrows(k1-1),1) + L.block(2,(k1-1)*m,numrows(k1-1),2)*(mu_tilde.block(0,k1,2,1)-mu_temp.block(0,k1-1,2,1));
       
       // update the covariance matrix LDLT decomposition
       struct TVBS_Matrix_two output_u = TVBS_LDLT_update_cpp(L.block(0,(k1-1)*m,numrows(k1-1),numrows(k1-1)),D.block(0,(k1-1)*m,numrows(k1-1),numrows(k1-1)),Omega.block(0,k1*2,2,2),2);
       Eigen::MatrixXd newL = output_u.mat1; 
       Eigen::MatrixXd newD = output_u.mat2; 
       numrows(k1)= newL.cols(); 
       L.block(0,k1*m,numrows(k1),numrows(k1)) = newL;
       D.block(0,k1*m,numrows(k1),numrows(k1)) = newD;
       Cor_mat_temp.block(0,k1*m,numrows(k1),numrows(k1)) = newL*newD*(newL.transpose());
       
       if(m>=2*k1+4)
       {
         // regular case
         
         // calculate next probability
         p_k1_nomi(k1) = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.block(0,k1,4,1), Cor_mat_temp.block(0,k1*m,4,4));
         gr_all.row(k1) = TVBS_grad_gen_qvn_cdf(x_temp.segment((2*k1+1)-1,4)- mu_temp.block(0,k1,4,1), Cor_mat_temp.block(0,m*k1,4,4));
         
         p_k1_denomi(k1) = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.block(0,k1,2,1), Cor_mat_temp.block(0,k1*m,2,2));
         grd_all.row(k1) = grad_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.block(0,k1,2,1), Cor_mat_temp.block(0,m*k1,2,2));
         
         p(k1) = p_k1_nomi(k1)/p_k1_denomi(k1);
         
       }
       else
       {
         // last iteration, if uneven dimension
         
         // calculate next probability
         p_k1_nomi(k1) = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.block(0,k1,3,1), Cor_mat_temp.block(0,k1*m,3,3));
         Eigen::VectorXd grn(9);
         grn.setZero();
         grn = grad_gen_tvn_cdf(x_temp.segment((2*k1+1)-1,3)- mu_temp.block(0,k1,3,1), Cor_mat_temp.block(0,m*k1,3,3)); 
         
         for (int jg=0;jg<9;jg++){
           gr_all(k1,jg)= grn(jg);
         }
         
         p_k1_denomi(k1) = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.block(0,k1,2,1), Cor_mat_temp.block(0,k1*m,2,2));
         grd_all.row(k1) = grad_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.block(0,k1,2,1), Cor_mat_temp.block(0,m*k1,2,2));
         
         p(k1) = p_k1_nomi(k1)/p_k1_denomi(k1);
         
       }
       
     }
     
     
     // npar 
     for (int np=0; np< npar; np++){
       // initialize the matrices of partial derivatives 
       dOmega.setZero();
       dL.setZero();
       dD.setZero();
       dCor_mat = Cor_mat; 
       dCor_mat.setZero();
       dmu_tilde.setZero();
       dx_temp.setZero(); 
       dmu_temp = dx_temp;
       
       if (ind_mat(np,0)>-1){
         dx_temp((int)ind_mat(np,0))=1; // parameter corresponds to b entry.
         dL = L.block(0,0,m,m)*0;
         dD = D.block(0,0,m,m)*0;
       } else { // parameter corresponds to entry in covariance matrix.
         int r = ind_mat(np,1);
         int c = ind_mat(np,2);
         dCor_mat(r,c)=1;
         dCor_mat(c,r)=1;
         
         Eigen::VectorXd ind_1(2);
         ind_1(0)=r; 
         ind_1(1)=c; 
         
         // derivative of initial LDLT
         struct TVBS_Matrix_two  output_dLDLT  = TVBS_grad_LDLT_decomp(Cor_mat_temp.block(0,0,m,m), ind_1);
         dL = output_dLDLT.mat1;
         dD = output_dLDLT.mat2;
       }
       
       
       // first four coordinates
       
       // calculate P1 using the first four variables
       //p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
       // Eigen::VectorXd gr = TVBS_grad_gen_qvn_cdf(x_temp.head(4) - mu_temp.block(0,0,4,1), Cor_mat_temp.block(0,0,4,4));
       gr = gr_all.row(0);
       
       Eigen::VectorXd dxmmu(4);
       dxmmu = dx_temp.head(4) - dmu_temp.head(4);
       
       int cur=0;
       for (int nm=0;nm<4;nm++){
         dp(0,np) += gr(nm) * dxmmu(nm);
         cur += 1;
       }
       for (int nm=0;nm<4;nm++){
         for (int nk=nm;nk<4;nk++){
           dp(0,np) += gr(cur) * dCor_mat(nm,nk);
           cur +=1;
         }
       }
       
       
       
       // start here from k1=1 to make indexing easier
       for(int k1=1; k1<k_tilde; k1++)
       {
         // calc gradient. 
         doutput_2 =   doutput_2_all.block(5*(k1-1),0,5,7);
         Eigen::VectorXd dmu_tilde(2);
         
         dmu_tilde(0) = doutput_2(0,0)*dx_temp((2*k1-1)-1) + doutput_2(0,1)*dx_temp((2*k1-1)) + doutput_2(0,2)*dmu_temp(0)+doutput_2(0,3)*dmu_temp(1);
         dmu_tilde(0) += doutput_2(0,4)* dD(0,0) + doutput_2(0,5)*dD(0,1) + doutput_2(0,6) * dD(1,1);
         dmu_tilde(1) = doutput_2(1,0)*dx_temp((2*k1-1)-1) + doutput_2(1,1)*dx_temp((2*k1-1)) + doutput_2(1,2)*dmu_temp(0)+doutput_2(1,3)*dmu_temp(1);
         dmu_tilde(1) += doutput_2(1,4)* dD(0,0) + doutput_2(1,5)*dD(0,1) + doutput_2(1,6) * dD(1,1);
         
         
         // Attention! doutput_2 ordered strangely as b0,b1,om00,om11,om01. 
         dOmega(0,0) = doutput_2(2,4)* dD(0,0) + doutput_2(2,5)*dD(0,1) + doutput_2(2,6) * dD(1,1);
         dOmega(0,0) += doutput_2(2,0)*dx_temp((2*k1-1)-1) + doutput_2(2,1)*dx_temp((2*k1-1)) + doutput_2(2,2)*dmu_temp(0)+doutput_2(2,3)*dmu_temp(1);
         dOmega(0,1) = doutput_2(3,4)* dD(0,0) + doutput_2(3,5)*dD(0,1) + doutput_2(3,6) * dD(1,1);
         dOmega(0,1) += doutput_2(3,0)*dx_temp((2*k1-1)-1) + doutput_2(3,1)*dx_temp((2*k1-1)) + doutput_2(3,2)*dmu_temp(0)+doutput_2(3,3)*dmu_temp(1);
         dOmega(1,1) = doutput_2(4,4)* dD(0,0) + doutput_2(4,5)*dD(0,1) + doutput_2(4,6) * dD(1,1);
         dOmega(1,1) += doutput_2(4,0)*dx_temp((2*k1-1)-1) + doutput_2(4,1)*dx_temp((2*k1-1)) + doutput_2(4,2)*dmu_temp(0)+doutput_2(4,3)*dmu_temp(1);
         dOmega(1,0) = dOmega(0,1);
         
         // update pi_kp1
         
         // gradient
         dmu_temp = dmu_temp.segment(2,m-2*k1) + dL.block(2,0,m-2*k1,2)*(mu_tilde.block(0,k1,2,1)-mu_temp.block(0,k1-1,2,1)) + L.block(2,(k1-1)*m,m-2*k1,2)*(dmu_tilde-dmu_temp.head(2));
         // update the covariance matrix LDLT decomposition
         struct TVBS_Matrix_two dout_2 = TVBS_LDLT_update_grad(L.block(0,(k1-1)*m,numrows(k1-1),numrows(k1-1)), D.block(0,(k1-1)*m,numrows(k1-1),numrows(k1-1)), Omega.block(0,(k1)*2,2,2), dL,  dD, dOmega, 2);
         
         // calc gradient. 
         
         Eigen::MatrixXd newL = L.block(0,(k1)*m,numrows(k1),numrows(k1));
         Eigen::MatrixXd newD = D.block(0,(k1)*m,numrows(k1),numrows(k1));
         dL = dout_2.mat1;
         dD = dout_2.mat2;
         
         dCor_mat.block(0,0,numrows(k1),numrows(k1)) = dL * newD*(newL.transpose()) + newL * dD * (newL.transpose()) + newL * newD *(dL.transpose());
         
         if(m>=2*k1+4)
         {
           // gradient 
           gr = gr_all.row(k1); // 
           Eigen::VectorXd dxmdmu = dx_temp.segment((2*k1+1)-1,4)- dmu_temp.head(4);
           dpn(k1,np) = gr(0)*dxmdmu(0) + gr(1)*dxmdmu(1) + gr(2)*dxmdmu(2) + gr(3)*dxmdmu(3);
           int cur = 3; 
           for (int nr=0;nr<4;nr++){
             for (int nc=nr;nc<4;nc++){
               cur += 1;
               dpn(k1,np) += gr(cur) * dCor_mat(nr,nc); 
             }
           } 
           // gradient
           grd = grd_all.row(k1); // 
           Eigen::VectorXd dxmdmud = dx_temp.segment((2*k1+1)-1,2)- dmu_temp.head(2);
           dpd(k1,np) = grd(0) * dxmdmud(0) +  grd(1) * dxmdmud(1);
           dpd(k1,np) += grd(2) * dCor_mat(0,0) + grd(3) * dCor_mat(1,1) + grd(4) * dCor_mat(0,1); 
           
           // quotient. 
           // gradient: quotient rule. 
           dp(k1,np) = (dpn(k1,np) * p_k1_denomi(k1) - p_k1_nomi(k1) *dpd(k1,np))/(p_k1_denomi(k1)* p_k1_denomi(k1)); 
         }
         else
         {
           // last iteration, if uneven dimension
           int ii = (2*k1+1)-1;
           
           Eigen::VectorXd gr(9);
           gr.setZero();
           gr = gr_all.row(k1);  
           Eigen::VectorXd dxmdmu = dx_temp.segment((2*k1+1)-1,3)- dmu_temp.head(3);
           
           dpmu(0,np) = dxmdmu(0);
           dpmu(1,np) = dxmdmu(1);
           dpmu(2,np) = dxmdmu(2);
           
           
           dpn(k1,np) = gr(0)*dxmdmu(0) + gr(1)*dxmdmu(1) + gr(2)*dxmdmu(2);
           int cur = 2; 
           for (int nr=0;nr<3;nr++){
             for (int nc=nr;nc<3;nc++){
               cur += 1;
               dpn(k1,np) += gr(cur) * dCor_mat(nr,nc); 
             }
           } 
           
           
           grd = grd_all.row(k1);  
           Eigen::VectorXd dxmdmud = dx_temp.segment((2*k1+1)-1,2)- dmu_temp.head(2);
           dpd(k1,np) = grd(0) * dxmdmud(0) +  grd(1) * dxmdmud(1);
           dpd(k1,np) += grd(2) * dCor_mat(0,0) + grd(3) * dCor_mat(1,1) + grd(4) * dCor_mat(0,1); 
           dp(k1,np) = (dpn(k1,np)* p_k1_denomi(k1) - p_k1_nomi(k1) *dpd(k1,np))/(p_k1_denomi(k1) * p_k1_denomi(k1)); 
           
         }
       }
       
       if(log_out == 0)
       {
         p_out(0) = p.prod();
       }
       else
       {
         p_out(0) = p.array().log().sum();
       }
       
       // gradient 
       Eigen::VectorXd ipv(k_tilde);
       for (int jj=0;jj<k_tilde;jj++){
         ipv(jj) = 1/p(jj);
       }
       
       grad = ipv.transpose() * dp;
     } 
   }
   
   Eigen::VectorXd output_final(1+grad.size()); 
   output_final << grad, p_out(0);
   
   return output_final;
 }


//////////////////////////////////////////
/// TVBS Hessian matrix calculation //////
//////////////////////////////////////////

//' TVBS approximation to multivariate Gaussian CDF: calculates Hessian of log CDF 
//' @description
//' The function computes the Hessian of the log of the CDF and the corresponding gradient for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS). 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @return 
//' matrix; Hessian of log probability. 
//' @keywords internal 
//'
// [[Rcpp::export]]
Eigen::MatrixXd TVBS_hess_new(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat)
{
   // determine dimension of the problem
   int m = x_norm.size();
   
   // calculate number of iterations
   int k_tilde = std::max((m-1)/2,1);
   int npar = m + m*(m-1)/2; // m entries in b, the rest in upper triangluar of Sigma. 
   Eigen::MatrixXd d1p(k_tilde,npar),d1pd(k_tilde,npar),d1pn(k_tilde,npar);
   d1p.setZero();
   
   Eigen::MatrixXd Hp(npar,npar),Hpd(npar,1),Hpn(npar,1);
   Hp.setZero();
   Hpd.setZero();
   Hpn.setZero();
   
   // initialize variables for the output
   Eigen::VectorXd p(k_tilde), p_out(1), grad(npar);
   p.setZero();
   p_out.setZero();
   grad.setZero();
   
   // // // // // // // // // // //
   // reorder the dimensions     //
   // // // // // // // // // // //
   
   struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
   
   Eigen::MatrixXd Cor_mat_temp = output_reorder.mat1;        // we'll work with this one
   Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
   Cor_mat = Cor_mat_temp; 
   
   Eigen::VectorXd x_temp = output_reorder.vec1;
   Eigen::VectorXd x_ordered = x_temp;            // we'll keep this as the original one
   x_norm= x_temp;
   
   Eigen::VectorXd temp1 = output_reorder.vec2;
   
   
   // // // // // // // // // // //
   // end of reordering          //
   // // // // // // // // // // //
   
   
   // initialize mu_zero variable
   Eigen::VectorXd mu_zero(m), mu_temp(m);
   mu_zero.setZero();
   mu_temp.setZero();
   
   Eigen::MatrixXd dpmu1(3,npar);
   
   
   // deal with low dimensional cases directly
   if(m==2)
   {
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2)));
     Eigen::VectorXd gr = grad_gen_cdf(x_temp - mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
     Eigen::MatrixXd Hess = Hess_gen_cdf(x_temp - mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
     
     // Hess_gen includes derivatives with respect to the diagonal entries. in the order: x0,x1, Sig00, Sig11, Sig01.
     Hp.block(0,0,2,2)=Hess.block(0,0,2,2);
     Hp(2,0)=Hess(4,0);
     Hp(2,1)=Hess(4,1);
     Hp(2,2)=Hess(4,4);
     Hp(1,2)=Hp(2,1);
     Hp(0,2)=Hp(2,0);
     
     Eigen::VectorXd grh(3);
     grh(0)=gr(0);
     grh(1)=gr(1);
     grh(2) =gr(4); 
     
     Hp = Hp/exp(p_out(0)) - grh*grh.transpose()/exp(2*p_out(0)); 
   }
   else if(m==3)
   {
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3)));
     Eigen::VectorXd gr = grad_gen_tvn_cdf(x_temp - mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
     Eigen::MatrixXd Hess = Hessian_gen_tvn_cdf(x_temp - mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
     Eigen::VectorXd grd(6);
     grd.head(3) = gr.head(3); 
     grd(3) = gr(4);
     grd(4)= gr(5);
     grd(5) = gr(7); 
     grd = grd/exp(p_out(0));
     
     Eigen::VectorXd indin(6);
     indin(0) = 0;
     indin(1)=1;
     indin(2)=2;
     indin(3)=4;
     indin(4)=5; 
     indin(5)=7;
     
     Hp = extract_submatrix(Hess,indin);
     
     Hp = Hp/exp(p_out(0));
     Hp += - grd* grd.transpose();
   }
   else if(m==4)
   {
     p_out(0) = log(TVBS_pmvnorm_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4)));
     Eigen::VectorXd gr = TVBS_grad_gen_qvn_cdf(x_temp - mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     Eigen::MatrixXd Hess = Hessian_gen_qvn_cdf(x_temp - mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     Eigen::VectorXd grd(10), indin(10); 
     
     indin(0)= 0;
     indin(1)= 1;
     indin(2)= 2;
     indin(3)= 3;
     indin(4)= 5;
     indin(5)= 6;
     indin(6)= 7;
     indin(7)= 9;
     indin(8)= 10;
     indin(9)=12;
     
     grd = extract_submatrix(gr,indin)/exp(p_out(0));
     Hp = extract_submatrix(Hess,indin);
     Hp = Hp/exp(p_out(0));
     Hp += - grd* grd.transpose();
   }
   else 
   {
     
     // calculate initial LDLT decomposition
     
     // initiate some variables
     struct TVBS_Matrix_three output_3;
     struct TVBS_Matrix_four output_4;
     Eigen::VectorXd mu_tilde;
     Eigen::VectorXd out_v(5);
     
     Eigen::MatrixXd Omega(2,2), HOmega(2,2), HCor_mat(m,m);
     Eigen::MatrixXd d1Omega(2*npar,2), d1L(m*npar,m), d1D(m*npar,m), d1Cor_mat(m*npar,m); 
     Eigen::MatrixXd d1mu_temp(m,npar), d1mu_tilde(m,npar), d1x_temp(m,npar), Hmu_temp(m*npar,npar); 
     Eigen::MatrixXd d1mu_temp_new(m,npar), Hmu_temp_new(m*npar,npar); 
     Eigen::VectorXd Hmu_tilde(m), Hx_temp(m), mu_temp_new(m); 
     
     d1mu_temp_new.setZero();
     Hmu_temp_new.setZero();
     
     Eigen::MatrixXd d1output_2, Houtput_2;
     struct Hessian_trunc_struct Hess_trunc;
     
     Omega.setZero(); 
     Hmu_temp.setZero();
     
     double p_k1_nomi, p_k1_denomi;
     
     
     //////////////////////////////////////////////////////////////////////////////////
     //  gradient calculations: iterate over the elements                          ///
     //  this allows to reuse matrices of the correct size and makes coding easier. //
     /////////////////////////////////////////////////////////////////////////////////
     
     d1p.setZero();
     d1pn.setZero();
     d1pd.setZero();
     Hp.setZero();
     Hpn.setZero();
     Hpd.setZero();
     
     
     // find the indices. 
     // ind_mat has three columns: first column: index in b, if applicable, -1 else; second column: row index; third column col index.
     Eigen::MatrixXd ind_mat (npar,3);
     ind_mat.setZero();
     
     int cur = -1;
     for (int nm=0; nm< m; nm++){
       cur+=1;
       ind_mat(nm,0)=nm;
     }
     for (int nm=0; nm< m; nm++){
       for (int nk=0;nk<m-nm-1;nk++){
         cur+=1;
         ind_mat(cur,0)=-1;
         ind_mat(cur,1)=nm;
         ind_mat(cur,2)=nm+nk+1;
       }
     }
     
     // calculations 
     // reset matrices that are changed in between
     Cor_mat_temp = Cor_mat_ordered;     
     mu_temp = mu_zero;
     x_temp = x_ordered;  
     
     struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
     Eigen::MatrixXd L = output_2.mat1;
     Eigen::MatrixXd D = output_2.mat2;
     
     // calculate P1 using the first four variables
     p(0) = TVBS_pmvnorm_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     // gradient of qvn w.r.t. its entries -> neede for all gradients. 
     Eigen::VectorXd gr= TVBS_grad_gen_qvn_cdf(x_temp.head(4) - mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
     
     d1L.setZero();
     d1D.setZero();
     d1Cor_mat.setZero();
     d1mu_tilde.setZero();
     d1x_temp.setZero(); 
     d1mu_temp.setZero();
     d1pd.setZero();
     d1pn.setZero();
     d1Omega.setZero();
     
     Eigen::VectorXd in1(npar);
     in1.setZero();
     for (int jj=0;jj<4;jj++){
       in1(jj)=1;
     }
     
     // gradient in first pass over 
     for (int np1=0; np1< npar; np1++){
       
       Eigen::VectorXd ind_1(2),ind_2(2);
       ind_1.setZero();
       ind_2.setZero();
       
       
       if (ind_mat(np1,0)>-1){
         d1x_temp((int)ind_mat(np1,0),np1)=1; // parameter corresponds to b entry.
         
       } else { // parameter corresponds to entry in covariance matrix.
         int r = ind_mat(np1,1);
         int c = ind_mat(np1,2);
         d1Cor_mat(np1*m+r,c)=1;
         d1Cor_mat(np1*m+c,r)=1;
         
         ind_1(0)=r; 
         ind_1(1)=c; 
         
         // derivative of initial LDLT
         struct TVBS_Matrix_two  output_dLDLT  = TVBS_grad_LDLT_decomp(Cor_mat_temp, ind_1);
         d1L.block(np1*m,0,m,m) = output_dLDLT.mat1;
         d1D.block(np1*m,0,m,m) = output_dLDLT.mat2;
       }
       
       // first four coordinates
       
       
       if ((ind_mat(np1,0)==-1)&(ind_mat(np1,1)<4)&(ind_mat(np1,2)<4)) {
         in1(np1)=1;
       }
       
       Eigen::VectorXd d1xmmu(4); //,d2xmmu(4);
       d1xmmu = d1x_temp.block(0,np1,4,1) - d1mu_temp.block(0,np1,4,1);
       
       d1p(0,np1)=0;
       
       int cur=0;
       if (in1(np1)>0){
         for (int nm=0;nm<4;nm++){
           d1p(0,np1) += gr(nm) * d1xmmu(nm);
           cur += 1;
         }
         for (int nm=0;nm<4;nm++){
           for (int nk=nm;nk<4;nk++){
             d1p(0,np1) += gr(cur) * d1Cor_mat(m*np1+nm,nk);
             cur +=1;
           }
         }
       }        
       
     } // loop over parameters for gradient calculations.
     
     
     
     Eigen::MatrixXd HL(m*npar,m*npar), HD(m*npar,m*npar);
     HL.setZero();
     HD.setZero();
     
     // start Hessian here 
     
     Eigen::MatrixXd Hess(13,13);
     Hess.setZero();
     Hess= Hessian_gen_qvn_cdf(x_temp.head(4) - mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
     
     
     for (int np1=0; np1< npar; np1++){
       Eigen::VectorXd d1vecin(14);
       d1vecin.setZero(); 
       
       if (in1(np1)>0){ // if one parameter is not affected, Hessian entry is zero!
         
         d1vecin.head(4) = d1x_temp.block(0,np1,4,1); 
         d1vecin.segment(4,4) = d1Cor_mat.block(np1*m+0,0,1,4);
         d1vecin.segment(8,3) = d1Cor_mat.block(np1*m+1,1,1,3);
         d1vecin.segment(11,2) = d1Cor_mat.block(np1*m+2,2,1,2);
         d1vecin(13) = d1Cor_mat(np1*m+3,3);
       } 
       
       
       for (int np2=np1; np2< npar; np2++){
         
         // for Hessian entry
         HOmega.setZero();
         HCor_mat = Cor_mat_ordered; 
         HCor_mat.setZero();
         
         Hmu_tilde.setZero();
         Hx_temp.setZero(); 
         Hpd.setZero();
         Hpn.setZero();
         
         if ((ind_mat(np1,0)==-1)&(ind_mat(np2,0)==-1)){
           Eigen::VectorXd ind_1(2),ind_2(2);
           ind_1.setZero();
           ind_2.setZero(); 
           
           int r = ind_mat(np1,1);
           int c = ind_mat(np1,2);
           
           ind_1(0)=r; 
           ind_1(1)=c;
           r = ind_mat(np2,1);
           c = ind_mat(np2,2);
           ind_2(0)=r; 
           ind_2(1)=c;
           
           struct TVBS_Matrix_two out_HLDLT = TVBS_Hess_LDLT_decomp(Cor_mat_temp, ind_1, ind_2);
           HL.block(np1*m,np2*m,m,m) = out_HLDLT.mat1;
           HD.block(np1*m,np2*m,m,m) = out_HLDLT.mat2; 
         }
         
         // 
         int inh=0;
         inh= in1(np1)*in1(np2);
         
         
         Eigen::VectorXd Hvecin(14), d2vecin(14); 
         Hvecin.setZero();
         d2vecin.setZero(); 
         
         
         if (in1(np2)>0){
           d2vecin.head(4) = d1x_temp.block(0,np2,4,1); 
           d2vecin.segment(4,4) = d1Cor_mat.block(np2*m+0,0,1,4);
           d2vecin.segment(8,3) = d1Cor_mat.block(np2*m+1,1,1,3);
           d2vecin.segment(11,2) = d1Cor_mat.block(np2*m+2,2,1,2);
           d2vecin(13) = d1Cor_mat(np2*m+3,3);
         }
         
         
         if (inh>0){
           Hvecin.segment(4,4) = HCor_mat.block(0,0,1,4);
           Hvecin.segment(8,3) = HCor_mat.block(1,1,1,3);
           Hvecin.segment(11,2) = HCor_mat.block(2,2,1,2);
           Hvecin(13) = HCor_mat(3,3);
           
           double Hpp = d1vecin.transpose() * Hess * d2vecin  + gr.dot(Hvecin); 
           Hp(np1,np2) = Hpp/p(0) - d1p(0,np1)*d1p(0,np2)/pow(p(0),2);
           Hp(np2,np1) = Hp(np1,np2);
         }        
       } // end loop over np2
     } // end loop over np1
     
     
     // start here from k1=1 to make indexing easier
     for(int k1=1; k1<k_tilde; k1++)
     {
       
       // calculate values 
       
       out_v = TVBS_biv_gen_trunc(x_temp.segment((2*k1-1)-1,2), mu_temp.head(2), D.block(0,0,2,2));
       mu_tilde = out_v.segment(0,2);
       
       Omega(0,0) = out_v(2);
       Omega(1,1) = out_v(4);
       Omega(0,1)= out_v(3);
       Omega(1,0)= out_v(3); 
       
       Eigen::MatrixXd Lo = L;
       Eigen::MatrixXd Do = D; 
       Eigen::MatrixXd d1Lo = d1L;
       Eigen::MatrixXd d1Do = d1D;
       
       output_2 = TVBS_LDLT_update_cpp(L,D,Omega,2);
       L = output_2.mat1;
       D = output_2.mat2;
       
       mu_temp_new = mu_temp.segment(2,m-2*k1) + Lo.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
       
       // deal with gradient 
       
       
       Eigen::MatrixXd d1mu_tilde(2,npar);
       d1mu_tilde.setZero();
       
       d1output_2 =  TVBS_Jacobian_biv_gen_trunc(x_temp.segment((2*k1-1)-1,2), mu_temp.head(2), Do.block(0,0,2,2));
       
       
       int he;
       he= 2*k1+4;
       in1.setZero();
       for (int np1=0; np1< npar; np1++){
         
         if (np1<he){
           in1(np1)=1;
         } else {
           if ((ind_mat(np1,1)<he)&(ind_mat(np1,2)<he)) {
             in1(np1)=1;
           }
         }
         
         
         
         if (in1(np1)>0){  
           
           d1mu_tilde(0,np1) = d1output_2(0,0)*d1x_temp((2*k1-1)-1,np1) + d1output_2(0,1)*d1x_temp((2*k1-1),np1) + d1output_2(0,2)*d1mu_temp(0,np1)+d1output_2(0,3)*d1mu_temp(1,np1);
           d1mu_tilde(0,np1) += d1output_2(0,4)* d1Do(np1*m+0,0) + d1output_2(0,5)*d1Do(np1*m+0,1) + d1output_2(0,6) * d1Do(np1*m+1,1);
           d1mu_tilde(1,np1) = d1output_2(1,0)*d1x_temp((2*k1-1)-1,np1) + d1output_2(1,1)*d1x_temp((2*k1-1),np1) + d1output_2(1,2)*d1mu_temp(0,np1)+d1output_2(1,3)*d1mu_temp(1,np1);
           d1mu_tilde(1,np1) += d1output_2(1,4)* d1Do(np1*m+0,0) + d1output_2(1,5)*d1Do(np1*m+0,1) + d1output_2(1,6) * d1Do(np1*m+1,1);
           
           // update mu_temp. 
           d1mu_temp_new.block(0,np1,m-2*k1,1) = d1mu_temp.block(2,np1,m-2*k1,1) + d1Lo.block(np1*m+2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2)) + Lo.block(2,0,m-2*k1,2)*(d1mu_tilde.block(0,np1,2,1)-d1mu_temp.block(0,np1,2,1));
           
           
           // Attention! doutput_2 ordered strangely as b0,b1,om00,om11,om01. 
           d1Omega(np1*2+0,0) = d1output_2(2,4)* d1Do(np1*m+0,0) + d1output_2(2,5)*d1Do(np1*m+0,1) + d1output_2(2,6) * d1Do(np1*m+1,1);
           d1Omega(np1*2+0,0) += d1output_2(2,0)*d1x_temp((2*k1-1)-1,np1) + d1output_2(2,1)*d1x_temp((2*k1-1),np1) + d1output_2(2,2)*d1mu_temp(0,np1)+d1output_2(2,3)*d1mu_temp(1,np1);
           d1Omega(np1*2+0,1) = d1output_2(3,4)* d1Do(np1*m+0,0) + d1output_2(3,5)*d1Do(np1*m+0,1) + d1output_2(3,6) * d1Do(np1*m+1,1);
           d1Omega(np1*2+0,1) += d1output_2(3,0)*d1x_temp((2*k1-1)-1,np1) + d1output_2(3,1)*d1x_temp((2*k1-1),np1) + d1output_2(3,2)*d1mu_temp(0,np1)+d1output_2(3,3)*d1mu_temp(1,np1);
           d1Omega(np1*2+1,1) = d1output_2(4,4)* d1Do(np1*m+0,0) + d1output_2(4,5)*d1Do(np1*m+0,1) + d1output_2(4,6) * d1Do(np1*m+1,1);
           d1Omega(np1*2+1,1) += d1output_2(4,0)*d1x_temp((2*k1-1)-1,np1) + d1output_2(4,1)*d1x_temp((2*k1-1),np1) + d1output_2(4,2)*d1mu_temp(0,np1)+d1output_2(4,3)*d1mu_temp(1,np1);
           d1Omega(np1*2+1,0) = d1Omega(np1*2+0,1);
           
         } else {
           // update mu_temp. 
           d1mu_temp_new.block(0,np1,m-2*k1,1) = d1mu_temp.block(2,np1,m-2*k1,1) + d1Lo.block(np1*m+2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2)) + Lo.block(2,0,m-2*k1,2)*(d1mu_tilde.block(0,np1,2,1)-d1mu_temp.block(0,np1,2,1));
         }
         
         
         // update LDLT
         struct TVBS_Matrix_two d1output_LD;
         int nL = Lo.rows(); 
         
         
         d1output_LD = TVBS_LDLT_update_grad(Lo, Do, Omega, d1Lo.block(np1*m,0,nL,nL),  d1Do.block(np1*m,0,nL,nL), d1Omega.block(np1*2,0,2,2), 2);
         
         d1L.block(np1*m,0,L.rows(),L.cols()) = d1output_LD.mat1;
         d1D.block(np1*m,0,L.rows(),L.cols()) = d1output_LD.mat2;
         
         
       } // end loop over np1 for gradient. 
       
       
       Eigen::MatrixXd HLc(m,m), HDc(m,m), Hmu_tilde(npar*2,npar), HOmega(npar*2,npar*2);
       HLc.setZero();
       HDc.setZero();
       Hmu_tilde.setZero();
       HOmega.setZero(); 
       
       Hess_trunc = TVBS_Hessian_biv_gen_trunc(x_temp.segment((2*k1-1)-1,2), mu_temp.head(2), Do.block(0,0,2,2));
       
       // Hessian 
       for (int np1=0; np1< npar; np1++){
         Eigen::VectorXd d1vH(7); 
         d1vH.setZero();
         
         if (in1(np1)>0){
           d1vH(0)= d1x_temp((2*k1-1)-1,np1);
           d1vH(1) =d1x_temp((2*k1-1),np1);
           d1vH(2) = d1mu_temp(0,np1);
           d1vH(3)= d1mu_temp(1,np1); 
           d1vH(4) = d1Do(np1*m+0,0);
           d1vH(5) = d1Do(np1*m+0,1);
           d1vH(6) = d1Do(np1*m+1,1); 
         }
         for (int np2=np1; np2< npar; np2++){
           int inh=0; 
           inh= in1(np1)*in1(np2);
           
           
           if (inh>0){
             HLc = HL.block(np1*m,np2*m,m,m);
             HDc = HD.block(np1*m,np2*m,m,m);
             
             Eigen::VectorXd d2vH(7); 
             d2vH.setZero();
             
             if (in1(np2)>0){
               d2vH(0)= d1x_temp((2*k1-1)-1,np2);
               d2vH(1) =d1x_temp((2*k1-1),np2);
               d2vH(2) = d1mu_temp(0,np2);
               d2vH(3)= d1mu_temp(1,np2); 
               d2vH(4) = d1Do(np2*m+0,0);
               d2vH(5) = d1Do(np2*m+0,1);
               d2vH(6) = d1Do(np2*m+1,1);
             }
             
             Hmu_tilde(0+np1*2,np2) = d1vH.transpose() * Hess_trunc.mat1 * d2vH; 
             Hmu_tilde(0+np1*2,np2) +=  d1output_2(0,2)*Hmu_temp(np1*m+0,np2)+d1output_2(0,3)*Hmu_temp(np1*m+1,np2);
             Hmu_tilde(0+np1*2,np2) += d1output_2(0,4)* HDc(0,0) + d1output_2(0,5)*HDc(0,1) + d1output_2(0,6) * HDc(1,1);
             
             Hmu_tilde(1+np1*2,np2) = d1vH.transpose() * Hess_trunc.mat2 * d2vH; 
             Hmu_tilde(1+np1*2,np2) +=  d1output_2(1,2)*Hmu_temp(np1*m+0,np2)+d1output_2(1,3)*Hmu_temp(np1*m+1,np2);
             Hmu_tilde(1+np1*2,np2) += d1output_2(1,4)* HDc(0,0) + d1output_2(1,5)*HDc(0,1) + d1output_2(1,6) * HDc(1,1);
             
             
             
             // update mu_temp 
             
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) = Hmu_temp.block(np1*m+2,np2,m-2*k1,1) + HL.block(np1*m+2,np2*m+0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) +=  d1Lo.block(np2*m+2,0,m-2*k1,2)*(d1mu_tilde.block(0,np1,2,1)-d1mu_temp.block(0,np1,2,1)) + d1Lo.block(np1*m+2,0,m-2*k1,2)*(d1mu_tilde.block(0,np2,2,1)-d1mu_temp.block(0,np2,2,1));
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) +=  Lo.block(2,0,m-2*k1,2)*(Hmu_tilde.block(np1*2,np2,2,1)-Hmu_temp.block(np1*m+0,np2,2,1));
             
             // update Omega
             HOmega(0+np1*2,np2*2+0) = d1vH.transpose() * Hess_trunc.mat3 * d2vH;
             int jjo=2;
             HOmega(0+np1*2,np2*2+0) += d1output_2(jjo,2)*Hmu_temp(np1*m+0,np2)+d1output_2(jjo,3)*Hmu_temp(np1*m+1,np2);
             HOmega(0+np1*2,np2*2+0) += d1output_2(jjo,4)* HDc(0,0) + d1output_2(jjo,5)*HDc(0,1) + d1output_2(jjo,6) * HDc(1,1);
             
             HOmega(0+np1*2,np2*2+1) = d1vH.transpose() * Hess_trunc.mat4 * d2vH;
             jjo=3;
             HOmega(0+np1*2,np2*2+1) += d1output_2(jjo,2)*Hmu_temp(np1*m+0,np2)+d1output_2(jjo,3)*Hmu_temp(np1*m+1,np2);
             //HOmega(0+np1*2,np2*2+1) += d1output_2(jjo,0)*Hx_temp((2*k1-1)-1) + d1output_2(jjo,1)*Hx_temp((2*k1-1)) + d1output_2(jjo,2)*Hmu_temp(0)+d1output_2(jjo,3)*Hmu_temp(1);
             HOmega(0+np1*2,np2*2+1) += d1output_2(jjo,4)* HDc(0,0) + d1output_2(jjo,5)*HDc(0,1) + d1output_2(jjo,6) * HDc(1,1);
             
             HOmega(1+np1*2,np2*2+1) = d1vH.transpose() * Hess_trunc.mat5 * d2vH;
             jjo=4;
             HOmega(1+np1*2,np2*2+1) += d1output_2(jjo,2)*Hmu_temp(np1*m+0,np2)+d1output_2(jjo,3)*Hmu_temp(np1*m+1,np2);
             //HOmega(1+np1*2,np2*2+1) += d1output_2(jjo,0)*Hx_temp((2*k1-1)-1) + d1output_2(jjo,1)*Hx_temp((2*k1-1)) + d1output_2(jjo,2)*Hmu_temp(0)+d1output_2(jjo,3)*Hmu_temp(1);
             HOmega(1+np1*2,np2*2+1) += d1output_2(jjo,4)* HDc(0,0) + d1output_2(jjo,5)*HDc(0,1) + d1output_2(jjo,6) * HDc(1,1);
             
             HOmega(1+np1*2,np2*2+0)= HOmega(0+np1*2,np2*2+1); 
             
           } else {
             // update mu_temp 
             
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) = Hmu_temp.block(np1*m+2,np2,m-2*k1,1) + HL.block(np1*m+2,np2*m+0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) +=  d1Lo.block(np2*m+2,0,m-2*k1,2)*(d1mu_tilde.block(0,np1,2,1)-d1mu_temp.block(0,np1,2,1)) + d1Lo.block(np1*m+2,0,m-2*k1,2)*(d1mu_tilde.block(0,np2,2,1)-d1mu_temp.block(0,np2,2,1));
             Hmu_temp_new.block(np1*m,np2,m-2*k1,1) +=  Lo.block(2,0,m-2*k1,2)*(Hmu_tilde.block(np1*2,np2,2,1)-Hmu_temp.block(np1*m+0,np2,2,1));
             
             
           }
           
           
           // update LDLT. 
           struct TVBS_Matrix_two Houtput_LD;
           int nL = Lo.rows(); 
           
           Houtput_LD = TVBS_LDLT_update_Hess(Lo,Do,Omega,d1Lo.block(np1*m,0,nL,nL),d1Do.block(np1*m,0,nL,nL),d1Omega.block(np1*2,0,2,2),d1Lo.block(np2*m,0,nL,nL),d1Do.block(np2*m,0,nL,nL),d1Omega.block(np2*2,0,2,2),HL.block(np1*m,np2*m,nL,nL),HD.block(np1*m,np2*m,nL,nL), HOmega.block(np1*2,np2*2,2,2),2);
           // Hessian 
           HL.block(np1*m,np2*m,L.rows(),L.cols()) = Houtput_LD.mat1;
           HD.block(np1*m,np2*m,L.rows(),L.cols()) = Houtput_LD.mat2; 
           
         } // end loop over np2
       } // end loop over np1
       
       
       mu_temp = mu_temp_new; 
       d1mu_temp = d1mu_temp_new;
       Hmu_temp = Hmu_temp_new;
       
       Cor_mat_temp = L*D*(L.transpose());
       
       int nL = L.rows();
       
       d1Cor_mat.setZero(); 
       
       for (int np1=0; np1< npar; np1++)
       {
         d1Cor_mat.block(np1*m,0,nL,nL) = d1L.block(np1*m,0,nL,nL) * D*(L.transpose()) + L * d1D.block(np1*m,0,nL,nL) * (L.transpose()) + L * D *((d1L.block(np1*m,0,nL,nL)).transpose());
       }
       
       // calculate new prob.
       if(m>=2*k1+4)
       {
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
         p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
         
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
         
         // quotient. 
         p(k1) = p_k1_nomi/p_k1_denomi;
         
         Eigen::VectorXd gr(13);
         gr.setZero();
         gr = TVBS_grad_gen_qvn_cdf(x_temp.segment((2*k1+1)-1,4)- mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
         
         Eigen::VectorXd grd(5);
         grd.setZero();
         grd = grad_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         
         
         for (int np1=0;np1<npar;np1++)
         {
           
           Eigen::VectorXd d1xmdmu = d1x_temp.block((2*k1+1)-1,np1,4,1)- d1mu_temp.block(0,np1,4,1);
           
           d1pn(k1,np1) = gr(0)*d1xmdmu(0) + gr(1)*d1xmdmu(1) + gr(2)*d1xmdmu(2) + gr(3)*d1xmdmu(3);
           int cur = 3; 
           for (int nr=0;nr<4;nr++){
             for (int nc=nr;nc<4;nc++){
               cur += 1;
               d1pn(k1,np1) += gr(cur) * d1Cor_mat(np1*m+nr,nc); 
             }
           } 
           
           
           d1pd(k1,np1) = grd(0) * d1xmdmu(0) +  grd(1) * d1xmdmu(1);
           d1pd(k1,np1) += grd(2) * d1Cor_mat(np1*m+0,0) + grd(3) * d1Cor_mat(np1*m+1,1) + grd(4) * d1Cor_mat(np1*m+0,1); 
           
           d1p(k1,np1) = (d1pn(k1,np1) * p_k1_denomi - p_k1_nomi *d1pd(k1,np1))/(p_k1_denomi * p_k1_denomi); 
           
           
         }
         
         
         Eigen::MatrixXd Hessn(5,5);
         Hessn.setZero(); 
         
         
         Hessn = Hess_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         
         Eigen::MatrixXd Hess(13,13);
         Hess.setZero(); 
         Hess = Hessian_gen_qvn_cdf(x_temp.segment((2*k1+1)-1,4)- mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
         
         
         for (int np1=0; np1< npar; np1++)
         {
           
           Eigen::VectorXd d1vecin(14); 
           d1vecin.setZero();
           
           Eigen::VectorXd gr1n(5); 
           gr1n.setZero();
           
           
           if (in1(np1)>0){
             d1vecin.head(4) = d1x_temp.block((2*k1+1)-1,np1,4,1)-d1mu_temp.block(0,np1,4,1); 
             d1vecin.segment(4,4) = d1Cor_mat.block(np1*m+0,0,1,4);
             d1vecin.segment(8,3) = d1Cor_mat.block(np1*m+1,1,1,3);
             d1vecin.segment(11,2) = d1Cor_mat.block(np1*m+2,2,1,2);
             d1vecin(13) = d1Cor_mat(np1*m+3,3);
             
             gr1n.head(2) = d1vecin.head(2);
             gr1n(2) = d1vecin(4);
             gr1n(3) = d1vecin(8); 
             gr1n(4) = d1vecin(5);
             
           }
           
           for (int np2=np1;np2<npar;np2++){
             // HCor_mat from HL, HD. 
             int nL = L.rows();
             Eigen::MatrixXd HCor_mat(nL,nL);
             HCor_mat.setZero(); 
             
             HCor_mat = HL.block(np1*m,np2*m,nL,nL) *D*(L.transpose()) + d1L.block(np1*m,0,nL,nL)* d1D.block(np2*m,0,nL,nL) *(L.transpose()) + d1L.block(np1*m,0,nL,nL) * D* (d1L.block(np2*m,0,nL,nL).transpose());
             HCor_mat += d1L.block(np2*m,0,nL,nL) * d1D.block(np1*m,0,nL,nL) * (L.transpose()) + L * HD.block(np1*m,np2*m,nL,nL) * (L.transpose()) + L * d1D.block(np1*m,0,nL,nL) *(d1L.block(np2*m,0,nL,nL).transpose());
             HCor_mat += d1L.block(np2*m,0,nL,nL) * D * (d1L.block(np1*m,0,nL,nL).transpose()) + L * d1D.block(np2*m,0,nL,nL) * (d1L.block(np1*m,0,nL,nL).transpose()) + L * D* (HL.block(np1*m,np2*m,nL,nL).transpose()); 
             
             // Hessian is 14 x 14: 4 b's, 4 first row, 3 sec. row, ...,1 4x4 entry. 
             Eigen::VectorXd Hvecin(14), d2vecin(14); 
             Hvecin.setZero();
             
             Eigen::VectorXd gr2n(5), Hn(5); 
             gr2n.setZero();
             Hn.setZero(); 
             
             
             if (in1(np2)>0){
               d2vecin.setZero(); 
               d2vecin.head(4) = d1x_temp.block((2*k1+1)-1,np2,4,1)-d1mu_temp.block(0,np2,4,1); 
               d2vecin.segment(4,4) = d1Cor_mat.block(np2*m+0,0,1,4);
               d2vecin.segment(8,3) = d1Cor_mat.block(np2*m+1,1,1,3);
               d2vecin.segment(11,2) = d1Cor_mat.block(np2*m+2,2,1,2);
               d2vecin(13) = d1Cor_mat(np2*m+3,3);
               
               gr2n.head(2) = d2vecin.head(2);
               gr2n(2) = d2vecin(4);
               gr2n(3) = d2vecin(8); 
               gr2n(4) = d2vecin(5);
               
             }  
             
             int inh=0;
             inh = in1(np1)*in1(np2);
             
             double Hpn = 0;
             double Hpd =0;
             
             if (inh>0){
               Hvecin.head(4) = -Hmu_temp.block(np1*m,np2,4,1); 
               Hvecin.segment(4,4) = HCor_mat.block(0,0,1,4);
               Hvecin.segment(8,3) = HCor_mat.block(1,1,1,3);
               Hvecin.segment(11,2) = HCor_mat.block(2,2,1,2);
               Hvecin(13) = HCor_mat(3,3);
               
               Hn.head(2) = Hvecin.head(2);
               Hn(2) = Hvecin(4);
               Hn(3) = Hvecin(8); 
               Hn(4) = Hvecin(5);
               
               Hpn = d1vecin.transpose() * Hess * d2vecin  + gr.dot(Hvecin); 
               
               // denominator 
               // Hessian 
               Hpd = gr1n.transpose() * Hessn * gr2n + grd.dot(Hn);
             }
             
             // Hessian. 
             double brack = d1pn(k1,np1)*p_k1_denomi - p_k1_nomi * d1pd(k1,np1); 
             double d2brack = Hpn*p_k1_denomi + d1pn(k1,np1)*d1pd(k1,np2) - d1pn(k1,np2)*d1pd(k1,np1) - p_k1_nomi * Hpd; 
             
             double Hpp = (d2brack * p_k1_denomi - 2*d1pd(k1,np2)*brack)/pow(p_k1_denomi,3);
             Hp(np1,np2) += Hpp/p(k1)-  d1p(k1,np1)*d1p(k1,np2)/(pow(p(k1),2)); 
             Hp(np2,np1) = Hp(np1,np2);
           } // end for np2
         } // end for np1
         
       } // end not last iteration 
       else
       {
         // last iteration, if uneven dimension
         
         
         // calculate next probability
         p_k1_nomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
         p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
         
         p_k1_denomi = TVBS_pmvnorm_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
         
         p(k1) = p_k1_nomi/p_k1_denomi;
         
         Eigen::VectorXd gr(9);
         gr.setZero();
         
         
         gr= grad_gen_tvn_cdf(x_temp.segment((2*k1+1)-1,3)- mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
         
         Eigen::VectorXd grd(5);
         grd.setZero();
         grd = grad_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         
         
         for (int np1=0;np1<npar;np1++)
         {
           
           
           Eigen::VectorXd d1xmdmu = d1x_temp.block((2*k1+1)-1,np1,3,1)- d1mu_temp.block(0,np1,3,1);
           
           
           if (in1(np1)){
             d1pn(k1,np1) = gr(0)*d1xmdmu(0) + gr(1)*d1xmdmu(1) + gr(2)*d1xmdmu(2);
             int cur = 2; 
             for (int nr=0;nr<3;nr++){
               for (int nc=nr;nc<3;nc++){
                 cur += 1;
                 d1pn(k1,np1) += gr(cur) * d1Cor_mat(np1*m+nr,nc); 
               }
             } 
             
             
             d1pd(k1,np1) = grd(0) * d1xmdmu(0) +  grd(1) * d1xmdmu(1);
             d1pd(k1,np1) += grd(2) * d1Cor_mat(np1*m+0,0) + grd(3) * d1Cor_mat(np1*m+1,1) + grd(4) * d1Cor_mat(np1*m+0,1); 
             
             
             d1p(k1,np1) = (d1pn(k1,np1) * p_k1_denomi - p_k1_nomi *d1pd(k1,np1))/(p_k1_denomi * p_k1_denomi); 
           }
           
         }
         
         
         Eigen::MatrixXd Hess(9,9);
         Hess.setZero(); 
         Hess = Hessian_gen_tvn_cdf(x_temp.segment((2*k1+1)-1,3)- mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
         
         Eigen::MatrixXd Hessn(5,5);
         Hessn.setZero(); 
         
         Hessn = Hess_gen_cdf(x_temp.segment((2*k1+1)-1,2)- mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
         
         
         
         for (int np1=0; np1< npar; np1++)
         {
           Eigen::VectorXd d1vecin(9); 
           d1vecin.setZero(); 
           
           Eigen::VectorXd gr1n(5);
           gr1n.setZero();
           
           if (in1(np1)>0){
             d1vecin.head(3) = d1x_temp.block((2*k1+1)-1,np1,3,1)-d1mu_temp.block(0,np1,3,1); 
             d1vecin.segment(3,3) = d1Cor_mat.block(np1*m+0,0,1,3);
             d1vecin.segment(6,2) = d1Cor_mat.block(np1*m+1,1,1,2);
             d1vecin(8) = d1Cor_mat(np1*m+2,2);
             
             
             
             gr1n.head(2) = d1vecin.head(2);
             gr1n(2) = d1vecin(3);
             gr1n(3) = d1vecin(6); 
             gr1n(4) = d1vecin(4);
           } 
           
           double brack = d1pn(k1,np1)*p_k1_denomi - p_k1_nomi * d1pd(k1,np1); 
           
           for (int np2=np1;np2<npar;np2++){
             // HCor_mat from HL, HD. 
             int nL = L.rows();
             Eigen::MatrixXd HCor_mat(nL,nL);
             HCor_mat.setZero(); 
             
             HCor_mat = HL.block(np1*m,np2*m,nL,nL) *D*(L.transpose()) + d1L.block(np1*m,0,nL,nL)* d1D.block(np2*m,0,nL,nL) *(L.transpose()) + d1L.block(np1*m,0,nL,nL) * D* (d1L.block(np2*m,0,nL,nL).transpose());
             HCor_mat += d1L.block(np2*m,0,nL,nL) * d1D.block(np1*m,0,nL,nL) * (L.transpose()) + L * HD.block(np1*m,np2*m,nL,nL) * (L.transpose()) + L * d1D.block(np1*m,0,nL,nL) *(d1L.block(np2*m,0,nL,nL).transpose());
             HCor_mat += d1L.block(np2*m,0,nL,nL) * D * (d1L.block(np1*m,0,nL,nL).transpose()) + L * d1D.block(np2*m,0,nL,nL) * (d1L.block(np1*m,0,nL,nL).transpose()) + L * D* (HL.block(np1*m,np2*m,nL,nL).transpose()); 
             
             
             
             
             // Hessian is 14 x 14: 4 b's, 4 first row, 3 sec. row, ...,1 4x4 entry. 
             Eigen::VectorXd Hvecin(9), d2vecin(9); 
             Hvecin.setZero();
             d2vecin.setZero(); 
             
             Eigen::VectorXd gr2n(5), Hn(5); 
             gr2n.setZero();
             Hn.setZero(); 
             
             if (in1(np2)>0){
               
               d2vecin.head(3) = d1x_temp.block((2*k1+1)-1,np2,3,1)-d1mu_temp.block(0,np2,3,1); 
               d2vecin.segment(3,3) = d1Cor_mat.block(np2*m+0,0,1,3);
               d2vecin.segment(6,2) = d1Cor_mat.block(np2*m+1,1,1,2);
               d2vecin(8) = d1Cor_mat(np2*m+2,2);
               
               
               gr2n.head(2) = d2vecin.head(2);
               gr2n(2) = d2vecin(3);
               gr2n(3) = d2vecin(6); 
               gr2n(4) = d2vecin(4);
               
               
             }
             
             int inh=0;
             inh = in1(np1)*in1(np2); 
             
             double Hpn=0;
             double Hpd = 0; 
             
             if (inh>0){
               Hvecin.head(3) = -Hmu_temp.block(np1*m,np2,3,1); 
               Hvecin.segment(3,3) = HCor_mat.block(0,0,1,3);
               Hvecin.segment(6,3) = HCor_mat.block(1,1,1,2);
               Hvecin(8) = HCor_mat(2,2);
               
               Hn.head(2) = Hvecin.head(2);
               Hn(2) = Hvecin(3);
               Hn(3) = Hvecin(6); 
               Hn(4) = Hvecin(4);
               
               
               Hpn = d1vecin.transpose() * Hess * d2vecin  + gr.dot(Hvecin); 
               
               // denominator 
               // Hessian 
               
               Hpd = gr1n.transpose() * Hessn * gr2n + grd.dot(Hn);
               
               // Hessian. 
             }
             
             double d2brack = Hpn*p_k1_denomi + d1pn(k1,np1)*d1pd(k1,np2) - d1pn(k1,np2)*d1pd(k1,np1) - p_k1_nomi * Hpd; 
             
             double Hpp = (d2brack * p_k1_denomi - 2*d1pd(k1,np2)*brack)/pow(p_k1_denomi,3);
             Hp(np1,np2) += Hpp/p(k1)-  d1p(k1,np1)*d1p(k1,np2)/(pow(p(k1),2)); 
             Hp(np2,np1) = Hp(np1,np2);
           } // end for np2
         } // end for np1
         
       } // end last iteration
       
     }   
   } // end if(m>4). 
   
   p_out(0) = p.array().log().sum();
   
   Eigen::MatrixXd Hp_reord(Hp.rows(),Hp.cols());
   Hp_reord.setZero();
   
   
   Eigen::VectorXd ind_ord = reorder_indices(temp1);
   
   Hp_reord = reorder_matrix(Hp, ind_ord);
   return Hp_reord;
 }

