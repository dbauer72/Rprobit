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
#include "TVBS_vdb_helper_2.h"

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


// // // // // // // // // // // // // // // // //
//                                              //
//    PDF AND CDF FUNCTIONS                     //
//                                              //
// // // // // // // // // // // // // // // // //

double cal_delta(double w0, double w1, double rho)
{
  double cdf_0 = std_normal_cdf((w1 - rho*w0)/sqrt(1-rho*rho));
  double std_0 = std_normal_pdf(w0);
  
  double delta0 = std_0*cdf_0;
  
  
  double cdf_1 = std_normal_cdf((w0 - rho*w1)/sqrt(1-rho*rho));
  double std_1 = std_normal_pdf(w1);
  
  double delta1 = std_1*cdf_1;
  
  double Phi2 = biv_normal_cdf(w0,w1,rho);
  double lambda_0 =  -(delta0 +rho*delta1)/Phi2;
  double lambda_1 =  -(delta1 +rho*delta0)/Phi2;
  
  // omega_00. 
  
  double phi2 = biv_normal_pdf(w0,w1,rho);
  
  // omega_01
  double exp_01 = rho*w1*delta1 + rho*w0*delta0  - (1-rho*rho)*phi2;
  double omega_01 = rho - exp_01/Phi2 - lambda_0*lambda_1; 
  
  return delta0;
}

Eigen::MatrixXd Hessian_delta(double w0, double w1, double rho)
{
  Eigen::MatrixXd Hess(3,3);
  
  // 
  double cdf_0 = std_normal_cdf((w1 - rho*w0)/sqrt(1-rho*rho));
  double std_0 = std_normal_pdf(w0);
  
  // gradient of std_0 
  Eigen::VectorXd grad_std_0(3);
  grad_std_0.setZero(); 
  grad_std_0(0)= std_0*(-w0);
  
  // gradient of tilde w0
  double tildew0 = (w1- rho*w0)/sqrt(1-rho*rho); 
  Eigen::VectorXd grad_tildew0(3);
  
  grad_tildew0.setZero();
  grad_tildew0(0) = -rho/sqrt(1-rho*rho); 
  grad_tildew0(1) = 1/sqrt(1-rho*rho);
  grad_tildew0(2) = (w1*rho-w0)/(1-rho*rho)/sqrt(1-rho*rho);
  
  // Hessian of tildew0 
  Eigen::MatrixXd H_tildew0(3,3); 
  H_tildew0.setZero();
  
  H_tildew0(1,2) = rho/(1-rho*rho)/sqrt(1-rho*rho);
  H_tildew0(0,2) = -rho*H_tildew0(1,2) - 1/sqrt(1-rho*rho);
  H_tildew0(2,0)= H_tildew0(0,2);
  H_tildew0(2,1)= H_tildew0(1,2);
  
  H_tildew0(2,2) = w1/(1-rho*rho)/sqrt(1-rho*rho)+ (w1*rho-w0)*rho*3/pow(sqrt(1-rho*rho),5);
  
  // Hessian of std_normal_pdf(w0).
  Eigen::MatrixXd H_std_0(3,3);
  H_std_0.setZero();
  H_std_0(0,0) = std_0*(w0*w0- 1);
  
  // calculate Hessian 
  Hess = std_0*std_normal_pdf(tildew0)*(-tildew0*grad_tildew0 * grad_tildew0.transpose() + H_tildew0)+H_std_0*cdf_0 + std_normal_pdf(tildew0)*(grad_std_0*grad_tildew0.transpose() + grad_tildew0*grad_std_0.transpose());
  
  return Hess; 
}




// // // // // // // // // // // // // // // // //
//                                              //
//    Seperation Of Variables Method (SOV)      //
//                                              //
// // // // // // // // // // // // // // // // //




double TVBS_std_normal_cdf_inv(Eigen::VectorXd w, Eigen::VectorXd b, Eigen::MatrixXd L)
{
  //double tol = 0.000001;
  int k = b.size();
  Eigen::VectorXd e(k), z_norm_inv(k-1);
  
  e(0) = std_normal_cdf(b(0)/L(0,0));
  
  for(int i = 1; i < k; i++ )
  {
    // THIS IS A FIX TO AVOID NUMERICAL PROBLEMS TAKING PHI^-1(0)
    z_norm_inv(i-1) = std_normal_cdf_inv_1(std::max(e(i-1)*w(i-1), tol));
    
    e(i) = std_normal_cdf(((b(i)-((L.row(i).head(i).array())*(z_norm_inv.transpose().head(i).array())).sum())/L(i,i)));
  }
  
  return e.prod();
}


// reorder variables for better numerical results
struct TVBS_Vec_Mat TVBS_reorder_chol_cpp(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  int k = b.size();
  Eigen::MatrixXd L_chol(k,k), Sigma_row_store(1,k), Sigma_col_store(k,1), Sigma_row_new(1,k), Sigma_col_new(k,1);
  L_chol.setZero();
  double b_store, s_ii, s_ij, s_ji, s_jj;
  
  for(int i_1=0; i_1<k; i_1++)
  {
    // check for smallest value
    int i_current = i_1;
    for(int i_2=i_1; i_2<k; i_2++){
      if( b(i_current)/sqrt(Sigma(i_current,i_current)) > b(i_2)/sqrt(Sigma(i_2,i_2)) )
      {
        i_current = i_2;
      }
    }
    
    // if you found a smaller value
    if((i_current > i_1) && (i_1<k-1))
    {
      // change variables
      b_store = b(i_1);
      b(i_1) = b(i_current);
      b(i_current) = b_store;
      
      s_ii = Sigma(i_1,i_1);
      s_ij = Sigma(i_1, i_current);
      s_ji = Sigma(i_current, i_1);
      s_jj = Sigma(i_current,i_current);
      Sigma_row_store = Sigma.block(i_1,0,1,k);
      
      Sigma_col_store = Sigma.block(0,i_1,k,1);
      
      Sigma_row_new = Sigma.block(i_current,0,1,k);
      
      Sigma_col_new = Sigma.block(0,i_current,k,1);
      
      Sigma.block(0,i_1,k,1) = Sigma_col_new;
      Sigma.block(i_1,0,1,k) = Sigma_row_new;
      Sigma.block(0,i_current,k,1) = Sigma_col_store;
      Sigma.block(i_current,0,1,k) = Sigma_row_store;
      
      Sigma(i_1,i_1) = s_jj;
      Sigma(i_1, i_current) = s_ji;
      Sigma(i_current, i_1) = s_ij;
      Sigma(i_current,i_current) = s_ii;
      
      Sigma_row_store = L_chol.block(i_1,0,1,k);
      L_chol.block(i_1,0,1,k) = L_chol.block(i_current,0,1,k);
      L_chol.block(i_current,0,1,k) = Sigma_row_store;
    }
    
    
    // calculate cholesky
    if(i_1==0)
    {
      L_chol.col(0) = Sigma.col(0)/sqrt(Sigma(0,0));
    }
    if(i_1 == 1)
    {
      L_chol(1,1) = sqrt(Sigma(1,1) - pow(L_chol(1,0),2));
      for(int i_3=2; i_3<k; i_3++)
      {
        L_chol(i_3,1) = (Sigma(i_3,1)-L_chol(1,0)*L_chol(i_3,0))/L_chol(1,1);
      }
    }
    if(i_1 > 1)
    {
      L_chol(i_1,i_1) = sqrt(Sigma(i_1,i_1) - ( (L_chol.row(i_1).array())*(L_chol.row(i_1).array()) ).sum() );
      for(int i_3=i_1+1; i_3<k; i_3++)
      {
        L_chol(i_3,i_1) = ( Sigma(i_3,i_1) - ( (L_chol.row(i_3).array())*(L_chol.row(i_1).array()) ).sum() )/L_chol(i_1,i_1);
      }
    }
  }
  
  struct TVBS_Vec_Mat output;
  output.vec1 = b;
  output.mat1 = L_chol;
  return output;
}

// reorder variables for better numerical results
struct TVBS_Vec_Mat_Vec TVBS_reorder_sig_cpp(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  int k = b.size();
  Eigen::VectorXd order_vec(k), b_tilde(k), y_i(k);
  y_i.setZero();
  order_vec = Eigen::VectorXd::LinSpaced(k,0,k-1);
  Eigen::MatrixXd L_chol(k,k), Sigma_row_store(1,k), Sigma_col_store(k,1), Sigma_row_new(1,k), Sigma_col_new(k,1);
  L_chol.setZero();
  double b_store, s_ii, s_ij, s_ji, s_jj, l_jj_current, l_jj_new, v_select_current, v_select_new;
  
  for(int i_1=0; i_1<k; i_1++)
  {
    
    // calculate candidate for b_tilde and l_jj
    if(i_1==0)
    {
      l_jj_current = sqrt(Sigma(i_1,i_1));
      b_tilde(0) = b(0)/l_jj_current;
    }
    if(i_1==1)
    {
      l_jj_current = sqrt(Sigma(1,1) - pow(L_chol(1,0),2));
      b_tilde(1) = (b(1) - L_chol(1,0)*y_i(0))/l_jj_current;
    }
    if(i_1>1)
    {
      l_jj_current = sqrt(Sigma(i_1,i_1) - ( (L_chol.row(i_1).array())*(L_chol.row(i_1).array()) ).sum() );
      b_tilde(i_1) = (b(i_1) - ( (L_chol.row(i_1).array())*(y_i.array()) ).sum())/l_jj_current;
    }
    
    
    // check for smallest value
    int i_current = i_1;
    v_select_current = 1 - b_tilde(i_1)*std_normal_pdf(b_tilde(i_1))/std_normal_cdf(b_tilde(i_1)) - pow( std_normal_pdf(b_tilde(i_1))/std_normal_cdf(b_tilde(i_1)),2 );
    for(int i_2=i_1; i_2<k; i_2++){
      
      // calculate candidate for b_tilde and l_jj
      if(i_1==0)
      {
        l_jj_new = sqrt(Sigma(i_2,i_2));
        b_tilde(i_2) = b(i_2)/l_jj_new;
      }
      if(i_1==1)
      {
        l_jj_new = sqrt(Sigma(i_2,i_2) - pow(L_chol(i_2,0),2));
        b_tilde(i_2) = (b(i_2) - L_chol(i_2,0)*y_i(0))/l_jj_new;
      }
      if(i_1>1)
      {
        l_jj_new = sqrt(Sigma(i_2,i_2) - ( (L_chol.row(i_2).array())*(L_chol.row(i_2).array()) ).sum() );
        b_tilde(i_2) = (b(i_2) - ( (L_chol.row(i_2).array())*(y_i.array()) ).sum())/l_jj_new;
      }
      
      
      v_select_new = 1 - b_tilde(i_2)*std_normal_pdf(b_tilde(i_2))/std_normal_cdf(b_tilde(i_2)) - pow( std_normal_pdf(b_tilde(i_2))/std_normal_cdf(b_tilde(i_2)),2 );
      if( v_select_new < v_select_current )
      {
        i_current = i_2;
        v_select_current = v_select_new;
        l_jj_current = l_jj_new;
      }
    }
    
    
    
    // calculate y value
    y_i(i_1) = -std_normal_pdf(b_tilde(i_current))/std_normal_cdf(b_tilde(i_current));
    
    
    
    // if you found a smaller value
    if((i_current > i_1) && (i_1<k-1))
    {
      
      // update order vector if a change was made
      b_store = order_vec(i_1);
      order_vec(i_1) = order_vec(i_current);
      order_vec(i_current) = b_store;
      
      // change variables
      b_store = b(i_1);
      b(i_1) = b(i_current);
      b(i_current) = b_store;
      
      s_ii = Sigma(i_1,i_1);
      s_ij = Sigma(i_1, i_current);
      s_ji = Sigma(i_current, i_1);
      s_jj = Sigma(i_current,i_current);
      Sigma_row_store = Sigma.block(i_1,0,1,k);
      
      Sigma_col_store = Sigma.block(0,i_1,k,1);
      
      Sigma_row_new = Sigma.block(i_current,0,1,k);
      
      Sigma_col_new = Sigma.block(0,i_current,k,1);
      
      Sigma.block(0,i_1,k,1) = Sigma_col_new;
      Sigma.block(i_1,0,1,k) = Sigma_row_new;
      Sigma.block(0,i_current,k,1) = Sigma_col_store;
      Sigma.block(i_current,0,1,k) = Sigma_row_store;
      
      Sigma(i_1,i_1) = s_jj;
      Sigma(i_1, i_current) = s_ji;
      Sigma(i_current, i_1) = s_ij;
      Sigma(i_current,i_current) = s_ii;
      
      Sigma_row_store = L_chol.block(i_1,0,1,k);
      L_chol.block(i_1,0,1,k) = L_chol.block(i_current,0,1,k);
      L_chol.block(i_current,0,1,k) = Sigma_row_store;
    }
    
    
    // calculate cholesky
    if(i_1==0)
    {
      L_chol.col(0) = Sigma.col(0)/sqrt(Sigma(0,0));
    }
    if(i_1 == 1)
    {
      //L_chol(1,1) = sqrt(Sigma(1,1) - pow(L_chol(1,0),2));
      L_chol(1,1) = l_jj_current;
      for(int i_3=2; i_3<k; i_3++)
      {
        L_chol(i_3,1) = (Sigma(i_3,1)-L_chol(1,0)*L_chol(i_3,0))/L_chol(1,1);
      }
    }
    if(i_1 > 1)
    {
      //L_chol(i_1,i_1) = sqrt(Sigma(i_1,i_1) - ( (L_chol.row(i_1).array())*(L_chol.row(i_1).array()) ).sum() );
      L_chol(i_1,i_1) = l_jj_current;
      for(int i_3=i_1+1; i_3<k; i_3++)
      {
        L_chol(i_3,i_1) = ( Sigma(i_3,i_1) - ( (L_chol.row(i_3).array())*(L_chol.row(i_1).array()) ).sum() )/L_chol(i_1,i_1);
      }
    }
    
  }
  
  struct TVBS_Vec_Mat_Vec output;
  output.vec1 = b;
  output.mat1 = Sigma;
  output.vec2 = order_vec;
  return output;
}

/// [[Rcpp::export]]
Eigen::VectorXd TVBS_biv_norm_trunc(Eigen::VectorXd w, double rho)
{
  Eigen::VectorXd lamb_ome(5); // lambda_{0,1}, Omega_{0,0},Omega_{0,1},Omega_{1,1} (indexation according to CPP to match code) 
  lamb_ome.setZero();
  
  // calculate deltas 
  double cdf_0 = std_normal_cdf((w(1) - rho*w(0))/sqrt(1-rho*rho));
  double delta_0 = std_normal_pdf(w(0))*cdf_0;
  double cdf_1 = std_normal_cdf((w(0) - rho*w(1))/sqrt(1-rho*rho));
  double delta_1 = std_normal_pdf(w(1))*cdf_1;
  
  
  // calculate lambdas
  double Phi2 = biv_normal_cdf(w(0),w(1),rho);
  double lambda_0 = -(delta_0+rho*delta_1)/Phi2;
  double lambda_1 = -(delta_1+rho*delta_0)/Phi2;
  
  // correct rounding errors for small values of w(0) and/or w(1)
  if (lambda_0> w(0)){ lambda_0 = w(0)-0.001;}
  if (lambda_1> w(1)){ lambda_1 = w(1)-0.001;}
  //Rcout << cdf_0 << cdf_1 << "del" << delta_0 << delta_1 << "Phi" << Phi2 << std::endl; 
  
  // corresponding mu_tildes 
  double mutilde0 = w(0) + lambda_0;
  
  
  // turn to the Omegas. 
  double phi2 = biv_normal_pdf(w(0),w(1),rho);
  
  // Omega_{00}
  double exp_0 = w(0)*delta_0 + rho*rho*w(1)*delta_1 - rho*(1-rho*rho)*phi2;
  double omega_00 = 1 - exp_0/Phi2 - lambda_0*lambda_0; 
  
  if (omega_00<0.0001) {omega_00=0.0001; }
  // Omega_{11}
  double exp_1 = w(1)*delta_1 + rho*rho*w(0)*delta_0 - rho*(1-rho*rho)*phi2;
  double omega_11 = 1 -  exp_1/Phi2 - lambda_1*lambda_1; 
  
  if (omega_11<0.0001) {omega_11=0.0001; }
  
  // Omega_{01}
  double exp_01 = rho*w(1)*delta_1 + rho*w(0)*delta_0 - (1-rho*rho)*phi2;
  double omega_01 = rho-  exp_01/Phi2 - lambda_0*lambda_1; 
  
  double rho_01;
  rho_01 = omega_01/std::sqrt(omega_00*omega_11);
  
  if (rho_01>0.999) { rho_01 = 0.999;}
  if (rho_01< -0.999) { rho_01 = 0.999;}
  omega_01 = rho_01 * std::sqrt(omega_00*omega_11);
  // fill in results 
  lamb_ome << lambda_0, lambda_1, omega_00, omega_01, omega_11; 
  
  // return matrix. 
  return lamb_ome;
}


/// [[Rcpp::export]]
Eigen::VectorXd TVBS_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd tildex = x-alpha; 
  Eigen::MatrixXd D(2,2);
  D.setZero();
  
  D(0,0)= sqrt(Sigma(0,0));
  D(1,1)= sqrt(Sigma(1,1));
  double rho = Sigma(0,1)/(D(0,0)*D(1,1));
  
  //Rcout << rho << tildex << D << std::endl; 
  
  tildex = D.inverse() * tildex; 
  Eigen::VectorXd lambda_Omega = TVBS_biv_norm_trunc(tildex, rho);
  
  //Rcout << lambda_Omega << std::endl; 
  lambda_Omega(0) = alpha(0) + D(0,0)*lambda_Omega(0);
  lambda_Omega(1) = alpha(1) + D(1,1)*lambda_Omega(1);
  lambda_Omega(2) = Sigma(0,0)*lambda_Omega(2);
  lambda_Omega(3) = D(0,0)*D(1,1)*lambda_Omega(3);
  lambda_Omega(4) = Sigma(1,1)*lambda_Omega(4);
  
  return lambda_Omega; 
}





double tildePhithree(Eigen::VectorXd x, Eigen::MatrixXd Sigma)
{
  double tilde_phi =0;
  
  // truncate the first two variables -> return new means 
  Eigen::VectorXd out = TVBS_biv_gen_trunc(x.segment(0,2), x.segment(0,2)*0, Sigma.block(0,0,2,2));
  
  Eigen::MatrixXd Omega(2,2);
  Omega(0,0)=out(2);
  Omega(1,1)=out(4);
  Omega(0,1)=out(3);
  Omega(1,0) = out(3);
  
  // update means and variances of three and four 
  double mu = 0; 
  double tildeSig = 0; 
  Eigen::MatrixXd L21(1,2);
  L21.setZero();
  
  L21 = Sigma.block(2,0,1,2) * Sigma.block(0,0,2,2).inverse();
  mu = x(2) - L21(0,0) * out(0)- L21(0,1)*out(1);
  tildeSig = Sigma(2,2) - Sigma(0,2) * L21(0,0) - Sigma(1,2) * L21(0,1) + L21(0,0) * L21(0,0)* Omega(0,0) + (2*L21(0,0) * L21(0,1)* Omega(0,1)) + L21(0,1) * L21(0,1)* Omega(1,1);
  
  // calculate probs: first two. 
  double Phi2 = biv_gen_cdf(x.segment(0,2), Sigma.block(0,0,2,2));
  
  // calculate fourth conditional on third (for first two truncated) as ratio of Phi_2 / Phi.
  double Phi3 =  std_normal_cdf(mu/std::sqrt(tildeSig)); 
  
  tilde_phi = Phi3 * Phi2; 
  
  if (NumericVector::is_na(tilde_phi)) {tilde_phi =tol;}
  if (tilde_phi<tol) {tilde_phi =tol;}
  return tilde_phi;
}

// // // // // // // // //
//  quadrovariate       //
// // // // // // // // //

double tildePhifour(Eigen::VectorXd x, Eigen::MatrixXd Sigma)
{
  double tilde_phi =0;
  
  // truncate the first two variables -> return new means 
  Eigen::VectorXd out = TVBS_biv_gen_trunc(x.segment(0,2), x.segment(0,2)*0, Sigma.block(0,0,2,2));
  
  Eigen::MatrixXd Omega(2,2);
  Omega(0,0)=out(2);
  Omega(1,1)=out(4);
  Omega(0,1)=out(3);
  Omega(1,0) = out(3);
  
  // update means and variances of three and four 
  Eigen::VectorXd mu(2);
  Eigen::MatrixXd L21(2,2);
  L21.setZero();
  mu.setZero();
  
  L21 = Sigma.block(2,0,2,2) * Sigma.block(0,0,2,2).inverse();
  mu = x.segment(2,2) - L21 * out.segment(0,2);
  Eigen::MatrixXd tildeSig(2,2);
  tildeSig.setZero();
  
  tildeSig = Sigma.block(2,2,2,2) + (L21*Omega - Sigma.block(2,0,2,2)) * L21.transpose();
  
  // calculate probs: first three. 
  double Phi3 = tildePhithree(x.segment(0,3), Sigma.block(0,0,3,3)); // TVBS_pmvnorm_cpp(x.segment(0,3),x.segment(0,3)*0 , Sigma.block(0,0,3,3));
  
  // calculate fourth conditional on third (for first two truncated) as ratio of Phi_2 / Phi.
  double Phi2 =  biv_gen_cdf(mu, tildeSig) / std_normal_cdf(mu(0)/std::sqrt(tildeSig(0,0))); 
  
  tilde_phi = Phi3 * Phi2;   
  //Rcout << Phi2 << Phi3 << std::endl; 
  
  if (NumericVector::is_na(tilde_phi)) {tilde_phi =tol;}
  if (tilde_phi<tol) {tilde_phi =tol;}
  
  return tilde_phi;
}

//' multivariate Gaussian CDF: calculates CDF up to four dimensions
//' @description
//' The function computes the  Gaussian CDF for 1, 2, 3 and four dimensional problems. 
//' @param upper
//' nx1 vector of points where to evaluate 
//' @param mu
//' nx1 vector; expectation
//' @param Sigma
//' nxn correlation matrix.
//' @return 
//' vector; gradient of log probability, double; (log of) probability. 
//' @keywords internal
//'
// [[Rcpp::export]]
double TVBS_pmvnorm_cpp(Eigen::VectorXd upper, Eigen::VectorXd mu, Eigen::MatrixXd Sigma)
{
   int k         = upper.rows();
   double output = 1;
   
   if( k == 1 ){
     output       = std_normal_cdf((upper(0) - mu(0)) / sqrt(Sigma(0,0)));
   }
   if( k == 2 ){
     output       = biv_normal_cdf((upper(0)-mu(0))/sqrt(Sigma(0,0)), (upper(1)-mu(1))/sqrt(Sigma(1,1)), Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1)));
   }
   if( k == 3 ) {
     output = tildePhithree(upper - mu, Sigma);    
   }
   if( k == 4 ) {
     output = tildePhifour(upper - mu, Sigma);    
   }
   if (NumericVector::is_na(output)) {output =tol;}
   if (output<tol) {output =tol;}
   
   return output;
 }



// // // // // // // // // // // // // // // // //
//                                              //
//    LDLT DECOMPOSITION AND UPDATING           //
//                                              //
// // // // // // // // // // // // // // // // //


struct TVBS_Matrix_two TVBS_LDLT_decomp_cpp(Eigen::MatrixXd Sigma)                        // create a function to calculate two-block LDL decomposition
{
  // Based on Algorithm 3.1 from the paper
  // Bivariate conditioning approximations for multivariate normal probabilities
  // By Giang Trinh and Alan Genz
  // DOI: 10.1007/s11222-014-9468-y
  int n = Sigma.rows();
  Eigen::MatrixXd L(n,n), D(n,n);
  L.setZero();
  D.setZero();
  L.diagonal().setOnes();
  int i_ext = 0;
  
  for(int i = 0; i < (n-2); i+=2)
  {
    D.block(i,i,2,2)                  = Sigma.block(i,i,2,2);
    L.block(i+2,i,n-i-2,2)            = Sigma.block(i+2,i,n-i-2,2) * ((D.block(i,i,2,2)).inverse());
    Sigma.block(i+2,i+2,n-i-2,n-i-2)  = Sigma.block(i+2,i+2,n-i-2,n-i-2) - L.block(i+2,i,n-i-2,2) * D.block(i,i,2,2) * (L.block(i+2,i,n-i-2,2)).transpose();
    i_ext = i;
  }
  D.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = Sigma.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  
  if(n<3)
  {
    D = Sigma;
  }
  
  
  struct TVBS_Matrix_two output_final = {L, D};
  
  return output_final;
}


struct TVBS_Matrix_two TVBS_grad_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind)                        // create a function to calculate two-block LDL decomposition
{
  // Based on Algorithm 3.1 from the paper
  // Bivariate conditioning approximations for multivariate normal probabilities
  // By Giang Trinh and Alan Genz
  // DOI: 10.1007/s11222-014-9468-y
  int n = Sigma.rows();
  Eigen::MatrixXd L(n,n), D(n,n), dL(n,n), dD(n,n), dSigma(n,n);
  L.setZero();
  D.setZero();
  dL.setZero();
  dD.setZero(); 
  dSigma.setZero();
  
  L.diagonal().setOnes();
  int i_ext = 0;
  
  int r = ind(0); // row index
  int c= ind(1);  // column index
  
  dSigma(r,c)=1; 
  dSigma(c,r)=1; 
  
  for(int i = 0; i < (n-2); i+=2)
  {
    D.block(i,i,2,2)                  = Sigma.block(i,i,2,2);
    L.block(i+2,i,n-i-2,2)            = Sigma.block(i+2,i,n-i-2,2) * ((D.block(i,i,2,2)).inverse());
    Sigma.block(i+2,i+2,n-i-2,n-i-2)  = Sigma.block(i+2,i+2,n-i-2,n-i-2) - L.block(i+2,i,n-i-2,2) * D.block(i,i,2,2) * (L.block(i+2,i,n-i-2,2)).transpose();
    i_ext = i;
    
    // adjust derivatives 
    dD.block(i,i,2,2)                  = dSigma.block(i,i,2,2);
    dL.block(i+2,i,n-i-2,2)            = (dSigma.block(i+2,i,n-i-2,2) - L.block(i+2,i,n-i-2,2) *dD.block(i,i,2,2)) * (D.block(i,i,2,2)).inverse();
    dSigma.block(i+2,i+2,n-i-2,n-i-2)  = dSigma.block(i+2,i+2,n-i-2,n-i-2) - dL.block(i+2,i,n-i-2,2) * Sigma.block(i+2,i,n-i-2,2).transpose(); 
    dSigma.block(i+2,i+2,n-i-2,n-i-2) -= L.block(i+2,i,n-i-2,2) * dSigma.block(i+2,i,n-i-2,2).transpose();
  }
  
  // final part. 
  D.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = Sigma.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  dD.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = dSigma.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  
  if(n<3)
  {
    D = Sigma;
    dD = dSigma;
  }
  
  
  struct TVBS_Matrix_two output_final = {dL, dD};
  
  return output_final;
}


// calculation of gradient w.r.t. entry (r_a,c_a).

struct TVBS_Matrix_two TVBS_Hess_LDLT_decomp(Eigen::MatrixXd Sigma, Eigen::VectorXd ind1, Eigen::VectorXd ind2)                        // create a function to calculate two-block LDL decomposition
{
  // Based on Algorithm 3.1 from the paper
  // Bivariate conditioning approximations for multivariate normal probabilities
  // By Giang Trinh and Alan Genz
  // DOI: 10.1007/s11222-014-9468-y
  int n = Sigma.rows();
  Eigen::MatrixXd L(n,n), D(n,n), dL1(n,n), dD1(n,n), dSigma1(n,n), dL2(n,n), dD2(n,n), dSigma2(n,n);
  Eigen::MatrixXd iD2(2,2);
  
  // for the Hessian 
  Eigen::MatrixXd hL12(n,n), hD12(n,n), hSigma12(n,n);
  hL12.setZero();
  hD12.setZero();
  hSigma12.setZero();
  
  iD2.setZero();
  L.setZero();
  D.setZero();
  dL1.setZero();
  dD1.setZero(); 
  dSigma1.setZero();
  dL2.setZero();
  dD2.setZero(); 
  dSigma2.setZero();
  
  L.diagonal().setOnes();
  
  int i_ext = 0;
  
  int r1 = ind1(0); // row index 1
  int c1= ind1(1);  // column index 1
  
  int r2 = ind2(0); // row index 2
  int c2= ind2(1);  // column index 2
  
  dSigma1(r1,c1)=1; 
  dSigma1(c1,r1)=1; 
  dSigma2(r2,c2)=1; 
  dSigma2(c2,r2)=1; 
  
  for(int i = 0; i < (n-2); i+=2)
  {
    D.block(i,i,2,2)                  = Sigma.block(i,i,2,2);
    iD2 = (D.block(i,i,2,2)).inverse();
    L.block(i+2,i,n-i-2,2)            = Sigma.block(i+2,i,n-i-2,2) * iD2;
    Sigma.block(i+2,i+2,n-i-2,n-i-2)  = Sigma.block(i+2,i+2,n-i-2,n-i-2) - L.block(i+2,i,n-i-2,2) * D.block(i,i,2,2) * (L.block(i+2,i,n-i-2,2)).transpose();
    i_ext = i;
    
    // adjust derivatives 1
    dD1.block(i,i,2,2)                  = dSigma1.block(i,i,2,2);
    dL1.block(i+2,i,n-i-2,2)            = (dSigma1.block(i+2,i,n-i-2,2) - L.block(i+2,i,n-i-2,2) *dD1.block(i,i,2,2)) * iD2;
    dSigma1.block(i+2,i+2,n-i-2,n-i-2)  = dSigma1.block(i+2,i+2,n-i-2,n-i-2) - dL1.block(i+2,i,n-i-2,2) * Sigma.block(i+2,i,n-i-2,2).transpose(); 
    dSigma1.block(i+2,i+2,n-i-2,n-i-2) -= L.block(i+2,i,n-i-2,2) * dSigma1.block(i+2,i,n-i-2,2).transpose();
    
    // adjust derivatives 2
    dD2.block(i,i,2,2)                  = dSigma2.block(i,i,2,2);
    dL2.block(i+2,i,n-i-2,2)            = (dSigma2.block(i+2,i,n-i-2,2) - L.block(i+2,i,n-i-2,2) *dD2.block(i,i,2,2)) * iD2;
    dSigma2.block(i+2,i+2,n-i-2,n-i-2)  = dSigma2.block(i+2,i+2,n-i-2,n-i-2) - dL2.block(i+2,i,n-i-2,2) * Sigma.block(i+2,i,n-i-2,2).transpose(); 
    dSigma2.block(i+2,i+2,n-i-2,n-i-2) -= L.block(i+2,i,n-i-2,2) * dSigma2.block(i+2,i,n-i-2,2).transpose();
    
    // recursion for Hessian 
    
    hD12.block(i,i,2,2)   = hSigma12.block(i,i,2,2);
    hL12.block(i+2,i,n-i-2,2) = hSigma12.block(i+2,i,n-i-2,2) * iD2;
    hL12.block(i+2,i,n-i-2,2) -= dSigma1.block(i+2,i,n-i-2,2) * iD2 * dD2.block(i,i,2,2) * iD2;
    hL12.block(i+2,i,n-i-2,2) -= dL2.block(i+2,i,n-i-2,2) * dD1.block(i,i,2,2) * iD2;
    hL12.block(i+2,i,n-i-2,2) -= L.block(i+2,i,n-i-2,2) * hD12.block(i,i,2,2) *  iD2;
    hL12.block(i+2,i,n-i-2,2) +=  L.block(i+2,i,n-i-2,2) * dD1.block(i,i,2,2) * iD2 * dD2.block(i,i,2,2) * iD2;
    
    hSigma12.block(i+2,i+2,n-i-2,n-i-2)  -=  hSigma12.block(i+2,i,n-i-2,2)*L.block(i+2,i,n-i-2,2).transpose();
    hSigma12.block(i+2,i+2,n-i-2,n-i-2)  -=  dSigma1.block(i+2,i,n-i-2,2)*dL2.block(i+2,i,n-i-2,2).transpose();
    hSigma12.block(i+2,i+2,n-i-2,n-i-2) += dL2.block(i+2,i,n-i-2,2) * dD1.block(i,i,2,2) * L.block(i+2,i,n-i-2,2).transpose(); 
    hSigma12.block(i+2,i+2,n-i-2,n-i-2) += L.block(i+2,i,n-i-2,2) * hD12.block(i,i,2,2) * L.block(i+2,i,n-i-2,2).transpose();   
    hSigma12.block(i+2,i+2,n-i-2,n-i-2) += L.block(i+2,i,n-i-2,2) * dD1.block(i,i,2,2) * dL2.block(i+2,i,n-i-2,2).transpose();   
    hSigma12.block(i+2,i+2,n-i-2,n-i-2) -= dL2.block(i+2,i,n-i-2,2) * dSigma1.block(i+2,i,n-i-2,2).transpose();   
    hSigma12.block(i+2,i+2,n-i-2,n-i-2) -= L.block(i+2,i,n-i-2,2) * hSigma12.block(i+2,i,n-i-2,2).transpose();   
  }
  
  // final part for gradient
  D.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = Sigma.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  dD1.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = dSigma1.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  dD2.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = dSigma2.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  
  hD12.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2) = hSigma12.block(i_ext+2,i_ext+2,n-i_ext-2,n-i_ext-2);
  
  
  
  // special case 
  if(n<3)
  {
    D = Sigma;
    dD1 = dSigma1;
  }
  
  
  struct TVBS_Matrix_two output_final = {hL12, hD12};
  
  return output_final;
}


struct TVBS_Matrix_two TVBS_LDLT_update_cpp(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, int m = 2)   
{
  // create a function to specifically do the rank-2 update of the Covariance Matrix as proposed in "Property 2" and calculate its LDLT-decomposition 
  int n     = L_k.rows();
  int n1    = n-m;
  Eigen::MatrixXd M1(m,m); 
  M1 = Omega;
  Eigen::MatrixXd B(n1,m), D_2(n1,n1), L_22(n1,n1), L_21(n1,m), z(n1,m), Cov_Y_BR(n1,n1);
  B.setZero();
  z.setZero();
  
  D_2       = D_k.block(m,m,n1,n1);
  L_22      = L_k.block(m,m,n1,n1);
  L_21      = L_k.block(m,0,n1,m);
  z = L_22.inverse() * L_21;
  
  Eigen::MatrixXd newL(n1,n1), newD(n1,n1), Dinv(m,m), t(m,m);
  Dinv.setZero();
  
  newL = Eigen::MatrixXd::Identity(n1, n1);
  newD.setZero(); 
  
  if (n1 <= m){
    newD = z* Omega* z.transpose() + D_2;
  }else {
    for (int jj=0;jj<n1; jj= jj +m){ // increment by m
      //Rcout << jj << std::endl;
      if (jj+m <= n1){
        t = M1* z.block(jj,0,m,m).transpose();
        newD.block(jj,jj,m,m) = D_2.block(jj,jj,m,m)+z.block(jj,0,m,m) * t;
        Dinv = newD.block(jj,jj,m,m).inverse();
        B.block(jj,0,m,m)= Dinv * t.transpose();
        M1 = M1 - t*B.block(jj,0,m,m);
        if (jj<(n1-m)){
          if (jj != (n1-m-1)){
            newL.block(jj+m,0,m,jj+m) = z.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose();
          } else if (jj == n1-m-1){
            newL.block(jj+m,0,1,jj+m) = z.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose();
          }
        }
      } else if (jj == n1-1){
        Eigen::VectorXd t(m);
        t = M1*z.row(jj).transpose();
        newD(jj,jj) = D_2(jj,jj) + z.row(jj)*t;
      }
    }
  }
  newL = L_22* newL;
  
  struct TVBS_Matrix_two LDLT_new = {newL, newD};
  
  return LDLT_new;
}

// gradient for updates
struct TVBS_Matrix_two TVBS_LDLT_update_grad(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd dL_k, Eigen::MatrixXd dD_k, Eigen::MatrixXd dOmega, int m = 2)   
{
  int n     = L_k.rows();
  int n1    = n-m;
  Eigen::MatrixXd M1(m,m), dM1(m,m); 
  M1 = Omega;
  dM1 = dOmega;
  
  Eigen::MatrixXd B(n1,m), D_2(n1,n1), L_22(n1,n1), L_21(n1,m), z(n1,m), dB(n1,m), dD_2(n1,n1), dL_22(n1,n1), dL_21(n1,m),dz(n1,m);
  B.setZero(); dB.setZero();
  z.setZero(); dz.setZero();
  
  D_2       = D_k.block(m,m,n1,n1);
  L_22      = L_k.block(m,m,n1,n1);
  L_21      = L_k.block(m,0,n1,m);
  
  dD_2       = dD_k.block(m,m,n1,n1);
  dL_22      = dL_k.block(m,m,n1,n1);
  dL_21      = dL_k.block(m,0,n1,m);
  
  Eigen::MatrixXd iL22(n1,n1);
  iL22.setZero();
  iL22 = L_22.inverse();
  
  z = iL22 * L_21;
  dz = iL22 * dL_21 - iL22 * dL_22 *z ;
  
  Eigen::MatrixXd newL(n1,n1), newD(n1,n1), Dinv(m,m), t(m,m), dnewL(n1,n1), dnewD(n1,n1),dDinv(m,m),dt(m,m);
  Dinv.setZero();dDinv.setZero();
  dnewL.setZero();
  dnewD.setZero();
  dDinv.setZero();
  dt.setZero();
  t.setZero(); 
  
  newL = Eigen::MatrixXd::Identity(n1, n1);
  newD.setZero(); 
  
  if (n1 <= m){
    newD = z* Omega* z.transpose() + D_2;
    dnewD = dz*Omega * z.transpose() + z * dOmega * z.transpose() + z*Omega * dz.transpose() + dD_2;
  }else {
    for (int jj=0;jj<n1; jj= jj +m){ // increment by m
      
      //Rcout << "jj" << jj << std::endl;
      if (jj+m <= n1){
        t = M1* z.block(jj,0,m,m).transpose();
        dt = dM1 * z.block(jj,0,m,m).transpose() + M1* dz.block(jj,0,m,m).transpose();
        
        newD.block(jj,jj,m,m) = D_2.block(jj,jj,m,m)+z.block(jj,0,m,m) * t;
        dnewD.block(jj,jj,m,m) = dD_2.block(jj,jj,m,m)+ dz.block(jj,0,m,m) * t+z.block(jj,0,m,m) * dt; 
        
        Dinv = newD.block(jj,jj,m,m).inverse();
        dDinv = - Dinv * dnewD.block(jj,jj,m,m) * Dinv;
        
        B.block(jj,0,m,m)= Dinv * t.transpose();
        dB.block(jj,0,m,m)= dDinv * t.transpose()+ Dinv * dt.transpose();
        
        M1 = M1 - t*B.block(jj,0,m,m);
        dM1 = dM1 - dt*B.block(jj,0,m,m)- t*dB.block(jj,0,m,m);
        
        //Rcout << "dnewL" << dnewL << std::endl; 
        
        if (jj<(n1-m)){
          if (jj != (n1-m-1)){
            newL.block(jj+m,0,m,jj+m) = z.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose();
            dnewL.block(jj+m,0,m,jj+m) = dz.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose()+z.block(jj+m,0,m,m)* dB.block(0,0,jj+m,m).transpose();
          } else if (jj == n1-m-1){
            newL.block(jj+m,0,1,jj+m) = z.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose();
            dnewL.block(jj+m,0,1,jj+m) = dz.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose()+ z.block(jj+m,0,1,m)* dB.block(0,0,jj+m,m).transpose();
          }
        }
      } else if (jj == n1-1){
        
        
        Eigen::VectorXd t(m), dt(m);
        t = M1*z.row(jj).transpose();
        
        dt = dM1*z.row(jj).transpose();  
        dt += M1*dz.row(jj).transpose();
        
        newD(jj,jj) = D_2(jj,jj) + z.row(jj)*t;
        dnewD(jj,jj) = dD_2(jj,jj) + dz.row(jj)*t + z.row(jj)*dt;
      }
    }
  }
  
  Eigen::MatrixXd nL(n1,n1);
  nL.setZero();
  nL = L_22* newL;
  dnewL = dL_22* newL + L_22 * dnewL;
  
  // return updated values 
  struct TVBS_Matrix_two LDLT_grad;
  LDLT_grad.mat1 = dnewL;
  LDLT_grad.mat2 = dnewD;
  
  return LDLT_grad;
}


// Hessian for updates
struct TVBS_Matrix_two TVBS_LDLT_update_Hess(Eigen::MatrixXd L_k, Eigen::MatrixXd D_k, Eigen::MatrixXd Omega, Eigen::MatrixXd d1L_k, Eigen::MatrixXd d1D_k, Eigen::MatrixXd d1Omega, Eigen::MatrixXd d2L_k, Eigen::MatrixXd d2D_k, Eigen::MatrixXd d2Omega, Eigen::MatrixXd hL_k, Eigen::MatrixXd hD_k, Eigen::MatrixXd hOmega, int m = 2)   
{
  
  int n     = L_k.rows();
  int n1    = n-m;
  Eigen::MatrixXd M1(m,m), d1M1(m,m), d2M1(m,m), hM1(m,m); 
  M1 = Omega;
  d1M1 = d1Omega;
  d2M1 = d2Omega;
  hM1 = hOmega;
  
  Eigen::MatrixXd B(n1,m), D_2(n1,n1), L_22(n1,n1), L_21(n1,m), z(n1,m);
  Eigen::MatrixXd d1B(n1,m), d1D_2(n1,n1), d1L_22(n1,n1), d1L_21(n1,m),d1z(n1,m);
  Eigen::MatrixXd d2B(n1,m), d2D_2(n1,n1), d2L_22(n1,n1), d2L_21(n1,m),d2z(n1,m);
  Eigen::MatrixXd hB(n1,m), hD_2(n1,n1), hL_22(n1,n1), hL_21(n1,m),hz(n1,m);
  
  B.setZero(); d1B.setZero(); d2B.setZero(); hB.setZero();
  z.setZero(); d1z.setZero();d2z.setZero(), hz.setZero();
  
  D_2       = D_k.block(m,m,n1,n1);
  L_22      = L_k.block(m,m,n1,n1);
  L_21      = L_k.block(m,0,n1,m);
  
  d1D_2       = d1D_k.block(m,m,n1,n1);
  d1L_22      = d1L_k.block(m,m,n1,n1);
  d1L_21      = d1L_k.block(m,0,n1,m);
  d2D_2       = d2D_k.block(m,m,n1,n1);
  d2L_22      = d2L_k.block(m,m,n1,n1);
  d2L_21      = d2L_k.block(m,0,n1,m);
  hD_2       = hD_k.block(m,m,n1,n1);
  hL_22      = hL_k.block(m,m,n1,n1);
  hL_21      = hL_k.block(m,0,n1,m);
  
  Eigen::MatrixXd L_22_i = L_22.inverse();
  z = L_22_i * L_21;
  d1z = L_22_i * d1L_21 - L_22_i * d1L_22 * z;
  d2z = L_22_i * d2L_21 - L_22_i * d2L_22 * z;
  hz = L_22_i * hL_21 - L_22_i * d1L_22 * L_22_i * d2L_21 - L_22_i * d2L_22 * d1z;
  hz = hz - L_22_i *hL_22 * z + L_22_i *d1L_22* L_22_i * d2L_22 * z;
  
  Eigen::MatrixXd newL(n1,n1), newD(n1,n1), t(m,m), Dinv(m,m);  
  Eigen::MatrixXd d1newL(n1,n1), d1newD(n1,n1),d1Dinv(m,m),d1t(m,m);
  Eigen::MatrixXd d2newL(n1,n1), d2newD(n1,n1),d2Dinv(m,m),d2t(m,m);
  Eigen::MatrixXd hnewL(n1,n1), hnewD(n1,n1),hDinv(m,m),ht(m,m);
  
  Dinv.setZero();   d1Dinv.setZero(); d2Dinv.setZero(); hDinv.setZero();
  d1newL.setZero();  d2newL.setZero();  hnewL.setZero();
  newD.setZero(); d1newD.setZero();d2newD.setZero(); hnewD.setZero();
  t.setZero(); d1t.setZero(); d2t.setZero(); ht.setZero();
  
  newL = Eigen::MatrixXd::Identity(n1, n1);
  
  if (n1 <= m){
    newD = z* Omega* z.transpose() + D_2;
    d1newD = d1z*Omega * z.transpose() + z * d1Omega * z.transpose() + z*Omega * d1z.transpose() + d1D_2;
    d2newD = d2z*Omega * z.transpose() + z * d2Omega * z.transpose() + z*Omega * d2z.transpose() + d2D_2;    
    hnewD = hz *Omega * z.transpose() + d2z * d1Omega * z.transpose() + d2z*Omega * d1z.transpose();
    hnewD += d1z * d2Omega * z.transpose() + z * hOmega * z.transpose() + z* d2Omega * d1z.transpose();
    hnewD += d1z * Omega * d2z.transpose() + z* d1Omega * d2z.transpose() + z *Omega * hz.transpose(); 
    hnewD += hD_2;    
    
  }else {
    for (int jj=0;jj<n1; jj= jj +m){ // increment by m
      //Rcout << jj << std::endl;
      if (jj+m <= n1){
        t = M1* z.block(jj,0,m,m).transpose();
        d1t = d1M1 * z.block(jj,0,m,m).transpose() + M1* d1z.block(jj,0,m,m).transpose();
        d2t = d2M1 * z.block(jj,0,m,m).transpose() + M1* d2z.block(jj,0,m,m).transpose();
        ht = hM1 * z.block(jj,0,m,m).transpose() + d2M1 * d1z.block(jj,0,m,m).transpose(); 
        ht += d1M1* d2z.block(jj,0,m,m).transpose() + M1* hz.block(jj,0,m,m).transpose();
        
        newD.block(jj,jj,m,m) = D_2.block(jj,jj,m,m)+z.block(jj,0,m,m) * t;
        d1newD.block(jj,jj,m,m) = d1D_2.block(jj,jj,m,m)+ d1z.block(jj,0,m,m) * t+z.block(jj,0,m,m) * d1t; 
        d2newD.block(jj,jj,m,m) = d2D_2.block(jj,jj,m,m)+ d2z.block(jj,0,m,m) * t+z.block(jj,0,m,m) * d2t; 
        hnewD.block(jj,jj,m,m) = hD_2.block(jj,jj,m,m)+ hz.block(jj,0,m,m) * t + d2z.block(jj,0,m,m) * d1t;
        hnewD.block(jj,jj,m,m) +=  d1z.block(jj,0,m,m) * d2t + z.block(jj,0,m,m) * ht; 
        
        Dinv = newD.block(jj,jj,m,m).inverse();
        d1Dinv = - Dinv * d1newD.block(jj,jj,m,m) * Dinv;
        d2Dinv = - Dinv * d2newD.block(jj,jj,m,m) * Dinv;
        hDinv = - d1Dinv * d2newD.block(jj,jj,m,m) * Dinv - Dinv* hnewD.block(jj,jj,m,m) * Dinv - Dinv * d2newD.block(jj,jj,m,m) * d1Dinv;
        
        B.block(jj,0,m,m)= Dinv * t.transpose();
        d1B.block(jj,0,m,m)= d1Dinv * t.transpose()+ Dinv * d1t.transpose();
        d2B.block(jj,0,m,m)= d2Dinv * t.transpose()+ Dinv * d2t.transpose();
        hB.block(jj,0,m,m)= hDinv * t.transpose()+ d2Dinv * d1t.transpose(); 
        hB.block(jj,0,m,m) += d1Dinv * d2t.transpose() + Dinv * ht.transpose();
        
        M1 = M1 - t*B.block(jj,0,m,m);
        d1M1 = d1M1 - d1t*B.block(jj,0,m,m)- t*d1B.block(jj,0,m,m);
        d2M1 = d2M1 - d2t*B.block(jj,0,m,m)- t*d2B.block(jj,0,m,m);
        hM1 = hM1 - ht*B.block(jj,0,m,m) - d2t*d1B.block(jj,0,m,m);
        hM1 = hM1 - d1t*d2B.block(jj,0,m,m)- t*hB.block(jj,0,m,m);
        
        if (jj<(n1-m)){
          if (jj != (n1-m-1)){
            newL.block(jj+m,0,m,jj+m) = z.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose();
            d1newL.block(jj+m,0,m,jj+m) = d1z.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose()+z.block(jj+m,0,m,m)* d1B.block(0,0,jj+m,m).transpose();
            d2newL.block(jj+m,0,m,jj+m) = d2z.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose()+z.block(jj+m,0,m,m)* d2B.block(0,0,jj+m,m).transpose();
            hnewL.block(jj+m,0,m,jj+m) = hz.block(jj+m,0,m,m)* B.block(0,0,jj+m,m).transpose() + d2z.block(jj+m,0,m,m)* d1B.block(0,0,jj+m,m).transpose();
            hnewL.block(jj+m,0,m,jj+m) +=  d1z.block(jj+m,0,m,m)* d2B.block(0,0,jj+m,m).transpose()+ z.block(jj+m,0,m,m)* hB.block(0,0,jj+m,m).transpose();
            
          } else if (jj == n1-m-1){
            newL.block(jj+m,0,1,jj+m) = z.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose();
            d1newL.block(jj+m,0,1,jj+m) = d1z.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose()+ z.block(jj+m,0,1,m)* d1B.block(0,0,jj+m,m).transpose();
            d2newL.block(jj+m,0,1,jj+m) = d2z.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose()+ z.block(jj+m,0,1,m)* d2B.block(0,0,jj+m,m).transpose();
            hnewL.block(jj+m,0,1,jj+m) = hz.block(jj+m,0,1,m)* B.block(0,0,jj+m,m).transpose()+d2z.block(jj+m,0,1,m)* d1B.block(0,0,jj+m,m).transpose(); 
            hnewL.block(jj+m,0,1,jj+m) += d1z.block(jj+m,0,1,m)* d2B.block(0,0,jj+m,m).transpose()+z.block(jj+m,0,1,m)* hB.block(0,0,jj+m,m).transpose();
            
          }
        }
      } else if (jj == n1-1){
        Eigen::VectorXd t(m), d1t(m), d2t(m), ht(m);
        t = M1*z.row(jj).transpose();
        d1t = d1M1*z.row(jj).transpose()+M1*d1z.row(jj).transpose();
        d2t = d2M1*z.row(jj).transpose()+M1*d2z.row(jj).transpose();
        ht = hM1*z.row(jj).transpose()+d2M1*d1z.row(jj).transpose();
        ht += d1M1*d2z.row(jj).transpose()+M1*hz.row(jj).transpose();
        
        newD(jj,jj) = D_2(jj,jj) + z.row(jj)*t;
        d1newD(jj,jj) = d1D_2(jj,jj) + d1z.row(jj)*t + z.row(jj)*d1t;
        d2newD(jj,jj) = d2D_2(jj,jj) + d2z.row(jj)*t + z.row(jj)*d2t;
        hnewD(jj,jj) = hD_2(jj,jj) + hz.row(jj)*t + d2z.row(jj)*d1t +  d1z.row(jj)*d2t + z.row(jj)*ht;
      }
    }
  }
  Eigen::MatrixXd nL(n1,n1), d1nL(n1,n1), d2nL(n1,n1);
  nL.setZero(); d1nL.setZero(); d2nL.setZero();
  nL = L_22* newL;
  d1nL = d1L_22* newL + L_22 * d1newL;
  d2nL = d2L_22* newL + L_22 * d2newL;
  hnewL = hL_22* newL + d2L_22 * d1newL + d1L_22 * d2newL+ L_22 * hnewL;
  
  // collect the pieces 
  struct TVBS_Matrix_two LDLT_Hess;
  LDLT_Hess.mat1 = hnewL;
  LDLT_Hess.mat2 = hnewD;
  
  return LDLT_Hess;
}


