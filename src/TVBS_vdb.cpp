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

int counter_check = 0;

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
//[[Rcpp::export]]
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


//Hessian_trunc_struct 
struct Hessian_trunc_struct TVBS_Hessian_biv_norm_trunc(Eigen::VectorXd w, double rho)
{
  //Hessian_trunc Hess; // Hessian of lambda_{0,1}, Omega_{0,0},Omega_{0,1},Omega_{1,1} as a function of w_0, w_1, rho. (indexation according to CPP to match code) 
  // structure contains five matrices, one Hessian for each element. 
  Eigen::MatrixXd H_lambda0(3,3);
  H_lambda0.setZero();
  
  Eigen::MatrixXd H_lambda1(3,3);
  H_lambda1.setZero();
  
  Eigen::MatrixXd H_Omega00(3,3);
  H_Omega00.setZero();
  
  Eigen::MatrixXd H_Omega01(3,3);
  H_Omega01.setZero();
  
  Eigen::MatrixXd H_Omega11(3,3);
  H_Omega11.setZero();
  
  
  // calculate deltas 
  double cdf_0 = std_normal_cdf((w(1) - rho*w(0))/sqrt(1-rho*rho));
  double delta_0 = std_normal_pdf(w(0))*cdf_0;
  double cdf_1 = std_normal_cdf((w(0) - rho*w(1))/sqrt(1-rho*rho));
  double delta_1 = std_normal_pdf(w(1))*cdf_1;
  
  // now for derivatives 
  // gradients of delta_0, 
  Eigen::MatrixXd grad_delta_0(3,1);
  grad_delta_0.setZero();
  
  grad_delta_0(0,0) = -std_normal_pdf(w(0))*cdf_0*w(0) - std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_0(1,0) = std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_0(2,0) = std_normal_pdf(w(0))*normal_pdf((w(1) - rho*w(0)),(1-rho*rho))*(-w(0)+(w(1)-rho*w(0))*rho/(1-rho*rho));
  
  // gradients of delta_1. 
  Eigen::MatrixXd grad_delta_1(3,1);
  grad_delta_1.setZero();
  grad_delta_1(1,0) = -std_normal_pdf(w(1))*cdf_1*w(1) - std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_1(0,0) = std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_1(2,0) = std_normal_pdf(w(1))*normal_pdf((w(0) - rho*w(1)),(1-rho*rho))*(-w(1)+(w(0)-rho*w(1))*rho/(1-rho*rho));  
  
  // Hessians 
  Eigen::MatrixXd H_delta0 = Hessian_delta(w(0),w(1),rho);
  Eigen::MatrixXd HH_delta1 = Hessian_delta(w(1),w(0),rho);
  Eigen::MatrixXd H_delta1 = HH_delta1;
  
  H_delta1(0,0)= HH_delta1(1,1);
  H_delta1(1,1)= HH_delta1(0,0);
  H_delta1(0,2) = HH_delta1(1,2);
  H_delta1(1,2) = HH_delta1(0,2);
  H_delta1(2,0) = HH_delta1(2,1);
  H_delta1(2,1) = HH_delta1(2,0);
  
  
  // now for the lambdas. 
  // calculate lambdas
  double Phi2 = biv_normal_cdf(w(0),w(1),rho);
  double lambda_0 = -(delta_0+rho*delta_1)/Phi2;
  double lambda_1 = -(delta_1+rho*delta_0)/Phi2;
  
  Eigen::VectorXd dPhi2 = grad_cdf(w(0),w(1),rho);
  Eigen::MatrixXd H_Phi2 = Hess_cdf(w(0),w(1),rho);
  
  Eigen::VectorXd drho(3);
  drho.setZero();
  drho(2)=1; 
  
  // Hessian for lambda_0 
  // three terms. 
  Eigen::VectorXd grad_lambda_0(3);
  grad_lambda_0.setZero();
  
  grad_lambda_0 = -(grad_delta_0 + rho*grad_delta_1)/Phi2 - drho*delta_1/Phi2 + (delta_0+rho*delta_1)*dPhi2/pow(Phi2,2);
  
  // first term
  H_lambda0 = -(H_delta0 + rho*H_delta1)/Phi2 + dPhi2*(grad_delta_0.transpose() + rho*grad_delta_1.transpose())/pow(Phi2,2); 
  
  H_lambda0 += -(drho*grad_delta_1.transpose())/Phi2;
  
  // second term.
  H_lambda0 += -(grad_delta_1* drho.transpose())/Phi2 + delta_1*(dPhi2 * drho.transpose())/pow(Phi2,2); 
  
  // third term 
  H_lambda0 +=  (grad_delta_0 + rho*grad_delta_1 + drho*delta_1) * dPhi2.transpose()/pow(Phi2,2);
  H_lambda0 +=  (delta_0+rho*delta_1) *(H_Phi2/pow(Phi2,2) - 2*dPhi2 * dPhi2.transpose()/pow(Phi2,3));
  
  
  
  // Hessian for lambda_1 
  // three terms. 
  Eigen::VectorXd grad_lambda_1(3);
  grad_lambda_1.setZero();
  
  grad_lambda_1 = -(grad_delta_1 + rho*grad_delta_0)/Phi2 - drho*delta_0/Phi2 + (delta_1+rho*delta_0)*dPhi2/pow(Phi2,2);
  
  
  // first term
  H_lambda1 = -(H_delta1 + rho*H_delta0)/Phi2 + dPhi2*(grad_delta_1.transpose() + rho*grad_delta_0.transpose())/pow(Phi2,2); 
  H_lambda1 += -(drho*grad_delta_0.transpose())/Phi2;
  
  // second term.
  H_lambda1 += -(grad_delta_0* drho.transpose())/Phi2 + delta_0*(dPhi2 * drho.transpose())/pow(Phi2,2); 
  
  // third term 
  H_lambda1 +=  (grad_delta_1 + rho*grad_delta_0 + drho*delta_0) * dPhi2.transpose()/pow(Phi2,2);
  H_lambda1 +=  (delta_1+rho*delta_0) *(H_Phi2/pow(Phi2,2) - 2*dPhi2 * dPhi2.transpose()/pow(Phi2,3));
  
  ///////////////////////////
  // Hessian for the Omegas 
  ///////////////////////////
  
  
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
  Eigen::MatrixXd H_phi2 = Hess_pdf(w(0),w(1),rho);
  Eigen::MatrixXd H_exp_0(3,3);
  H_exp_0.setZero();
  
  Eigen::VectorXd dw0(3);
  dw0.setZero();
  dw0(0)=1; 
  
  Eigen::VectorXd dw1(3);
  dw1.setZero();
  dw1(1)=1; 
  
  // first term 
  H_exp_0 = w(0)*H_delta0 + dw0 *grad_delta_0.transpose() + grad_delta_0 * dw0.transpose(); 
  
  // second term 
  H_exp_0 += rho*rho*(w(1)*H_delta1 + dw1 * grad_delta_1.transpose() + grad_delta_1 * dw1.transpose()); // derive delta1 + w1; twice delta1
  H_exp_0 += w(1)*2*rho*(drho * grad_delta_1.transpose() + grad_delta_1 *drho.transpose()); // derive delta1 and rho 
  H_exp_0 += w(1)*delta_1*2* drho *drho.transpose(); // twice rho 
  H_exp_0 += delta_1*2*rho*(drho * dw1.transpose() + dw1 *drho.transpose()); // derive w1 and rho
  
  // third term 
  H_exp_0 += -rho*(1-rho*rho)*H_phi2 - (dphi2*drho.transpose() + drho*dphi2.transpose())*(1-3*rho*rho); // twice phi2; phi2 and rho. 
  H_exp_0(2,2) += 3*2*rho*phi2;
  
  // gradient of exp_0. 
  Eigen::VectorXd grad_exp_0(3);
  grad_exp_0.setZero();
  for (int k=0;k<3;k++){
    grad_exp_0(k) = w(0)*grad_delta_0(k) + rho*rho*w(1)*grad_delta_1(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_0(0) += delta_0; 
  grad_exp_0(1) += rho*rho*delta_1; 
  grad_exp_0(2) += 2*rho*w(1)*delta_1 - phi2*(1-3*rho*rho);
  
  
  //double omega_00 = 1 - exp_0/Phi2 - lambda_0*lambda_0; 
  
  H_Omega00 = -H_exp_0/Phi2 + (dPhi2 * grad_exp_0.transpose() +grad_exp_0 * dPhi2.transpose() )/pow(Phi2,2) + exp_0*H_Phi2/pow(Phi2,2) -2*exp_0*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega00 += - 2*(H_lambda0 *lambda_0 + grad_lambda_0* grad_lambda_0.transpose());
  
  // Omega_{11}
  double exp_1 = w(1)*delta_1 + rho*rho*w(0)*delta_0 - rho*(1-rho*rho)*phi2;
  Eigen::MatrixXd H_exp_1(3,3);
  H_exp_1.setZero();
  
  // first term 
  H_exp_1 = w(1)*H_delta1 + dw1 *grad_delta_1.transpose() + grad_delta_1 * dw1.transpose(); 
  
  // second term 
  H_exp_1 += rho*rho*(w(0)*H_delta0 + dw0 * grad_delta_0.transpose() + grad_delta_0 * dw0.transpose()); // derive delta1 + w1; twice delta1
  H_exp_1 += w(0)*2*rho*(drho * grad_delta_0.transpose() + grad_delta_0 *drho.transpose()); // derive delta1 and rho 
  H_exp_1 += w(0)*delta_0*2* drho *drho.transpose(); // twice rho 
  H_exp_1 += delta_0*2*rho*(drho * dw0.transpose() + dw0 *drho.transpose()); // derive w1 and rho
  
  // third term 
  H_exp_1 += -rho*(1-rho*rho)*H_phi2 - (dphi2*drho.transpose() + drho*dphi2.transpose())*(1-3*rho*rho); // twice phi2; phi2 and rho. 
  H_exp_1(2,2) += 3*2*rho*phi2;
  
  // gradient of exp_1. 
  Eigen::VectorXd grad_exp_1(3);
  grad_exp_1.setZero();
  for (int k=0;k<3;k++){
    grad_exp_1(k) = w(1)*grad_delta_1(k) + rho*rho*w(0)*grad_delta_0(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_1(0) += rho*rho*delta_0; 
  grad_exp_1(1) += delta_1; 
  grad_exp_1(2) += 2*rho*w(0)*delta_0 - phi2*(1-3*rho*rho);
  
  //double omega_11 = 1 - exp_1/Phi2 - lambda_1*lambda_1; 
  
  H_Omega11 = -H_exp_1/Phi2 + (dPhi2 * grad_exp_1.transpose() +grad_exp_1 * dPhi2.transpose() )/pow(Phi2,2) + exp_1*H_Phi2/pow(Phi2,2) -2*exp_1*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega11 += - 2*(H_lambda1 *lambda_1 + grad_lambda_1* grad_lambda_1.transpose());
  
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
  
  
  Eigen::MatrixXd H_exp_01(3,3);
  H_exp_01.setZero();
  
  // first term 
  H_exp_01 = w(1)*(drho * grad_delta_1.transpose() + grad_delta_1 * drho.transpose()) +  delta_1*(drho*dw1.transpose() + dw1 * drho.transpose());
  H_exp_01 += rho*(dw1 * grad_delta_1.transpose() + grad_delta_1 * dw1.transpose());
  H_exp_01 += w(1)*rho*H_delta1; 
  
  // second term 
  H_exp_01 += w(0)*(drho * grad_delta_0.transpose() + grad_delta_0 * drho.transpose()) +  delta_0*(drho*dw0.transpose() + dw0 * drho.transpose());
  H_exp_01 += rho*(dw0 * grad_delta_0.transpose() + grad_delta_0 * dw0.transpose());
  H_exp_01 += w(0)*rho*H_delta0; 
  
  // third term 
  H_exp_01 += 2*rho*(dphi2*drho.transpose() + drho*dphi2.transpose()) + (rho*rho-1)*H_phi2; 
  H_exp_01(2,2) += 2*phi2;
  
  //double omega_01 = rho - exp_01/Phi2 - lambda_0*lambda_1; 
  
  // finally the last Hessian
  H_Omega01 = -H_exp_01/Phi2 + (dPhi2 * grad_exp_01.transpose() +grad_exp_01 * dPhi2.transpose() )/pow(Phi2,2) + exp_01*H_Phi2/pow(Phi2,2) -2*exp_01*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega01 += - H_lambda0 *lambda_1 - grad_lambda_0* grad_lambda_1.transpose() - grad_lambda_1* grad_lambda_0.transpose() - lambda_0 * H_lambda1;
  
  
  // return matrices. 
  struct Hessian_trunc_struct output_final = {H_lambda0,H_lambda1, H_Omega00,H_Omega01,H_Omega11};
  return output_final;
}


//Hessian_trunc_struct 
Eigen::MatrixXd TVBS_Hessian_biv_norm_trunc_R(Eigen::VectorXd w, double rho, int out)
{
  Eigen::MatrixXd mat_out(3,3);
  
  struct Hessian_trunc_struct output_final = TVBS_Hessian_biv_norm_trunc(w, rho);
  
  mat_out = output_final.mat1;
  if (out==2){
    mat_out = output_final.mat2;
  }
  if (out==3){
    mat_out = output_final.mat3;
  }
  if (out==4){
    mat_out = output_final.mat4;
  }
  if (out==5){
    mat_out = output_final.mat5;
  }
  
  return mat_out;
}


Eigen::MatrixXd TVBS_Jacobian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma)
{
  Eigen::MatrixXd Jacob(5,7);
  Eigen::MatrixXd Jacob_norm(3,7); // Jacobian of normalizing transformation. 
  Jacob_norm.setZero(); 
  
  Eigen::MatrixXd D(2,2);
  D.setZero();
  
  D(0,0)= sqrt(Sigma(0,0));
  D(1,1)= sqrt(Sigma(1,1));
  double rho = Sigma(0,1)/(D(0,0)*D(1,1));
  
  Eigen::VectorXd tildex = D.inverse() *(x-alpha); 
  
  // w(0) = sigma(0,0)^(-1)*(x(0)-alpha(0))
  Jacob_norm(0,0) = 1/D(0,0); 
  Jacob_norm(0,2) = -Jacob_norm(0,0);
  Jacob_norm(0,4) =- tildex(0)/(2*Sigma(0,0)); 
  
  // w(1) = sigma(1,1)^(-1)*(x(1)-alpha(1))
  Jacob_norm(1,1) = 1/D(1,1); 
  Jacob_norm(1,3) = -Jacob_norm(1,1);
  Jacob_norm(1,6) = -tildex(1)/(2*Sigma(1,1)); 
  
  // rho = Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1))
  Jacob_norm(2,4) = -rho/(2*Sigma(0,0));
  Jacob_norm(2,5) = 1/(D(0,0)*D(1,1));
  Jacob_norm(2,6) = -rho/(2*Sigma(1,1));
  
  // calculate using the normalized values 
  Eigen::VectorXd lambda_Omega = TVBS_biv_norm_trunc(tildex,rho); 
  Eigen::MatrixXd Jacob_lambda_Omega = TVBS_Jacobian_biv_norm_trunc(tildex,rho);
  
  // rescale 
  Eigen::MatrixXd Scale(5,5);
  Scale.setZero(); 
  Scale(0,0) = D(0,0);
  Scale(1,1)= D(1,1);
  Scale(2,2)= Sigma(0,0);
  Scale(3,3) = D(0,0)*D(1,1);
  Scale(4,4) = Sigma(1,1); 
  
  Jacob = Scale * Jacob_lambda_Omega * Jacob_norm; 
  Jacob(0,4) += lambda_Omega(0)/(2*D(0,0));
  Jacob(1,6) += lambda_Omega(1)/(2*D(1,1));
  Jacob(2,4) += lambda_Omega(2);
  Jacob(3,4) += lambda_Omega(3)*(D(1,1)/(2*D(0,0)));
  Jacob(3,6) += lambda_Omega(3)*(D(0,0)/(2*D(1,1)));
  Jacob(4,6) += lambda_Omega(4);
  
  Jacob(0,2) += 1;
  Jacob(1,3) += 1; 
  
  return Jacob;
  
}

struct Hessian_trunc_struct TVBS_Hessian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma)
{
  Eigen::MatrixXd nH_lambda0(5,5), nH_lambda1(5,5), nH_Omega00(5,5), nH_Omega01(5,5), nH_Omega11(5,5); // normalized matrices 
  Eigen::MatrixXd H_lambda0(5,5), H_lambda1(5,5), H_Omega00(5,5), H_Omega01(5,5), H_Omega11(5,5); // general including scaling.
  
  Eigen::MatrixXd Jacob(5,7);
  Eigen::MatrixXd Jacob_norm(3,7); // Jacobian of normalizing transformation. 
  Jacob_norm.setZero(); 
  
  Eigen::MatrixXd D(2,2);
  D.setZero();
  
  D(0,0)= sqrt(Sigma(0,0));
  D(1,1)= sqrt(Sigma(1,1));
  double rho = Sigma(0,1)/(D(0,0)*D(1,1));
  
  Eigen::VectorXd tildex = D.inverse() *(x-alpha); 
  
  // w(0) = sigma(0,0)^(-1)*(x(0)-alpha(0))
  Jacob_norm(0,0) = 1/D(0,0); 
  Jacob_norm(0,2) = -Jacob_norm(0,0);
  Jacob_norm(0,4) =- tildex(0)/(2*Sigma(0,0)); 
  
  // w(1) = sigma(1,1)^(-1)*(x(1)-alpha(1))
  Jacob_norm(1,1) = 1/D(1,1); 
  Jacob_norm(1,3) = -Jacob_norm(1,1);
  Jacob_norm(1,6) = -tildex(1)/(2*Sigma(1,1)); 
  
  // rho = Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1))
  Jacob_norm(2,4) = -rho/(2*Sigma(0,0));
  Jacob_norm(2,5) = 1/(D(0,0)*D(1,1));
  Jacob_norm(2,6) = -rho/(2*Sigma(1,1));
  
  // calculate using the normalized values 
  Eigen::VectorXd lambda_Omega = TVBS_biv_norm_trunc(tildex,rho); 
  // Jacobian of lambdaj,Omegaij w.r.t normalized w, rho.: 5x3. 
  Eigen::MatrixXd Jacob_lambda_Omega = TVBS_Jacobian_biv_norm_trunc(tildex,rho);
  
  // Jacobian without rescaling for easier Hessian calculation: 5x 7.  
  Jacob = Jacob_lambda_Omega * Jacob_norm; 
  
  // H_fmat: Hessian for normalized quantities: 5 matrices of size 3x3 each. 
  struct Hessian_trunc_struct H_fmat = TVBS_Hessian_biv_norm_trunc(tildex, rho);
  
  // Hessian for normalisation of each entry separately. 
  Eigen::MatrixXd HT_1(7,7), HT_2(7,7), HT_3(7,7);
  HT_1.setZero();
  HT_2.setZero();
  HT_3.setZero();
  
  // HT_1 Hessian of w0 w.r.t. x,alpha, Sigma. 
  HT_1(0,4)= -1/(2*sqrt(Sigma(0,0))*Sigma(0,0));
  HT_1(2,4)= 1/(2*sqrt(Sigma(0,0))*Sigma(0,0));
  HT_1(4,4)= 3*tildex(0)/(4*Sigma(0,0)*Sigma(0,0));
  HT_1(4,0) = -0.5/(Sigma(0,0)*sqrt(Sigma(0,0)));
  HT_1(4,2) = 0.5/(Sigma(0,0)*sqrt(Sigma(0,0)));
  
  // HT_2 Hessian of w1 w.r.t. x,alpha, Sigma. 
  HT_2(1,6)= -1/(2*sqrt(Sigma(1,1))*Sigma(1,1));
  HT_2(3,6)= 1/(2*sqrt(Sigma(1,1))*Sigma(1,1));
  HT_2(6,6)= 3*tildex(1)/(4*Sigma(1,1)*Sigma(1,1));
  HT_2(6,1) = -0.5/(Sigma(1,1)*sqrt(Sigma(1,1)));
  HT_2(6,3) = 0.5/(Sigma(1,1)*sqrt(Sigma(1,1)));
  
  // HT_3 Hessian of rho w.r.t. x,alpha, Sigma. 
  HT_3(4,4)= 3*rho/(4*Sigma(0,0)*Sigma(0,0));
  HT_3(5,4) = -0.5/(Sigma(0,0)*sqrt(Sigma(0,0)*Sigma(1,1)));
  HT_3(6,4)= rho/(4*Sigma(0,0)*Sigma(1,1));
  HT_3(4,5) = HT_3(5,4);
  
  HT_3(6,6)= 3*rho/(4*Sigma(1,1)*Sigma(1,1));
  HT_3(5,6) = -0.5/(Sigma(1,1)*sqrt(Sigma(0,0)*Sigma(1,1)));
  HT_3(4,6)= rho/(4*Sigma(0,0)*Sigma(1,1));
  HT_3(6,5) = HT_3(5,6);

  // Hessian without rescaling:
  // Hessian for lambda_0
  nH_lambda0 = Jacob_norm.transpose() * H_fmat.mat1 * Jacob_norm;
  nH_lambda0 += HT_1* Jacob_lambda_Omega(0,0) +  HT_2* Jacob_lambda_Omega(0,1) +  HT_3* Jacob_lambda_Omega(0,2);
  
  // Hessian for lambda_1
  nH_lambda1 = Jacob_norm.transpose() * H_fmat.mat2 * Jacob_norm;
  nH_lambda1 += HT_1* Jacob_lambda_Omega(1,0) +  HT_2* Jacob_lambda_Omega(1,1) +  HT_3* Jacob_lambda_Omega(1,2);
  
  // Hessian for Omega_00
  nH_Omega00 = Jacob_norm.transpose() * H_fmat.mat3 * Jacob_norm;
  nH_Omega00 += HT_1* Jacob_lambda_Omega(2,0) +  HT_2* Jacob_lambda_Omega(2,1) +  HT_3* Jacob_lambda_Omega(2,2);
  
  // Hessian for Omega_01
  nH_Omega01 = Jacob_norm.transpose() * H_fmat.mat4 * Jacob_norm;
  nH_Omega01 += HT_1* Jacob_lambda_Omega(3,0) +  HT_2* Jacob_lambda_Omega(3,1) +  HT_3* Jacob_lambda_Omega(3,2);
  
  // Hessian for Omega_11
  nH_Omega11 = Jacob_norm.transpose() * H_fmat.mat5 * Jacob_norm;
  nH_Omega11 += HT_1* Jacob_lambda_Omega(4,0) +  HT_2* Jacob_lambda_Omega(4,1) +  HT_3* Jacob_lambda_Omega(4,2);
  
  // rescale 
  Eigen::MatrixXd DS_1(7,1),DS_2(7,1),DS_3(7,1),DS_4(7,1),DS_5(7,1);
  DS_1.setZero();
  DS_2.setZero();
  DS_3.setZero();
  DS_4.setZero();
  DS_5.setZero();
  
  DS_1(4,0) = 0.5/sqrt(Sigma(0,0));
  DS_2(6,0)= 0.5/sqrt(Sigma(1,1));
  DS_3(4,0) = 1;
  DS_4(4,0) = 0.5*sqrt(Sigma(1,1)/Sigma(0,0));
  DS_4(6,0) = 0.5*sqrt(Sigma(0,0)/Sigma(1,1));
  DS_5(6,0)= 1;
  
  Eigen::MatrixXd HS_1(7,7),HS_2(7,7),HS_3(7,7),HS_4(7,7),HS_5(7,7);
  HS_1.setZero();
  HS_2.setZero();
  HS_3.setZero();
  HS_4.setZero();
  HS_5.setZero();
  
  HS_1(4,4) = -DS_1(4,0)/(2*Sigma(0,0));
  HS_2(6,6) = -DS_2(6,0)/(2*Sigma(1,1));
  
  HS_4(4,4) = -DS_4(4,0)/(2*Sigma(0,0));
  HS_4(6,6) = -DS_4(6,0)/(2*Sigma(1,1));
  HS_4(4,6) = 0.25/sqrt(Sigma(0,0)*Sigma(1,1));
  HS_4(6,4) = 0.25/sqrt(Sigma(0,0)*Sigma(1,1));
  
  // HS_5 and HS_3 are zero. 
  H_lambda0 = lambda_Omega(0)*HS_1 + DS_1 * Jacob.block(0,0,1,7)  + Jacob.block(0,0,1,7).transpose() * DS_1.transpose() + nH_lambda0 *sqrt(Sigma(0,0));
  H_lambda1 = lambda_Omega(1)*HS_2 + DS_2 * Jacob.block(1,0,1,7)  + Jacob.block(1,0,1,7).transpose() * DS_2.transpose() + nH_lambda1 *sqrt(Sigma(1,1));
  
  H_Omega00 = lambda_Omega(2)*HS_3 + DS_3 * Jacob.block(2,0,1,7)  + Jacob.block(2,0,1,7).transpose() * DS_3.transpose() + nH_Omega00 * (Sigma(0,0));
  H_Omega01 = lambda_Omega(3)*HS_4 + DS_4 * Jacob.block(3,0,1,7)  + Jacob.block(3,0,1,7).transpose() * DS_4.transpose() + nH_Omega01 * sqrt(Sigma(0,0)*Sigma(1,1));
  H_Omega11 = lambda_Omega(4)*HS_5 + DS_5 * Jacob.block(4,0,1,7)  + Jacob.block(4,0,1,7).transpose() * DS_5.transpose() + nH_Omega11 * (Sigma(1,1));
  
  
  // return matrices. 
  struct Hessian_trunc_struct output_final = {H_lambda0,H_lambda1, H_Omega00,H_Omega01,H_Omega11};
  return output_final;
  
}

Eigen::MatrixXd TVBS_Hessian_trunc_gen_R(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma, int out)
{
  struct Hessian_trunc_struct  output = TVBS_Hessian_biv_gen_trunc(x, alpha, Sigma);
  
  Eigen::MatrixXd outp = output.mat1;
  if (out==2){ outp = output.mat2;}
  if (out==3){ outp = output.mat3;}
  if (out==4){ outp = output.mat4;}
  if (out==5){ outp = output.mat5;}
  
  return outp;
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
    Eigen::VectorXd temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k));
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



Eigen::MatrixXd TVBS_g_a_omega_b_cpp(Eigen::MatrixXd x1, Eigen::MatrixXd x2, int omsymmetric = 0, int omdiagonal = 0) 
{
  int k = x1.cols();
  
  Eigen::MatrixXd g = kroneckerProduct(x1.transpose(), x2);
  int nrows = g.rows();
  int ncols = g.cols();
  
  if(omsymmetric_global==1)
  {
    Eigen::MatrixXd temp1 = TVBS_vec_symmetry_cpp(k);
    g = temp1*g;
    nrows = g.rows();
    ncols = g.cols();
    
    if(omdiagonal_global==1)
    {
      Eigen::VectorXd temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k));
      Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix);
      nrows = row_select.size();
      for(int j = 0; j<nrows; j++)
      {
        g.row(j) = g.row(row_select(j));
      }
    }
  }
  
  Eigen::MatrixXd g_out = g.block(0,0,nrows,ncols);
  return g_out;
}



Eigen::MatrixXd TVBS_g_cholesky_cov_cpp(Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd sigma_diag_sqrt = Sigma.diagonal().array().sqrt();
  
  int parm_dim = sigma_diag_sqrt.size();
  int parm_cov = parm_dim*(parm_dim+1)/2;
  
  // cholesky decomposition t(L)%*%L=Sigma
  Eigen::MatrixXd Sigma_star_chol = Sigma.llt().matrixU();
  
  Eigen::VectorXd Sigma_star_chol_diag = Sigma_star_chol.diagonal();
  
  Eigen::MatrixXd dg(parm_cov,parm_cov);
  dg.setZero();
  
  Eigen::VectorXd Sigma_star_chol_vech = TVBS_lower_tri_entries_cpp(Sigma_star_chol.transpose());
  
  int j = 1;
  int mm = 1;
  
  if(parm_dim>1)
  {
    while(j<=parm_cov)
    {
      int gss_c = j;
      int l = 1;
      while(gss_c<=parm_cov)
      {
        dg.block(j+l-1-1,gss_c-1,(j+parm_dim-mm) - (j+l-1) + 1, (gss_c+parm_dim-l+1-mm) - gss_c + 1) = (Eigen::MatrixXd::Identity(parm_dim-l+2-mm, parm_dim-l+2-mm).array()*Sigma_star_chol_vech(j+l-1-1)).matrix();
        dg.block(j+l-1-1,gss_c-1,1,(gss_c+parm_dim-l+1-mm) - gss_c + 1) += (Sigma_star_chol_vech.segment(j+l-1-1, (gss_c+parm_dim-l+1-mm) - gss_c + 1)).transpose();
        gss_c += parm_dim-l+1-mm+1;
        l++;
      }
      j += parm_dim-mm+1;
      mm++;
    }
  }
  if(parm_dim==1)
  {
    dg.setOnes();
    dg = dg.block(0,0,1,1);
  }
  
  return dg;
}




Eigen::MatrixXd TVBS_g_cholesky_cor_cpp(Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd sigma_diag_sqrt = Sigma.diagonal().array().sqrt();
  
  int parm_dim = sigma_diag_sqrt.size();
  int parm_cor = parm_dim*(parm_dim-1)/2;
  
  // cholesky decomposition t(L)%*%L=Sigma
  Eigen::MatrixXd Sigma_star_chol = Sigma.llt().matrixU();
  Eigen::VectorXd Sigma_star_chol_diag = Sigma_star_chol.diagonal();
  
  Eigen::MatrixXd dg(parm_cor,parm_cor);
  dg.setZero();
  
  
  int j = 1;
  int mm = 1;
  
  if(parm_dim>2)
  {
    for(mm; mm<=parm_dim; mm++)
    {
      Eigen::MatrixXd ss, ss1;
      int gss_c = j;
      int l = 1;
      if(j<parm_cor)
      {
        if(parm_dim-mm-1>1)
        {
          Eigen::MatrixXd ss_lhs(parm_dim-mm-1,parm_dim-mm-1), ss_rhs(parm_dim-mm-1,parm_dim-mm-1), temp_matrix_ones(parm_dim-mm-1,parm_dim-mm-1);
          ss_lhs.setZero();
          ss_rhs.setZero();
          temp_matrix_ones.setOnes();
          ss_lhs.triangularView<Eigen::Upper>() = temp_matrix_ones.triangularView<Eigen::Upper>();
          Eigen::MatrixXd ss_rhs_diag = (Sigma_star_chol.block(mm-1, mm+2-1, 1, parm_dim - (mm+2)+1));
          ss_rhs.diagonal() = ss_rhs_diag;
          ss = ss_lhs*ss_rhs;
          
          
          Eigen::MatrixXd ss1_1((parm_dim-1)-(mm+1)+1, (parm_dim-1)-(mm+1)+1);
          ss1_1.setZero();
          ss1_1.diagonal() = (1.0/Sigma_star_chol_diag.segment(mm+1-1, (parm_dim-1)-(mm+1)+1).array());
          
          Eigen::MatrixXd ss1_2_vec = Sigma_star_chol.block(mm-1, mm+1-1, 1, parm_dim-1 - (mm+1)+1);
          Eigen::MatrixXd ss1_2(parm_dim-1 - (mm+1)+1, parm_dim-1 - (mm+1)+1);
          ss1_2.setZero();
          ss1_2.diagonal() = ss1_2_vec;
          
          Eigen::MatrixXd ss1_3_full = Sigma_star_chol.block(mm+1-1, mm+2-1, parm_dim-1-(mm+1)+1, parm_dim - (mm+2)+1);
          Eigen::MatrixXd ss1_3(parm_dim-1-(mm+1)+1, parm_dim-1-(mm+1)+1);
          ss1_3.setZero();
          ss1_3.triangularView<Eigen::Upper>() = ss1_3_full.triangularView<Eigen::Upper>();
          
          ss1 = ss - ss1_1*ss1_2*ss1_3;
        }
        if(parm_dim-mm-1==1)
        {
          ss = Sigma_star_chol.block(mm-1, (mm+2)-1,1,1);
          ss1 = (ss.array() - (1.0/Sigma_star_chol_diag((mm+1)-1))*Sigma_star_chol(mm-1, (mm+1)-1)*Sigma_star_chol((mm+1)-1, (mm+2)-1)).matrix();
        }
      }
      while(gss_c<=parm_cor)
      {
        dg.block((j+l-1)-1, gss_c-1, (j+parm_dim-mm-1)-(j+l-1)+1, (gss_c+parm_dim-l-mm)-gss_c+1) = Eigen::MatrixXd::Identity(parm_dim-l+1-mm, parm_dim-l+1-mm)*Sigma_star_chol(mm-1,(mm+l-1)-1);
        
        if(gss_c<parm_cor)
        {
          dg.block((j+l-1)-1, (gss_c+parm_dim-l-mm+1)-1, 1, (gss_c+2*parm_dim-2*l-2*mm)-(gss_c+parm_dim-l-mm+1)+1) = ss1.block(l-1, l-1, 1, ss.cols()-l+1);
        }
        gss_c = gss_c+parm_dim-l+1-mm;
        l++;
        
      }
      j += parm_dim-mm;
    }
  }
  if(parm_dim==2)
  {
    dg.setOnes();
    dg = dg.block(0,0,1,1);
  }
  
  return dg;
}



// input: x_norm, Sigma_diag_sqrt, Cor_mat
struct TVBS_Matrix_two TVBS_grad_cor_cov_cpp(Eigen::VectorXd x_norm, Eigen::VectorXd Sigma_diag_sqrt, Eigen::MatrixXd Cor_mat)
{
  // dimension of the matrix/problem
  int dim_lcl = Sigma_diag_sqrt.size();
  
  // number of lower diagonal entries (with diagonal)
  int row_lcl = dim_lcl*(dim_lcl+1)/2;
  
  // number of lower diagonal entries (without diagonal)
  int col_lcl = dim_lcl*(dim_lcl-1)/2;
  
  Eigen::VectorXd nu_1_lcl = Sigma_diag_sqrt.array()*Sigma_diag_sqrt.array();
  Eigen::VectorXd nu_2_lcl = Sigma_diag_sqrt;
  
  Eigen::MatrixXd temp(dim_lcl, dim_lcl);
  temp.setZero();
  temp.diagonal() = Sigma_diag_sqrt;
  Eigen::MatrixXd Cov_mat = temp*Cor_mat*temp;
  
  Eigen::MatrixXd grm_lcl(row_lcl, dim_lcl);
  grm_lcl.setZero();
  
  Eigen::VectorXd cc1_lcl = 0.5*(x_norm.array()/nu_1_lcl.array());
  int l_lcl = 0;
  
  for(int idimvn=1; idimvn<=dim_lcl; idimvn++)
  {
    grm_lcl((l_lcl+1)-1,idimvn-1) = cc1_lcl(idimvn-1);
    l_lcl += dim_lcl+1-idimvn;
  }
  
  Eigen::MatrixXd grc_lcl(row_lcl, col_lcl);
  grc_lcl.setZero();
  
  int j = 1;
  int mm = 1;
  int c1 = 1;
  
  while(j != row_lcl)
  {
    Eigen::VectorXd n_diag_1_vec = 1.0/(nu_2_lcl(mm-1)*nu_2_lcl.segment((mm+1)-1, dim_lcl-(mm+1)+1).array());
    Eigen::MatrixXd n_diag_1(dim_lcl-mm, dim_lcl-mm);
    n_diag_1.setZero();
    n_diag_1.diagonal() = n_diag_1_vec;
    
    Eigen::MatrixXd Cor_mat_block = Cor_mat.block(mm-1, (mm+1)-1, 1, dim_lcl-(mm+1)+1);
    Eigen::MatrixXd n_hor_1 = (-0.5*Cor_mat_block.array()/Cov_mat(mm-1,mm-1));
    
    Cor_mat_block = Cor_mat.block(mm-1, (mm+1)-1, 1, dim_lcl-(mm+1)+1);
    Eigen::MatrixXd n_hor_2 = -0.5*Cor_mat_block.array()/(nu_1_lcl.segment((mm+1)-1, dim_lcl-(mm+1)+1)).array();
    
    grc_lcl.block((j+1)-1, c1-1, (j+dim_lcl-mm)-(j+1)+1, (c1+dim_lcl-mm-1)-c1+1) = n_diag_1;
    grc_lcl.block(j-1, c1-1, 1, (c1+dim_lcl-1-mm)-c1+1) = n_hor_1;
    
    int k = j;
    int l = 1;
    while(l!=(dim_lcl+1-mm))
    {
      grc_lcl((k+dim_lcl+2-l-mm)-1, (l+c1-1)-1) = n_hor_2(l-1,0);
      k += dim_lcl+2-l-mm;
      l++;
    }
    c1 += dim_lcl-mm;
    j += dim_lcl+1-mm;
    mm++;
  }
  
  grm_lcl = -1.0*grm_lcl;
  struct TVBS_Matrix_two output = {grm_lcl, grc_lcl};
  
  return output;
}



// issues arrise for some cases due to divergence to behaviour compared to R code (when tempselmat length is not equal to number of rows/columns to select from)
struct TVBS_Matrix_two TVBS_g_both_x_omega_x_cpp(Eigen::MatrixXd Mat_1, Eigen::MatrixXd Mat_2)
{
  // input dimensions
  int n = Mat_1.rows();
  int k = Mat_1.rows();
  
  Eigen::MatrixXd temp = Mat_2.transpose()*Mat_1.transpose();
  
  Eigen::MatrixXd t1(n*temp.rows(), n*n);
  t1.setZero();
  
  for(int i1 = 0; i1 < n; i1++)
  {
    t1.block(0,i1*n,n*temp.rows(),n) = kroneckerProduct(Eigen::MatrixXd::Identity(n, n), temp.col(i1));
  }
  
  Eigen::MatrixXd g_Mat1_cov = kroneckerProduct(Eigen::MatrixXd::Identity(n, n), Mat_2*Mat_1.transpose()) + t1;
  Eigen::MatrixXd g_Mat2_cov = kroneckerProduct(Mat_1.transpose(), Mat_1.transpose());
  
  if(x2symmetric_global==1)
  {
    Eigen::MatrixXd temp_sel_matrix(n,n);
    temp_sel_matrix.setZero();
    Eigen::MatrixXd temp_sel_matrix_ones(n,n);
    temp_sel_matrix_ones.setOnes();
    temp_sel_matrix.triangularView<Eigen::Lower>() = temp_sel_matrix_ones.triangularView<Eigen::Lower>();
    temp_sel_matrix.resize(n*n,1);
    Eigen::VectorXi col_select = TVBS_logic2position_cpp(temp_sel_matrix);
    
    // sub-matrix by columns
    {
      Eigen::MatrixXd g_Mat1_cov_new(g_Mat1_cov.rows(), col_select.size());
      g_Mat1_cov_new.setZero();
      for(int i1=0; i1<col_select.size(); i1++ )
      {
        g_Mat1_cov_new.col(i1) = g_Mat1_cov.col(col_select(i1));
      }
      g_Mat1_cov = g_Mat1_cov_new;
    }
    
    g_Mat2_cov = TVBS_g_asym_to_sym_cpp(g_Mat2_cov);
    
    temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k));
    Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix);
    
    if(x1symmetric_global==1)
    {
      Eigen::MatrixXd temp1 = TVBS_vec_symmetry_cpp(n);;
      g_Mat1_cov = temp1*g_Mat1_cov;
    }
    if(x2diagonal_global==0)
    {
      if(x1diagonal_global==1)
      {
        // sub-matrix by rows
         {
          Eigen::MatrixXd g_Mat1_cov_new(row_select.size(), g_Mat1_cov.cols());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.row(i1) = g_Mat1_cov.row(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
      }
    }
    if(x2diagonal_global==1)
    {
      // sub-matrix by row
      {
        Eigen::MatrixXd g_Mat2_cov_new(row_select.size(), g_Mat2_cov.cols());
        g_Mat2_cov_new.setZero();
        for(int i1=0; i1<row_select.size(); i1++)
        {
          g_Mat2_cov_new.row(i1) = g_Mat2_cov.row(row_select(i1));
        }
        g_Mat2_cov = g_Mat2_cov_new;
      }
      
      if(x1diagonal_global==1)
      {
        // sub-matrix by row
       {
          Eigen::MatrixXd g_Mat1_cov_new(row_select.size(), g_Mat1_cov.cols());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.row(i1) = g_Mat1_cov.row(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
        
        // sub-matrix by column
        {
          Eigen::MatrixXd g_Mat1_cov_new(g_Mat1_cov.rows(), row_select.size());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.col(i1) = g_Mat1_cov.col(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
        
        // sub-matrix by column
        {
          Eigen::MatrixXd g_Mat2_cov_new(g_Mat2_cov.rows(), row_select.size());
          g_Mat2_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat2_cov_new.col(i1) = g_Mat2_cov.col(row_select(i1));
          }
          g_Mat2_cov = g_Mat2_cov_new;
        }
        
      }
    }
    
    if( (x2correlation_global==1) && (x2diagonal_global==1) )
    {
      Eigen::MatrixXd g_Mat2_cov_new = g_Mat2_cov.block(0,0,1,1);
      g_Mat2_cov_new.setZero();
      g_Mat2_cov = g_Mat2_cov_new;
    }
    if( (x2correlation_global==1) && (x2diagonal_global==0) )
    {
      Eigen::MatrixXd temp_sel_matrix_inv = 1.0-temp_sel_matrix.array();
      Eigen::VectorXi row_select_inv = TVBS_logic2position_cpp(temp_sel_matrix_inv);
      
      // sub-matrix by not-row
      {
        Eigen::MatrixXd g_Mat2_cov_new(row_select_inv.size(), g_Mat2_cov.cols());
        g_Mat2_cov_new.setZero();
        for(int i1=0; i1<row_select_inv.size(); i1++)
        {
          g_Mat2_cov_new.row(i1) = g_Mat2_cov.row(row_select_inv(i1));
        }
        g_Mat2_cov = g_Mat2_cov_new;
      }
      
    }
    
  }
  struct TVBS_Matrix_two output = {g_Mat1_cov, g_Mat2_cov};
  
  return output;
}



struct TVBS_Matrix_two TVBS_g_cond_cov_cpp(Eigen::MatrixXd Id_mat, Eigen::MatrixXd Sigma)
{
  // initialize output
  Eigen::MatrixXd g_g_y, g_g_x;
  
  int dim1 = Id_mat.rows();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  Eigen::MatrixXd Sigma_11 = Sigma.block(0,0,dimdiff,dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, dim2-dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_22 = Sigma.block(dimdiff, dimdiff, dim2-dimdiff, dim2-dimdiff);
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  Eigen::MatrixXd Sigma_22_condcov = Sigma_22 - Sigma_12*Sigma_11_inv*Sigma_12.transpose();
  
  // update global variables
  x2symmetric_global     = 1;
  x2diagonal_global      = 0;
  x1symmetric_global     = 1;
  x1diagonal_global      = 1; 
  x2correlation_global   = 0;
  
  struct TVBS_Matrix_two output = TVBS_g_both_x_omega_x_cpp(Id_mat, Sigma_22_condcov);
  g_g_y = output.mat1;
  Eigen::MatrixXd g_g_2 = output.mat2;
  
  // update global variables
  x2symmetric_global     = 1;
  x2diagonal_global      = 0;
  x1symmetric_global     = 0;
  x1diagonal_global      = 0;
  x2correlation_global   = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_12, Sigma_11_inv);
  Eigen::MatrixXd g_g_22 = output.mat1;
  Eigen::MatrixXd g_g_23 = output.mat2;
  g_g_22 = -1.0*g_g_22*g_g_2;
  
  
  Eigen::MatrixXd indic_1(dim1*dimdiff,1);
  indic_1 = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
  indic_1.resize(dimdiff,dim1);
  Eigen::MatrixXd indic_1_t = indic_1.transpose();
  indic_1_t.resize(dim1*dimdiff,1);
  
  Eigen::MatrixXd g_g_22_new(indic_1_t.rows(), g_g_22.cols());
  for(int i1 = 0; i1<indic_1_t.rows(); i1++)
  {
    g_g_22_new.row(i1) = g_g_22.row(indic_1_t(i1));
  }
  
  
  if(condcov_global==1)
  {
    // update global variables
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
    
    Eigen::MatrixXd mat_n_seq(ddcov,1);
    mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq);
    
    Eigen::MatrixXd indic_2_2_temp = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcov, dimdiff*dim1+ddcov-1);
    indic_2_2_temp.resize(dim1,dimdiff);
    Eigen::MatrixXd indic_2_2 = indic_2_2_temp.transpose();
    
    Eigen::MatrixXd indic_2(dimdiff, dim1+dimdiff);
    indic_2.block(0, 0, dimdiff, dimdiff) = mat_d_up;
    indic_2.block(0, dimdiff, dimdiff, dim1) = indic_2_2;
    
    mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, ddcov+dimdiff*dim1+dd1cov-1);
    
    mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq);
    
    
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up.cols());
    indic_3.setZero();
    indic_3.block(0,dimdiff,dim1,mat_d_up.cols()) = mat_d_up;
    
    Eigen::MatrixXd indic_4(dimdiff+dim1, dim1+dimdiff);
    indic_4.block(0, 0, dimdiff, dim1+dimdiff) = indic_2;
    indic_4.block(dimdiff, 0, dim1, dim1+dimdiff) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose());
    
    g_g_23 = -1.0*TVBS_g_inverse_cpp(Sigma_11)*g_g_23*g_g_2;
    Eigen::MatrixXd g_g_x_n1(g_g_23.rows()+g_g_22_new.rows()+g_g_2.rows(), g_g_23.cols());
    g_g_x_n1.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
    g_g_x_n1.block(g_g_23.rows(), 0, g_g_22_new.rows(), g_g_23.cols()) = g_g_22_new;
    g_g_x_n1.block(g_g_23.rows()+g_g_22_new.rows(), 0, g_g_2.rows(), g_g_23.cols()) = g_g_2;
    
    Eigen::MatrixXd g_g_x_n2(g_g_23.rows()+g_g_22_new.rows()+g_g_2.rows(), g_g_23.cols());
    for(int i_row = 0; i_row < indic.size(); i_row++)
    {
      g_g_x_n2.row(i_row) = g_g_x_n1.row(indic(i_row));
    }
    g_g_x = g_g_x_n2;
    
  }
  
  // INDEX MATRITZEN MUESSEN VORHER DEFINIERT WERDEN UND DANN IN GROESSE UND INHALT ANGEPASST WERDEN!!!
  
  
  if(condcov_global==0)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
    
    // define indic matrices
    Eigen::MatrixXd indic_2, indic_3;
    
    if(ddcor == 0)
    {
      Eigen::MatrixXd indic_2_rhs = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_rhs.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, dim1+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, dimdiff, dim1) = indic_2_rhs.transpose();
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      //int indic_2_lhs_dim = (1+sqrt(1+8*ddcor))/2;
      Eigen::MatrixXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor, 0, ddcor-1);
      Eigen::MatrixXd indic_2_lhs = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_2_rhs = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_rhs.resize(dim1, dimdiff);
      
      Eigen::MatrixXd indic_2_int(dimdiff, indic_2_lhs.cols()+dim1);
      indic_2_int.block(0, 0, indic_2_lhs.rows(), indic_2_lhs.cols()) = indic_2_lhs;
      indic_2_int.block(0, indic_2_lhs.cols(), dimdiff, dim1) = indic_2_rhs.transpose();
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1,dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      //int indic_2_lhs_dim = (1+sqrt(1+8*ddcor))/2;
      Eigen::MatrixXd mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+dd1cor-1);
      Eigen::MatrixXd indic_3_rhs = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+indic_3_rhs.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, dim1, indic_3_rhs.cols()) = indic_3_rhs;
      indic_3 = indic_3_int;
      
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0,0,indic_2.rows(),indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(),0,indic_3.rows(),indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(), 0);
    
    
    Eigen::VectorXd g_g_2_new_indic_logic = 1.0-TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(dim1, dim1)).array();
    Eigen::VectorXi g_g_2_new_indic = TVBS_logic2position_cpp(g_g_2_new_indic_logic);
    
    Eigen::MatrixXd g_g_2_new(g_g_2_new_indic.size(), g_g_2.cols());
    for(int i_rows=0; i_rows<g_g_2_new_indic.size(); i_rows++ )
    {
      g_g_2_new.row(i_rows) = g_g_2.row(g_g_2_new_indic(i_rows));
    }
    
    g_g_23 = -1*TVBS_g_inverse_cpp(Sigma_11)*g_g_23*g_g_2;
    
    if((ddcor==0)&&(dd1cor==0))
    {
      g_g_x = g_g_22_new;
    }
    if((ddcor==0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_22_new.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_22_new.rows(), g_g_22_new.cols()) = g_g_22_new;
      g_g_x_new.block(g_g_22_new.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      g_g_x = g_g_x_new;
    }
    if((ddcor!=0)&&(dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_23.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
      g_g_x_new.block(g_g_23.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      
      Eigen::MatrixXd g_g_x_new_select(indic.size(), g_g_x_new.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_new_select.row(i_rows) = g_g_x_new.row(indic(i_rows));
      }
      g_g_x = g_g_x_new_select;
    }
    if((ddcor!=0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_23.rows()+g_g_22_new.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
      g_g_x_new.block(g_g_23.rows(), 0, g_g_22_new.rows(), g_g_2_new.cols()) = g_g_22_new;
      g_g_x_new.block(g_g_23.rows()+g_g_22_new.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      
      Eigen::MatrixXd g_g_x_new_select(indic.size(), g_g_x_new.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_new_select.row(i_rows) = g_g_x_new.row(indic(i_rows));
      }
      g_g_x = g_g_x_new_select;
    }
  }
  
  if((cholesky_global==1)&&(condcov_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1)&&(condcov_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  
  struct TVBS_Matrix_two output_final = {g_g_y, g_g_x};
  return output_final;
  
}





struct TVBS_Matrix_three TVBS_g_cond_cov_trunc_cpp(Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x)
{
  // initiate outputs
  Eigen::MatrixXd g_g_x;
  
  int dim1 = Sigma.rows()-x.size();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  Eigen::MatrixXd Sigma_11 = Sigma.block(0, 0, dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, Sigma.rows()-dimdiff, dimdiff);
  
  Eigen::VectorXd mu_trunc;
  Eigen::MatrixXd sigma_trunc, g_mu_trunc, g_sigma_trunc;
  if(dimdiff==1)
  {
    // THIS CASE DOES NOT OCCUR WITH TVBS
  }
  if(dimdiff==2)
  {
    struct TVBS_Matrix_two output = TVBS_truncate_bi_normal_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    mu_trunc = output.mat1;
    sigma_trunc = output.mat2;
    
    output = TVBS_truncate_bi_normal_gradient_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    g_mu_trunc = output.mat1;
    g_sigma_trunc = output.mat2;
  }
  
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  Eigen::MatrixXd b = Sigma_11_inv*sigma_trunc*Sigma_11_inv;
  
  // update some global variables
  cholesky_global = 0;
  if(condcovsigtrunc_global==0)
  {
    condcov_global = 0;
    x2correlation_global = 1;
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
  }
  if(condcovsigtrunc_global!=0)
  {
    condcov_global = 1;
    x2correlation_global = 0;
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
  }
  
  struct TVBS_Matrix_two output = TVBS_g_cond_cov_cpp(Eigen::MatrixXd::Identity(dim1, dim1), Sigma);
  Eigen::MatrixXd g_y = output.mat1;
  Eigen::MatrixXd g_d_1 = output.mat2;
  
  // update global variables
  x2symmetric_global = 1;
  x2diagonal_global = 0;
  x1symmetric_global = 1;
  x1diagonal_global = 0;
  x2correlation_global = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_11_inv, sigma_trunc);
  Eigen::MatrixXd g_Sigma_11_inv = output.mat1;
  Eigen::MatrixXd g_b_omega = output.mat2;
  
  Eigen::MatrixXd g_bpsi_11 = TVBS_g_inverse_cpp(Sigma_11)*g_Sigma_11_inv;
  
  // update global variables
  x2symmetric_global = 1;
  x2diagonal_global = 0;
  x1symmetric_global = 0;
  x1diagonal_global = 0;
  x2correlation_global = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_12, b);
  Eigen::MatrixXd g_psi_12 = output.mat1;
  Eigen::MatrixXd g_b = output.mat2;
  Eigen::MatrixXd g_psi_11 = g_bpsi_11*g_b;
  Eigen::MatrixXd g_omega = g_b_omega*g_b;
  Eigen::MatrixXd g_sigma_tilde_new = g_sigma_trunc*g_omega;
  
  Eigen::MatrixXd g_mu(dimdiff+dim1, g_sigma_tilde_new.cols());
  g_mu.setZero();
  g_mu.block(0, 0, dimdiff, g_sigma_tilde_new.cols()) = g_sigma_tilde_new.block(0, 0, dimdiff, g_sigma_tilde_new.cols());
  
  Eigen::MatrixXd g_c = g_sigma_tilde_new.block(g_sigma_tilde_new.rows()-dimdiff, 0, dimdiff, g_sigma_tilde_new.cols());
  
  Eigen::MatrixXd indic_1_r(dim1*dimdiff,1);
  indic_1_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
  indic_1_r.resize(dimdiff, dim1);
  Eigen::MatrixXd indic_1 = indic_1_r.transpose();
  indic_1.resize(dim1*dimdiff, 1);
  
  {
    Eigen::MatrixXd g_psi_12_foo(indic_1.size(), g_psi_12.cols());
    for(int i_rows=0; i_rows<indic_1.size(); i_rows++)
    {
      g_psi_12_foo.row(i_rows)=g_psi_12.row(indic_1(i_rows));
    }
    g_psi_12 = g_psi_12_foo;
  }
  
  Eigen::MatrixXd g_psi_22(g_d_1.rows(), g_d_1.cols());
  g_psi_22.setZero();
  
  
  if(condcov_global==1)
  {
    g_psi_11 += g_sigma_tilde_new.block(dimdiff, 0, g_sigma_tilde_new.rows()-dimdiff-dimdiff, g_sigma_tilde_new.cols());
    
    Eigen::MatrixXd mat_n_seq(ddcov,1);
    mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq);
    
    Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
    indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcov, ddcov+dim1*dimdiff-1);
    indic_2_2_r.resize(dim1, dimdiff);
    Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
    
    Eigen::MatrixXd indic_2(mat_d_up.rows(), mat_d_up.cols()+indic_2_2.cols());
    indic_2.block(0, 0, mat_d_up.rows(), mat_d_up.cols()) = mat_d_up;
    indic_2.block(0, mat_d_up.cols(), mat_d_up.rows(), indic_2_2.cols()) = indic_2_2;
    
    
    Eigen::MatrixXd mat_n_seq_2(dd1cov,1);
    mat_n_seq_2 = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, dd1cov+ddcov+dimdiff*dim1-1);
    Eigen::MatrixXd mat_d_up_2 = TVBS_vec_2_upper_diag_cpp(mat_n_seq_2);
    
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up_2.cols());
    indic_3.setZero();
    indic_3.block(0, dimdiff, dim1, mat_d_up_2.cols()) = mat_d_up_2;
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose());
    
    Eigen::MatrixXd g_g_x_1(g_psi_11.rows()+g_psi_12.rows()+g_psi_22.rows(), g_psi_22.cols());
    g_g_x_1.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
    g_g_x_1.block(g_psi_11.rows(), 0, g_psi_12.rows(), g_psi_11.cols()) = g_psi_12;
    g_g_x_1.block(g_psi_11.rows()+g_psi_12.rows(), 0, g_psi_22.rows(), g_psi_11.cols()) = g_psi_22;
    Eigen::MatrixXd g_g_x_2(indic.size(), g_g_x_1.cols());
    for(int i_rows=0; i_rows<indic.size(); i_rows++)
    {
      g_g_x_2.row(i_rows) = g_g_x_1.row(indic(i_rows));
    }
    
    g_g_x = g_d_1 + g_g_x_2;
  }
  if(condcov_global==0)
  {
    Eigen::MatrixXd indic_2, indic_3;
    
    if(dimdiff==2)
    {
      g_psi_11 += g_sigma_tilde_new.row(dimdiff+1);
    }
    
    if(ddcor==0)
    {
      Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
      indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_2_r.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
      
      Eigen::MatrixXd indic_2_int(indic_2_2.rows(), indic_2_2.cols()+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, indic_2_2.rows(), indic_2_2.cols()) = indic_2_2;
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      Eigen::MatrixXd mat_n_seq(ddcor,1);
      mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor,0,ddcor-1);
      Eigen::MatrixXd mat_nd_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq,0);
      
      Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
      indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_2_r.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
      
      Eigen::MatrixXd indic_2_int(indic_2_2.rows(), mat_nd_up.cols()+indic_2_2.cols());
      indic_2_int.block(0, 0, mat_nd_up.rows(), mat_nd_up.cols()) = mat_nd_up;
      indic_2_int.block(0, mat_nd_up.cols(), indic_2_2.rows(), indic_2_2.cols()) = indic_2_2;
      
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      Eigen::MatrixXd mat_n_seq(dd1cor,1);
      mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+dd1cor-1);
      Eigen::MatrixXd mat_nd_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+mat_nd_up.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, mat_nd_up.rows(), mat_nd_up.cols()) = mat_nd_up;
      indic_3 = indic_3_int;
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(),0);
    
    if((ddcor==0)&&(dd1cor==0))
    {
      g_g_x = g_d_1 + g_psi_12;
    }
    if((ddcor==0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_12.rows()+dd1cor, g_psi_12.cols());
      g_g_x_rhs.setZero();
      g_g_x_rhs.block(0, 0, g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      g_g_x = g_d_1 + g_g_x_rhs;
    }
    if((ddcor!=0)&&(dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_11.rows()+g_psi_12.rows(), g_psi_11.cols());
      g_g_x_rhs.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
      g_g_x_rhs.block(0, g_psi_11.rows(), g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      Eigen::MatrixXd g_g_x_rhs_ordered(indic.size(), g_g_x_rhs.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_rhs_ordered.row(i_rows) = g_g_x_rhs.row(indic(i_rows));
      }
      g_g_x = g_d_1+g_g_x_rhs_ordered;
    }
    if((ddcor!=0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_11.rows()+g_psi_12.rows()+g_psi_22.rows(), g_psi_11.cols());
      g_g_x_rhs.setZero();
      g_g_x_rhs.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
      g_g_x_rhs.block(g_psi_11.rows(), 0, g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      g_g_x_rhs.block(g_psi_11.rows()+g_psi_12.rows(), 0, g_psi_22.rows(), g_psi_22.cols()) = g_psi_22;
      Eigen::MatrixXd g_g_x_rhs_ordered(indic.size(), g_g_x_rhs.cols());
      g_g_x_rhs_ordered.setZero();
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_rhs_ordered.row(i_rows) = g_g_x_rhs.row(indic(i_rows));
      }
      g_g_x = g_d_1+g_g_x_rhs_ordered;
      
    }
  }
  
  if((cholesky_global==1)&&(condcov_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1)&&(condcov_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_g_x, g_c};
  return output_final;
}




struct TVBS_Matrix_four TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd Id_mat, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x)
{
  // initialize output
  Eigen::MatrixXd g_g_x, g_c;
  
  // extract dimensions of the inputs
  int dim1 = Id_mat.rows();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  // for correct sequencing of the entries in Sigma
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  // extract submatrices of the Covariance matrix
  Eigen::MatrixXd Sigma_11 = Sigma.block(0,0,dimdiff,dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, Sigma.rows()-dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  
  
  Eigen::MatrixXd g_mu_tilde = (Id_mat*Sigma_12*Sigma_11_inv).transpose();
  Eigen::MatrixXd g_mu(g_mu_tilde.rows()+dim1, dim1);
  g_mu.block(0, 0, g_mu_tilde.rows(), g_mu_tilde.cols()) = -1.0*g_mu_tilde;
  g_mu.block(g_mu_tilde.rows(), 0, dim1, dim1) = Eigen::MatrixXd::Identity(dim1, dim1);
  
  
  // truncate
  Eigen::VectorXd mu_trunc;
  Eigen::MatrixXd sigma_trunc, g_mu_trunc, g_sigma_trunc;
  if(dimdiff==1)
  {
    // not relevant for TVBS 
  }
  if(dimdiff==2)
  {
    struct TVBS_Matrix_two output = TVBS_truncate_bi_normal_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    mu_trunc = output.mat1;
    sigma_trunc = output.mat2;
    output = TVBS_truncate_bi_normal_gradient_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    g_mu_trunc = output.mat1;
    g_sigma_trunc = output.mat2;
  }
  
  Eigen::MatrixXd g_y_1 = Sigma_12*Sigma_11_inv*(mu_trunc-mu.head(dimdiff));
  Eigen::MatrixXd g_g_y(g_y_1.rows(), g_y_1.rows());
  g_g_y.setZero();
  g_g_y.diagonal() = g_y_1;
  
  
  Eigen::MatrixXd g_mu_tilde_new = g_mu_trunc*g_mu_tilde;
  g_mu.block(0, 0, dimdiff, g_mu.cols()) += g_mu_tilde_new.block(0, 0, dimdiff, g_mu_tilde_new.cols());
  
  Eigen::MatrixXd g_x_12 = kroneckerProduct(Id_mat, Sigma_11_inv*(mu_trunc-mu.head(dimdiff)));
  
  // reordering indices for matrix transposion
  {
    Eigen::MatrixXd indic_1_t = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
    indic_1_t.resize(dimdiff, dim1);
    Eigen::MatrixXd indic_1 = indic_1_t.transpose();
    indic_1.resize(dim1*dimdiff,1);
    {
      Eigen::MatrixXd g_x_12_ordered(indic_1.rows(), g_x_12.cols());
      for(int i_rows = 0; i_rows<indic_1.rows(); i_rows++)
      {
        g_x_12_ordered.row(i_rows) = g_x_12.row(indic_1(i_rows));
      }
      g_x_12 = g_x_12_ordered;
    }
  }
  
  // set global variables
  omsymmetric_global = 1;
  omdiagonal_global = 0;
  
  if(condcovmeantrunc_global==1)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
    
    Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq);
    
    Eigen::MatrixXd indic_2_1_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcov, dimdiff*dim1+ddcov-1);
    indic_2_1_t.resize(dim1, dimdiff);
    Eigen::MatrixXd indic_2_1 = indic_2_1_t.transpose();
    
    Eigen::MatrixXd indic_2(mat_d_up.rows(), mat_d_up.cols()+indic_2_1.cols() );
    indic_2.block(0, 0, mat_d_up.rows(), mat_d_up.cols()) = mat_d_up;
    indic_2.block(0, mat_d_up.cols(), indic_2_1.rows(), indic_2_1.cols()) = indic_2_1;
    
    mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, ddcov+dimdiff*dim1+dd1cov-1);
    mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq);
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up.cols());
    indic_3.setZero();
    indic_3.block(0, dimdiff, dim1, mat_d_up.cols()) = mat_d_up;
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose());
    
    Eigen::MatrixXd g_x_11 = TVBS_g_inverse_cpp(Sigma_11)*TVBS_g_a_omega_b_cpp(Id_mat*Sigma_12, mu_trunc-mu.head(dimdiff));
    g_x_11 += g_mu_tilde_new.block(dimdiff, 0, (dimdiff+dimdiff*(dimdiff+1)/2)-(dimdiff+1)+1, g_mu_tilde_new.cols());
    
    Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows()+dd1cov, dim1);
    g_g_x_int.setZero();
    g_g_x_int.block(0, 0, g_x_11.rows(), dim1) = g_x_11;
    g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), dim1) = g_x_12;
    {
      Eigen::MatrixXd g_g_x_ordered(indic.size(), dim1);
      for(int i_rows=0; i_rows<g_g_x_int.rows(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
  }
  
  if(condcovmeantrunc_global==0)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
    
    Eigen::MatrixXd g_x_11 = TVBS_g_inverse_cpp(Sigma_11)*TVBS_g_a_omega_b_cpp(Id_mat*Sigma_12, mu_trunc-mu.head(dimdiff));
    
    if(dimdiff==2)
    {
      g_x_11 += g_mu_tilde_new.row(dimdiff+1);
    }
    
    // for reordering purposes:
    Eigen::MatrixXd indic_2, indic_3;
    if(ddcor==0)
    {
      Eigen::MatrixXd indic_2_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_t.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, dim1+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, dimdiff, dim1) = indic_2_t.transpose();
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor, 0, ddcor-1);
      Eigen::MatrixXd mat_nd_up_diag_zero = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      Eigen::MatrixXd indic_2_2_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_2_t.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, mat_nd_up_diag_zero.cols()+dim1);
      indic_2_int.block(0, 0, dimdiff, mat_nd_up_diag_zero.cols()) = mat_nd_up_diag_zero;
      indic_2_int.block(0, mat_nd_up_diag_zero.cols(), dimdiff, dim1) = indic_2_2_t.transpose();
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+ddcor-1);
      Eigen::MatrixXd mat_nd_up_diag_zero = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+mat_nd_up_diag_zero.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, mat_nd_up_diag_zero.rows(), mat_nd_up_diag_zero.cols()) = mat_nd_up_diag_zero;
      indic_3 = indic_3_int;
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_3.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(),0);
    
    
    if((ddcor==0) && (dd1cor==0))
    {
      g_g_x = g_x_12;
    }
    if((ddcor==0) && (dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_12.rows()+dd1cor, dim1);
      g_g_x_int.setZero();
      g_g_x_int.block(0, 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      g_g_x = g_g_x_int;
    }
    if((ddcor!=0) && (dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows(), g_x_12.cols());
      g_g_x_int.block(0, 0, g_x_11.rows(), g_x_11.cols()) = g_x_11;
      g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      Eigen::MatrixXd g_g_x_ordered(indic.size(), g_x_12.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
    if((ddcor!=0) && (dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows()+dd1cor, g_x_12.cols());
      g_g_x_int.setZero();
      g_g_x_int.block(0, 0, g_x_11.rows(), g_x_11.cols()) = g_x_11;
      g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      Eigen::MatrixXd g_g_x_ordered(indic.size(), g_x_12.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
    
  }
  
  if((cholesky_global==1) && (condcovmeantrunc_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1) && (condcovmeantrunc_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  g_c = g_mu_tilde_new.block(dimdiff+(dimdiff*(dimdiff+1))/2, 0, g_mu_tilde_new.rows()-(dimdiff+dimdiff*((dimdiff+1))/2), g_mu_tilde_new.cols());
  
  
  struct TVBS_Matrix_four output_final = {g_g_y, g_mu, g_g_x, g_c};
  return output_final;
}



// // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//                                                                                  //
//    Functions to evaluate bi, tri- and quadro-variate NCDF and the gradients      //
//                                                                                  //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // //


// // // // // // // // //
//  bivariate           //
// // // // // // // // //

struct TVBS_double_three TVBS_grad_cdf_bvn_cpp(Eigen::VectorXd x_temp, Eigen::MatrixXd Rho)
{
  double rhotilde, tr1, tr2, pdf2, g_w1, g_w2;
  
  
  rhotilde = sqrt(1.0-pow(Rho(0,1), 2));
  
  //Rcout << "rhotilde" << rhotilde << "x_t" << x_temp << std::endl; 
  
  tr1 = (x_temp(1) - Rho(0,1)*x_temp(0))/rhotilde;
  tr2 = (x_temp(0) - Rho(0,1)*x_temp(1))/rhotilde;
  pdf2 = std_normal_pdf(x_temp(0))*std_normal_pdf(tr1)/rhotilde;
  g_w1 = std_normal_pdf(x_temp(0))*std_normal_cdf(tr1);
  g_w2 = std_normal_pdf(x_temp(1))*std_normal_cdf(tr2);
  
  struct TVBS_double_three output_final = {g_w1, g_w2, pdf2};

  return output_final;
}


struct TVBS_double_three TVBS_grad_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd w, Eigen::MatrixXd Rho)
{
  double bivar_cdf, univar_cdf, g_w1, g_w2, grho, gw1new;
  Eigen::VectorXd mu(2);
  mu.setZero();

  bivar_cdf = TVBS_pmvnorm_cpp(w, mu, Rho);
  univar_cdf = std_normal_cdf(w(0));
  
  if (univar_cdf<tol) { univar_cdf = tol; }
  
  struct TVBS_double_three output = TVBS_grad_cdf_bvn_cpp(w, Rho);
  g_w1 = output.d1;
  g_w2 = output.d2;
  grho = output.d3;
  
  gw1new = (univar_cdf*g_w1 - bivar_cdf*std_normal_pdf(w(0)))/pow(univar_cdf, 2);
  
  struct TVBS_double_three output_final = {gw1new, g_w2/univar_cdf, grho/univar_cdf};
  return output_final;
}




struct TVBS_Matrix_three TVBS_grad_non_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  Eigen::VectorXd Sigma_diag_sqrt, x_norm, g_w(2);
  Eigen::MatrixXd Cor_mat, g_b_corcov, g_omega_corcov, g_cov, g_x, g_mu;
  double rho, g_w1, g_w2, grho;
  
  Sigma_diag_sqrt = cov.diagonal().array().sqrt();

  x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  
  Eigen::MatrixXd Sigma_diag_sqrt_inv(2,2);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();

  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_by_cdfn_cpp(x_norm.head(2), Cor_mat);
  g_w1 = output_3.d1;
  g_w2 = output_3.d2;
  grho = output_3.d3;
  g_w << g_w1, g_w2;
  

  struct TVBS_Matrix_two output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  g_b_corcov = output_2.mat1;
  g_omega_corcov = output_2.mat2;
  
  g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  g_cov = g_b_corcov*g_w + g_omega_corcov*grho;
  g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  return output_final;
}




// // // // // // // // //
//  trivariate          //
// // // // // // // // //

struct TVBS_Matrix_two TVBS_grad_cdf_tvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr)
{
  // initialize variables
  int enteredloop1, enteredloop2;
  double epst, pt, h1, h2, h3, r12, r13, r23, tvn;
  Eigen::VectorXd gg(6), d_tvn(6);
  
  // safeguard threshold against numerical errors 
  epst = pow(10, -7);
  
  pt = pi_global/2;
  h1 = x_norm(0);
  h2 = x_norm(1);
  h3 = x_norm(2);
  r12 = Corr(0,1);
  r13 = Corr(0,2);
  r23 = Corr(1,2);
  
  // initialize vector to store gradient output
  gg.setZero();
  
  // following are safeguards against numerical errors (implemented by a reordering)
  enteredloop1 = 0;     // control variables to check which safeguard was used
  enteredloop2 = 0;
  
  // check if one of the correlations is bigger than the other
  // (kind of optimal ordering)
  if( (abs(r12)-abs(r13)) > epst)
  {
    // 1
    h2 = h3;
    h3 = x_norm(1);
    r12 = r13;
    r13 = Corr(0,1);
    enteredloop1 = 1;
  }
  if( (abs(r13)-abs(r23)) > epst)
  {
    // 2
    h1 = h2;
    h2 = x_norm(0);
    r23 = r13;
    r13 = Corr(1,2);
    enteredloop2 = 1;
  }
  
  // take care of some cases with small values seperately
  if( (abs(h1) + abs(h2) + abs(h3)) < epst )
  {
    // 3
    gg(3) = ( 2* 1.0/sqrt(1.0-pow(r12,2))/pi_global )/8;
    gg(4) = ( 2* 1.0/sqrt(1.0-pow(r13,2))/pi_global )/8;
    gg(5) = ( 2* 1.0/sqrt(1.0-pow(r23,2))/pi_global )/8;
    d_tvn = gg;
  } 
  else if( (abs(r12) + abs(r13))<epst )
  {
    // 4
    double tempval1 = std_normal_cdf(h1);
    
    Eigen::VectorXd h23(2), mu_h23(2);
    mu_h23.setZero();
    h23 << h2, h3;
    Eigen::MatrixXd corr_h23(2,2);
    corr_h23.setOnes();
    corr_h23(0,1) = r23;
    corr_h23(1,0) = r23;
    double tempval2 = TVBS_pmvnorm_cpp(h23, mu_h23, corr_h23);
    
    gg(0) = std_normal_pdf(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h23, corr_h23);
    gg(1) = tempval1*temp_bvn_grad.d1;
    gg(2) = tempval1*temp_bvn_grad.d2;
    gg(5) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r13) + abs(r23))<epst )
  {
    // 5
    double tempval1 = std_normal_cdf(h3);
    
    Eigen::VectorXd h12(2), mu_h12(2);
    mu_h12.setZero();
    h12 << h1, h2;
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    double tempval2 = TVBS_pmvnorm_cpp(h12, mu_h12, corr_h12);
    
    gg(2) = std_normal_pdf(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h12, corr_h12);
    gg(0) = tempval1*temp_bvn_grad.d1;
    gg(1) = tempval1*temp_bvn_grad.d2;
    gg(3) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r12) + abs(r23))<epst )
  {
    // 6
    double tempval1 = std_normal_cdf(h2);
    
    Eigen::VectorXd h13(2), mu_h13(2);
    mu_h13.setZero();
    h13 << h1, h3;
    Eigen::MatrixXd corr_h13(2,2);
    corr_h13.setOnes();
    corr_h13(0,1) = r13;
    corr_h13(1,0) = r13;
    double tempval2 = TVBS_pmvnorm_cpp(h13, mu_h13, corr_h13);
    
    // tvn should have no influence on this function
    tvn = std_normal_cdf(h2)*tempval2;
    
    
    gg(1) = std_normal_pdf(h2)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h13, corr_h13);
    gg(0) = tempval1*temp_bvn_grad.d1;
    gg(2) = tempval1*temp_bvn_grad.d2;
    gg(4) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (1.0-r23)<epst )
  {
    // 7
    double h_23_min = std::min(h2,h3);
    
    Eigen::VectorXd h1_23_min(2), mu_h13(2);
    mu_h13.setZero();
    h1_23_min << h1, h_23_min;
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    
    // tvn should have no influence on this function
    tvn = TVBS_pmvnorm_cpp(h1_23_min, mu_h13, corr_h12);
    
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h1_23_min, corr_h12);
    
    gg(0) = temp_bvn_grad.d1;
    if(h2==h3)
    {
      gg(1) = temp_bvn_grad.d2;
      gg(2) = temp_bvn_grad.d2;
    }
    else if(h2<h3)
    {
      gg(1) = temp_bvn_grad.d2;
    }
    else
    {
      gg(2) = temp_bvn_grad.d2;
    }
    gg(3) = temp_bvn_grad.d3;
    d_tvn = gg;
    
  }
  else if( (r23+1.0<epst) && (h2>-1.0*h3) )
  {
    // 8
    Eigen::VectorXd h_12(2), h_1m3(2);
    h_12 << h1, h2;
    h_1m3 << h1, -1.0*h3;
    
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h_12, corr_h12);
    struct TVBS_double_three temp_bvn_grad2 = TVBS_grad_cdf_bvn_cpp(h_1m3, corr_h12);
    gg(0) = temp_bvn_grad.d1 - temp_bvn_grad2.d1;
    gg(1) = temp_bvn_grad.d2;
    gg(2) = temp_bvn_grad2.d2;
    gg(3) = temp_bvn_grad.d3 - temp_bvn_grad2.d3;
    d_tvn = gg;
    
  }
  else
  {
    // now the general case....
    Eigen::Vector2d h_23, mu_h23;
    mu_h23.setZero();
    h_23 << h2, h3;
    Eigen::Matrix2d corr_h23;
    corr_h23 << 1, r23, r23, 1;
    
    double cdfbvn_h2_h3 = TVBS_pmvnorm_cpp(h_23, mu_h23, corr_h23);
    
    double cdfn_h1 = std_normal_cdf(h1);
    tvn = cdfn_h1*cdfbvn_h2_h3;
    
    double d_a = std_normal_pdf(h2)*std_normal_cdf((h3-r23*h2)/sqrt(1.0-pow(r23,2)));
    double d_b = std_normal_pdf(h3)*std_normal_cdf((h2-r23*h3)/sqrt(1.0-pow(r23,2)));
    double d_corr = (exp(-0.5*(pow(h2,2) + pow(h3,2) - 2*r23 * h2 * h3   ) / (1.0-r23*r23) )) / sqrt(1.0-r23 * r23)/(2*pi_global);
    
    d_tvn << (std_normal_pdf(h1)*cdfbvn_h2_h3), cdfn_h1*d_a, (cdfn_h1*d_b), 0, 0, (cdfn_h1*d_corr);
    
    
    // bunch of auxillary variables
    Eigen::VectorXd wg(12), xg(12);
    wg << 0.1279381953467518, 0.1258374563468280, 0.1216704729278031, 0.1155056680537265, 0.1074442701159659, 0.09761865210411358, 0.08619016153195296, 0.07334648141108081, 0.05929858491543594, 0.04427743881742087, 0.02853138862893389, 0.01234122979998693;
    xg << 0.06405689286260559, 0.1911188674736164, 0.3150426796961635, 0.4337935076260450, 0.5454214713888396, 0.6480936519369754, 0.7401241915785546, 0.8200019859739028, 0.8864155270044012, 0.9382745520027329, 0.9747285559713096, 0.9951872199970214;
    double rua = asin(r12);
    double rub = asin(r13);
    double res = 0.0;
    double d_rua_r12 = 1.0/sqrt(1.0-pow(r12,2));
    double d_rub_r13 = 1.0/sqrt(1.0-pow(r13,2));
    double d_res_h1 = 0;
    double d_res_h2 = 0;
    double d_res_h3 = 0;
    double d_res_r12 = 0;
    double d_res_r13 = 0;
    double d_res_r23 = 0;
    
    // here comes the (numerical) magic
    for(int j_int = 0; j_int<12; j_int++)
    {
      double fc = 0;
      double d_fc_h1 = 0;
      double d_fc_h2 = 0;
      double d_fc_h3 = 0;
      double d_fc_r12 = 0;
      double d_fc_r13 = 0;
      double d_fc_r23 = 0;
      
      Eigen::VectorXd temp_sincs_grad = TVBS_sincs_grad_cpp(rua*(1.0-xg(j_int))/2);
      double r12t = temp_sincs_grad(2);
      double rr2 = temp_sincs_grad(3);
      double d_r12t_r12 = temp_sincs_grad(0)*( 1.0 - xg(j_int) )/2*d_rua_r12;
      double d_rr2_r12 = temp_sincs_grad(1)*( 1.0 - xg(j_int) )/2*d_rua_r12;
      
      temp_sincs_grad = TVBS_sincs_grad_cpp(rub*(1.0-xg(j_int))/2);
      double r13t = temp_sincs_grad(2);
      double rr3 = temp_sincs_grad(3);
      double d_r13t_r13 = temp_sincs_grad(0)*( 1.0 - xg(j_int) )/2*d_rub_r13;
      double d_rr3_r13 = temp_sincs_grad(1)*( 1.0 - xg(j_int) )/2*d_rub_r13;
      
      if( abs(rua)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h2, h3, r13t, r23, r12t, rr2 );
        fc += rua*temp_pntgnd_grad(7);
        d_fc_h1 += rua*temp_pntgnd_grad(0);
        d_fc_h2 += rua*temp_pntgnd_grad(1);
        d_fc_h3 += rua*temp_pntgnd_grad(2);
        d_fc_r12 += d_rua_r12*temp_pntgnd_grad(7) + rua*(temp_pntgnd_grad(5)*d_r12t_r12 + temp_pntgnd_grad(6)*d_rr2_r12 );
        d_fc_r13 += rua*(temp_pntgnd_grad(3)*d_r13t_r13);
        d_fc_r23 += rua*temp_pntgnd_grad(4);
        
      }
      if( abs(rub)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h3, h2, r12t, r23, r13t, rr3 );
        fc += rub*temp_pntgnd_grad(7);
        d_fc_h1 += rub*temp_pntgnd_grad(0);
        d_fc_h2 += rub*temp_pntgnd_grad(2);
        d_fc_h3 += rub*temp_pntgnd_grad(1);
        d_fc_r12 += rub*temp_pntgnd_grad(3)*d_r12t_r12;
        d_fc_r13 += d_rub_r13*temp_pntgnd_grad(7) + rub*(temp_pntgnd_grad(5)*d_r13t_r13 + temp_pntgnd_grad(6)*d_rr3_r13);
        d_fc_r23 += rub*temp_pntgnd_grad(4);
        
      }
      
      temp_sincs_grad = TVBS_sincs_grad_cpp( rua*( 1.0 + xg(j_int) )/2 );
      r12t = temp_sincs_grad(2);
      rr2 = temp_sincs_grad(3);
      d_r12t_r12 = temp_sincs_grad(0)*( 1.0 + xg(j_int) )/2*d_rua_r12;
      d_rr2_r12 = temp_sincs_grad(1)*( 1.0 + xg(j_int) )/2*d_rua_r12;
      
      temp_sincs_grad = TVBS_sincs_grad_cpp( rub*( 1.0 + xg(j_int) )/2 );
      r13t = temp_sincs_grad(2);
      rr3 = temp_sincs_grad(3);
      d_r13t_r13 = temp_sincs_grad(0)*( 1.0 + xg(j_int) )/2*d_rub_r13;
      d_rr3_r13 = temp_sincs_grad(1)*( 1.0 + xg(j_int) )/2*d_rub_r13;
      
      if( abs(rua)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h2, h3, r13t, r23, r12t, rr2 );
        fc += rua*temp_pntgnd_grad(7);
        d_fc_h1 += rua*temp_pntgnd_grad(0);
        d_fc_h2 += rua*temp_pntgnd_grad(1);
        d_fc_h3 += rua*temp_pntgnd_grad(2);
        d_fc_r12 += d_rua_r12*temp_pntgnd_grad(7) + rua*(temp_pntgnd_grad(5)*d_r12t_r12 + temp_pntgnd_grad(6)*d_rr2_r12 );
        d_fc_r13 += rua*(temp_pntgnd_grad(3)*d_r13t_r13);
        d_fc_r23 += rua*temp_pntgnd_grad(4);
      }
      if( abs(rub)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h3, h2, r12t, r23, r13t, rr3 );
        fc += rub*temp_pntgnd_grad(7);
        d_fc_h1 += rub*temp_pntgnd_grad(0);
        d_fc_h2 += rub*temp_pntgnd_grad(2);
        d_fc_h3 += rub*temp_pntgnd_grad(1);
        d_fc_r12 += rub*temp_pntgnd_grad(3)*d_r12t_r12;
        d_fc_r13 += d_rub_r13*temp_pntgnd_grad(7) + rub*(temp_pntgnd_grad(5)*d_r13t_r13 + temp_pntgnd_grad(6)*d_rr3_r13);
        d_fc_r23 += rub*temp_pntgnd_grad(4);
      }
      
      res += wg(j_int)*fc;
      d_res_h1 += wg(j_int)*d_fc_h1;
      d_res_h2 += wg(j_int)*d_fc_h2;
      d_res_h3 += wg(j_int)*d_fc_h3;
      d_res_r12 += wg(j_int)*d_fc_r12;
      d_res_r13 += wg(j_int)*d_fc_r13;
      d_res_r23 += wg(j_int)*d_fc_r23;
    }
    
    tvn += res/(4*pi_global);
    Eigen::VectorXd d_tvn_add(6);
    d_tvn_add << d_res_h1, d_res_h2, d_res_h3, d_res_r12, d_res_r13, d_res_r23;
    d_tvn += d_tvn_add/(4*pi_global);
    
  }
  
  // reorder according to which case was present
  if(enteredloop2==1)
  {
    Eigen::VectorXd d_tvn_unsorted = d_tvn;
    d_tvn_unsorted(0) = d_tvn(1);
    d_tvn_unsorted(1) = d_tvn(0);
    d_tvn_unsorted(4) = d_tvn(5);
    d_tvn_unsorted(5) = d_tvn(4);
    d_tvn = d_tvn_unsorted;
  }
  if(enteredloop1==1)
  {
    Eigen::VectorXd d_tvn_unsorted = d_tvn;
    d_tvn_unsorted(1) = d_tvn(2);
    d_tvn_unsorted(2) = d_tvn(1);
    d_tvn_unsorted(3) = d_tvn(4);
    d_tvn_unsorted(4) = d_tvn(3);
    d_tvn = d_tvn_unsorted;
  }
  
  struct TVBS_Matrix_two output_final = {d_tvn.head(3), d_tvn.tail(3)};
  return output_final;
}




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
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2));
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
    g_c_new(temp1(i_row)) = g_c(i_row);
    for(int i_col=0; i_col<m; i_col++)
    {
      g_rho_rho_mat(temp1(i_row),temp1(i_col)) = g_rho_rho_mat_new(i_row,i_col);
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
      // claculate the truncated mean and variance
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2));
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
    g_c_new(temp1(i_row)) = g_c(i_row);
    for(int i_col=0; i_col<m; i_col++)
    {
      g_rho_rho_mat(temp1(i_row),temp1(i_col)) = g_rho_rho_mat_new(i_row,i_col);
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
      // claculate the truncated mean and variance
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2));
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
      struct TVBS_Matrix_two output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.block(0,k1-1,2,1), D.block(0,(k1-1)*m,2,2), x_temp.segment((2*k1-1)-1,2));
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
        dx_temp(ind_mat(np,0))=1; // parameter corresponds to b entry.
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
        d1x_temp(ind_mat(np1,0),np1)=1; // parameter corresponds to b entry.
        
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
