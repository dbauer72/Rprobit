#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Matrix3d;    // one of the eigenvalue solvers


//Mendel-Elston Approximation
//INPUTS
//- x: a vector containing the upper integration limits
//- r: a correaltion matrix
//OUTPUT
//-logp: First the gradient infos (linear coeffs and than the off-diagonal correlations) and then the log-probability


//[[Rcpp::depends(RcppEigen)]]
//using namespace Eigen; 
using Eigen::Map;                                               // 'maps' rather than copies
using Eigen::MatrixXd;                                          // variable size matrix, double precision
using Eigen::VectorXd;          

int ind_diag(int r, int c, int m);                                // index of element (r,c) in vech of upper diagonal including the diagonal 
int ind_nodiag(int r, int c, int m);                              // index of element (r,c) in vech of upper diagonal excluding the diagonal 
double ME(Eigen::VectorXd x, Eigen::MatrixXd r);                  // log of (approximated) probability
Eigen::VectorXd dlcond_ME(Eigen::VectorXd x,Eigen::MatrixXd r);   // gradient of log probability
Eigen::MatrixXd ME_hess_new(Eigen::VectorXd x,Eigen::MatrixXd r); // Hessian of log probability
