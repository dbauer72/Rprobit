

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Matrix3d;    // one of the eigenvalue solvers


//Solow-Joe Approximation
//INPUTS
//- x: a vector containing the upper integration limits
//- r: a correaltion matrix
//OUTPUT
//-sj: First the gradient infos (linear coeffs and than the off-diagonal correlations) and then the log-probability



//[[Rcpp::depends(RcppEigen)]]
//using namespace Eigen; 
using Eigen::Map;                                               // 'maps' rather than copies
using Eigen::MatrixXd;                                          // variable size matrix, double precision
using Eigen::VectorXd;          

double std_normal_cdf(double x);                          // create a function for the standard normal cumulative distribution function
double normal_cdf(double x, double v);                          // create a function for the standard normal cumulative distribution function
double std_normal_pdf(double x);                                 // create a function for the standard normal probability distribution function
double normal_pdf(double x, double v);                                 // create a function for the standard normal probability distribution function
double biv_normal_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
double biv_normal_pdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
Eigen::VectorXd grad_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
Eigen::MatrixXd Hess_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
Eigen::VectorXd trilCols(int M);
Eigen::MatrixXd Omega(Eigen::VectorXd x,Eigen::MatrixXd r);
double SJ(Eigen::VectorXd x, Eigen::MatrixXd r);
List dOmega(Eigen::VectorXd x,Eigen::MatrixXd r); // calculates the gradient of Omega.
Eigen::VectorXd dlcond(Eigen::VectorXd x,Eigen::MatrixXd r); // calculates the gradient of Omega.
List d2Omega(Eigen::VectorXd x,Eigen::MatrixXd r); // calculates the Hessian of Omega.
Eigen::MatrixXd SJ_hess_new(Eigen::VectorXd x,Eigen::MatrixXd r);

