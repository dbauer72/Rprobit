using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Matrix3d;    // one of the eigenvalue solvers

using Eigen::Map;                                               // 'maps' rather than copies

double std_normal_cdf(double x);                          // create a function for the standard normal cumulative distribution function
double normal_cdf(double x, double v);                          // create a function for the standard normal cumulative distribution function
double std_normal_pdf(double x);                                 // create a function for the standard normal probability distribution function
double normal_pdf(double x, double v);                                 // create a function for the standard normal probability distribution function
double std_normal_cdf_inv_1(double y);
Eigen::VectorXd std_normal_cdf_inv(Eigen::VectorXd y);  
double biv_normal_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
double biv_normal_pdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
double TVBS_biv_std_norm_pdf(Eigen::VectorXd x, double rho); 
Eigen::VectorXd grad_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
Eigen::MatrixXd Hess_cdf(double x0, double x1, double r12);                                 // create a function for the standard normal probability distribution function
Eigen::MatrixXd Hess_pdf(double w0, double w1, double rho);
double biv_gen_pdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
double biv_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::VectorXd grad_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::MatrixXd Hess_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::MatrixXd Hessian_gen_tvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::VectorXd grad_gen_tvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
double tvn_gen_pdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
Eigen::MatrixXd J_wR_bSig(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
