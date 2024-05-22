#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;


//[[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies


// [[Rcpp::plugins(cpp11)]]

double maxEV(Eigen::MatrixXf mat);
struct state_space_system param_to_system(Eigen::VectorXd param, int s, int n, int grad_bool, bool stationary);
Rcpp::List  param_to_system_R(Eigen::VectorXd param, int s, int n, int grad_bool, bool stationary); 

struct state_space_system param_to_system_grad(int coord, int s, int n);
Eigen::MatrixXd solve_Lyapunov_equation( Eigen::MatrixXd A, Eigen::VectorXd vQ);
Eigen::MatrixXd calculate_M(struct state_space_system syst , Eigen::MatrixXd Sigma);
Eigen::MatrixXd calculate_observability(Eigen::MatrixXd A,Eigen::MatrixXd C,Eigen::VectorXd dtime);
Eigen::MatrixXd calculate_deriv_observability(Eigen::MatrixXd A,Eigen::MatrixXd C,Eigen::MatrixXd dA,Eigen::MatrixXd dC,Eigen::VectorXd dtime);

Eigen::MatrixXd calculate_Cov_Seq(struct state_space_system syst, Eigen::MatrixXd P0, Eigen::MatrixXd Sigma, Eigen::VectorXd time);
Eigen::MatrixXd calculate_dCov_Seq(struct state_space_system syst, struct state_space_system dsyst,Eigen::MatrixXd P0, Eigen::MatrixXd dP0, Eigen::MatrixXd Sigma, Eigen::MatrixXd dSigma, Eigen::VectorXd time);

//struct system build_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time);
//Rcpp::List build_system_from_model_R(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time);
//Eigen::MatrixXd build_derived_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time, int coord);

