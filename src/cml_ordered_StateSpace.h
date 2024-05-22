#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;


//[[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies


// [[Rcpp::plugins(cpp11)]]

struct system build_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time);
Eigen::MatrixXd build_derived_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time, int coord);
int find_ind(double x,Eigen::VectorXd vec);
Eigen::VectorXd get_observation_indices(Eigen::VectorXd time,int quest,Eigen::VectorXd time_n, Eigen::VectorXd quest_n);
Eigen::MatrixXd subset(Eigen::MatrixXd Gam, Eigen::VectorXd vec);
Eigen::MatrixXd elim_ind(int n, Eigen::VectorXd ind);
NumericVector ll_macml_o_StSp(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control);
Eigen::MatrixXd pred_probit_ordered_approx_StSp(Eigen::VectorXd theta, Eigen::MatrixXd Xn, Eigen::VectorXd yn, Rcpp::List mod, Eigen::VectorXd time, int quest, Eigen::VectorXd timen, Eigen::VectorXd questn);
