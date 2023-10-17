
#include <RcppEigen.h>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd; 
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int precision

#include "SJ_vdb.h"
#include "TVBS_vdb.h"
#include "TVBS.h"
#include "lin_alg.h"
#include "distrib.h"

#include "toms462.h"

Eigen::MatrixXd J_diff(int mm, int kk);
Eigen::MatrixXd J_2dM(int mm1, int mm2, int k1, int k2);
Eigen::MatrixXd J_theta(Eigen::MatrixXd X, Eigen::VectorXd b, Eigen::MatrixXd sigma);
Eigen::MatrixXd H_theta(Eigen::MatrixXd X, Eigen::VectorXd b, Eigen::MatrixXd sigma);
Eigen::MatrixXd J_chol(Eigen::VectorXd x);  
Eigen::MatrixXd J_chol1(Eigen::VectorXd x, Eigen::MatrixXd L, Eigen::MatrixXd K, Eigen::MatrixXd D);
Eigen::MatrixXd J_bvOvSig_bLOLSig(Eigen::VectorXd tho, Eigen::VectorXd thl, int M, Eigen::MatrixXd Lo, Eigen::MatrixXd Ko, Eigen::MatrixXd Do, Eigen::MatrixXd Ll, Eigen::MatrixXd Kl, Eigen::MatrixXd Dl);
Eigen::MatrixXd J_bvSig_bLSig(Eigen::VectorXd thl, int M, Eigen::MatrixXd Ll, Eigen::MatrixXd Kl, Eigen::MatrixXd Dl);
Eigen::MatrixXd H_chol(Eigen::MatrixXd J);
Eigen::MatrixXd H_bvOvSig_bLOLSig(Eigen::MatrixXd J, int M, int O, int L);
Eigen::MatrixXd J_Omega(Eigen::MatrixXd X);
Eigen::MatrixXd J_Omega1(Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D);
Eigen::MatrixXd makeM(int m, int y);
Eigen::MatrixXd J1(int mm1, int mm2, int k1, int k2, int M, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D);
Eigen::MatrixXd J_TpdM(int alt, int Tp, Eigen::VectorXd y);
Eigen::MatrixXd J_TpdM_red(int alt, Eigen::VectorXd y);
Eigen::MatrixXd J_M_bOmSig(int alt, int Tp, Eigen::VectorXd y, int lthbb, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D);
Eigen::MatrixXd J_M_bOmSig_red(int alt, IntegerVector ind_BI, Eigen::VectorXd y, int lthbb, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D);
Eigen::MatrixXd J_bvOv(Eigen::VectorXd tho, Eigen::VectorXd thl, int M, Eigen::MatrixXd Lo, Eigen::MatrixXd Ko, Eigen::MatrixXd Do);
double prob_ordered_2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk);
Eigen::VectorXd  grad_po2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk);
Eigen::MatrixXd  hess_po2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk);
double prob_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type);
Eigen::VectorXd grad_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type);
Eigen::MatrixXd hess_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type);
Eigen::MatrixXd hess_ordered_CML_approx(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type);
double prob_ordered_1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk);
Eigen::VectorXd grad_po1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk);
Eigen::MatrixXd hess_po1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk);
