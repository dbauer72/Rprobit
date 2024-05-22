#include <RcppEigen.h>
#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;

Eigen::VectorXd int2bin(int num,int noBits);
Eigen::VectorXd vectorize( Eigen::MatrixXd Q);
Eigen::MatrixXd identity( int n);
Eigen::MatrixXd duplmat( int n);
Eigen::MatrixXd vechor(int k);
Eigen::VectorXd fdiag(int k);
Eigen::MatrixXd vechor_diag(int k);
Eigen::MatrixXd commmat( int r, int c);
Eigen::MatrixXd elimmat( int n);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
IntegerVector NA_ind(Eigen::MatrixXd Xc);

Eigen::MatrixXd reorder_matrix(Eigen::MatrixXd M, Eigen::VectorXd indord);
Eigen::MatrixXd extract_submatrix(Eigen::MatrixXd Mold, Eigen::VectorXd indin);
Eigen::VectorXd reorder_indices(Eigen::VectorXd ind_ord);
///Eigen::MatrixXd J_wR_bSig(Eigen::VectorXd b, Eigen::MatrixXd Sigma);
