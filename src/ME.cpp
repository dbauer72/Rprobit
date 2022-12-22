#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>


using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Matrix3d;    // one of the eigenvalue solvers

using namespace Eigen;
using Eigen::Map;                                               // 'maps' rather than copies
#include "SJ_vdb.h"
#include "toms462.h"                                          // allows bivariate normal cdf calculations
#include "distrib.h"

extern double bound_max;
extern double tol;

// [[Rcpp::plugins(cpp11)]]


// returns the index of the element in the upper diagonal of a mxm square matrix corresponding to row r and col c including the diagonals
int ind_diag(int r, int c, int m)
{
  int ind=0;
  int curr = 0;
  for (int jj=0;jj<r;jj++){
    curr += m-jj;
  }

  ind = curr+c-r;
  return ind;
}

// returns the index of the element in the upper diagonal of a mxm square matrix corresponding to row r and col c excluding the diagonals
int ind_nodiag(int r, int c, int m)
{
  int ind=0;
  int curr = 0;
  for (int jj=0;jj<r;jj++){
    curr += m-jj-1;
  }

  ind = curr+c-r-1;
  return ind;
}


//////////////////////////////////////////////////////////////
//' Mendell-Elston approximation to multivariate Gaussian CDF.
//' @description
//' The function computes the  CDF for a multivariate Gaussian distribution approximated according to the method of Mendell-Elston.
//' @param x
//' nx1 vector of point to evaluate the CDF at.
//' @param r
//' nxn correlation matrix.
//' @return
//' double; log of probability.
//' @keywords internal
//'
// [[Rcpp::export]]
double ME(Eigen::VectorXd x, Eigen::MatrixXd r){

  double lcond=0;

  int m = x.rows();

  for (int i = 0; i < m; i++)
  {                                                                 // for i in 1:m               essentially limit all x(i) to [-6,6]
    if (x(i) < -bound_max)
    {                                                                   // if x(i) < -6
      x(i) = -bound_max;                                                            // set x(i) = -6
    }
    if (x(i) > bound_max)
    {                                                                   // if x(i) > 6
      x(i) = bound_max;                                                             // set x(i) = 6
    }
  }

  // initialize vectors and matrices to work with
  Eigen::VectorXd Z(m);
  Z.setZero();
  Z=x;

  Eigen::MatrixXd Rij(m,m);
  Rij.setZero();
  Rij = r;

  Eigen::VectorXd sig(m);
  sig.setZero();

  double Ujjm1 =0;
  double U =0;
  double  ssquare = 0;

  for (int jj=0; jj<m-1 ; jj++)
    { // cycle over entries
    lcond += std::log(std_normal_cdf(Z[jj]));
    Ujjm1 = std_normal_pdf(Z[jj])/std_normal_cdf(Z[jj]);
    U = Ujjm1*(Ujjm1+Z[jj]);
    for (int ij=jj+1;ij<m;ij++){ // update all elements upstream
      ssquare = 1- Rij(jj,ij)*Rij(jj,ij)*U;
      sig[ij] = std::sqrt(1- Rij(jj,ij)*Rij(jj,ij)*U);
      Z[ij] = (Z[ij]+Ujjm1*Rij(jj,ij))/sig[ij];
    }
    for (int mj=jj+1;mj<m;mj++){
      for (int nj=mj;nj<m;nj++){
        Rij(mj,nj) = (Rij(mj,nj) - Rij(jj,mj)*Rij(jj,nj)*U)/(sig[mj]*sig[nj]);
      }
    }
  }

  lcond += std::log(std_normal_cdf(Z[m-1]));

  return lcond;
}

//////////////////////////////////////////////////////////////
//' Gradient calculation for Mendell-Elston approximation to multivariate Gaussian CDF.
//' @description
//' The function computes the gradient of the CDF for a multivariate Gaussian distribution approximated according to the method of Mendell-Elston.
//' @param x
//' nx1 vector of point to evaluate the CDF at.
//' @param r
//' nxn correlation matrix.
//' @return
//' vector; gradient of log of probability.
//' @keywords internal
//'
// [[Rcpp::export]]
Eigen::VectorXd dlcond_ME(Eigen::VectorXd x,Eigen::MatrixXd r){
  int m = x.rows();                                                 // save the number of rows of x in m


  for (int i = 0; i < m; i++)
  {                                                                 // for i in 1:m               essentially limit all x(i) to [-6,6]
    if (x(i) < -bound_max)
    {                                                                   // if x(i) < -6
      x(i) = -bound_max;                                                            // set x(i) = -6
    }
    if (x(i) > bound_max)
    {                                                                   // if x(i) > 6
      x(i) = bound_max;                                                             // set x(i) = 6
    }
  }

  // define gradient and start calculations.
  double lcond =0.0;

  int nd;
  nd = m + m*(m-1)/2;

  Eigen::VectorXd dU(nd), dUUZ(nd), dlcond(nd);
  dlcond.setZero();
  dU.setZero();
  dUUZ.setZero();

  Eigen::MatrixXd dsig(m,nd),dZ(m,nd), dRij(nd,nd);
  dsig.setZero();
  dZ.setZero();
  dRij.setZero();

  // initialize vectors and matrices to work with
  Eigen::VectorXd Z(m);
  Z.setZero();
  Z=x;

  for (int jj=0;jj<m;jj++){
    dZ(jj,jj)=1;
  }

  // gradient for entries of correlation matrix.
  Eigen::MatrixXd Rij(m,m);
  Rij.setZero();
  Rij = r;

  for (int mj=0;mj<m;mj++){
    for (int nj=mj+1;nj<m;nj++){
      dRij(ind_diag(mj,nj,m),ind_nodiag(mj,nj,m)+m)=1;
    }
  }


  Eigen::VectorXd sig(m), dsigsquare(nd);
  sig.setZero();
  dsigsquare.setZero();

  // now gradients are initialized.
  double Ujjm1 =0;
  double logpj = 0;
  Ujjm1 = std_normal_pdf(Z[0])/std_normal_cdf(Z[0]);

  int curjj = 0;
  for (int jj=0; jj<m-1 ; jj++)
  { // cycle over entries
    logpj = std::log(std_normal_cdf(Z[jj]));

    lcond += logpj;
    dlcond += std_normal_pdf(Z[jj])*dZ.row(jj).transpose()/std::exp(logpj);

    Ujjm1 = std_normal_pdf(Z[jj])/std_normal_cdf(Z[jj]);

    dU = -dZ.row(jj).transpose()*(Ujjm1*Z[jj]+Ujjm1*Ujjm1);
    dUUZ = 2*dU*Ujjm1+ Ujjm1*dZ.row(jj).transpose()+dU*Z[jj];


    for (int ij=jj+1;ij<m;ij++){ // update all elements upstream
      sig[ij] = std::sqrt(1- Rij(jj,ij)*Rij(jj,ij)*Ujjm1*(Ujjm1+Z[jj]));
      curjj = ind_diag(jj,ij,m);

      dsigsquare = 2*Rij(jj,ij)*Ujjm1*(Ujjm1+Z[jj])*dRij.row(curjj).transpose()+Rij(jj,ij)*Rij(jj,ij)*dUUZ;
      dsig.row(ij) = -dsigsquare.transpose()/(2*sig[ij]);

      Z[ij] = (Z[ij]+Ujjm1*Rij(jj,ij))/sig[ij];
      dZ.row(ij)= (dZ.row(ij) + dU.transpose()*Rij(jj,ij)+Ujjm1*dRij.row(curjj))/sig[ij] -Z[ij]/sig[ij]*dsig.row(ij);
    }
    for (int mj=jj+1;mj<m;mj++){
      for (int nj=mj;nj<m;nj++){
        int currmn = ind_diag(mj,nj,m);
        int currjm = ind_diag(jj,mj,m);
        int currjn = ind_diag(jj,nj,m);

        Rij(mj,nj) = (Rij(mj,nj) - Rij(jj,mj)*Rij(jj,nj)*Ujjm1*(Ujjm1+Z[jj]))/(sig[mj]*sig[nj]);

        dRij.row(currmn) = (dRij.row(currmn) - (dRij.row(currjm)*Rij(jj,nj)+dRij.row(currjn)*Rij(jj,mj))*Ujjm1*(Ujjm1+Z[jj])-Rij(jj,mj)*Rij(jj,nj)*dUUZ.transpose())/(sig[mj]*sig[nj]);
        dRij.row(currmn) -= Rij(mj,nj)/sig[mj]*dsig.row(mj) + Rij(mj,nj)/sig[nj]*dsig.row(nj);
      }
    }

  }

  logpj = std::log(std_normal_cdf(Z[m-1]));
  lcond += logpj;
  dlcond += std_normal_pdf(Z[m-1])*dZ.row(m-1).transpose()/std::exp(logpj);

  Eigen::VectorXd output_final(1+dlcond.size());
  output_final.setZero();

  output_final << dlcond, lcond;
  return output_final;

}

//////////////////////////////////////////////////////////////
//' Hessian for Mendell-Elston approximation to multivariate Gaussian CDF.
//' @description
//' The function computes the Hessian of the CDF for a multivariate Gaussian distribution approximated according to the method of Mendell-Elston.
//' @param x
//' nx1 vector of point to evaluate the CDF at.
//' @param r
//' nxn correlation matrix.
//' @return
//' Matrix; Hessian of log of probability.
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd ME_hess_new(Eigen::VectorXd x,Eigen::MatrixXd r)
{                                                               // returns a column vector of doubles of variable length, input are a vector x and a quadratic matrix r

 int m = x.rows();                                                 // save the number of rows of x in m


  for (int i = 0; i < m; i++)
  {                                                                 // for i in 1:m               essentially limit all x(i) to [-6,6]
    if (x(i) < -bound_max)
    {                                                                   // if x(i) < -6
      x(i) = -bound_max;                                                            // set x(i) = -6
    }
    if (x(i) > bound_max)
    {                                                                   // if x(i) > 6
      x(i) = bound_max;                                                             // set x(i) = 6
    }
  }

  // define gradient and start calculations.
  double lcond;

  int nd;
  nd = m + m*(m-1)/2;

  Eigen::MatrixXd Hess(nd,nd);
  Hess.setZero();

  Eigen::VectorXd dU(nd), dUUZ(nd), dlcond(nd);
  dlcond.setZero();
  dU.setZero();
  dUUZ.setZero();

  Eigen::MatrixXd dsig(m,nd),dZ(m,nd), dRij(nd,nd), HU(nd,nd), HUUZ(nd,nd), Hsig(m*nd,nd),HZ(m*nd,nd), HRij(nd*nd,nd);
  dsig.setZero();
  dZ.setZero();
  dRij.setZero();
  HU.setZero();
  HUUZ.setZero();
  Hsig.setZero();
  HZ.setZero();
  HRij.setZero();

  // initialize vectors and matrices to work with
  Eigen::VectorXd Z(m);
  Z.setZero();
  Z=x;

  for (int jj=0;jj<m;jj++){
    dZ(jj,jj)=1;
  }

  // gradient for entries of correlation matrix.
  Eigen::MatrixXd Rij(m,m);
  Rij.setZero();
  Rij = r;

  for (int mj=0;mj<m;mj++){
    for (int nj=mj+1;nj<m;nj++){
      dRij(ind_diag(mj,nj,m),ind_nodiag(mj,nj,m)+m)=1;
    }
  }

  //Rcout << dRij << std::endl;

  Eigen::VectorXd sig(m), dsigsquare(nd);
  sig.setZero();
  dsigsquare.setZero();

  Eigen::MatrixXd Hsigsquare(nd,nd);
  Hsigsquare.setZero();

  Eigen::MatrixXd HQ(nd,nd),HS(nd,nd);
  HQ.setZero();
  HS.setZero();

  Eigen::VectorXd dQ(nd);
  dQ.setZero();

    // now gradients are initialized.
  double Ujjm1 =0;
  double logpj = 0;
  double UUZ =0;

  int curjj = 0;
  for (int jj=0; jj<m-1 ; jj++)
  { // cycle over entries
    logpj = std::log(std_normal_cdf(Z[jj]));

    lcond += logpj;
    dlcond += std_normal_pdf(Z[jj])*dZ.row(jj).transpose()/std::exp(logpj);
    Ujjm1 = std_normal_pdf(Z[jj])/std_normal_cdf(Z[jj]);
    Hess += Ujjm1*HZ.block(jj*nd,0,nd,nd)- dZ.row(jj).transpose() * dZ.row(jj)*(Ujjm1*Z[jj]+Ujjm1*Ujjm1);


    dU = -dZ.row(jj)*(Ujjm1*Z[jj]+Ujjm1*Ujjm1);
    HU = -HZ.block(jj*nd,0,nd,nd)*(Ujjm1*Z[jj]+Ujjm1*Ujjm1) -dZ.row(jj).transpose()*(dU.transpose()*Z[jj]+Ujjm1*dZ.row(jj)+2*Ujjm1*dU.transpose());


    UUZ = (Ujjm1*Z[jj]+Ujjm1*Ujjm1);
    dUUZ = 2*dU*Ujjm1+ Ujjm1*dZ.row(jj).transpose()+dU*Z[jj];
    HUUZ = 2*HU*Ujjm1 + 2*dU*dU.transpose()+ dU*dZ.row(jj) + Ujjm1*HZ.block(jj*nd,0,nd,nd)+HU*Z[jj]+dU*dZ.row(jj);

    for (int ij=jj+1;ij<m;ij++){ // update all elements upstream
      sig[ij] = std::sqrt(1- Rij(jj,ij)*Rij(jj,ij)*Ujjm1*(Ujjm1+Z[jj]));
      curjj = ind_diag(jj,ij,m);

      dsigsquare = 2*Rij(jj,ij)*Ujjm1*(Ujjm1+Z[jj])*dRij.row(curjj).transpose()+Rij(jj,ij)*Rij(jj,ij)*dUUZ;
      Hsigsquare = 2*(dRij.row(curjj).transpose()*Ujjm1*(Ujjm1+Z[jj]) +  Rij(jj,ij)*dUUZ)*dRij.row(curjj) +  2*Rij(jj,ij)*Ujjm1*(Ujjm1+Z[jj])*HRij.block(curjj*nd,0,nd,nd);
      Hsigsquare += 2*dRij.row(curjj).transpose()*Rij(jj,ij)*dUUZ.transpose() + Rij(jj,ij)*Rij(jj,ij)*HUUZ;

      dsig.row(ij) = -dsigsquare.transpose()/(2*sig[ij]);
      Hsig.block(ij*nd,0,nd,nd) = -Hsigsquare/(2*sig[ij]) - dsigsquare*dsigsquare.transpose()/(4*pow(sig[ij],3));

      Z[ij] = (Z[ij]+Ujjm1*Rij(jj,ij))/sig[ij];
      dQ = dZ.row(ij) + dU.transpose()*Rij(jj,ij)+Ujjm1*dRij.row(curjj);
      dZ.row(ij)= (dZ.row(ij) + dU.transpose()*Rij(jj,ij)+Ujjm1*dRij.row(curjj))/sig[ij] -Z[ij]/sig[ij]*dsig.row(ij);
      HQ = HZ.block(ij*nd,0,nd,nd)+HU*Rij(jj,ij)+dU*dRij.row(curjj)+dRij.row(curjj).transpose()*dU.transpose()+Ujjm1*HRij.block(curjj*nd,0,nd,nd);

      HZ.block(ij*nd,0,nd,nd) = HQ/sig[ij]-(dQ*dsig.row(ij)+dsig.row(ij).transpose()*dQ.transpose())/(sig[ij]*sig[ij]);
      HZ.block(ij*nd,0,nd,nd) += 2*Z[ij]/(sig[ij]*sig[ij])*dsig.row(ij).transpose()*dsig.row(ij)-Z[ij]/sig[ij]*Hsig.block(ij*nd,0,nd,nd);
    }
    for (int mj=jj+1;mj<m;mj++){
      for (int nj=mj;nj<m;nj++){
        int currmn = ind_diag(mj,nj,m);
        int currjm = ind_diag(jj,mj,m);
        int currjn = ind_diag(jj,nj,m);

        double  Q = Rij(mj,nj) - Rij(jj,mj)*Rij(jj,nj)*Ujjm1*(Ujjm1+Z[jj]);
        Rij(mj,nj) = (Rij(mj,nj) - Rij(jj,mj)*Rij(jj,nj)*Ujjm1*(Ujjm1+Z[jj]))/(sig[mj]*sig[nj]);
        double S = sig[mj]*sig[nj];

        Eigen::VectorXd dQ(nd), dS(nd);
        dQ.setZero();
        dS.setZero();
        dQ = (dRij.row(currmn) - (dRij.row(currjm)*Rij(jj,nj)+dRij.row(currjn)*Rij(jj,mj))*Ujjm1*(Ujjm1+Z[jj])-Rij(jj,mj)*Rij(jj,nj)*dUUZ.transpose());
        dS = sig[mj]*dsig.row(nj) + sig[nj]*dsig.row(mj);

        dRij.row(currmn) = dQ/S - Rij(mj,nj)/S*dS;

        // numerator: Rij(mj,nj) - Rij(jj,mj)*Rij(jj,nj)*UUZ
        HQ.setZero();
        HQ = HRij.block(currmn*nd,0,nd,nd) - HRij.block(currjm*nd,0,nd,nd)*Rij(jj,nj)*UUZ - dRij.row(currjm).transpose()*(dRij.row(currjn)*UUZ+ Rij(jj,nj)*dUUZ.transpose());
        HQ -= dRij.row(currjn).transpose()*(dRij.row(currjm)*UUZ+ Rij(jj,mj)*dUUZ.transpose()) + HRij.block(currjn*nd,0,nd,nd)*Rij(jj,mj)*UUZ;
        HQ -= dUUZ*(dRij.row(currjm)*Rij(jj,nj)+Rij(jj,mj)*dRij.row(currjn)) + Rij(jj,mj)*Rij(jj,nj)*HUUZ;
        // denominator: sigma_m sigma_n.
        HS = Hsig.block(mj*nd,0,nd,nd)*sig[nj]+dsig.row(mj).transpose()*dsig.row(nj)+dsig.row(nj).transpose()*dsig.row(mj)+ sig[mj]*Hsig.block(nj*nd,0,nd,nd);

        // combine numerator and denominator
        HRij.block(currmn*nd,0,nd,nd) = HQ/S- dQ*dS.transpose()/(S*S) - dS*dQ.transpose()/(S*S)- Q/(S*S)*HS + 2*Q*dS*dS.transpose()/pow(S,3);

      }
    }

  }

  logpj = std::log(std_normal_cdf(Z[m-1]));
  lcond += logpj;
  dlcond +=  std_normal_pdf(Z[m-1])*dZ.row(m-1).transpose()/std::exp(logpj);

  Ujjm1 = std_normal_pdf(Z[m-1])/std_normal_cdf(Z[m-1]);
  Hess += HZ.block((m-1)*nd,0,nd,nd)*Ujjm1 -dZ.row(m-1).transpose() * dZ.row(m-1)*(Ujjm1*Z[m-1]+Ujjm1*Ujjm1);

  return Hess;
}

