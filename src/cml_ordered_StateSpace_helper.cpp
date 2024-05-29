#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
//#include <omp.h>

#include "SJ_vdb.h"
#include "TVBS.h"
#include "ME.h"
#include "macml_grad_Hessian.h"
#include "lin_alg.h"
#include "toms462.h"                                         // allows bivariate normal cdf calculations


using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
//using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
//using Eigen::EigenSolver;    // one of the eigenvalue solvers


// calculate maximal eigenvalue
// [[Rcpp::export]]
double maxEV(Eigen::MatrixXf mat){
  //  // compute eigenvalues 
  Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
  eigensolver.compute(mat);
  Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real(); 
  Eigen::VectorXf eigen_values_imag = eigensolver.eigenvalues().imag(); 
  
  //Rcout << "eigen_values" << eigen_values << std::endl;
  
  //  // calc modulus
  double maxMod = std::pow(eigen_values[0],2) + std::pow(eigen_values_imag[0],2);
  //  // calc max
  int n = eigen_values.size();
  
  double max_help=0.0;

  if (n>1){
    for (int jj=1;jj<n;jj++){
      max_help = std::pow(eigen_values[jj],2) + std::pow(eigen_values_imag[jj],2);
      if (max_help > maxMod){
        maxMod = max_help;
      }
    }
  }
  
  return maxMod;
}



//////////////////////////////////////////////////////////////////////
////// helper functions to calculate the state space system       ////
////// and corresponding derivatives of the covariances with      ////
////// respect to the parameters (using echelon canonical form)   ////
//////////////////////////////////////////////////////////////////////

struct system
{
  Eigen::MatrixXd Omega, Sigma, GammaT, Sigma_chol_vech, Omega_chol_vech;
  Eigen::VectorXd beta, tauk;

  double factor;

};

struct state_space_system
{
  Eigen::MatrixXd A, C, K;

  double factor;
};


///' converts parameters to state space system matrices. 
///' @description
///' Computes the state space system with output dimension s and state dimension n, 
///' corresponding to the parameters param.  
///' @param param 
///' parameter vector
///' @param s
///' integer; output dimension 
///' @param n
///' integer; state dimension 
///' @param grad_bool
///' integer; if 0 the system is calculated, if non-zero the returned system is the derivative. The difference lies in entries fixed to one in the system which are zero for the derivative. 
///' @return
///' state_space_system structure containing the system. 
///' @keywords internal 
///' 


struct state_space_system param_to_system(Eigen::VectorXd param, int s, int n, int grad_bool, bool stationary) 
{
  // check if number of params is correct 
  int npar = 2*s*n;
  int parn = param.size();
  
  Eigen::MatrixXd A(n,n), K(n,s), C(s,n);
  A.setZero();
  K.setZero();
  C.setZero();
  
  int curp = 0; // current number of parameter.
  
  // does the vector contain the right number of parameters?
  if (param.size() == npar){
    // ok, then start with calculating the Kronecker indices 
    Eigen::VectorXd kron_index(s);
    kron_index.setZero(); 
    
    std::div_t dv{};
    dv = std::div(n, s); // division with remainder. The first dv.rem components have KI dv.quot+1, the remaining one dv.quot. 
    
    if (n >=s){ // all rows of C are filled with exactly one 1. 
      int cur = 0; // current component
      for (int co = 0; co<s;co++){ // cycle over components 
        if (grad_bool == 0) {  C(co,cur) = 1;  }
        if (co<dv.rem){
          if (grad_bool == 0) { A.block(cur,cur+1,dv.quot,dv.quot) = identity(  dv.quot);} // diagonal blocks of have shift structure. 
          cur = cur + dv.quot+1; 
        } else {
          if (grad_bool == 0) { A.block(cur,cur+1,dv.quot-1,dv.quot-1) = identity( dv.quot-1 );} // diagonal blocks of have shift structure. 
          cur = cur + dv.quot; 
        }
      }
      A.row(cur-1) = param.segment(curp,n); // last rows in each block contain parameters. 
      curp = curp + n; 
      
      // parameters for K.
      for (int co = 0; co<s;co++){ // cycle over columns
        K.col(co) = param.segment(curp,n);
        curp = curp + n; 
      }
    }
    if (n< s){ // less states than outputs 
      //int cur = 0;
      for (int co = 0; co<n;co++){ // cycle over components, but only up to n, not s!
        if (grad_bool == 0) { C(co,co) = 1;}
      }
      for (int co = n; co<s; co++){ // fill up remaining rows of C with parameters
        C.row(co)= param.segment(curp,n);
        curp = curp + n; 
      }
      
      for (int con= 0; con<n;con++){ // for A all entries are parameters. 
        A.row(con) = param.segment(curp,n);
        curp = curp + n; 
      }
      // parameters for K.
      for (int co = 0; co<s;co++){ // cycle over columns
        K.col(co) = param.segment(curp,n);
        curp = curp + n; 
      }
    } 
  }
  
  double maxev = 0; 
  
  if ((!grad_bool)&(stationary)){
    
    // correct, if system unstable::
    Eigen::MatrixXf Af = A.cast <float> ();   // Matrix of floats.

    maxev = std::sqrt(maxEV(Af));
  
    if (maxev >0.99){ // bigger than 0.99 implies instability. 
      A = A/maxev * 0.99;
    }
  }
  
  struct state_space_system syst;
  syst.A = A;
  syst.K = K;
  syst.C = C;

  syst.factor = maxev;

  return syst; 
};

//' converts parameters to state space system matrices and returns a list to R. 
//' @description
//' Computes the state space system with output dimension s and state dimension n, 
//' corresponding to the parameters param. Passes the parameters to the C++ code. 
//' @param param 
//' parameter vector
//' @param s
//' integer; output dimension 
//' @param n
//' integer; state dimension 
//' @param grad_bool
//' integer; if 0 the system is calculated, if non-zero the returned system is the derivative. The difference lies in entries fixed to one in the system which are zero for the derivative. 
//' @param stationary 
//' boolean; if true, the state is started at its stationary distribution, else at zero. 
//' @return
//' state_space_system as a list containing the system. 
//' @export
//'
// [[Rcpp::export]]
Rcpp::List  param_to_system_R(Eigen::VectorXd param, int s, int n, int grad_bool, bool stationary) 
{
   struct state_space_system syst; 
   syst = param_to_system(param, s, n, grad_bool, stationary); 
   
   Rcpp::List out; 
   out["A"] = syst.A; 
   out["C"] = syst.C; 
   out["K"] = syst.K; 
   
   return out; 

};
 
 
///' derivative of matrices w.r.t. parameters
///' @description
///' Computes the derivative of a state space system with output dimension s and state dimension n, 
///' corresponding to the parameter with number coord-1.  
///' @param coord
///' integer; index number in C terms (first entry being 0)
///' @param s
///' integer; output dimension 
///' @param n
///' integer; state dimension 
///' @return
///' state_space_system structure containing the derived system. 
///' @keywords internal 
///' 
/// [[Rcpp::export]]
 
struct state_space_system param_to_system_grad(int coord, int s, int n){
   
   int npar = 2*s*n; 
   Eigen::VectorXd grad(npar);
   grad.setZero(); 
   
   if (coord< npar){ grad(coord) = 1.0; }
   
   struct state_space_system syst;
   int grad_bool = 1;
   bool stationary = false;

   
   syst = param_to_system(grad, s, n, grad_bool, stationary ); // calculate the system corresponding to the derivative. 
   
   return syst; 
}
 
 
//' Solution to the Lyapunov equation (assumes A is stable)
//' @description
//' Computes the solution to the Lyapunov equation P = APA' + Q   
//' @param A
//' Eigen::MatrixXd; matrix describing the dynamics. 
//' @param vQ
//' Eigen::VectorXd; vectorized version of the matrix Q.  
//' @return
//' Eigen::MatrixXd; solution matrix.  
//' @keywords internal 
//' 
// [[Rcpp::export]]
Eigen::MatrixXd solve_Lyapunov_equation( Eigen::MatrixXd A, Eigen::VectorXd vQ)
{
   // initialize 
   Eigen::MatrixXd hP = A;
   hP.setZero(); 
   
   int n = A.rows(); // state dimension. 
  
   // correct A, if it is unstable 
   
   Eigen::MatrixXf Af = A.cast <float> ();   // Matrix of floats.
   
   double maxev = std::sqrt(maxEV(Af));
   
   if (maxev >0.99){ // bigger than 0.99 implies instability. 
     A = A/maxev * 0.99;
   }

   // calculate kronecker product A kron A. 
   Eigen::MatrixXd kA = kroneckerProduct(A,A);
   Eigen::MatrixXd IminuskA = identity(n*n)- kA; 
   
   Eigen::VectorXd vP = IminuskA.inverse() * vQ; // calculate the vectorized version. 
   
   //Rcout << vP << std::endl; 
   // convert to matrix
   for (int co = 0; co<n;co++){
     hP.row(co) = vP.segment(co*n,n);
   }
   
   // add to transpose to make sure it is symmetric. 
   Eigen::MatrixXd P = (hP + hP.transpose())/2; 
   
   // return value
   return P; 

}


 
 
///' Provides the matrix M = Cov(x_{t+1},y_t). 
///' @description
///' Computation of the matrix M involved in the covariance sequence.   
///' @param syst
///' state_space_system structure; contains the transfer function.  
///' @param Sigma
///' Eigen::MatrixXd; innovation variance sequence. 
///' @return
///' Eigen::MatrixXd; solution matrix.  
///' @keywords internal 
///' 
/// [[Rcpp::export]]
 
Eigen::MatrixXd calculate_M(struct state_space_system syst , Eigen::MatrixXd Sigma)
{
   // dimensions 
   int n = syst.A.rows();
   int s= syst.C.rows(); 
   
   // calculate vQ
   Eigen::MatrixXd KSig = syst.K * Sigma;
   Eigen::MatrixXd Q = KSig * syst.K.transpose(); 
   
   Eigen::VectorXd vQ(n*n);
   vQ = vectorize(Q); 
   
   // calculate P 
   Eigen::MatrixXd P = solve_Lyapunov_equation(syst.A, vQ); 
   
   // M = APC' + Sigma 
   Eigen::MatrixXd M(n,s);
   M.setZero();
   
   M = syst.A * P * syst.C.transpose() + syst.K * Sigma;
   return (M);   

}
 
//' observability matrix. 
//' @description
//' Computes the observability matrix.
//' @param A
//' Eigen::MatrixXd; matrix describing the dynamics. 
//' @param C
//' Eigen::MatrixXd; matrix describing the relation between state and output. 
//' @param dtime
//' Eigen::VectorXd; vector of integers, contains the time lags.   
//' @return
//' Eigen::MatrixXd; observability matrix.  
//' @keywords internal 
//' 
// [[Rcpp::export]]
Eigen::MatrixXd calculate_observability(Eigen::MatrixXd A,Eigen::MatrixXd C,Eigen::VectorXd dtime)
{

   // dimensions 
   int n = A.rows();
   int s = C.rows(); 
   int T = dtime.size(); 
   
   // initialize matrix
   Eigen::MatrixXd Ot(T*s,n), Aj(n,n);
   Ot.setZero(); 
   Aj = identity(n);
   
   // iterate over time: first time step minus 1!
   dtime(0)=dtime(0)-1; 
   for (int t=0;t<T;t++){
     for (int dt=0;dt<dtime(t); dt++){
       Aj = Aj * A; 
     }
     Ot.block(t*s,0,s,n) = C * Aj; 
   }
   
   return (Ot); 

}
 
//' derivative of the observability matrix. 
//' @description
//' Computes the derivative of the observability matrix. 
//' @param A
//' Eigen::MatrixXd; matrix describing the dynamics. 
//' @param C
//' Eigen::MatrixXd; matrix describing the relation between state and output. 
//' @param dA
//' Eigen::MatrixXd; derivative of the matrix describing the dynamics. 
//' @param dC
//' Eigen::MatrixXd; derivative of the matrix describing the relation between state and output. 
//' @param dtime
//' Eigen::VectorXd; vector of integers, contains the time lags between observations. 
//' @return
//' Eigen::MatrixXd; derivative of the observability matrix.  
//' @keywords internal 
//' 
// [[Rcpp::export]]
Eigen::MatrixXd calculate_deriv_observability(Eigen::MatrixXd A,Eigen::MatrixXd C,Eigen::MatrixXd dA,Eigen::MatrixXd dC,Eigen::VectorXd dtime)
{
   // dimensions 
   int n = A.rows();
   int s = C.rows(); 
   int T = dtime.size(); 
   
   // initialize matrix
   Eigen::MatrixXd Ot(T*s,n), Aj(n,n), dOt(T*s,n), dAj(n,n);

   //Ot.setZero(); 

   dOt.setZero(); 
   Aj = identity(n);
   dAj = 0*Aj; 
   // iterate over time: first time step minus 1!
   dtime(0)=dtime(0)-1; 
   for (int t=0;t<T;t++){
     for (int dt=0;dt<dtime(t); dt++){
       dAj = dAj * A + Aj * dA;
       Aj = Aj * A; 
     }

     Ot.block(t*s,0,s,n) = C * Aj;

     dOt.block(t*s,0,s,n) = dC * Aj + C * dAj; 
   }
   
   return (dOt); 
}


 ///' Calculates the covariance sequence  
 ///' @description
 ///' Computes the covariance sequence corresponding to a transfer function and a noise covariance matrix. 
 ///' @param syst
 ///' state_space_system structure; contains the transfer function.  
 ///' @param P0
 ///' Eigen::MatrixXd; matrix of initial state variance. If the (1,1) entry is negative, stationary initialization is assumed. 
 ///' @param Sigma
 ///' Eigen::MatrixXd; innovation variance sequence. 
 ///' @param time
 ///' Eigen::VectorXd; vector of integers, contains the observation times.    
 ///' @return
 ///' Eigen::MatrixXd; covariance matrix of all observations.   
 ///' @keywords internal 
 /// [[Rcpp::export]]
 
 Eigen::MatrixXd calculate_Cov_Seq(struct state_space_system syst, Eigen::MatrixXd P0, Eigen::MatrixXd Sigma, Eigen::VectorXd time)
 {
   // dimensions 
   int n = syst.A.rows();
   int s = syst.C.rows(); 
   
   int T = time.size(); 
   
   // initialize matrices
   Eigen::MatrixXd A, C, K; 
   A = syst.A;
   C = syst.C;
   K = syst.K; 
   
   //Rcout << "A:" << A << "K:" << K << "C:" << C << std::endl; 
   
   Eigen::MatrixXd KSK(n,n), KI(n+s,s), AC(n+s,n), KISKI(n+s,n+s);
   KSK = K * Sigma * K.transpose(); 
   KI.block(0,0,n,s) = K;
   KI.block(n,0,s,s) = identity(s); 
   AC.block(0,0,n,n) = A;
   AC.block(n,0,s,n) = C; 
   KISKI = KI * Sigma * KI.transpose(); 
   
   // return matrix 
   Eigen::MatrixXd GammaT(s*T,s*T);
   GammaT.setZero(); 
   
   // find time increments 
   Eigen::VectorXd dtime(T-1);
   dtime.setZero();
   for (int j=0;j<T-1;j++){
     dtime(j) = time(j+1)-time(j);
   }
   
   //Rcout << dtime << std::endl; 
   
   // find out, if P0 contains a start covariance matrix or no -> initialize with stationary distribution. 
   if (P0(0,0)<0){ 
     // stationary start. 
     Eigen::MatrixXd P0(n,n);
     P0.setZero(); 
     Eigen::VectorXd vQ = vectorize(KSK); 
     P0 = solve_Lyapunov_equation(A, vQ); 
   }
   // calculate via recursion over time, use covariance of x_{t+1},y_t. 
   // start matrix. 
   Eigen::MatrixXd Put = AC * P0 * AC.transpose() + KISKI; 
   GammaT.block(0,0,s,s) = Put.block(n,n,s,s)/2;
   Eigen::MatrixXd Mt = Put.block(0,n,n,s); 
   Eigen::MatrixXd Ot = calculate_observability(A,C,dtime);
   GammaT.block(s,0,(T-1)*s,s) = Ot * Mt; 
   
   //Rcout << Ot << std::endl; 
   
   for (int t=1; t<T; t++){ // iterate over time steps. 
     // update covariance matrix.
     Eigen::MatrixXd Putn; 
     for (int dt=0;dt<dtime(t-1); dt++){
       Putn = AC * Put.block(0,0,n,n) * AC.transpose() + KISKI;  
       Put = Putn; 
     }
     
     // extract Gamma matrix
     GammaT.block(s*t,s*t,s,s) = Put.block(n,n,s,s)/2;
     // extract Mt 
     Mt = Put.block(0,n,n,s); 
     if (T-t-1>0){
       Eigen::MatrixXd Otn = calculate_observability(A,C,dtime.segment(t,T-t-1));
       GammaT.block(s*(t+1),s*t,(T-t-1)*s,s) = Otn * Mt;      
     }
   }
   
   // fill in upper half
   Eigen::MatrixXd Gam(s*T,s*T);
   Gam = GammaT + GammaT.transpose(); 
   GammaT = Gam; 
   
   return (GammaT);  
 }


///' Calculates the derivative of the covariance sequence  
///' @description
///' Computes the derivative of the covariance sequence corresponding to a transfer function and a noise covariance matrix. 
///' @param syst
///' state_space_system structure; contains the transfer function.  
///' @param dsyst
///' state_space_system structure; contains the derivative of the transfer function.  
///' @param P0
///' Eigen::MatrixXd; matrix of initial state variance. If the (1,1) entry is negative, stationary initialization is assumed. 
///' @param dP0
///' Eigen::MatrixXd; matrix of the derivative of the initial state variance. If the (1,1) entry is negative, stationary initialization is assumed. 
///' @param Sigma
///' Eigen::MatrixXd; innovation variance sequence. 
///' @param dSigma
///' Eigen::MatrixXd; derivative of the innovation variance sequence. 
///' @param time
///' Eigen::VectorXd; vector of integers, contains the observation times.    
///' @return
///' Eigen::MatrixXd; derivative of the covariance matrix of all observations.   
///' @keywords internal
/// [[Rcpp::export]] 

Eigen::MatrixXd calculate_dCov_Seq(struct state_space_system syst, struct state_space_system dsyst,Eigen::MatrixXd P0, Eigen::MatrixXd dP0, Eigen::MatrixXd Sigma, Eigen::MatrixXd dSigma, Eigen::VectorXd time)
{
  // dimensions 
  int n = syst.A.rows();
  int s = syst.C.rows(); 
  
  int T = time.size(); 
  
  // initialize matrices
  Eigen::MatrixXd A, C, K; 
  A = syst.A;
  C = syst.C;
  K = syst.K; 
  
  //Rcout << "A:" << A << "K:" << K << "C:" << C << std::endl; 
  
  Eigen::MatrixXd dA, dC, dK; 
  dA = dsyst.A;
  dC = dsyst.C;
  dK = dsyst.K; 
  
  //Rcout << "dA:" << dA << "dK:" << dK << "dC:" << dC << std::endl; 
  
  
  Eigen::MatrixXd KSK(n,n), KI(n+s,s), AC(n+s,n), KISKI(n+s,n+s);
  KSK = K * Sigma * K.transpose(); 
  KI.block(0,0,n,s) = K;
  KI.block(n,0,s,s) = identity(s); 
  AC.block(0,0,n,n) = A;
  AC.block(n,0,s,n) = C; 
  KISKI = KI * Sigma * KI.transpose(); 
  
  Eigen::MatrixXd dKSK(n,n), dKI(n+s,s), dAC(n+s,n), dKISKI(n+s,n+s);
  dKSK = dK * Sigma * K.transpose() + K * dSigma * K.transpose() + K* Sigma * dK.transpose(); 
  dKI.block(0,0,n,s) = dK;
  dKI.block(n,0,s,s) = 0*identity(s); 
  dAC.block(0,0,n,n) = dA;
  dAC.block(n,0,s,n) = dC; 
  dKISKI = dKI * Sigma * KI.transpose() + KI * dSigma * KI.transpose() +KI * Sigma * dKI.transpose(); 
  
  
  // return matrix 
  Eigen::MatrixXd GammaT(s*T,s*T), dGammaT(s*T,s*T);
  dGammaT.setZero(); 
  GammaT.setZero();
  
  // find time increments 
  Eigen::VectorXd dtime(T-1);
  dtime.setZero();
  for (int j=0;j<T-1;j++){
    dtime(j) = time(j+1)-time(j);
  }
  
  
  // find out, if P0 contains a start covariance matrix or no -> initialize with stationary distribution. 
  //if (P0(0,0)<0){ 
  // stationary start. 
  //   Eigen::VectorXd vQ = vectorize(KSK); 
  //   P0 = solve_Lyapunov_equation(A, vQ); 
  
  // derivative 
  //   Eigen::MatrixXd dP0(n,n);
  //   dP0.setZero(); 
  
  //       dQ =  0*dKSK + dA * P0 * A.transpose() + 0*A * P0 * dA.transpose();
  
  //   vQ = vectorize(dKSK);
  //   dP0 = solve_Lyapunov_equation(A, vQ); 
  //}
  // 
  // 
  // // calculate via recursion over time, use covariance of x_{t+1},y_t. 
  // // start matrix. 
  Eigen::MatrixXd dPut = dAC * P0 * AC.transpose() + AC * dP0 * AC.transpose() + AC * P0 * dAC.transpose() +dKISKI; 
  Eigen::MatrixXd Put = AC * P0 * AC.transpose() + KISKI; 
  
  
  GammaT.block(0,0,s,s) = Put.block(n,n,s,s)/2;
  dGammaT.block(0,0,s,s) = dPut.block(n,n,s,s)/2;
  
  Eigen::MatrixXd Mt = Put.block(0,n,n,s); 
  Eigen::MatrixXd dMt = dPut.block(0,n,n,s); 
  
  // covariances 
  Eigen::MatrixXd Ot = calculate_observability(A,C,dtime);
  GammaT.block(s,0,(T-1)*s,s) = Ot * Mt; 
  
  // derivatives 
  Eigen::MatrixXd dOt = calculate_deriv_observability(A,C,dA,dC,dtime);
  dGammaT.block(s,0,(T-1)*s,s) = dOt * Mt + Ot * dMt; 
  
  for (int t=1; t<T; t++){ // iterate over time steps. 
    // update covariance matrix.
    Eigen::MatrixXd Putn, dPutn; 
    for (int dt=0;dt<dtime(t-1); dt++){
      dPutn = dAC * Put.block(0,0,n,n) * AC.transpose() + AC * dPut.block(0,0,n,n) * AC.transpose() + AC * Put.block(0,0,n,n) * dAC.transpose() + dKISKI;
      Putn = AC * Put.block(0,0,n,n) * AC.transpose() + KISKI;  
      Put = Putn; 
      dPut = dPutn; 
    }
    
    // extract Gamma matrix
    GammaT.block(s*t,s*t,s,s) = Put.block(n,n,s,s)/2;
    dGammaT.block(s*t,s*t,s,s) = dPut.block(n,n,s,s)/2;
    // extract Mt 
    Mt = Put.block(0,n,n,s); 
    dMt = dPut.block(0,n,n,s);
    
    if (T-t-1>0){
      Eigen::MatrixXd Otn = calculate_observability(A,C,dtime.segment(t,T-t-1));
      GammaT.block(s*(t+1),s*t,(T-t-1)*s,s) = Otn * Mt;  
      // deriv 
      Eigen::MatrixXd dOtn = calculate_deriv_observability(A,C,dA,dC,dtime.segment(t,T-t-1));
      dGammaT.block(s*(t+1),s*t,(T-t-1)*s,s) = dOtn * Mt + Otn * dMt; 
    }
  }
  
  
  // fill in upper half
  Eigen::MatrixXd Gam(s*T,s*T);
  Gam = GammaT + GammaT.transpose(); 
  GammaT = Gam; 
  
  // fill in upper half
  Eigen::MatrixXd dGam(s*T,s*T);
  dGam =  dGammaT + dGammaT.transpose(); 
  dGammaT = dGam; 
  
  return (dGammaT);  
}

  

