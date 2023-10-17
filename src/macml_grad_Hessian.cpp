
#include <RcppEigen.h>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd; 
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int precision

#include "SJ_vdb.h"
#include "TVBS_vdb.h"
#include "distrib.h"
#include "lin_alg.h"

#include "toms462.h"

// This file contains helper functions, that are needed to 
// calculate the gradient and the Hessian of the probit and MacML likelihoods.
//
// J_XY denotes functions for Jacobians needed for gradients and Hessian
// H_XY denotes functions of second derivatives. 




Eigen::MatrixXd J_diff(int mm, int kk){
  //Computes the first derivative (Jacobian) of  vech(M_k Sig M_k') w.r.t. vech(Sig), where Sig is a mm X mm matrix
  //Note that vech(M_k Sig M_k) has dimension "smaller by one"
  //Note that the Hessian is zero
  //Diff w.r.t kk (the function uses k = kk-1 to comply with cpp indexing)!
  //mm is the dimension of the matrix
  
  int m =  (((mm*mm)-mm)/2 + mm);
  int nn = mm - 1;
  int n =  (((nn*nn)-nn)/2 + nn);
  
  Eigen::MatrixXd J(n,m);
  J.setZero();
  
  //compute entries of Sig
  
  
  Eigen::MatrixXd indSig = vechor_diag(mm);
  
  int k = kk-1;
  
  //compute entries of M*Sig*t(M)
  
  VectorXi row_selected = (indSig.col(0).array() != k).cast<int>();
  VectorXi col_selected = (indSig.col(1).array() != k).cast<int>();
  VectorXi selected = (row_selected.array() * col_selected.array());
  
  Eigen::MatrixXd indMSigM (selected.sum(), 2);
  int rownew = 0;
  for (int i = 0; i < indSig.rows(); ++i) {
    if (selected[i]) {
      indMSigM.row(rownew) = indSig.row(i);
      rownew++;
    }
  }
  
  
  
  //j is the rows of J
  //i is the cols of J
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < n; ++j){
      
      //Every entry gets diffed by the diag entry
      if((indSig(i,0)==k)&(indSig(i,1)==k)){
        J(j,i) = 1;
      }
      
      //Derivative for each entry with itself
      if(indSig.row(i) == indMSigM.row(j)){
        J(j,i) = 1;
      }
      
      //Off-Diagnoal entry in the diffed matrix I
      if((indSig(i,1)==k)&(indSig(i,0)==indMSigM(j,0))){
        J(j,i) = -1;
      }
      
      //Off-Diagnoal entry in the diffed matrix II
      if((indSig(i,1)==k)&(indSig(i,0)==indMSigM(j,1))){
        J(j,i) = -1;
      }
      
      //Off-Diagnoal entry in the diffed matrix III
      if((indSig(i,0)==k)&(indSig(i,1)==indMSigM(j,1))){
        J(j,i) = -1;
      }
      
      //Off-Diagnoal entry in the diffed matrix IV
      if((indSig(i,0)==k)&(indSig(i,1)==indMSigM(j,0))){
        J(j,i) = -1;
      }
      
      //Diagnoal entry in the diffed matrix I
      if((indSig(i,0)==k)&(indSig(i,1)==indMSigM(j,0))&(indMSigM(j,0)==indMSigM(j,1))){
        J(j,i) = -2;
      }
      
      //Diagnoal entry in the diffed matrix II
      if((indSig(i,1)==k)&(indSig(i,0)==indMSigM(j,0))&(indMSigM(j,0)==indMSigM(j,1))){
        J(j,i) = -2;
      }
    }
    
  }
  
  
  return J;
}



Eigen::MatrixXd J_2dM(int mm1, int mm2, int k1, int k2){
  //Computes the first derivative (Jacobian) of  M w.r.t. LL',
  // M = diag( M_k1 Sigma M_k1', M_k2 Sigma M_k2'). 
  //
  //AT THE MOMENT ONLY SITUATIONS WHERE THE PERSON HAS THE SAME NUMBER OF CHOICE ALTERNATIVES
  // where Sig is a LL matrix which appears as a blockdiagonal with different differencing
  // mm1, mm2 is dimnesion of LL (so, equal to the number of choice alternatives at the time point), where M is of dimension (mm1-1) + (mm2- 1)
  Eigen::MatrixXd J1 = J_diff(mm1, k1);
  Eigen::MatrixXd  J2 = J_diff(mm2, k2);
  
  Eigen::MatrixXd J((2*(mm1-1)*2*(mm2-1)-2*(mm2-1))/2 + 2*(mm2-1), J1.cols());
  J.setZero();
  
  Eigen::MatrixXd fill((mm2- 1), J1.cols());
  fill.setZero();
  
  Eigen::VectorXd ind = vechor_diag(mm1-1).col(1).array();
  
  int k = ind[0];
  int g = 0;
  
  //Add the Jacobian for the first block and tamper with zeros after each column
  for(int i = 0; i < ind.rows(); ++i){
    if(k == ind[i]){
      J.row(g) += J1.row(i);
      g++;
    }else{
      g = g + (mm2-1);
      J.row(g) += J1.row(i);
      g++;
      k = ind[i];
    }
  }
  
  g = g + (mm2-1);
  
  //Add the Jacobian for the second block
  
  for(int i = 0; i < J2.rows(); ++i){
    J.row(g) += J2.row(i);
    g++;
  }
  
  return J;
  
}



Eigen::MatrixXd J_theta(Eigen::MatrixXd X, Eigen::VectorXd b, Eigen::MatrixXd sigma){
  
  //x: covariate matrix
  //b: linear coefficients (k)
  //sigma: covariance matrix XOOX + LL (k x k)
  
  //returns the jacobian and w.r.t. b and vech(M), M = XOOX + LL (k x k)
  //Hence, moving vom theta = [b/diag(sqrt(s_11)), ..., s*_12] to b, vech(M)
  
  int n = X.cols(); //picking the length of b from X because problems with t(b) are avoided
  
  
  Eigen::VectorXd sig = sigma.diagonal();
  int k = sigma.rows();
  
  Eigen::VectorXd kk = (-1*(X*b).array());
  Eigen::VectorXd kk1 = kk.array()/sig.array().sqrt();
  Eigen::MatrixXd J(k+((k*k)-k)/2, n+((k*k)-k)/2+k);
  J.setZero();
  
  //Derivative of the lin-utilities w.r.t. b
  Eigen::MatrixXd XX = X;
  XX.array().colwise() /= sig.array().sqrt();
  J.block(0,0,X.rows(), X.cols()) = -1*XX;
  
  
  int m = n; //cpp indexing (-1)-R
  
  
  Eigen::VectorXd id(k);
  id.setZero();
  
  //Derivative of the lin-utilities w.r.t. diag(sqrt(Sig))
  
  Eigen::VectorXd ssig = 1/sig.array().pow(1.5);
  Eigen::VectorXd sqsig = sig.array().sqrt();
  
  
  for(int i = 0; i < k; ++i){
    id(i) = m;
    
    J(i,m) = -0.5*kk[i]*ssig[i]; //using square bracets for clarification
    int mm = m + (k-i);
    m = mm;
  }
  
  Eigen::MatrixXd ind = vechor(k);
  
  //Derivative of the off-diagonals w.r.t. diag(sqrt(Sig)) 
  //Could be done in one if-clause but it works and hence will not be changed
  for(int i = 0; i < k; ++i){ //diagonal element
    for(int j = k; j < (k+((k*k)-k)/2); ++j){ //off-diagonal elements
      
      if((ind((j-k),0)==i) & (ind((j-k),1)!=i)){   
        
        J(j,id[i]) = (sigma(ind((j-k),0),ind((j-k),1))/sqsig[ind((j-k),1)])*-0.5*ssig[i];
      }
      
      if((ind((j-k),1)==i) & (ind((j-k),0)!=i)){
        J(j,id[i]) = (sigma(ind((j-k),0),ind((j-k),1))/sqsig[ind((j-k),0)])*(-0.5*ssig[i]);
      }
      
    }
  }
  
  Eigen::MatrixXd ind1 = vechor_diag(k);
  
  //Derivative of the off-diagonals w.r.t. the off-diagonals
  
  for(int i = k; i < J.rows(); ++i){ //loop over rows
    for(int j = n; j < J.cols(); ++j){ //loop over cols
      
      if((ind(i-k,0) == ind1(j-n,0))&(ind(i-k,1) == ind1(j-n,1))){
        J(i,j) = 1/(sqsig[ind((i-k),0)]*sqsig[ind((i-k),1)]);
      }
      
    }
  }
  
  
  return J;
}



Eigen::MatrixXd H_theta(Eigen::MatrixXd X, Eigen::VectorXd b, Eigen::MatrixXd sigma){
  
  //x: covariate matrix
  //b: linear coefficients (k)
  //sigma: covariance matrix XOOX + LL (k x k)
  
  //returns the Hessian and w.r.t. b and vech(sigma), sigma = XOOX + LL (k x k)
  //Hence, moving vom theta = [b/diag(sqrt(s_11)), ..., s*_12] to b, vech(sigma)
  
  int n = X.cols(); //picking the length of b from X because problems with t(b) are avoided
  Eigen::VectorXd sig = sigma.diagonal();
  int k = sigma.rows();
  
  //Set up the Hessian
  int h1 = k+((k*k)-k)/2;
  int h2 = n+((k*k)-k)/2+k;
  Eigen::MatrixXd H(h2*h1,h2);
  H.setZero();
  
  //identify indeices and diags
  Eigen::MatrixXd ind = vechor(k);
  Eigen::MatrixXd ind1 = vechor_diag(k);
  Eigen::VectorXd dd = fdiag(k);
  
  //precompute relevant values
  Eigen::VectorXd kk = (-1*(X*b).array());
  Eigen::VectorXd ssig15 = 1/sig.array().pow(1.5);
  Eigen::VectorXd ssig25 = 1/sig.array().pow(2.5);
  Eigen::MatrixXd ssigma15 = 1/sigma.array().pow(1.5);
  Eigen::MatrixXd ssigma25 = 1/sigma.array().pow(2.5);
  Eigen::MatrixXd ssigma = 1/sigma.array().sqrt();
  
  
  int bb = 0;
  for(int i = 0; i < h1; ++i){ //loop over rows
    
    //Hessian for d b / d diag and d diag / d diag
    //dbdb = 0
    //dbdOD = 0
    
    if(i < (k)){
      
      H.block(bb,dd[i]+n, n, 1) = ((-1)*X.row(i)*(-0.5*ssig15[i])).transpose();
      H.block(dd[i]+n+bb,0, 1, n) = ((-1)*X.row(i)*(-0.5*ssig15[i]));
      H((dd[i]+n+bb),(dd[i]+n)) = (kk[i])*((0.75)*ssig25[i]);
      
    } else {
      
      int od = (i -k) + 1 + ind((i-k),1); //off-diagnal entry (shifted from the original index by the diags, hence 1 plus the number of cols) 
      int d1 = dd(ind((i-k),1)); //diag entries corresponding to the off-diagnal entry 
      int d2 = dd(ind((i-k),0));
      
      //d od/d d1 (d1 and d2 are diags)
      H(bb+n+od, n+d1) = ssigma(ind1(d2,0), ind1(d2,0)) *(-0.5)*ssigma15(ind1(d1,0), ind1(d1,0));
      H(bb+n+d1, n+od) = H(bb+n+od, n+d1);
      
      //d od/d d2
      H(bb+n+od, n+d2) = ssigma(ind1(d1,0), ind1(d1,0)) *(-0.5)*ssigma15(ind1(d2,0), ind1(d2,0));
      H(bb+n+d2, n+od) = H(bb+n+od, n+d2);
      
      
      //d d2/d d1
      
      H(bb+n+d1, n+d2) = sigma(ind1(od,0), ind1(od,1))*(-0.5)*ssigma15(ind1(d2,0), ind1(d2,0))*(-0.5)*ssigma15(ind1(d1,0), ind1(d1,0));
      H(bb+n+d2, n+d1) = H(bb+n+d1, n+d2);
      
      
      //d d1/d d1
      H(bb+n+d1, n+d1) = sigma(ind1(od,0), ind1(od,1))*ssigma(ind1(d2,0), ind1(d2,0))*(0.75)*ssigma25(ind1(d1,0), ind1(d1,0));
      
      //d d2/d d2
      H(bb+n+d2, n+d2) = sigma(ind1(od,0), ind1(od,1))*ssigma(ind1(d1,0), ind1(d1,0))*(0.75)*ssigma25(ind1(d2,0), ind1(d2,0));
      
      
      
    }
    bb += h2;
  }
  
  
  
  return H;
  
}


//' returns the dvech(LL')/dvech(L)
//' @description
//' Returns the Jacobian of the matrix-vectorization with respect to the vech of the Cholesky factor  
//' @param x
//' real vector of vectorization of lower triangular of Cholesky factor 
//' @return 
//' Jacobian matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd J_chol(Eigen::VectorXd x){
  
  //function that only takes a vector of Cholesky factors (lower triag - column major!) 
  //and returns the dvech(LL')/dvech(L) Jacobian and the correspoding Hessian
  
  
  int p = x.rows();
  int m =  -0.5 + sqrt(0.25 + 2*p); //compute the dimension of the matrix from the length of the vector
  
  Eigen::MatrixXd C(m,m);
  C.setZero();
  
  int ii = 0;
  int jj = 0;
  
  for(int i = 0; i < p; ++i){ 
    C(ii,jj) = x[i];
    ii++;
    if(ii == m){
      ii = jj +1;
      jj = ii;
    }
  }
  
  Eigen::MatrixXd L = elimmat(m);
  Eigen::MatrixXd K = commmat(m,m);
  Eigen::MatrixXd D = duplmat(m);
  
  Eigen::MatrixXd Lt = L.transpose();
  
  Eigen::MatrixXd I(m,m);
  Eigen::MatrixXd II(m*m,m*m);
  I.setIdentity();
  II.setIdentity();
  
  Eigen::MatrixXd  J = L*(II + K)*kroneckerProduct(C,I)*Lt; //colum-major dvech(LL')/dvech(L) Jacobian
  
  
  return J;
}



Eigen::MatrixXd J_chol1(Eigen::VectorXd x, Eigen::MatrixXd L, Eigen::MatrixXd K, Eigen::MatrixXd D){
  
  //function that only takes a vector of Cholesky factors (lower triag - column major!) 
  //and returns the dvech(LL')/dvech(L) Jacobian and the correspoding Hessian
  //
  // same as J_chol, but uses precomputed values. 
  
  
  int p = x.rows();
  int m =  -0.5 + sqrt(0.25 + 2*p); //compute the dimension of the matrix from the length of the vector
  
  Eigen::MatrixXd C(m,m);
  C.setZero();
  
  int ii = 0;
  int jj = 0;
  
  for(int i = 0; i < p; ++i){ 
    C(ii,jj) = x[i];
    ii++;
    if(ii == m){
      ii = jj +1;
      jj = ii;
    }
  }
  
  Eigen::MatrixXd Lt = L.transpose();
  
  Eigen::MatrixXd I(m,m);
  Eigen::MatrixXd II(m*m,m*m);
  I.setIdentity();
  II.setIdentity();
  
  Eigen::MatrixXd  J = L*(II + K)*kroneckerProduct(C,I)*Lt; //colum-major dvech(LL')/dvech(L) Jacobian
  
  
  return J;
}




Eigen::MatrixXd J_bvOvSig_bLOLSig(Eigen::VectorXd tho, Eigen::VectorXd thl, int M, Eigen::MatrixXd Lo, Eigen::MatrixXd Ko, Eigen::MatrixXd Do, Eigen::MatrixXd Ll, Eigen::MatrixXd Kl, Eigen::MatrixXd Dl){
  
  //function that only takes a two vector of Cholesky factors (lower triag - column major!) for O and L
  // and the number of linear parameters M
  //
  // partial (b,vech(Omega),vech(Sigma)) / partial (b, vech(L_Omega),vech(L_Sigma))
  
  int O = tho.rows();
  int L = thl.rows();
  
  
  
  
  Eigen::MatrixXd  J1 = J_chol1(tho, Lo, Ko, Do); //colum-major dvech(LL')/dvech(L) Jacobian
  Eigen::MatrixXd  J2 = J_chol1(thl, Ll, Kl, Dl); //colum-major dvech(LL')/dvech(L) Jacobian
  
  Eigen::MatrixXd JJ(M+O+L,M+O+L), Ind(M,M);
  JJ.setZero();
  Ind.setIdentity();
  
  JJ.block(0,0,M,M) = Ind;
  JJ.block(M,M,O,O) = J1;
  JJ.block(M+O,M+O,L,L) = J2;
  
  return JJ;
}

Eigen::MatrixXd J_bvSig_bLSig(Eigen::VectorXd thl, int M, Eigen::MatrixXd Ll, Eigen::MatrixXd Kl, Eigen::MatrixXd Dl){
  
  //function that only takes a two vector of Cholesky factors (lower triag - column major!) for O and L
  // and the number of linear parameters M
  //
  // partial (b,vech(Sigma)) / partial (b,vech(L_Sigma))
  
  int L = thl.rows();
  
  Eigen::MatrixXd  J2 = J_chol1(thl, Ll, Kl, Dl); //colum-major dvech(LL')/dvech(L) Jacobian
  
  Eigen::MatrixXd JJ(M+L,M+L), Ind(M,M);
  JJ.setZero();
  Ind.setIdentity();
  
  JJ.block(0,0,M,M) = Ind;
  JJ.block(M,M,L,L) = J2;
  
  return JJ;
}


Eigen::MatrixXd H_chol(Eigen::MatrixXd J){
  
  //function that only takes the Jacobi matrix and computes the Hessian
  //Compute the Hessian as a stack version 
  // It is always ones for each entry but 2 for the diagonal entries => Check for the positive entries in the Jacobian and for diagonals
  int h1 = J.rows();
  int h2 = J.cols();
  
  Eigen::MatrixXd H(h1*h2,h2);
  H.setZero();
  Eigen::MatrixXd JJ(h1,h2);
  JJ.setZero();
  JJ = (J.array() != 0).select(1, JJ);
  
  int m =  -0.5 + sqrt(0.25 + 2*h1);
  
  int bb = 0;
  
  int k = 0;
  int jj = 0;
  
  for(int i = 0; i < h1; ++i){ //loop over rows
    int e1 = -1;
    
    //Identify the diag entries in the final matrix because those are squares whenever J is positive
    if( k ==i){
      for(int j = 0; j < h2; ++j){ //loop over rows
        if(JJ(i,j)!=0){
          H(bb+j,j) = 2;
        }
      }
      k = k +(m-jj);
      jj++;
    } else {
      for(int j = 0; j < h2; ++j){ //loop over rows
        
        if(JJ(i,j)!=0){
          if(e1==-1){
            e1 = j;
          } else {
            H(bb+j,(e1)) = 1;
            H(bb+(e1),(j)) = 1; 
            e1 = -1;      
          }
        }
      }
    }
    bb += h2;
  }
  
  
  return H;
}




Eigen::MatrixXd H_bvOvSig_bLOLSig(Eigen::MatrixXd J, int M, int O, int L){
  
  
  //function that only takes the Jacobian matrix and computes the Hessian
  //Compute the Hessian as a stack version 
  // It is always ones for each entry but 2 for the diagonal entries => Check for the positive entries in the Jacobian and for diagonals
  int h1 = J.rows();
  int h2 = J.cols();
  
  Eigen::MatrixXd H(h1*h2,h2);
  H.setZero();
  
  
  Eigen::MatrixXd JJ(h1,h2);
  JJ.setZero();
  JJ = (J.array() != 0).select(1, JJ);
  
  int LL =  -0.5 + sqrt(0.25 + 2*L);
  int OO =  -0.5 + sqrt(0.25 + 2*O);
  
  int bb = 0; //number oof parameters
  
  int kl = (M+O);
  int ko = M;
  
  int jjl = 0;
  int jjo = 0;
  
  
  for(int i = 0; i < h1; ++i){ //loop over rows
    
    if((i > (M-1))& (i < (M+O) )){ //Skip all the linear parameters because their Hessian is 0
      //Move through the O paras
      
      int e1 = -1;
      //Identify the diag entries in the final matrix because those are squares whenever J is positive
      if( ko ==i){
        for(int j = M; j < h2; ++j){ //loop over cols and start at the relevant section of parameters
          if(JJ(i,j)!=0){
            H(bb+j,j) = 2;
          }
        }
        ko = ko +(OO-jjo);
        jjo++;
      } else {
        for(int j = (M); j < h2; ++j){ //loop over cols and start at the relevant section of parameters
          
          if(JJ(i,j)!=0){
            if(e1==-1){
              e1 = j;
            } else {
              H(bb+j,(e1)) = 1;
              H(bb+(e1),(j)) = 1; 
              e1 = -1;      
            }
          }
          
        }
      }
    }
    
    
    
    if(i > (M+O-1)){
      
      //Move through the L paras
      
      int e1 = -1;
      //Identify the diag entries in the final matrix because those are squares whenever J is positive
      if( kl ==i){
        for(int j = (M+O); j < h2; ++j){ //loop over cols
          if(JJ(i,j)!=0){
            H(bb+j,j) = 2;
          }
        }
        kl = kl +(LL-jjl);
        jjl++;
      } else {
        for(int j = (M+O); j < h2; ++j){ //loop over cols
          
          if(JJ(i,j)!=0){
            if(e1==-1){
              e1 = j;
            } else {
              H(bb+j,(e1)) = 1;
              H(bb+(e1),(j)) = 1; 
              e1 = -1;
            }
          }
        }
      }
    }
    
    bb += h2;
    
  }
  
  
  
  return H;
}



Eigen::MatrixXd J_Omega(Eigen::MatrixXd X){
  //Computes the Jacobian for vech(Omega)
  //X: Coefficient matrix or (better) submatrix relevant to RE
  
  //Basically just Theorem 2, page 35 in Magnus/Neudecker + Elim-Mat for vec() to vech()
  // partial vech( X Omega X') / partial vech(Omega)
  int p = X.rows();
  int n = X.cols();
  
  Eigen::MatrixXd L = elimmat(p);
  Eigen::MatrixXd D = duplmat(n);
  
  Eigen::MatrixXd J = L*kroneckerProduct(X,X)*D;
  
  return J;
}



Eigen::MatrixXd J_Omega1(Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D){
  //Precomputed dupl and elm matrices for better performance
  
  //Computes the JAcobian for vech(Omega)
  //X: Coefficient matrix or (better) submatrix relevant to RE
  
  //Basically just Theorem 2, page 35 in Magnus/Neudecker + Elim-Mat for vec() to vech()
  // same as J_Omega with precomputed matrices.
  Eigen::MatrixXd J = L*kroneckerProduct(X,X)*D;
  
  return J;
  
}



Eigen::MatrixXd makeM(int m, int y){
  //m: size of the unrestricted/undiffed Matrix
  //y: row to remove, col to set to -1 (usually this will be choice - 1 due to cpp indexing)
  Eigen::MatrixXd M(m,m);
  M.setIdentity();
  M.col(y).setConstant(-1);
  removeRow(M, y);
  return M;
}



Eigen::MatrixXd J1(int mm1, int mm2, int k1, int k2, int M, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D){
  // partial (b, vech( X Omega X'), vech(M) / partial (b, vech(Omega), vech(Sigma))
  //Precomputed dupl and elm matrices for better performance
  
  //Computes the first derivative (Jacobian) for the first step
  //b is ident
  //O is D(XOOX)/D(OO)
  //L is D(block(diff(LL,k1), diff(LL, k2)))/D(LL)
  //AT THE MOMENT ONLY SITUATIONS WHERE THE PERSON HAS THE SAME NUMBER OF CHOICE ALTERNATIVES
  
  //INPUTS:
  
  // mm1, mm2 is dimnesion of LL (so, equal to the number of choice alternatives at the time point)
  //k1, k2 choices at the timepoints
  //M number of linear parameters
  //X covariate (sub-)matrix of the REs
  //L elimination matrix for J_Omega1
  //D duplication matrix for J_Omega1
  
  Eigen::MatrixXd J = J_2dM(mm1, mm2, k1, k2);
  Eigen::MatrixXd JO = J_Omega1(X, L, D);
  Eigen::MatrixXd JJ(J.rows() + M, J.cols()+ JO.cols()+M), ind(M,M);
  JJ.setZero();
  ind.setIdentity();
  
  JJ.block(0,0,M,M) = ind;
  
  JJ.block(M,M,JO.rows(),JO.cols()) = JO;
  JJ.block(M,M+JO.cols(),J.rows(),J.cols()) = J;
  
  
  return JJ;
  
}



Eigen::MatrixXd J_TpdM(int alt, int Tp, Eigen::VectorXd y){
  //Computes the first derivative (Jacobian) of  M_Sig w.r.t. Sigma,
  // M_sig = diag( M_k1 Sigma M_k1', M_k2 Sigma M_k2',...). 
  //
  // where Sig is a LL matrix which appears as a blockdiagonal with different differencing
  // mm1, mm2 is dimnesion of LL (so, equal to the number of choice alternatives at the time point), where M is of dimension (mm1-1) + (mm2- 1)
  int dim_w = Tp*(alt-1);
  
  Eigen::MatrixXd J(dim_w*(dim_w+1)/2, (alt+1)*alt/2);
  J.setZero();
  
  int cur = 0;
  //int cur_col = 0;
  
  for (int jj=0; jj< Tp;jj++){
    Eigen::MatrixXd J1 = J_diff(alt, y[jj]); // choice.
    int curJ = 0;
    for (int a=0;a<alt-1;a++){ // a is column index
      for (int b=0;b<alt-1-a;b++){ // b is row index
        if (cur< J.rows()){
          J.row(cur) = J1.row(curJ);
          curJ++;
          cur++;
        }
      } // end loop b
      cur += dim_w-(jj+1)*(alt-1);
    } // end loop a
  }
  
  return J;
}



Eigen::MatrixXd J_TpdM_red(int alt, Eigen::VectorXd y){
  //Computes the first derivative (Jacobian) of  M_Sig w.r.t. Sigma,
  // M_sig = diag( M_k1 Sigma M_k1', M_k2 Sigma M_k2',...).
  //
  // where Sig is a LL matrix which appears as a blockdiagonal with different differencing
  // alt is dimnesion of LL (so, equal to the number of choice alternatives at the time point), 
  int Tp = y.rows();
  int dim_w = Tp*(alt-1);
  
  Eigen::MatrixXd J(dim_w*(dim_w+1)/2, (alt+1)*alt/2);
  J.setZero();
  
  int cur = 0;
  //int cur_col = 0;
  
  for (int jj=0; jj< Tp;jj++){
    Eigen::MatrixXd J1 = J_diff(alt, y[jj]); // choice.
    int curJ = 0;
    for (int a=0;a<alt-1;a++){ // a is column index
      for (int b=0;b<alt-1-a;b++){ // b is row index
        if (cur< J.rows()){
          J.row(cur) = J1.row(curJ);
          curJ++;
          cur++;
        }
      } // end loop b
      cur += dim_w-(jj+1)*(alt-1);
    } // end loop a
  }
  
  return J;
}



Eigen::MatrixXd J_M_bOmSig(int alt, int Tp, Eigen::VectorXd y, int lthbb, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D){
  // partial (b, vech(M) / partial (b, vech(Omega), vech(Sigma))
  //Precomputed dupl and elm matrices for better performance
  
  //Computes the first derivative (Jacobian) for the first step
  //b is ident
  //O is D(XOOX)/D(OO)
  //L is D(block(diff(LL,k1), diff(LL, k2)))/D(LL)
  //AT THE MOMENT ONLY SITUATIONS WHERE THE PERSON HAS THE SAME NUMBER OF CHOICE ALTERNATIVES
  
  //INPUTS:
  
  // mm1, mm2 is dimnesion of LL (so, equal to the number of choice alternatives at the time point)
  //k1, k2 choices at the timepoints
  //M number of linear parameters
  //X covariate (sub-)matrix of the REs
  //L elimination matrix for J_Omega1
  //D duplication matrix for J_Omega1
  
  Eigen::MatrixXd J = J_TpdM(alt, Tp, y);
  Eigen::MatrixXd JO = J_Omega1(X, L, D);
  Eigen::MatrixXd JJ(J.rows() + lthbb, J.cols()+ JO.cols()+lthbb), ind(lthbb,lthbb);
  JJ.setZero();
  ind.setIdentity();
  
  JJ.block(0,0,lthbb,lthbb) = ind;
  
  JJ.block(lthbb,lthbb,JO.rows(),JO.cols()) = JO;
  JJ.block(lthbb,lthbb+JO.cols(),J.rows(),J.cols()) = J;
  
  
  return JJ;
  
}



Eigen::MatrixXd J_M_bOmSig_red(int alt, IntegerVector ind_BI, Eigen::VectorXd y, int lthbb, Eigen::MatrixXd X, Eigen::MatrixXd L, Eigen::MatrixXd D){
  // partial (b, vech(M) / partial (b, vech(Omega), vech(Sigma))
  //Precomputed dupl and elm matrices for better performance
  
  //Computes the first derivative (Jacobian) for the first step
  //b is ident
  //O is D(XOOX)/D(OO)
  //L is D(block(diff(LL,k1), diff(LL, k2)))/D(LL)
  //AT THE MOMENT ONLY SITUATIONS WHERE THE PERSON HAS THE SAME NUMBER OF CHOICE ALTERNATIVES
  
  //INPUTS:
  
  // mm1, mm2 is dimnesion of LL (so, equal to the number of choice alternatives at the time point)
  //k1, k2 choices at the timepoints
  //M number of linear parameters
  //X covariate (sub-)matrix of the REs
  //L elimination matrix for J_Omega1
  //D duplication matrix for J_Omega1
  
  int cur = X.rows();
  
  
  Eigen::MatrixXd J = J_TpdM_red(alt,  y);
  int Tp = y.rows();
  int rows_full = Tp*(alt-1);
  VectorXi ind( rows_full);
  ind.setZero();
  
  for (int jj=0;jj<cur;jj++){
    ind[ind_BI[jj]]=1;
  }
  // reduce due to missing values
  int cur_red = 0;
  int cur_full = 0;
  
  for (int a=0;a<rows_full;a++){
    if (ind[a]==1){ // column is not reduced
      for (int b=a;b<rows_full;b++){
        if(ind[b]==1){ // row is not reduced
          J.row(cur_red) = J.row(cur_full);
          cur_red++;
        }
        cur_full++;
      } 
    }
    if (ind[a]==0) { // column is reduced -> increase counter accordingly.
      cur_full = cur_full + rows_full - (a);
    }
  }
  
  
  Eigen::MatrixXd J_red = J.block(0,0,cur_red,J.cols());
  // derviative of vech(X*Omega*X') with respect to entries in lower diagonal of Omega (including diagonal) 
  Eigen::MatrixXd JO = J_Omega1(X, L, D);
  
  
  // initialize output matrix.
  Eigen::MatrixXd JJ(J_red.rows() + lthbb, J_red.cols()+ JO.cols()+lthbb), BdB(lthbb,lthbb);
  JJ.setZero();
  BdB.setIdentity();
  
  JJ.block(0,0,lthbb,lthbb) = BdB;
  
  JJ.block(lthbb,lthbb,JO.rows(),JO.cols()) = JO;
  JJ.block(lthbb,lthbb+JO.cols(),J_red.rows(),J_red.cols()) = J_red;
  
  return JJ;
}

Eigen::MatrixXd J_bvOv(Eigen::VectorXd tho, Eigen::VectorXd thl, int M, Eigen::MatrixXd Lo, Eigen::MatrixXd Ko, Eigen::MatrixXd Do){
  
  //function that only takes a vector of Cholesky factors (lower triag - column major!) for O 
  // and the number of linear parameters M
  //
  // partial (b,vech(Omega)) / partial (b, vech(L_Omega))
  
  int O = tho.rows();
  
  Eigen::MatrixXd  J1 = J_chol1(tho, Lo, Ko, Do); //colum-major dvech(LL')/dvech(L) Jacobian
  
  Eigen::MatrixXd JJ(M+O,M+O), Ind(M,M);
  JJ.setZero();
  Ind.setIdentity();
  
  JJ.block(0,0,M,M) = Ind;
  JJ.block(M,M,O,O) = J1;
  
  return JJ;
}

//' calculates the probability of pair of choices in an ordered probit model
//' @description
//' calculates the probability of pair of choices in an ordered probit model
//' @param xb
//' vector of systematic utilities
//' @param y
//' vector of choices
//' @param Lambda 
//' correlation matrix
//' @param alt
//' integer; number of alternatives
//' @param dtauk
//' vector of category limits. 
//' @return 
//' matrix
//' @keywords internal
//'
//[[Rcpp::export]]
double prob_ordered_2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk)
{
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  

  double y1 = y(0);
  double y2 = y(1);
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  if (y2<1){ y2=1;}
  if (y2>alt){ y2=alt;}
  
  // upper and lower bounds 
  double lower_bound1,lower_bound2;
  lower_bound1 = tauk(0)-10.0;
  lower_bound2 = tauk(0)-10.0;
  
  double upper_bound1, upper_bound2;
  upper_bound1= tauk(alt-2)+10.0;
  upper_bound2= tauk(alt-2)+10.0;
  
  // first values 
  if (y1==1){
    upper_bound1 =tauk(0);
  }
  if (y1==alt){
    lower_bound1 = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound1 = tauk(y1-2);
      upper_bound1 = tauk(y1-1);
    }
  }
  // second value
  if (y2==1){
    upper_bound2 =tauk(0);
  }
  if (y2==alt){
    lower_bound2 = tauk(alt-2);
  }
  
  if (y2>1){
    if(y2<alt){
      lower_bound2 = tauk(y2-2);
      upper_bound2 = tauk(y2-1);
    }
  }
  
  // normalized bounds.
  double ub1 = upper_bound1/lambda(0,0)-x_norm(0);
  double lb1 = lower_bound1/lambda(0,0)-x_norm(0);
  double ub2 = upper_bound2/lambda(1,0)-x_norm(1);
  double lb2 = lower_bound2/lambda(1,0)-x_norm(1);
  
  // probabilities 
  double rho = Lambda_cor(0,1);
  
  double pru1u2 = biv_normal_cdf(ub1,ub2,rho);
  double pru1l2 = biv_normal_cdf(ub1,lb2,rho);
  double prl1u2 = biv_normal_cdf(lb1,ub2,rho);
  double prl1l2 = biv_normal_cdf(lb1,lb2,rho);
  
  pr = pru1u2 - pru1l2 - prl1u2 + prl1l2;
  if (pr < 0.0000000001){ pr = 0.0000000001; } // regularize, if percentage is too small. 
  double lpr = log(pr);
  
  // return value.
  return lpr;
}

Eigen::VectorXd  grad_po2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk)
{
  Eigen::VectorXd grad(5+dtauk.size());
  grad.setZero(); 
  
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  

  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
 
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  double y1 = y(0);
  double y2 = y(1);
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  if (y2<1){ y2=1;}
  if (y2>alt){ y2=alt;}
  
  // upper and lower bounds 
  double lower_bound1,lower_bound2;
  lower_bound1 = tauk(0)-10.0;
  lower_bound2 = tauk(0)-10.0;
  
  double upper_bound1, upper_bound2;
  upper_bound1= tauk(alt-2)+10.0;
  upper_bound2= tauk(alt-2)+10.0;
  
  // first values 
  if (y1==1){
    upper_bound1 =tauk(0);
  }
  if (y1==alt){
    lower_bound1 = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound1 = tauk(y1-2);
      upper_bound1 = tauk(y1-1);
    }
  }
  // second value
  if (y2==1){
    upper_bound2 =tauk(0);
  }
  if (y2==alt){
    lower_bound2 = tauk(alt-2);
  }
  
  if (y2>1){
    if(y2<alt){
      lower_bound2 = tauk(y2-2);
      upper_bound2 = tauk(y2-1);
    }
  }
  
  // normalized bounds.
  double ub1 = upper_bound1/lambda(0,0)-x_norm(0);
  double lb1 = lower_bound1/lambda(0,0)-x_norm(0);
  double ub2 = upper_bound2/lambda(1,0)-x_norm(1);
  double lb2 = lower_bound2/lambda(1,0)-x_norm(1);
  
  // probabilities 
  double rho = Lambda_cor(0,1);
  
  double pru1u2 = biv_normal_cdf(ub1,ub2,rho);
  double pru1l2 = biv_normal_cdf(ub1,lb2,rho);
  double prl1u2 = biv_normal_cdf(lb1,ub2,rho);
  double prl1l2 = biv_normal_cdf(lb1,lb2,rho);
  
  pr = pru1u2 - pru1l2 - prl1u2 + prl1l2;
  
  if (pr < 0.0000000001){ pr = 0.0000000001; } // regularize, if percentage is too small. 
  
  
  // derivative of pru1u2
  Eigen::VectorXd grad_u1u2(3);
  grad_u1u2 = grad_cdf(ub1,ub2,rho);
  
  grad(0) += grad_u1u2(0)*(-1/lambda(0,0));
  grad(1) += grad_u1u2(1)*(-1/lambda(1,0));
  grad(2) += grad_u1u2(0)*(-ub1/(2*Lambda(0,0))) - grad_u1u2(2)*(rho/(2*Lambda(0,0)));
  grad(3) += grad_u1u2(1)*(-ub2/(2*Lambda(1,1))) - grad_u1u2(2)*(rho/(2*Lambda(1,1)));
  grad(4) += grad_u1u2(2)/(lambda(0,0)*lambda(1,0));
  
  Eigen::VectorXd iota(alt-1);
  //iota.setOnes();
  iota = edtauk; 
  iota(0)=1.0; 
  
  if (y1<alt){
    grad.segment(5,y1) +=  grad_u1u2(0)*iota.head(y1)/lambda(0,0);   
  }
  if (y2<alt){
    grad.segment(5,y2) +=  grad_u1u2(1)*iota.head(y2)/lambda(1,0);   
  }
  
  // derivative of pru1l2
  Eigen::VectorXd grad_u1l2(3);
  grad_u1l2 = grad_cdf(ub1,lb2,rho);
  
  grad(0) -= grad_u1l2(0)*(-1/lambda(0,0));
  grad(1) -= grad_u1l2(1)*(-1/lambda(1,0));
  grad(2) -= grad_u1l2(0)*(-ub1/(2*Lambda(0,0))) - grad_u1l2(2)*(rho/(2*Lambda(0,0)));
  grad(3) -= grad_u1l2(1)*(-lb2/(2*Lambda(1,1))) - grad_u1l2(2)*(rho/(2*Lambda(1,1)));
  grad(4) -= grad_u1l2(2)/(lambda(0,0)*lambda(1,0));
  
  if (y1<alt){
    grad.segment(5,y1) -=  grad_u1l2(0)*iota.head(y1)/lambda(0,0);   
  }
  if (y2>1){
    grad.segment(5,y2-1) -=  grad_u1l2(1)*iota.head(y2-1)/lambda(1,0);   
  }
  
  // derivative of prl1u2
  Eigen::VectorXd grad_l1u2(3);
  grad_l1u2 = grad_cdf(lb1,ub2,rho);
  
  grad(0) -= grad_l1u2(0)*(-1/lambda(0,0));
  grad(1) -= grad_l1u2(1)*(-1/lambda(1,0));
  grad(2) -= grad_l1u2(0)*(-lb1/(2*Lambda(0,0))) - grad_l1u2(2)*(rho/(2*Lambda(0,0)));
  grad(3) -= grad_l1u2(1)*(-ub2/(2*Lambda(1,1))) - grad_l1u2(2)*(rho/(2*Lambda(1,1)));
  grad(4) -= grad_l1u2(2)/(lambda(0,0)*lambda(1,0));
  
  if (y1>1){
    grad.segment(5,y1-1) -=  grad_l1u2(0)*iota.head(y1-1)/lambda(0,0);   
  }
  if (y2<alt){
    grad.segment(5,y2) -=  grad_l1u2(1)*iota.head(y2)/lambda(1,0);   
  }
  
  // derivative of prl1l2
  Eigen::VectorXd grad_l1l2(3);
  grad_l1l2 = grad_cdf(lb1,lb2,rho);
  
  grad(0) += grad_l1l2(0)*(-1/lambda(0,0));
  grad(1) += grad_l1l2(1)*(-1/lambda(1,0));
  grad(2) += grad_l1l2(0)*(-lb1/(2*Lambda(0,0))) - grad_l1l2(2)*(rho/(2*Lambda(0,0)));
  grad(3) += grad_l1l2(1)*(-lb2/(2*Lambda(1,1))) - grad_l1l2(2)*(rho/(2*Lambda(1,1)));
  grad(4) += grad_l1l2(2)/(lambda(0,0)*lambda(1,0));
  
  if (y1>1){
    grad.segment(5,y1-1) +=  grad_l1l2(0)*iota.head(y1-1)/lambda(0,0);   
  }
  if (y2>1){
    grad.segment(5,y2-1) +=  grad_l1l2(1)*iota.head(y2-1)/lambda(1,0);   
  }
  
  // take into account logarithm. 
  grad = grad / pr; 
  // return value.
  return grad;
}

Eigen::MatrixXd  hess_po2(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk)
{
  Eigen::MatrixXd Hess(5+dtauk.size(),5+dtauk.size()); 
  Hess.setZero(); 
  
  Eigen::VectorXd grad(5+dtauk.size());
  grad.setZero(); 
  
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }

  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  

  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  double y1 = y(0);
  double y2 = y(1);
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  if (y2<1){ y2=1;}
  if (y2>alt){ y2=alt;}
  
  // upper and lower bounds 
  double lower_bound1,lower_bound2;
  lower_bound1 = tauk(0)-10.0;
  lower_bound2 = tauk(0)-10.0;
  
  double upper_bound1, upper_bound2;
  upper_bound1= tauk(alt-2)+10.0;
  upper_bound2= tauk(alt-2)+10.0;
  
  // first values 
  if (y1==1){
    upper_bound1 =tauk(0);
  }
  if (y1==alt){
    lower_bound1 = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound1 = tauk(y1-2);
      upper_bound1 = tauk(y1-1);
    }
  }
  // second value
  if (y2==1){
    upper_bound2 =tauk(0);
  }
  if (y2==alt){
    lower_bound2 = tauk(alt-2);
  }
  
  if (y2>1){
    if(y2<alt){
      lower_bound2 = tauk(y2-2);
      upper_bound2 = tauk(y2-1);
    }
  }
  
  // normalized bounds.
  double ub1 = upper_bound1/lambda(0,0)-x_norm(0);
  double lb1 = lower_bound1/lambda(0,0)-x_norm(0);
  double ub2 = upper_bound2/lambda(1,0)-x_norm(1);
  double lb2 = lower_bound2/lambda(1,0)-x_norm(1);
  
  // probabilities 
  double rho = Lambda_cor(0,1);
  
  double pru1u2 = biv_normal_cdf(ub1,ub2,rho);
  double pru1l2 = biv_normal_cdf(ub1,lb2,rho);
  double prl1u2 = biv_normal_cdf(lb1,ub2,rho);
  double prl1l2 = biv_normal_cdf(lb1,lb2,rho);
  
  pr = pru1u2 - pru1l2 - prl1u2 + prl1l2;
  
  if (pr < 0.0000000001){ pr = 0.0000000001; } // regularize, if percentage is too small. 

  // derivative of pru1u2
  Eigen::VectorXd grad_u1u2(3);
  grad_u1u2 = grad_cdf(ub1,ub2,rho);
  
  Eigen::MatrixXd Ju1u2(3,grad.size());
  Ju1u2.setZero(); 
  
  Ju1u2(0,0) = -1/lambda(0,0);
  Ju1u2(1,1) = -1/lambda(1,0);
  Ju1u2(0,2) = -ub1/(2*Lambda(0,0));
  Ju1u2(2,2) = -rho/(2*Lambda(0,0));
  Ju1u2(1,3) = -ub2/(2*Lambda(1,1));
  Ju1u2(2,3) = -rho/(2*Lambda(1,1));
  Ju1u2(2,4) = 1/(lambda(0,0)*lambda(1,0));
  
  Eigen::VectorXd iota(alt-1);
  iota = edtauk;
  iota(0)= 1.0; 
  
  if (y1<alt){
    Ju1u2.block(0,5,1,y1) =   iota.head(y1).transpose()/lambda(0,0);
  }
  if (y2<alt){
    Ju1u2.block(1,5,1,y2) =   iota.head(y2).transpose()/lambda(1,0);
  }
  
  grad = grad_u1u2.transpose() * Ju1u2;
  
  // derivative of pru1l2
  Eigen::VectorXd grad_u1l2(3);
  grad_u1l2 = grad_cdf(ub1,lb2,rho);
  
  Eigen::MatrixXd Ju1l2(3,grad.size());
  Ju1l2.setZero(); 
  
  Ju1l2(0,0) = -1/lambda(0,0);
  Ju1l2(1,1) = -1/lambda(1,0);
  Ju1l2(0,2) = -ub1/(2*Lambda(0,0));
  Ju1l2(2,2) = -rho/(2*Lambda(0,0));
  Ju1l2(1,3) = -lb2/(2*Lambda(1,1));
  Ju1l2(2,3) = -rho/(2*Lambda(1,1));
  Ju1l2(2,4) = 1/(lambda(0,0)*lambda(1,0));
  
  if (y1<alt){
    Ju1l2.block(0,5,1,y1) =   iota.head(y1).transpose()/lambda(0,0);
  }
  if (y2>1){
    Ju1l2.block(1,5,1,y2-1) =   iota.head(y2-1).transpose()/lambda(1,0);
  }
  
  grad -= grad_u1l2.transpose() * Ju1l2;
  
  // derivative of prl1u2
  Eigen::VectorXd grad_l1u2(3);
  grad_l1u2 = grad_cdf(lb1,ub2,rho);
  
  Eigen::MatrixXd Jl1u2(3,grad.size());
  Jl1u2.setZero(); 
  
  Jl1u2(0,0) = -1/lambda(0,0);
  Jl1u2(1,1) = -1/lambda(1,0);
  Jl1u2(0,2) = -lb1/(2*Lambda(0,0));
  Jl1u2(2,2) = -rho/(2*Lambda(0,0));
  Jl1u2(1,3) = -ub2/(2*Lambda(1,1));
  Jl1u2(2,3) = -rho/(2*Lambda(1,1));
  Jl1u2(2,4) = 1/(lambda(0,0)*lambda(1,0));
  
  if (y2<alt){
    Jl1u2.block(1,5,1,y2) =   iota.head(y2).transpose()/lambda(1,0);
  }
  if (y1>1){
    Jl1u2.block(0,5,1,y1-1) =   iota.head(y1-1).transpose()/lambda(0,0);
  }
  
  grad -= grad_l1u2.transpose() * Jl1u2;
  
  
  // derivative of prl1l2
  Eigen::VectorXd grad_l1l2(3);
  grad_l1l2 = grad_cdf(lb1,lb2,rho);
  
  Eigen::MatrixXd Jl1l2(3,grad.size());
  Jl1l2.setZero(); 
  
  Jl1l2(0,0) = -1/lambda(0,0);
  Jl1l2(1,1) = -1/lambda(1,0);
  Jl1l2(0,2) = -lb1/(2*Lambda(0,0));
  Jl1l2(2,2) = -rho/(2*Lambda(0,0));
  Jl1l2(1,3) = -lb2/(2*Lambda(1,1));
  Jl1l2(2,3) = -rho/(2*Lambda(1,1));
  Jl1l2(2,4) = 1/(lambda(0,0)*lambda(1,0));
  
  if (y2>1){
    Jl1l2.block(1,5,1,y2-1) =   iota.head(y2-1).transpose()/lambda(1,0);
  }
  if (y1>1){
    Jl1l2.block(0,5,1,y1-1) =   iota.head(y1-1).transpose()/lambda(0,0);
  }
  
  grad += grad_l1l2.transpose() * Jl1l2;
  
  
  
  /////////////////////////////
  //  Hessian calculation 
  /////////////////////////////
  
  /////////////////////////////////////////////////////
  // u1u2 derivative of grad.
  ////////////////////////////////////////////////////
  
  Eigen::MatrixXd Hess_u1u2(5,5);
  Hess_u1u2 = Hess_cdf(ub1,ub2,rho);
  
  Hess = Ju1u2.transpose() * Hess_u1u2 * Ju1u2;
  
  // deriv w.r.t. x(0). 
  Hess(0,2) += grad_u1u2(0)/(2*Lambda(0,0)*lambda(0,0));
  // deriv w.r.t. x(1). 
  Hess(1,3) += grad_u1u2(1)/(2*Lambda(1,1)*lambda(1,0));
  // deriv w.r.t. Lambda(0,0).
  Hess(2,0) += grad_u1u2(0)/(2*Lambda(0,0)*lambda(0,0));
  Hess(2,2) += grad_u1u2(0)*3*ub1/(4*Lambda(0,0)*Lambda(0,0))+ 3*grad_u1u2(2)*rho/(4*Lambda(0,0)*Lambda(0,0));
  Hess(2,3) += grad_u1u2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(2,4) -= grad_u1u2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(0,0));
  if (y1<alt){
    Hess.block(2,5,1,y1) -= grad_u1u2(0)*iota.head(y1).transpose()/(2*lambda(0,0)*Lambda(0,0));
  }
  // deriv w.r.t. Lambda(1,1).
  Hess(3,1) += grad_u1u2(1)/(2*Lambda(1,1)*lambda(1,0));
  Hess(3,2) += grad_u1u2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(3,3) += grad_u1u2(1)*3*ub2/(4*Lambda(1,1)*Lambda(1,1))+ 3*grad_u1u2(2)*rho/(4*Lambda(1,1)*Lambda(1,1));
  Hess(3,4) -= grad_u1u2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(1,1));
  if (y2<alt){
    Hess.block(3,5,1,y2) -= grad_u1u2(1)*iota.head(y2).transpose()/(2*lambda(1,0)*Lambda(1,1));
  }
  // deriv w.r.t. Lambda(0,1).
  Hess(4,2) -=  grad_u1u2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(0,0));
  Hess(4,3) -=  grad_u1u2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(1,1));
  // deriv w.r.t. tauk.
  if (y1<alt){
    Hess.block(5,2,y1,1) -= grad_u1u2(0)*iota.head(y1)/(2*lambda(0,0)*Lambda(0,0));
  }
  if (y2<alt){
    Hess.block(5,3,y2,1) -= grad_u1u2(1)*iota.head(y2)/(2*lambda(1,0)*Lambda(1,1));
  }
  
  
  /////////////////////////////////////////////////////
  // u1l2 derivative of grad.
  ////////////////////////////////////////////////////
  
  Eigen::MatrixXd Hess_u1l2(5,5);
  Hess_u1l2 = Hess_cdf(ub1,lb2,rho);
  

  Hess += -Ju1l2.transpose() * Hess_u1l2 * Ju1l2;
  
  // deriv w.r.t. x(0). 
  Hess(0,2) -= grad_u1l2(0)/(2*Lambda(0,0)*lambda(0,0));
  // deriv w.r.t. x(1). 
  Hess(1,3) -= grad_u1l2(1)/(2*Lambda(1,1)*lambda(1,0));
  // deriv w.r.t. Lambda(0,0).
  Hess(2,0) -= grad_u1l2(0)/(2*Lambda(0,0)*lambda(0,0));
  Hess(2,2) -= grad_u1l2(0)*3*ub1/(4*Lambda(0,0)*Lambda(0,0))+ 3*grad_u1l2(2)*rho/(4*Lambda(0,0)*Lambda(0,0));
  Hess(2,3) -= grad_u1l2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(2,4) += grad_u1l2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(0,0));
  if (y1<alt){
    Hess.block(2,5,1,y1) += grad_u1l2(0)*iota.head(y1).transpose()/(2*lambda(0,0)*Lambda(0,0));
  }
  // deriv w.r.t. Lambda(1,1).
  Hess(3,1) -= grad_u1l2(1)/(2*Lambda(1,1)*lambda(1,0));
  Hess(3,2) -= grad_u1l2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(3,3) -= grad_u1l2(1)*3*lb2/(4*Lambda(1,1)*Lambda(1,1))+ 3*grad_u1l2(2)*rho/(4*Lambda(1,1)*Lambda(1,1));
  Hess(3,4) += grad_u1l2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(1,1));
  if (y2>1){
    Hess.block(3,5,1,y2-1) += grad_u1l2(1)*iota.head(y2-1).transpose()/(2*lambda(1,0)*Lambda(1,1));
  }
  // deriv w.r.t. Lambda(0,1).
  Hess(4,2) +=  grad_u1l2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(0,0));
  Hess(4,3) +=  grad_u1l2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(1,1));
  // deriv w.r.t. tauk.
  if (y1<alt){
    Hess.block(5,2,y1,1) += grad_u1l2(0)*iota.head(y1)/(2*lambda(0,0)*Lambda(0,0));
  }
  if (y2>1){
    Hess.block(5,3,y2-1,1) += grad_u1l2(1)*iota.head(y2-1)/(2*lambda(1,0)*Lambda(1,1));
  }
  
  
  /////////////////////////////////////////////////////
  // l1u2 derivative of grad.
  ////////////////////////////////////////////////////
  
  Eigen::MatrixXd Hess_l1u2(5,5);
  Hess_l1u2 = Hess_cdf(lb1,ub2,rho);
  
  Hess += -Jl1u2.transpose() * Hess_l1u2 * Jl1u2;
  
  // deriv w.r.t. x(0). 
  Hess(0,2) -= grad_l1u2(0)/(2*Lambda(0,0)*lambda(0,0));
  // deriv w.r.t. x(1). 
  Hess(1,3) -= grad_l1u2(1)/(2*Lambda(1,1)*lambda(1,0));
  // deriv w.r.t. Lambda(0,0).
  Hess(2,0) -= grad_l1u2(0)/(2*Lambda(0,0)*lambda(0,0));
  Hess(2,2) -= grad_l1u2(0)*3*lb1/(4*Lambda(0,0)*Lambda(0,0))+ 3*grad_l1u2(2)*rho/(4*Lambda(0,0)*Lambda(0,0));
  Hess(2,3) -= grad_l1u2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(2,4) += grad_l1u2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(0,0));
  if (y1>1){
    Hess.block(2,5,1,y1-1) += grad_l1u2(0)*iota.head(y1-1).transpose()/(2*lambda(0,0)*Lambda(0,0));
  }
  // deriv w.r.t. Lambda(1,1).
  Hess(3,1) -= grad_l1u2(1)/(2*Lambda(1,1)*lambda(1,0));
  Hess(3,2) -= grad_l1u2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(3,3) -= grad_l1u2(1)*3*ub2/(4*Lambda(1,1)*Lambda(1,1))+ 3*grad_l1u2(2)*rho/(4*Lambda(1,1)*Lambda(1,1));
  Hess(3,4) += grad_l1u2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(1,1));
  if (y2<alt){
    Hess.block(3,5,1,y2) += grad_l1u2(1)*iota.head(y2).transpose()/(2*lambda(1,0)*Lambda(1,1));
  }
  // deriv w.r.t. Lambda(0,1).
  Hess(4,2) +=  grad_l1u2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(0,0));
  Hess(4,3) +=  grad_l1u2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(1,1));
  // deriv w.r.t. tauk.
  if (y1>1){
    Hess.block(5,2,y1-1,1) += grad_l1u2(0)*iota.head(y1-1)/(2*lambda(0,0)*Lambda(0,0));
  }
  if (y2<alt){
    Hess.block(5,3,y2,1) += grad_l1u2(1)*iota.head(y2)/(2*lambda(1,0)*Lambda(1,1));
  }
  
  /////////////////////////////////////////////////////
  // l1l2 derivative of grad.
  ////////////////////////////////////////////////////
  
  Eigen::MatrixXd Hess_l1l2(5,5);
  Hess_l1l2 = Hess_cdf(lb1,lb2,rho);

  Hess += Jl1l2.transpose() * Hess_l1l2 * Jl1l2;
  
  // deriv w.r.t. x(0). 
  Hess(0,2) += grad_l1l2(0)/(2*Lambda(0,0)*lambda(0,0));
  // deriv w.r.t. x(1). 
  Hess(1,3) += grad_l1l2(1)/(2*Lambda(1,1)*lambda(1,0));
  // deriv w.r.t. Lambda(0,0).
  Hess(2,0) += grad_l1l2(0)/(2*Lambda(0,0)*lambda(0,0));
  Hess(2,2) += grad_l1l2(0)*3*lb1/(4*Lambda(0,0)*Lambda(0,0))+ 3*grad_l1l2(2)*rho/(4*Lambda(0,0)*Lambda(0,0));
  Hess(2,3) += grad_l1l2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(2,4) -= grad_l1l2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(0,0));
  if (y1>1){
    Hess.block(2,5,1,y1-1) -= grad_l1l2(0)*iota.head(y1-1).transpose()/(2*lambda(0,0)*Lambda(0,0));
  }
  // deriv w.r.t. Lambda(1,1).
  Hess(3,1) += grad_l1l2(1)/(2*Lambda(1,1)*lambda(1,0));
  Hess(3,2) += grad_l1l2(2)*rho/(4*Lambda(0,0)*Lambda(1,1));
  Hess(3,3) += grad_l1l2(1)*3*lb2/(4*Lambda(1,1)*Lambda(1,1))+ 3*grad_l1l2(2)*rho/(4*Lambda(1,1)*Lambda(1,1));
  Hess(3,4) -= grad_l1l2(2)/(2*lambda(1,0)*lambda(0,0)*Lambda(1,1));
  if (y2>1){
    Hess.block(3,5,1,y2-1) -= grad_l1l2(1)*iota.head(y2-1).transpose()/(2*lambda(1,0)*Lambda(1,1));
  }
  // deriv w.r.t. Lambda(0,1).
  Hess(4,2) -=  grad_l1l2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(0,0));
  Hess(4,3) -=  grad_l1l2(2)/(2*lambda(0,0)*lambda(1,0)*Lambda(1,1));
  // deriv w.r.t. tauk.
  if (y1>1){
    Hess.block(5,2,y1-1,1) -= grad_l1l2(0)*iota.head(y1-1)/(2*lambda(0,0)*Lambda(0,0));
  }
  if (y2>1){
    Hess.block(5,3,y2-1,1) -= grad_l1l2(1)*iota.head(y2-1)/(2*lambda(1,0)*Lambda(1,1));
  }
  
  // take into account logarithm.
  Hess = Hess/pr; 
  Hess -= grad* grad.transpose()/(pr*pr); 
  
  // return value.
  return Hess;
}


//' calculates the CML of a number of choices in an ordered probit model
//' @description
//' calculates the CML of a number of choices in an ordered probit model
//' @param xb
//' vector of systematic utilities
//' @param y
//' vector of choices
//' @param Lambda 
//' correlation matrix
//' @param alt
//' integer; number of alternatives
//' @param dtauk
//' vector of category limits. 
//' @param cml_pair_type
//' integer; indicating pair type structure
//' @return 
//' matrix
//' @keywords internal
//'
//[[Rcpp::export]]
double prob_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type){
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  int Tp_n= xb.size();
  int lth_pairs_n = (Tp_n-1)*Tp_n/2;
  
  // Define a matrix that allows to specify the pairs (possible extension to higher order CMLs)
  Eigen::MatrixXd pairs_n;
  int pairs_check = 0;
  
  // if there were no valid pairs for the decision maker, construct them
  if(pairs_check==0){
    Eigen::MatrixXd pairs_temp(lth_pairs_n,3);
    pairs_temp.setOnes();
    // if the pairs are not predefined, build them
    int i_pair = 0;
    for(int i_pair_0 = 0; i_pair_0 < (Tp_n-1); i_pair_0 ++){
      for(int i_pair_1 = (i_pair_0+1);i_pair_1 < (Tp_n); i_pair_1 ++){
        // first choice occasion of the pair
        pairs_temp(i_pair,0) = i_pair_0;
        // second choice occasion of the pair
        pairs_temp(i_pair,1) = i_pair_1;
        
        // define the according weights to exclude undesired pairs
        if(cml_pair_type == 1){
          // for adjacent pairs looped:
          if( !(((i_pair_0+1) == i_pair_1) || ((i_pair_0 == 0) && (i_pair_1 == (Tp_n-1)))) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        if(cml_pair_type == 2){
          // for adjacent pairs chained
          if( !((i_pair_0+1) == i_pair_1) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        i_pair++;
      }
    }
    pairs_n = pairs_temp;
    pairs_check = 1;
  }
  
  for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){
    
    // only calculate the pair if the weight is !=0
    if(pairs_n(i_pair,2)!=0){
      // get the two coordinates 
      int i_1= pairs_n(i_pair,0);
      int i_2 = pairs_n(i_pair,1);
      
      Eigen::VectorXd yp(2);
      yp(0)= y(i_1);
      yp(1) = y(i_2);
      
      Eigen::VectorXd xbp(2);
      xbp(0)= xb(i_1);
      xbp(1)= xb(i_2);
      
      Eigen::MatrixXd Lambdap(2,2);
      Lambdap.setZero();
      Lambdap(0,0)= Lambda(i_1,i_1);
      Lambdap(0,1)= Lambda(i_1,i_2);
      Lambdap(1,0)= Lambda(i_2,i_1);
      Lambdap(1,1)= Lambda(i_2,i_2);
      
      double prp = prob_ordered_2(xbp,yp,Lambdap,alt,dtauk);

      pr += prp;
      
    }
  }
  
  // return value. 
  return pr;
}

Eigen::VectorXd grad_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type){
  
  // calculate the number of parameters. 
  int Tp_n= xb.size();
  
  int n_par = Tp_n + Tp_n*(Tp_n+1)/2 + alt-1;
  int np = n_par - alt+1;
  
  Eigen::VectorXd grad(n_par);
  grad.setZero();
  
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  int lth_pairs_n = (Tp_n-1)*Tp_n/2;
  
  // Define a matrix that allows to specify the pairs (possible extension to higher order CMLs)
  Eigen::MatrixXd pairs_n;
  int pairs_check = 0;
  
  // if there were no valid pairs for the decision maker, construct them
  if(pairs_check==0){
    Eigen::MatrixXd pairs_temp(lth_pairs_n,3);
    pairs_temp.setOnes();
    // if the pairs are not predefined, build them
    int i_pair = 0;
    for(int i_pair_0 = 0; i_pair_0 < (Tp_n-1); i_pair_0 ++){
      for(int i_pair_1 = (i_pair_0+1);i_pair_1 < (Tp_n); i_pair_1 ++){
        // first choice occasion of the pair
        pairs_temp(i_pair,0) = i_pair_0;
        // second choice occasion of the pair
        pairs_temp(i_pair,1) = i_pair_1;
        
        // define the according weights to exclude undesired pairs
        if(cml_pair_type == 1){
          // for adjacent pairs looped:
          if( !(((i_pair_0+1) == i_pair_1) || ((i_pair_0 == 0) && (i_pair_1 == (Tp_n-1)))) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        if(cml_pair_type == 2){
          // for adjacent pairs chained
          if( !((i_pair_0+1) == i_pair_1) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        i_pair++;
      }
    }
    pairs_n = pairs_temp;
    pairs_check = 1;
  }
  
  for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){
    
    // only calculate the pair if the weight is !=0
    if(pairs_n(i_pair,2)!=0){
      // get the two coordinates 
      int i_1= pairs_n(i_pair,0);
      int i_2 = pairs_n(i_pair,1);
      
      Eigen::VectorXd yp(2);
      yp(0)= y(i_1);
      yp(1) = y(i_2);
      
      Eigen::VectorXd xbp(2);
      xbp(0)= xb(i_1);
      xbp(1)= xb(i_2);
      
      Eigen::MatrixXd Lambdap(2,2);
      Lambdap.setZero();
      Lambdap(0,0)= Lambda(i_1,i_1);
      Lambdap(0,1)= Lambda(i_1,i_2);
      Lambdap(1,0)= Lambda(i_2,i_1);
      Lambdap(1,1)= Lambda(i_2,i_2);
      
      double prp = prob_ordered_2(xbp,yp,Lambdap,alt,dtauk);

      pr += prp;
      
      // gradient calculation 
      Eigen::VectorXd gradp(5+alt-1);
      gradp.setZero(); 
      gradp = grad_po2(xbp,yp,Lambdap,alt,dtauk);
      
      // deriv due to xb. 
      grad(i_1) += gradp(0);
      grad(i_2) += gradp(1);
      // deriv w.r.t. Lambda entries. 
      int ind_i1i1 = i_1*Tp_n - i_1*(i_1-1)/2;
      int ind_i2i2 = i_2*Tp_n - i_2*(i_2-1)/2;
      int ind_i1i2 = ind_i1i1 + i_2-i_1;
      
      grad(Tp_n+ind_i1i1) += gradp(2);
      grad(Tp_n+ind_i2i2) +=gradp(3);
      grad(Tp_n+ind_i1i2) +=gradp(4);
      grad.segment(np,alt-1) += gradp.segment(5,alt-1);
    }
  }
  
  // return value. 
  return grad;
}

Eigen::MatrixXd hess_ordered_CML_approx(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type){
  
  // calculate the number of parameters. 
  int Tp_n= xb.size();
  
  int n_par = Tp_n + Tp_n*(Tp_n+1)/2 + alt-1;
  int np = n_par - alt+1;
  
  Eigen::VectorXd grad(n_par);
  grad.setZero();
  
  Eigen::MatrixXd hess(n_par,n_par);
  hess.setZero();
  
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  int lth_pairs_n = (Tp_n-1)*Tp_n/2;
  
  // Define a matrix that allows to specify the pairs (possible extension to higher order CMLs)
  Eigen::MatrixXd pairs_n;
  int pairs_check = 0;
  
  // if there were no valid pairs for the decision maker, construct them
  if(pairs_check==0){
    Eigen::MatrixXd pairs_temp(lth_pairs_n,3);
    pairs_temp.setOnes();
    // if the pairs are not predefined, build them
    int i_pair = 0;
    for(int i_pair_0 = 0; i_pair_0 < (Tp_n-1); i_pair_0 ++){
      for(int i_pair_1 = (i_pair_0+1);i_pair_1 < (Tp_n); i_pair_1 ++){
        // first choice occasion of the pair
        pairs_temp(i_pair,0) = i_pair_0;
        // second choice occasion of the pair
        pairs_temp(i_pair,1) = i_pair_1;
        
        // define the according weights to exclude undesired pairs
        if(cml_pair_type == 1){
          // for adjacent pairs looped:
          if( !(((i_pair_0+1) == i_pair_1) || ((i_pair_0 == 0) && (i_pair_1 == (Tp_n-1)))) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        if(cml_pair_type == 2){
          // for adjacent pairs chained
          if( !((i_pair_0+1) == i_pair_1) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        i_pair++;
      }
    }
    pairs_n = pairs_temp;
    pairs_check = 1;
  }
  
  for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){
    
    // only calculate the pair if the weight is !=0
    if(pairs_n(i_pair,2)!=0){
      grad.setZero();
      
      // get the two coordinates 
      int i_1= pairs_n(i_pair,0);
      int i_2 = pairs_n(i_pair,1);
      
      Eigen::VectorXd yp(2);
      yp(0)= y(i_1);
      yp(1) = y(i_2);
      
      Eigen::VectorXd xbp(2);
      xbp(0)= xb(i_1);
      xbp(1)= xb(i_2);
      
      Eigen::MatrixXd Lambdap(2,2);
      Lambdap.setZero();
      Lambdap(0,0)= Lambda(i_1,i_1);
      Lambdap(0,1)= Lambda(i_1,i_2);
      Lambdap(1,0)= Lambda(i_2,i_1);
      Lambdap(1,1)= Lambda(i_2,i_2);
      
      double prp = prob_ordered_2(xbp,yp,Lambdap,alt,dtauk);
      
      pr += prp;
      
      // gradient calculation 
      Eigen::VectorXd gradp(5+alt-1);
      gradp.setZero(); 
      gradp = grad_po2(xbp,yp,Lambdap,alt,dtauk);
      
      // deriv due to xb. 
      grad(i_1) += gradp(0);
      grad(i_2) += gradp(1);
      // deriv w.r.t. Lambda entries. 
      int ind_i1i1 = i_1*Tp_n - i_1*(i_1-1)/2;
      int ind_i2i2 = i_2*Tp_n - i_2*(i_2-1)/2;
      int ind_i1i2 = ind_i1i1 + i_2-i_1;
      
      grad(Tp_n+ind_i1i1) += gradp(2);
      grad(Tp_n+ind_i2i2) +=gradp(3);
      grad(Tp_n+ind_i1i2) +=gradp(4);
      grad.segment(np,alt-1) += gradp.segment(5,alt-1);
      
      // update the approximation of the Hessian
      hess += grad *grad.transpose();
    }
  }
  
  // return value. 
  return hess;
}

Eigen::MatrixXd hess_ordered_CML(Eigen::VectorXd xb, Eigen::VectorXd y, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, int cml_pair_type){
  
  // calculate the number of parameters. 
  int Tp_n= xb.size();
  
  int n_par = Tp_n + Tp_n*(Tp_n+1)/2 + alt-1;
  int np = n_par - alt+1;
  
  Eigen::VectorXd grad(n_par);
  grad.setZero();
  
  Eigen::MatrixXd Hess(n_par,n_par); 
  Hess.setZero(); 
  double pr = 0;
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  // normalize variance matrix 
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();
  
  //............................................................//
  // Normalize upper limits and covariance matrix (to correlation matrix)
  //............................................................//
  
  Eigen::VectorXd x_norm = xb.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
  
  int lth_pairs_n = (Tp_n-1)*Tp_n/2;
  
  // Define a matrix that allows to specify the pairs (possible extension to higher order CMLs)
  Eigen::MatrixXd pairs_n;
  int pairs_check = 0;
  
  // if there were no valid pairs for the decision maker, construct them
  if(pairs_check==0){
    Eigen::MatrixXd pairs_temp(lth_pairs_n,3);
    pairs_temp.setOnes();
    // if the pairs are not predefined, build them
    int i_pair = 0;
    for(int i_pair_0 = 0; i_pair_0 < (Tp_n-1); i_pair_0 ++){
      for(int i_pair_1 = (i_pair_0+1);i_pair_1 < (Tp_n); i_pair_1 ++){
        // first choice occasion of the pair
        pairs_temp(i_pair,0) = i_pair_0;
        // second choice occasion of the pair
        pairs_temp(i_pair,1) = i_pair_1;
        
        // define the according weights to exclude undesired pairs
        if(cml_pair_type == 1){
          // for adjacent pairs looped:
          if( !(((i_pair_0+1) == i_pair_1) || ((i_pair_0 == 0) && (i_pair_1 == (Tp_n-1)))) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        if(cml_pair_type == 2){
          // for adjacent pairs chained
          if( !((i_pair_0+1) == i_pair_1) ){
            pairs_temp(i_pair,2) = 0;
          }
        }
        i_pair++;
      }
    }
    pairs_n = pairs_temp;
    pairs_check = 1;
  }
  
  for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){
    
    // only calculate the pair if the weight is !=0
    if(pairs_n(i_pair,2)!=0){
      // get the two coordinates 
      int i_1= pairs_n(i_pair,0);
      int i_2 = pairs_n(i_pair,1);
      
      Eigen::VectorXd yp(2);
      yp(0)= y(i_1);
      yp(1) = y(i_2);
      
      Eigen::VectorXd xbp(2);
      xbp(0)= xb(i_1);
      xbp(1)= xb(i_2);
      
      Eigen::MatrixXd Lambdap(2,2);
      Lambdap.setZero();
      Lambdap(0,0)= Lambda(i_1,i_1);
      Lambdap(0,1)= Lambda(i_1,i_2);
      Lambdap(1,0)= Lambda(i_2,i_1);
      Lambdap(1,1)= Lambda(i_2,i_2);
      
      double prp = prob_ordered_2(xbp,yp,Lambdap,alt,dtauk);
 
      pr += prp;
      
      // gradient calculation 
      Eigen::VectorXd gradp(5+alt-1);
      gradp.setZero(); 
      gradp = grad_po2(xbp,yp,Lambdap,alt,dtauk);
      
      // match indices 
      Eigen::VectorXd ind(4+alt);
      ind.setZero();
      
      // deriv due to xb. 
      grad(i_1) += gradp(0);
      grad(i_2) += gradp(1);
      ind(0)= i_1;
      ind(1)= i_2; 
      // deriv w.r.t. Lambda entries. 
      int ind_i1i1 = i_1*Tp_n - i_1*(i_1-1)/2;
      int ind_i2i2 = i_2*Tp_n - i_2*(i_2-1)/2;
      int ind_i1i2 = ind_i1i1 + i_2-i_1;
      
      ind(2) = Tp_n+ind_i1i1;
      ind(3) = Tp_n+ind_i2i2;
      ind(4) = Tp_n+ind_i1i2;
      for (int j=0;j<alt-1;j++){
        ind(5+j)=np+j;
      }
      
      grad(Tp_n+ind_i1i1) += gradp(2);
      grad(Tp_n+ind_i2i2) +=gradp(3);
      grad(Tp_n+ind_i1i2) +=gradp(4);
      grad.segment(np,alt-1) += gradp.segment(5,alt-1);
      

      //////////////////
      // Hessian 
      //////////////////
      Eigen::MatrixXd Hessp(4+alt,4+alt);
      Hessp.setZero(); 
      Hessp = hess_po2(xbp,yp,Lambdap,alt,dtauk);
      
      for (int ia=0;ia<ind.size();ia++){
        for (int ib=0;ib<ind.size();ib++){
          Hess(ind(ia),ind(ib)) += Hessp(ia,ib);
        }
      }     
      
      
    }
  }
  
  // return value. 
  return Hess;
}

//' calculates the probability of one choice in an ordered probit model
//' @description
//' calculates the probability of one choice in an ordered probit model
//' @param xb
//' real; systematic utility
//' @param y1
//' integer; choice
//' @param Lambda 
//' real; variance
//' @param alt
//' integer; number of alternatives
//' @param dtauk
//' vector of category limits. 
//' @return 
//' matrix
//' @keywords internal
//'
//[[Rcpp::export]]
double prob_ordered_1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk)
{
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  double lam = sqrt(Lambda);
  
  double x_norm = xb/lam;
  
  double lower_bound;
  lower_bound = tauk(0)-10.0;
  
  double upper_bound;
  upper_bound= tauk(alt-2)+10.0;
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  
  if (y1==1){
    upper_bound =tauk(0);
  }
  if (y1==alt){
    lower_bound = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound = tauk(y1-2);
      upper_bound = tauk(y1-1);
    }
  }
  
  //
  double ub = upper_bound/lam-x_norm;
  double lb = lower_bound/lam-x_norm;
  
  double pu = std_normal_cdf(ub);
  double pl = std_normal_cdf(lb);
  
  double pr = pu-pl;
  double lpr = log(pr);
  return lpr;
}



Eigen::VectorXd grad_po1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk)
{
  // derivative with respect to xb, Lambda, tauk
  // gradu equals derivative of upper bound prob, gradl der. of prob lower bound. 
  
  
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  Eigen::VectorXd grad(tauk.size()+2), gradu(tauk.size()+2), gradl(tauk.size()+2);
  grad.setZero();
  gradu.setZero();
  gradl.setZero();
  
  // calculation for probability
  double lam = sqrt(Lambda);
  double x_norm = xb/lam;
  
  double lower_bound;
  lower_bound = tauk(0)-10.0;
  double upper_bound;
  upper_bound= tauk(alt-2)+10.0;
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  
  if (y1==1){
    upper_bound =tauk(0);
  }
  if (y1==alt){
    lower_bound = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound = tauk(y1-2);
      upper_bound = tauk(y1-1);
    }
  }
  
  //
  double ub = upper_bound/lam-x_norm;
  double lb = lower_bound/lam-x_norm;
  

  double pu = std_normal_cdf(ub);
  double pl = std_normal_cdf(lb);
  
  double pr = pu-pl;
  // now start calculations for gradient 
  // at upper bound 
  double gr_u = std_normal_pdf(ub);
  
  // deriv w.r.t. tau_k.
  Eigen::VectorXd iota(alt-1);
  //iota.setOnes(); 
  iota = edtauk; 
  iota(0)=1; 
  
  gradu(0) = -1/lam;
  gradu(1) = -ub/(2*Lambda);
  if (y1<alt){
    gradu.segment(2,y1) = iota.head(y1)/lam;
  }

  gradu = gradu*gr_u;
  
  // at lower bound 
  double gr_l = std_normal_pdf(lb);
  
  gradl(0) = -1/lam;
  gradl(1) = -lb/(2*Lambda);
  if (y1>1){
    gradl.segment(2,y1-1) = iota.head(y1-1)/lam;
  }
  
  gradl = gradl*gr_l;

  // put pieces together
  grad = gradu -gradl; 
  
  // logarithms.
  grad = grad/pr;
  
  return grad;
}



Eigen::MatrixXd hess_po1(double xb, int y1, double Lambda, int alt, Eigen::VectorXd dtauk)
{
  
  Eigen::MatrixXd Hess(2+dtauk.size(),2+dtauk.size());
  Hess.setZero();
  
  // calculate difference 
  // tauk's. 
  Eigen::VectorXd tauk(alt-1),edtauk(alt-1);
  tauk.setZero();
  edtauk.setZero();
  edtauk(0)=1.0; 
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    edtauk(j) = std::exp(dtauk(j));
    tauk(j)=tauk(j-1)+edtauk(j);
  }
  
  Eigen::VectorXd grad(tauk.size()+2), gradu(tauk.size()+2), gradl(tauk.size()+2);
  grad.setZero();
  gradu.setZero();
  gradl.setZero();
  
  /////////////////////////////////////
  // calculation for probability
  /////////////////////////////////////
  double lam = sqrt(Lambda);
  double x_norm = xb/lam;
  
  double lower_bound;
  lower_bound = tauk(0)-10.0;
  double upper_bound;
  upper_bound= tauk(alt-2)+10.0;
  
  if (y1<1){ y1=1;}
  if (y1>alt){ y1=alt;}
  
  if (y1==1){
    upper_bound =tauk(0);
  }
  if (y1==alt){
    lower_bound = tauk(alt-2);
  }
  
  if (y1>1){
    if(y1<alt){
      lower_bound = tauk(y1-2);
      upper_bound = tauk(y1-1);
    }
  }
  
  //
  double ub = upper_bound/lam-x_norm;
  double lb = lower_bound/lam-x_norm;
  
  //Rcout << "ub: " << ub << ", lb:  " << lb << std::endl; 
  
  double pu = std_normal_cdf(ub);
  double pl = std_normal_cdf(lb);
  
  double pr = pu-pl;
  
  /////////////////////////////////
  // now start calculations for gradient 
  // at upper bound 
  double gr_u = std_normal_pdf(ub);
  
  // deriv w.r.t. tau_k.
  Eigen::VectorXd iota(alt-1);
  iota = edtauk; 
  iota(0) = 1.0;
  
  gradu(0) = -1/lam;
  gradu(1) = -ub/(2*Lambda);
  if (y1<alt){
    gradu.segment(2,y1) = iota.head(y1)/lam;
  }
  // at lower bound 
  double gr_l = std_normal_pdf(lb);
  
  gradl(0) = -1/lam;
  gradl(1) = -lb/(2*Lambda);
  if (y1>1){
    gradl.segment(2,y1-1) = iota.head(y1-1)/lam;
  }
  
  // put pieces together
  grad = gradu *gr_u -gradl*gr_l; 
  
  /////////////////////////////
  // Hessian calculation 
  /////////////////////////////
  
  Eigen::MatrixXd Hessu(2+tauk.size(),2+tauk.size()), Hessl(2+tauk.size(),2+tauk.size()); 
  Hessu.setZero();
  Hessl.setZero();
  
  // Hessian of upper bound prob. 
  Hessu += gradu * gradu.transpose()* (-ub);
  Hessu(0,1) += 1/(2*Lambda*lam);
  Hessu(1,0) += 1/(2*Lambda*lam);
  Hessu(1,1) += 3*ub/(4*Lambda*Lambda);
  
  if (y1<alt){
    Hessu.block(2,1,y1,1) += -iota.head(y1)/(2*Lambda*lam);
    Hessu.block(1,2,1,y1) += -iota.head(y1).transpose()/(2*Lambda*lam);
  }
  
  // same for lower bound.
  Hessl += gradl * gradl.transpose()* (-lb);
  Hessl(0,1) += 1/(2*Lambda*lam);
  Hessl(1,0) += 1/(2*Lambda*lam);
  Hessl(1,1) += 3*lb/(4*Lambda*Lambda);
  
  if (y1>1){
    Hessl.block(2,1,y1-1,1) += -iota.head(y1-1)/(2*Lambda*lam);
    Hessl.block(1,2,1,y1-1) += -iota.head(y1-1).transpose()/(2*Lambda*lam);
  }
  
  // put pieces together 
  Hess = Hessu * gr_u - Hessl * gr_l;
  
  // respect taking logarithms. 
  Hess = Hess/pr - grad * grad.transpose()/(pr*pr); 
  
  // return values
  
  return Hess;
}

