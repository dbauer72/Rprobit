#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

//' duplication matrix
//' @description
//' Computes the duplication matrix converting from vech to vec..
//' @param n
//' integer; dimension of matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd duplmat( int n){
  // This function returns a matrix sith n * n rows and n * (n + 1 ) / 2
  //columns that transforms vech( A ) to vec( A ) where A is a symmetric n by n matrix
  
  // n = a positive integer value for the order of the matrix
  
  int p = n * ( n + 1 ) / 2;
  int     nsq = n * n;
  Eigen::MatrixXd Dt(p,nsq);
  Dt.setZero();
  
  Eigen::MatrixXd II(n,n);
  II.setIdentity();
  
  Eigen::MatrixXd I(p,p);
  I.setIdentity();
  
  Eigen::MatrixXd k(n,n);
  k.setZero();
  
  int ii = 0;
  int jj = 0;
  
  for(int i = 0; i < p; ++i){ 
    k(ii,jj) = i;
    ii++;
    if(ii == n){
      ii = jj +1;
      jj = ii;
    }
  }
  
  for(int j = 0; j < n; ++j){ 
    
    for(int i = j; i < n; ++i){ 
      
      Eigen::MatrixXd Iij = I.col(k(i,j));
      Eigen::VectorXd cc = II.row(j);
      Eigen::VectorXd rr = II.row(i);
      Eigen::MatrixXd Eij = (rr*cc.transpose());
      Eigen::MatrixXd Eji = (cc*rr.transpose());
      Eigen::MatrixXd Tij(Eij.rows(), Eij.cols());
      if(i == j){
        Tij = Eij;
      }else{
        Tij = Eij + Eji;
      }
      Tij.resize(1,n*n);
      
      Dt += (Iij*Tij);
    }
  }
  Eigen::MatrixXd DD = Dt.transpose();
  
  return DD;
}


//' vech indices
//' @description
//' Returns indices of vech(X) w.r.t to the entries of X WITHOUT THE DIAGONAL
//' @param k
//' integer; dimension of matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd vechor(int k) {
  //Returns indices of vech(X) w.r.t to the entries of X WITHOUT THE DIAGONAL
  int m = ((k*k)-k)/2 + k; //number of paramters in vech(X)
  
  Eigen::MatrixXd ind((m-k),2);
  ind.setZero();
  
  int gg = 1;//indices for cpp
  int hh = 0;//indices for cpp
  
  for(int i = 0; i < (m-k); ++i){
    ind(i,0) = gg;
    ind(i,1) = hh;
    gg++;
    if(gg == (k)){//indices for cpp
      hh = hh + 1;
      gg = hh + 1;
    }
    
  }
  
  
  return ind;
}

//' indices of diagonal elements in vectorized kxk square matrix.
//' @description
//' Returns indices of diagonal elements in vectorized kxk square matrix.
//' @param k
//' integer; dimension of matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::VectorXd fdiag(int k) {
  Eigen::VectorXd ind(k);
  ind.setZero();
  int kk = k;
  ind[0] = 0;//indices for cpp
  for(int i = 1;i < (k); i ++){
    ind[i] = ind[i-1]+kk;
    kk = (kk-1);
  }
  
  return ind;
}

//' Returns indices of vech(X) w.r.t to the entries of X
//' @description
//' Returns indices of vech(X) w.r.t to the entries of X
//' @param k
//' integer; dimension of matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd vechor_diag(int k) {
  //Returns indices of vech(X) w.r.t to the entries of X
  int m = ((k*k)-k)/2 + k; //number of paramters in vech(X)
  
  Eigen::MatrixXd ind(m,2);
  ind.setZero();
  
  int gg = 0;//indices for cpp
  int hh = 0;//indices for cpp
  
  for(int i = 0; i < (m); ++i){
    ind(i,0) = gg;
    ind(i,1) = hh;
    gg++;
    if(gg == (k)){//indices for cpp
      hh = hh + 1;
      gg = hh ;
    }
    
  }
  
  
  return ind;
}



//' Returns commutation matrix such that C vec(X) = vec(X')
//' @description
//' Returns commutation matrix such that C vec(X) = vec(X')
//' @param r
//' integer; number of rows in the matrix.
//' @param c
//' integer; number of columns in matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd commmat( int r, int c){
  // This function constructs the rc by rc commutation matrix
  
  // r = a positive integer value for the number of rows in the H matrices
  // c = a positive integer value for the number of columns in the H matrices
  
  int p = r*c;
  Eigen::MatrixXd K(p,p);
  K.setZero();
  Eigen::MatrixXd Ir(r,r);
  Eigen::MatrixXd Ic(c,c);
  Ir.setIdentity();
  Ic.setIdentity();
  
  for(int i = 0; i < r; ++i){ 
    for(int j = 0; j < c; ++j){ 
      Eigen::VectorXd cc = Ic.row(j);
      Eigen::VectorXd rr = Ir.row(i);
      Eigen::MatrixXd Hij = rr*cc.transpose();
      Eigen::MatrixXd tHij = Hij.transpose();
      K += kroneckerProduct(Hij,tHij);
    }
  }
  
  
  return K;
}


//' Returns elimination matrix such that C vec(X) = vech(X)
//' @description
//' Returns elimination matrix such that C vec(X) = vech(X)
//' @param n
//' integer; dimension of matrix.
//' @return 
//' matrix
//' @keywords internal 
//'
//[[Rcpp::export]]
Eigen::MatrixXd elimmat( int n){
  // This function constructs the p by n*n elemination matrix
  //Which maps from vec(A) to vech(A)
  
  // n = a positive integer value for the order of the matrix
  
  int p = n * ( n + 1 ) / 2;
  int     nsq = n * n;
  Eigen::MatrixXd L(p,nsq);
  L.setZero();
  
  Eigen::MatrixXd I(p,p);
  I.setIdentity();
  
  Eigen::MatrixXd II(n,n);
  II.setIdentity();
  
  Eigen::MatrixXd k(n,n);
  k.setZero();
  
  int ii = 0;
  int jj = 0;
  
  for(int i = 0; i < p; ++i){ 
    k(ii,jj) = i;
    ii++;
    if(ii == n){
      ii = jj +1;
      jj = ii;
      
    }
  }
  
  for(int j = 0; j < n; ++j){ 
    
    for(int i = j; i < n; ++i){ 
      
      Eigen::MatrixXd Iij = I.col(k(i,j));
      Eigen::VectorXd cc = II.row(j);
      Eigen::VectorXd rr = II.row(i);
      Eigen::MatrixXd IIij = (rr*cc.transpose());
      IIij.resize(1,n*n);
      L += (Iij*IIij);
    }
  }
  
  return L;
}

//https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}


// extracts the indices of the elements of xc that are not 'NA'.
IntegerVector NA_ind(Eigen::MatrixXd Xc){
  
  IntegerVector vecNA(Xc.rows());
  
  int j = 0;   //counter for vx
  
  for (int i=0;i<vecNA.size();i++) {
    if (!(R_IsNA(Xc(i,0)))) {
      //replace and update j
      vecNA[j] = i; //Xc[j];
      j++;
    }
  }
  return vecNA[Rcpp::Range(0, j-1)];
}

////////////////////////////////////////////7
/// helpers for TVBS     /////
//////////////////////////////


Eigen::MatrixXd reorder_matrix(Eigen::MatrixXd M, Eigen::VectorXd indord)
{
  Eigen::MatrixXd M_ord(M.rows(),M.cols()); 
  
  int m = M.rows(); 
  
  Eigen::MatrixXd T_ord(m,m);
  T_ord.setZero(); 
  
  for (int r=0;r<m;r++){
    T_ord(r,indord(r))=1;
  }
  
  M_ord = T_ord.transpose() * M * T_ord; 
  
  return M_ord; 
}



/// extract submatrix. 
Eigen::MatrixXd extract_submatrix(Eigen::MatrixXd Mold, Eigen::VectorXd indin)
{
  int n = indin.size();
  int m = Mold.cols();
  if (m>1){
    m = n;
  }
  
  
  Eigen::MatrixXd MM(n,m);
  MM.setZero(); 
  if (m>1){
    for (int a=0;a<n;a++){
      for (int b=0;b<n;b++){
        MM(a,b) = Mold(indin(a),indin(b));
      }
    }
  }
  if (m==1){
    // reduce a vector.
    
    for (int a=0;a<n;a++){
      MM(a,0) = Mold(indin(a),0);
    }
  }
  return MM;
}


Eigen::VectorXd reorder_indices(Eigen::VectorXd ind_ord)
{
  int m= ind_ord.size();
  int cur = 0;
  int npar = m+ m*(m-1)/2;
  
  Eigen::MatrixXd T_reord(m,m);
  T_reord.setZero();
  for (int r=0;r<m;r++){
    T_reord(r,ind_ord(r))=1;
  }
  
  Eigen::VectorXd ind_reord(npar);
  ind_reord.setZero(); 
  ind_reord.head(m) = ind_ord;
  
  Eigen::MatrixXd mat_ind(m,m);
  mat_ind.setZero(); 
  cur = m;
  for (int r=0;r<(m-1);r++){
    for (int c=(r+1);c<m;c++){
      mat_ind(r,c)= cur;
      mat_ind(c,r)= cur;
      cur++;
    }
  }
  
  
  mat_ind = T_reord * mat_ind * T_reord.transpose();
  
  cur = m;
  for (int r=0;r<(m-1);r++){
    for (int c=(r+1);c<m;c++){
      ind_reord(cur) = mat_ind(r,c);
      cur++;
    }
  }
  
  return ind_reord;
}

