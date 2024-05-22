#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies


//#include "toms462.h"  // allows bivariate normal cdf calculations
//#include "distrib.h"
//#include "lin_alg.h"
#include "TVBS_vdb.h"
#include "TVBS_vdb_helper.h"
#include "TVBS_vdb_helper_2.h"


#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
// #include <random>       // for random numbers

extern double pi_global;

extern double tol;


struct TVBS_Matrix_two
{
  Eigen::MatrixXd mat1, mat2;
};

struct TVBS_Matrix_three
{
  Eigen::MatrixXd mat1, mat2, mat3;
};

struct TVBS_Matrix_four
{
  Eigen::MatrixXd mat1, mat2, mat3, mat4;
};

struct TVBS_double_three
{
  double d1, d2, d3;
};

struct TVBS_double_Vector_two
{
  // not used
  double d1;
  Eigen::VectorXd v1, v2;
};

struct TVBS_Vector_three
{
  Eigen::VectorXd v1, v2, v3;
};

struct TVBS_Vector_four
{
  Eigen::VectorXd v1, v2, v3, v4;
};

struct TVBS_Vec_Mat
{
  Eigen::VectorXd vec1;
  Eigen::MatrixXd mat1;
};

struct TVBS_Vec_Mat_Vec
{
  Eigen::VectorXd vec1;
  Eigen::MatrixXd mat1;
  Eigen::VectorXd vec2;
};


extern int x1symmetric_global;
extern int x2symmetric_global;
extern int x1diagonal_global;
extern int x2diagonal_global;    // = 0 always
extern int x2correlation_global;
extern int xinvsymmetric_global;
extern int xinvcorrelation_global;
extern int xinvdiagonal_global;    // = 0 always
extern int omsymmetric_global;
extern int omdiagonal_global;    // = 0 always
extern int condcov_global;
extern int cholesky_global;    // = 0 always
extern int condcovmeantrunc_global;
extern int condcovsigtrunc_global;
extern int optimal_global;
extern int covarr_global;    // = 1 always

extern int counter_check;

extern double pi_global;
extern double tol;


// ex from vdb 





Eigen::MatrixXd TVBS_g_a_omega_b_cpp(Eigen::MatrixXd x1, Eigen::MatrixXd x2) //, int omsymmetric = 0, int omdiagonal = 0) 
{
  int k = x1.cols();
  
  Eigen::MatrixXd g = kroneckerProduct(x1.transpose(), x2);
  int nrows = g.rows();
  int ncols = g.cols();
  
  if(omsymmetric_global==1)
  {
    Eigen::MatrixXd temp1 = TVBS_vec_symmetry_cpp(k);
    g = temp1*g;
    nrows = g.rows();
    ncols = g.cols();
    
    if(omdiagonal_global==1)
    {
      int diag_inc = 1;
      Eigen::VectorXd temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k), diag_inc);
      Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix);
      nrows = row_select.size();
      for(int j = 0; j<nrows; j++)
      {
        g.row(j) = g.row(row_select(j));
      }
    }
  }
  
  Eigen::MatrixXd g_out = g.block(0,0,nrows,ncols);
  return g_out;
}




Eigen::MatrixXd TVBS_g_cholesky_cov_cpp(Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd sigma_diag_sqrt = Sigma.diagonal().array().sqrt();
  
  int parm_dim = sigma_diag_sqrt.size();
  int parm_cov = parm_dim*(parm_dim+1)/2;
  
  // cholesky decomposition t(L)%*%L=Sigma
  Eigen::MatrixXd Sigma_star_chol = Sigma.llt().matrixU();
  
  Eigen::VectorXd Sigma_star_chol_diag = Sigma_star_chol.diagonal();
  
  Eigen::MatrixXd dg(parm_cov,parm_cov);
  dg.setZero();
  
  int diag_inc =1; 
  Eigen::VectorXd Sigma_star_chol_vech = TVBS_lower_tri_entries_cpp(Sigma_star_chol.transpose(), diag_inc);
  
  int j = 1;
  int mm = 1;
  
  if(parm_dim>1)
  {
    while(j<=parm_cov)
    {
      int gss_c = j;
      int l = 1;
      while(gss_c<=parm_cov)
      {
        dg.block(j+l-1-1,gss_c-1,(j+parm_dim-mm) - (j+l-1) + 1, (gss_c+parm_dim-l+1-mm) - gss_c + 1) = (Eigen::MatrixXd::Identity(parm_dim-l+2-mm, parm_dim-l+2-mm).array()*Sigma_star_chol_vech(j+l-1-1)).matrix();
        dg.block(j+l-1-1,gss_c-1,1,(gss_c+parm_dim-l+1-mm) - gss_c + 1) += (Sigma_star_chol_vech.segment(j+l-1-1, (gss_c+parm_dim-l+1-mm) - gss_c + 1)).transpose();
        gss_c += parm_dim-l+1-mm+1;
        l++;
      }
      j += parm_dim-mm+1;
      mm++;
    }
  }
  if(parm_dim==1)
  {
    dg.setOnes();
    dg = dg.block(0,0,1,1);
  }
  
  return dg;
}




Eigen::MatrixXd TVBS_g_cholesky_cor_cpp(Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd sigma_diag_sqrt = Sigma.diagonal().array().sqrt();
  
  int parm_dim = sigma_diag_sqrt.size();
  int parm_cor = parm_dim*(parm_dim-1)/2;
  
  // cholesky decomposition t(L)%*%L=Sigma
  Eigen::MatrixXd Sigma_star_chol = Sigma.llt().matrixU();
  Eigen::VectorXd Sigma_star_chol_diag = Sigma_star_chol.diagonal();
  
  Eigen::MatrixXd dg(parm_cor,parm_cor);
  dg.setZero();
  
  
  int j = 1;
  int mm = 1;
  
  if(parm_dim>2)
  {
    for(mm; mm<=parm_dim; mm++)
    {
      Eigen::MatrixXd ss, ss1;
      int gss_c = j;
      int l = 1;
      if(j<parm_cor)
      {
        if(parm_dim-mm-1>1)
        {
          Eigen::MatrixXd ss_lhs(parm_dim-mm-1,parm_dim-mm-1), ss_rhs(parm_dim-mm-1,parm_dim-mm-1), temp_matrix_ones(parm_dim-mm-1,parm_dim-mm-1);
          ss_lhs.setZero();
          ss_rhs.setZero();
          temp_matrix_ones.setOnes();
          ss_lhs.triangularView<Eigen::Upper>() = temp_matrix_ones.triangularView<Eigen::Upper>();
          Eigen::MatrixXd ss_rhs_diag = (Sigma_star_chol.block(mm-1, mm+2-1, 1, parm_dim - (mm+2)+1));
          ss_rhs.diagonal() = ss_rhs_diag;
          ss = ss_lhs*ss_rhs;
          
          
          Eigen::MatrixXd ss1_1((parm_dim-1)-(mm+1)+1, (parm_dim-1)-(mm+1)+1);
          ss1_1.setZero();
          ss1_1.diagonal() = (1.0/Sigma_star_chol_diag.segment(mm+1-1, (parm_dim-1)-(mm+1)+1).array());
          
          Eigen::MatrixXd ss1_2_vec = Sigma_star_chol.block(mm-1, mm+1-1, 1, parm_dim-1 - (mm+1)+1);
          Eigen::MatrixXd ss1_2(parm_dim-1 - (mm+1)+1, parm_dim-1 - (mm+1)+1);
          ss1_2.setZero();
          ss1_2.diagonal() = ss1_2_vec;
          
          Eigen::MatrixXd ss1_3_full = Sigma_star_chol.block(mm+1-1, mm+2-1, parm_dim-1-(mm+1)+1, parm_dim - (mm+2)+1);
          Eigen::MatrixXd ss1_3(parm_dim-1-(mm+1)+1, parm_dim-1-(mm+1)+1);
          ss1_3.setZero();
          ss1_3.triangularView<Eigen::Upper>() = ss1_3_full.triangularView<Eigen::Upper>();
          
          ss1 = ss - ss1_1*ss1_2*ss1_3;
        }
        if(parm_dim-mm-1==1)
        {
          ss = Sigma_star_chol.block(mm-1, (mm+2)-1,1,1);
          ss1 = (ss.array() - (1.0/Sigma_star_chol_diag((mm+1)-1))*Sigma_star_chol(mm-1, (mm+1)-1)*Sigma_star_chol((mm+1)-1, (mm+2)-1)).matrix();
        }
      }
      while(gss_c<=parm_cor)
      {
        dg.block((j+l-1)-1, gss_c-1, (j+parm_dim-mm-1)-(j+l-1)+1, (gss_c+parm_dim-l-mm)-gss_c+1) = Eigen::MatrixXd::Identity(parm_dim-l+1-mm, parm_dim-l+1-mm)*Sigma_star_chol(mm-1,(mm+l-1)-1);
        
        if(gss_c<parm_cor)
        {
          dg.block((j+l-1)-1, (gss_c+parm_dim-l-mm+1)-1, 1, (gss_c+2*parm_dim-2*l-2*mm)-(gss_c+parm_dim-l-mm+1)+1) = ss1.block(l-1, l-1, 1, ss.cols()-l+1);
        }
        gss_c = gss_c+parm_dim-l+1-mm;
        l++;
        
      }
      j += parm_dim-mm;
    }
  }
  if(parm_dim==2)
  {
    dg.setOnes();
    dg = dg.block(0,0,1,1);
  }
  
  return dg;
}





// input: x_norm, Sigma_diag_sqrt, Cor_mat
struct TVBS_Matrix_two TVBS_grad_cor_cov_cpp(Eigen::VectorXd x_norm, Eigen::VectorXd Sigma_diag_sqrt, Eigen::MatrixXd Cor_mat)
{
  // dimension of the matrix/problem
  int dim_lcl = Sigma_diag_sqrt.size();
  
  // number of lower diagonal entries (with diagonal)
  int row_lcl = dim_lcl*(dim_lcl+1)/2;
  
  // number of lower diagonal entries (without diagonal)
  int col_lcl = dim_lcl*(dim_lcl-1)/2;
  
  Eigen::VectorXd nu_1_lcl = Sigma_diag_sqrt.array()*Sigma_diag_sqrt.array();
  Eigen::VectorXd nu_2_lcl = Sigma_diag_sqrt;
  
  Eigen::MatrixXd temp(dim_lcl, dim_lcl);
  temp.setZero();
  temp.diagonal() = Sigma_diag_sqrt;
  Eigen::MatrixXd Cov_mat = temp*Cor_mat*temp;
  
  Eigen::MatrixXd grm_lcl(row_lcl, dim_lcl);
  grm_lcl.setZero();
  
  Eigen::VectorXd cc1_lcl = 0.5*(x_norm.array()/nu_1_lcl.array());
  int l_lcl = 0;
  
  for(int idimvn=1; idimvn<=dim_lcl; idimvn++)
  {
    grm_lcl((l_lcl+1)-1,idimvn-1) = cc1_lcl(idimvn-1);
    l_lcl += dim_lcl+1-idimvn;
  }
  
  Eigen::MatrixXd grc_lcl(row_lcl, col_lcl);
  grc_lcl.setZero();
  
  int j = 1;
  int mm = 1;
  int c1 = 1;
  
  while(j != row_lcl)
  {
    Eigen::VectorXd n_diag_1_vec = 1.0/(nu_2_lcl(mm-1)*nu_2_lcl.segment((mm+1)-1, dim_lcl-(mm+1)+1).array());
    Eigen::MatrixXd n_diag_1(dim_lcl-mm, dim_lcl-mm);
    n_diag_1.setZero();
    n_diag_1.diagonal() = n_diag_1_vec;
    
    Eigen::MatrixXd Cor_mat_block = Cor_mat.block(mm-1, (mm+1)-1, 1, dim_lcl-(mm+1)+1);
    Eigen::MatrixXd n_hor_1 = (-0.5*Cor_mat_block.array()/Cov_mat(mm-1,mm-1));
    
    Cor_mat_block = Cor_mat.block(mm-1, (mm+1)-1, 1, dim_lcl-(mm+1)+1);
    Eigen::MatrixXd n_hor_2 = -0.5*Cor_mat_block.array()/(nu_1_lcl.segment((mm+1)-1, dim_lcl-(mm+1)+1)).array();
    
    grc_lcl.block((j+1)-1, c1-1, (j+dim_lcl-mm)-(j+1)+1, (c1+dim_lcl-mm-1)-c1+1) = n_diag_1;
    grc_lcl.block(j-1, c1-1, 1, (c1+dim_lcl-1-mm)-c1+1) = n_hor_1;
    
    int k = j;
    int l = 1;
    while(l!=(dim_lcl+1-mm))
    {
      grc_lcl((k+dim_lcl+2-l-mm)-1, (l+c1-1)-1) = n_hor_2(l-1,0);
      k += dim_lcl+2-l-mm;
      l++;
    }
    c1 += dim_lcl-mm;
    j += dim_lcl+1-mm;
    mm++;
  }
  
  grm_lcl = -1.0*grm_lcl;
  struct TVBS_Matrix_two output = {grm_lcl, grc_lcl};
  
  return output;
}



// issues arrise for some cases due to divergence to behaviour compared to R code (when tempselmat length is not equal to number of rows/columns to select from)
struct TVBS_Matrix_two TVBS_g_both_x_omega_x_cpp(Eigen::MatrixXd Mat_1, Eigen::MatrixXd Mat_2)
{
  // input dimensions
  int n = Mat_1.rows();
  int k = Mat_1.rows();
  
  Eigen::MatrixXd temp = Mat_2.transpose()*Mat_1.transpose();
  
  Eigen::MatrixXd t1(n*temp.rows(), n*n);
  t1.setZero();
  
  for(int i1 = 0; i1 < n; i1++)
  {
    t1.block(0,i1*n,n*temp.rows(),n) = kroneckerProduct(Eigen::MatrixXd::Identity(n, n), temp.col(i1));
  }
  
  Eigen::MatrixXd g_Mat1_cov = kroneckerProduct(Eigen::MatrixXd::Identity(n, n), Mat_2*Mat_1.transpose()) + t1;
  Eigen::MatrixXd g_Mat2_cov = kroneckerProduct(Mat_1.transpose(), Mat_1.transpose());
  
  if(x2symmetric_global==1)
  {
    Eigen::MatrixXd temp_sel_matrix(n,n);
    temp_sel_matrix.setZero();
    Eigen::MatrixXd temp_sel_matrix_ones(n,n);
    temp_sel_matrix_ones.setOnes();
    temp_sel_matrix.triangularView<Eigen::Lower>() = temp_sel_matrix_ones.triangularView<Eigen::Lower>();
    temp_sel_matrix.resize(n*n,1);
    Eigen::VectorXi col_select = TVBS_logic2position_cpp(temp_sel_matrix);
    
    // sub-matrix by columns
    {
      Eigen::MatrixXd g_Mat1_cov_new(g_Mat1_cov.rows(), col_select.size());
      g_Mat1_cov_new.setZero();
      for(int i1=0; i1<col_select.size(); i1++ )
      {
        g_Mat1_cov_new.col(i1) = g_Mat1_cov.col(col_select(i1));
      }
      g_Mat1_cov = g_Mat1_cov_new;
    }
    
    g_Mat2_cov = TVBS_g_asym_to_sym_cpp(g_Mat2_cov);
    
    int diag_inc = 1;
    temp_sel_matrix = TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(k, k), diag_inc);
    Eigen::VectorXi row_select = TVBS_logic2position_cpp(temp_sel_matrix);
    
    if(x1symmetric_global==1)
    {
      Eigen::MatrixXd temp1 = TVBS_vec_symmetry_cpp(n);;
      g_Mat1_cov = temp1*g_Mat1_cov;
    }
    if(x2diagonal_global==0)
    {
      if(x1diagonal_global==1)
      {
        // sub-matrix by rows
        {
          Eigen::MatrixXd g_Mat1_cov_new(row_select.size(), g_Mat1_cov.cols());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.row(i1) = g_Mat1_cov.row(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
      }
    }
    if(x2diagonal_global==1)
    {
      // sub-matrix by row
      {
        Eigen::MatrixXd g_Mat2_cov_new(row_select.size(), g_Mat2_cov.cols());
        g_Mat2_cov_new.setZero();
        for(int i1=0; i1<row_select.size(); i1++)
        {
          g_Mat2_cov_new.row(i1) = g_Mat2_cov.row(row_select(i1));
        }
        g_Mat2_cov = g_Mat2_cov_new;
      }
      
      if(x1diagonal_global==1)
      {
        // sub-matrix by row
        {
          Eigen::MatrixXd g_Mat1_cov_new(row_select.size(), g_Mat1_cov.cols());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.row(i1) = g_Mat1_cov.row(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
        
        // sub-matrix by column
        {
          Eigen::MatrixXd g_Mat1_cov_new(g_Mat1_cov.rows(), row_select.size());
          g_Mat1_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat1_cov_new.col(i1) = g_Mat1_cov.col(row_select(i1));
          }
          g_Mat1_cov = g_Mat1_cov_new;
        }
        
        // sub-matrix by column
        {
          Eigen::MatrixXd g_Mat2_cov_new(g_Mat2_cov.rows(), row_select.size());
          g_Mat2_cov_new.setZero();
          for(int i1=0; i1<row_select.size(); i1++)
          {
            g_Mat2_cov_new.col(i1) = g_Mat2_cov.col(row_select(i1));
          }
          g_Mat2_cov = g_Mat2_cov_new;
        }
        
      }
    }
    
    if( (x2correlation_global==1) && (x2diagonal_global==1) )
    {
      Eigen::MatrixXd g_Mat2_cov_new = g_Mat2_cov.block(0,0,1,1);
      g_Mat2_cov_new.setZero();
      g_Mat2_cov = g_Mat2_cov_new;
    }
    if( (x2correlation_global==1) && (x2diagonal_global==0) )
    {
      Eigen::MatrixXd temp_sel_matrix_inv = 1.0-temp_sel_matrix.array();
      Eigen::VectorXi row_select_inv = TVBS_logic2position_cpp(temp_sel_matrix_inv);
      
      // sub-matrix by not-row
      {
        Eigen::MatrixXd g_Mat2_cov_new(row_select_inv.size(), g_Mat2_cov.cols());
        g_Mat2_cov_new.setZero();
        for(int i1=0; i1<row_select_inv.size(); i1++)
        {
          g_Mat2_cov_new.row(i1) = g_Mat2_cov.row(row_select_inv(i1));
        }
        g_Mat2_cov = g_Mat2_cov_new;
      }
      
    }
    
  }
  struct TVBS_Matrix_two output = {g_Mat1_cov, g_Mat2_cov};
  
  return output;
}



struct TVBS_Matrix_two TVBS_g_cond_cov_cpp(Eigen::MatrixXd Id_mat, Eigen::MatrixXd Sigma)
{
  // initialize output
  Eigen::MatrixXd g_g_y, g_g_x;
  
  int dim1 = Id_mat.rows();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  Eigen::MatrixXd Sigma_11 = Sigma.block(0,0,dimdiff,dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, dim2-dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_22 = Sigma.block(dimdiff, dimdiff, dim2-dimdiff, dim2-dimdiff);
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  Eigen::MatrixXd Sigma_22_condcov = Sigma_22 - Sigma_12*Sigma_11_inv*Sigma_12.transpose();
  
  // update global variables
  x2symmetric_global     = 1;
  x2diagonal_global      = 0;
  x1symmetric_global     = 1;
  x1diagonal_global      = 1; 
  x2correlation_global   = 0;
  
  struct TVBS_Matrix_two output = TVBS_g_both_x_omega_x_cpp(Id_mat, Sigma_22_condcov);
  g_g_y = output.mat1;
  Eigen::MatrixXd g_g_2 = output.mat2;
  
  // update global variables
  x2symmetric_global     = 1;
  x2diagonal_global      = 0;
  x1symmetric_global     = 0;
  x1diagonal_global      = 0;
  x2correlation_global   = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_12, Sigma_11_inv);
  Eigen::MatrixXd g_g_22 = output.mat1;
  Eigen::MatrixXd g_g_23 = output.mat2;
  g_g_22 = -1.0*g_g_22*g_g_2;
  
  
  Eigen::MatrixXd indic_1(dim1*dimdiff,1);
  indic_1 = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
  indic_1.resize(dimdiff,dim1);
  Eigen::MatrixXd indic_1_t = indic_1.transpose();
  indic_1_t.resize(dim1*dimdiff,1);
  
  Eigen::MatrixXd g_g_22_new(indic_1_t.rows(), g_g_22.cols());
  for(int i1 = 0; i1<indic_1_t.rows(); i1++)
  {
    g_g_22_new.row(i1) = g_g_22.row(indic_1_t(i1));
  }
  
  
  if(condcov_global==1)
  {
    // update global variables
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
    
    Eigen::MatrixXd mat_n_seq(ddcov,1);
    mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    int diag_inc = 1;
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq, diag_inc);
    
    Eigen::MatrixXd indic_2_2_temp = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcov, dimdiff*dim1+ddcov-1);
    indic_2_2_temp.resize(dim1,dimdiff);
    Eigen::MatrixXd indic_2_2 = indic_2_2_temp.transpose();
    
    Eigen::MatrixXd indic_2(dimdiff, dim1+dimdiff);
    indic_2.block(0, 0, dimdiff, dimdiff) = mat_d_up;
    indic_2.block(0, dimdiff, dimdiff, dim1) = indic_2_2;
    
    mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, ddcov+dimdiff*dim1+dd1cov-1);
    
    mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq, diag_inc);
    
    
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up.cols());
    indic_3.setZero();
    indic_3.block(0,dimdiff,dim1,mat_d_up.cols()) = mat_d_up;
    
    Eigen::MatrixXd indic_4(dimdiff+dim1, dim1+dimdiff);
    indic_4.block(0, 0, dimdiff, dim1+dimdiff) = indic_2;
    indic_4.block(dimdiff, 0, dim1, dim1+dimdiff) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(), diag_inc);
    
    g_g_23 = -1.0*TVBS_g_inverse_cpp(Sigma_11)*g_g_23*g_g_2;
    Eigen::MatrixXd g_g_x_n1(g_g_23.rows()+g_g_22_new.rows()+g_g_2.rows(), g_g_23.cols());
    g_g_x_n1.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
    g_g_x_n1.block(g_g_23.rows(), 0, g_g_22_new.rows(), g_g_23.cols()) = g_g_22_new;
    g_g_x_n1.block(g_g_23.rows()+g_g_22_new.rows(), 0, g_g_2.rows(), g_g_23.cols()) = g_g_2;
    
    Eigen::MatrixXd g_g_x_n2(g_g_23.rows()+g_g_22_new.rows()+g_g_2.rows(), g_g_23.cols());
    for(int i_row = 0; i_row < indic.size(); i_row++)
    {
      g_g_x_n2.row(i_row) = g_g_x_n1.row(indic(i_row));
    }
    g_g_x = g_g_x_n2;
    
  }
  
  // INDEX MATRITZEN MUESSEN VORHER DEFINIERT WERDEN UND DANN IN GROESSE UND INHALT ANGEPASST WERDEN!!!
  
  
  if(condcov_global==0)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
    
    // define indic matrices
    Eigen::MatrixXd indic_2, indic_3;
    
    if(ddcor == 0)
    {
      Eigen::MatrixXd indic_2_rhs = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_rhs.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, dim1+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, dimdiff, dim1) = indic_2_rhs.transpose();
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      //int indic_2_lhs_dim = (1+sqrt(1+8*ddcor))/2;
      Eigen::MatrixXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor, 0, ddcor-1);
      Eigen::MatrixXd indic_2_lhs = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_2_rhs = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_rhs.resize(dim1, dimdiff);
      
      Eigen::MatrixXd indic_2_int(dimdiff, indic_2_lhs.cols()+dim1);
      indic_2_int.block(0, 0, indic_2_lhs.rows(), indic_2_lhs.cols()) = indic_2_lhs;
      indic_2_int.block(0, indic_2_lhs.cols(), dimdiff, dim1) = indic_2_rhs.transpose();
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1,dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      //int indic_2_lhs_dim = (1+sqrt(1+8*ddcor))/2;
      Eigen::MatrixXd mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+dd1cor-1);
      Eigen::MatrixXd indic_3_rhs = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+indic_3_rhs.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, dim1, indic_3_rhs.cols()) = indic_3_rhs;
      indic_3 = indic_3_int;
      
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0,0,indic_2.rows(),indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(),0,indic_3.rows(),indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(), 0);
    
    int diag_inc = 1;
    Eigen::VectorXd g_g_2_new_indic_logic = 1.0-TVBS_lower_tri_entries_cpp(Eigen::MatrixXd::Identity(dim1, dim1), diag_inc).array();
    Eigen::VectorXi g_g_2_new_indic = TVBS_logic2position_cpp(g_g_2_new_indic_logic);
    
    Eigen::MatrixXd g_g_2_new(g_g_2_new_indic.size(), g_g_2.cols());
    for(int i_rows=0; i_rows<g_g_2_new_indic.size(); i_rows++ )
    {
      g_g_2_new.row(i_rows) = g_g_2.row(g_g_2_new_indic(i_rows));
    }
    
    g_g_23 = -1*TVBS_g_inverse_cpp(Sigma_11)*g_g_23*g_g_2;
    
    if((ddcor==0)&&(dd1cor==0))
    {
      g_g_x = g_g_22_new;
    }
    if((ddcor==0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_22_new.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_22_new.rows(), g_g_22_new.cols()) = g_g_22_new;
      g_g_x_new.block(g_g_22_new.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      g_g_x = g_g_x_new;
    }
    if((ddcor!=0)&&(dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_23.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
      g_g_x_new.block(g_g_23.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      
      Eigen::MatrixXd g_g_x_new_select(indic.size(), g_g_x_new.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_new_select.row(i_rows) = g_g_x_new.row(indic(i_rows));
      }
      g_g_x = g_g_x_new_select;
    }
    if((ddcor!=0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_new(g_g_23.rows()+g_g_22_new.rows()+g_g_2_new.rows(), g_g_2_new.cols());
      g_g_x_new.block(0, 0, g_g_23.rows(), g_g_23.cols()) = g_g_23;
      g_g_x_new.block(g_g_23.rows(), 0, g_g_22_new.rows(), g_g_2_new.cols()) = g_g_22_new;
      g_g_x_new.block(g_g_23.rows()+g_g_22_new.rows(), 0, g_g_2_new.rows(), g_g_2_new.cols()) = g_g_2_new;
      
      Eigen::MatrixXd g_g_x_new_select(indic.size(), g_g_x_new.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_new_select.row(i_rows) = g_g_x_new.row(indic(i_rows));
      }
      g_g_x = g_g_x_new_select;
    }
  }
  
  if((cholesky_global==1)&&(condcov_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1)&&(condcov_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  
  struct TVBS_Matrix_two output_final = {g_g_y, g_g_x};
  return output_final;
  
}





struct TVBS_Matrix_three TVBS_g_cond_cov_trunc_cpp(Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x)
{
  // initiate outputs
  Eigen::MatrixXd g_g_x;
  
  int dim1 = Sigma.rows()-x.size();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  Eigen::MatrixXd Sigma_11 = Sigma.block(0, 0, dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, Sigma.rows()-dimdiff, dimdiff);
  
  Eigen::VectorXd mu_trunc;
  Eigen::MatrixXd sigma_trunc, g_mu_trunc, g_sigma_trunc;
  if(dimdiff==1)
  {
    // THIS CASE DOES NOT OCCUR WITH TVBS
  }
  if(dimdiff==2)
  {
    double tol = 0.000001;
    struct TVBS_Matrix_two output = TVBS_truncate_bi_normal_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x,tol);
    mu_trunc = output.mat1;
    sigma_trunc = output.mat2;
    
    output = TVBS_truncate_bi_normal_gradient_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    g_mu_trunc = output.mat1;
    g_sigma_trunc = output.mat2;
  }
  
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  Eigen::MatrixXd b = Sigma_11_inv*sigma_trunc*Sigma_11_inv;
  
  // update some global variables
  cholesky_global = 0;
  if(condcovsigtrunc_global==0)
  {
    condcov_global = 0;
    x2correlation_global = 1;
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
  }
  if(condcovsigtrunc_global!=0)
  {
    condcov_global = 1;
    x2correlation_global = 0;
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
  }
  
  struct TVBS_Matrix_two output = TVBS_g_cond_cov_cpp(Eigen::MatrixXd::Identity(dim1, dim1), Sigma);
  Eigen::MatrixXd g_y = output.mat1;
  Eigen::MatrixXd g_d_1 = output.mat2;
  
  // update global variables
  x2symmetric_global = 1;
  x2diagonal_global = 0;
  x1symmetric_global = 1;
  x1diagonal_global = 0;
  x2correlation_global = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_11_inv, sigma_trunc);
  Eigen::MatrixXd g_Sigma_11_inv = output.mat1;
  Eigen::MatrixXd g_b_omega = output.mat2;
  
  Eigen::MatrixXd g_bpsi_11 = TVBS_g_inverse_cpp(Sigma_11)*g_Sigma_11_inv;
  
  // update global variables
  x2symmetric_global = 1;
  x2diagonal_global = 0;
  x1symmetric_global = 0;
  x1diagonal_global = 0;
  x2correlation_global = 0;
  
  output = TVBS_g_both_x_omega_x_cpp(Sigma_12, b);
  Eigen::MatrixXd g_psi_12 = output.mat1;
  Eigen::MatrixXd g_b = output.mat2;
  Eigen::MatrixXd g_psi_11 = g_bpsi_11*g_b;
  Eigen::MatrixXd g_omega = g_b_omega*g_b;
  Eigen::MatrixXd g_sigma_tilde_new = g_sigma_trunc*g_omega;
  
  Eigen::MatrixXd g_mu(dimdiff+dim1, g_sigma_tilde_new.cols());
  g_mu.setZero();
  g_mu.block(0, 0, dimdiff, g_sigma_tilde_new.cols()) = g_sigma_tilde_new.block(0, 0, dimdiff, g_sigma_tilde_new.cols());
  
  Eigen::MatrixXd g_c = g_sigma_tilde_new.block(g_sigma_tilde_new.rows()-dimdiff, 0, dimdiff, g_sigma_tilde_new.cols());
  
  Eigen::MatrixXd indic_1_r(dim1*dimdiff,1);
  indic_1_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
  indic_1_r.resize(dimdiff, dim1);
  Eigen::MatrixXd indic_1 = indic_1_r.transpose();
  indic_1.resize(dim1*dimdiff, 1);
  
  {
    Eigen::MatrixXd g_psi_12_foo(indic_1.size(), g_psi_12.cols());
    for(int i_rows=0; i_rows<indic_1.size(); i_rows++)
    {
      g_psi_12_foo.row(i_rows)=g_psi_12.row(indic_1(i_rows));
    }
    g_psi_12 = g_psi_12_foo;
  }
  
  Eigen::MatrixXd g_psi_22(g_d_1.rows(), g_d_1.cols());
  g_psi_22.setZero();
  
  
  if(condcov_global==1)
  {
    g_psi_11 += g_sigma_tilde_new.block(dimdiff, 0, g_sigma_tilde_new.rows()-dimdiff-dimdiff, g_sigma_tilde_new.cols());
    
    Eigen::MatrixXd mat_n_seq(ddcov,1);
    mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 1);
    
    Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
    indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcov, ddcov+dim1*dimdiff-1);
    indic_2_2_r.resize(dim1, dimdiff);
    Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
    
    Eigen::MatrixXd indic_2(mat_d_up.rows(), mat_d_up.cols()+indic_2_2.cols());
    indic_2.block(0, 0, mat_d_up.rows(), mat_d_up.cols()) = mat_d_up;
    indic_2.block(0, mat_d_up.cols(), mat_d_up.rows(), indic_2_2.cols()) = indic_2_2;
    
    
    Eigen::MatrixXd mat_n_seq_2(dd1cov,1);
    mat_n_seq_2 = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, dd1cov+ddcov+dimdiff*dim1-1);
    Eigen::MatrixXd mat_d_up_2 = TVBS_vec_2_upper_diag_cpp(mat_n_seq_2,1);
    
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up_2.cols());
    indic_3.setZero();
    indic_3.block(0, dimdiff, dim1, mat_d_up_2.cols()) = mat_d_up_2;
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(), 1);
    
    Eigen::MatrixXd g_g_x_1(g_psi_11.rows()+g_psi_12.rows()+g_psi_22.rows(), g_psi_22.cols());
    g_g_x_1.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
    g_g_x_1.block(g_psi_11.rows(), 0, g_psi_12.rows(), g_psi_11.cols()) = g_psi_12;
    g_g_x_1.block(g_psi_11.rows()+g_psi_12.rows(), 0, g_psi_22.rows(), g_psi_11.cols()) = g_psi_22;
    Eigen::MatrixXd g_g_x_2(indic.size(), g_g_x_1.cols());
    for(int i_rows=0; i_rows<indic.size(); i_rows++)
    {
      g_g_x_2.row(i_rows) = g_g_x_1.row(indic(i_rows));
    }
    
    g_g_x = g_d_1 + g_g_x_2;
  }
  if(condcov_global==0)
  {
    Eigen::MatrixXd indic_2, indic_3;
    
    if(dimdiff==2)
    {
      g_psi_11 += g_sigma_tilde_new.row(dimdiff+1);
    }
    
    if(ddcor==0)
    {
      Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
      indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_2_r.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
      
      Eigen::MatrixXd indic_2_int(indic_2_2.rows(), indic_2_2.cols()+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, indic_2_2.rows(), indic_2_2.cols()) = indic_2_2;
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      Eigen::MatrixXd mat_n_seq(ddcor,1);
      mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor,0,ddcor-1);
      Eigen::MatrixXd mat_nd_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq,0);
      
      Eigen::MatrixXd indic_2_2_r(dimdiff*dim1,1);
      indic_2_2_r = Eigen::VectorXd::LinSpaced(dim1*dimdiff, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_2_r.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_2 = indic_2_2_r.transpose();
      
      Eigen::MatrixXd indic_2_int(indic_2_2.rows(), mat_nd_up.cols()+indic_2_2.cols());
      indic_2_int.block(0, 0, mat_nd_up.rows(), mat_nd_up.cols()) = mat_nd_up;
      indic_2_int.block(0, mat_nd_up.cols(), indic_2_2.rows(), indic_2_2.cols()) = indic_2_2;
      
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      Eigen::MatrixXd mat_n_seq(dd1cor,1);
      mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+dd1cor-1);
      Eigen::MatrixXd mat_nd_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+mat_nd_up.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, mat_nd_up.rows(), mat_nd_up.cols()) = mat_nd_up;
      indic_3 = indic_3_int;
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(),0);
    
    if((ddcor==0)&&(dd1cor==0))
    {
      g_g_x = g_d_1 + g_psi_12;
    }
    if((ddcor==0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_12.rows()+dd1cor, g_psi_12.cols());
      g_g_x_rhs.setZero();
      g_g_x_rhs.block(0, 0, g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      g_g_x = g_d_1 + g_g_x_rhs;
    }
    if((ddcor!=0)&&(dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_11.rows()+g_psi_12.rows(), g_psi_11.cols());
      g_g_x_rhs.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
      g_g_x_rhs.block(0, g_psi_11.rows(), g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      Eigen::MatrixXd g_g_x_rhs_ordered(indic.size(), g_g_x_rhs.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_rhs_ordered.row(i_rows) = g_g_x_rhs.row(indic(i_rows));
      }
      g_g_x = g_d_1+g_g_x_rhs_ordered;
    }
    if((ddcor!=0)&&(dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_rhs(g_psi_11.rows()+g_psi_12.rows()+g_psi_22.rows(), g_psi_11.cols());
      g_g_x_rhs.setZero();
      g_g_x_rhs.block(0, 0, g_psi_11.rows(), g_psi_11.cols()) = g_psi_11;
      g_g_x_rhs.block(g_psi_11.rows(), 0, g_psi_12.rows(), g_psi_12.cols()) = g_psi_12;
      g_g_x_rhs.block(g_psi_11.rows()+g_psi_12.rows(), 0, g_psi_22.rows(), g_psi_22.cols()) = g_psi_22;
      Eigen::MatrixXd g_g_x_rhs_ordered(indic.size(), g_g_x_rhs.cols());
      g_g_x_rhs_ordered.setZero();
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_rhs_ordered.row(i_rows) = g_g_x_rhs.row(indic(i_rows));
      }
      g_g_x = g_d_1+g_g_x_rhs_ordered;
      
    }
  }
  
  if((cholesky_global==1)&&(condcov_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1)&&(condcov_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_g_x, g_c};
  return output_final;
}




struct TVBS_Matrix_four TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd Id_mat, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, Eigen::VectorXd x)
{
  // initialize output
  Eigen::MatrixXd g_g_x, g_c;
  
  // extract dimensions of the inputs
  int dim1 = Id_mat.rows();
  int dim2 = Sigma.rows();
  int dimdiff = dim2-dim1;
  
  // for correct sequencing of the entries in Sigma
  int ddcov = dimdiff*(dimdiff+1)/2;
  int ddcor = dimdiff*(dimdiff-1)/2;
  int dd1cov = dim1*(dim1+1)/2;
  int dd1cor = dim1*(dim1-1)/2;
  
  // extract submatrices of the Covariance matrix
  Eigen::MatrixXd Sigma_11 = Sigma.block(0,0,dimdiff,dimdiff);
  Eigen::MatrixXd Sigma_12 = Sigma.block(dimdiff, 0, Sigma.rows()-dimdiff, dimdiff);
  Eigen::MatrixXd Sigma_11_inv = Sigma_11.inverse();
  
  
  Eigen::MatrixXd g_mu_tilde = (Id_mat*Sigma_12*Sigma_11_inv).transpose();
  Eigen::MatrixXd g_mu(g_mu_tilde.rows()+dim1, dim1);
  g_mu.block(0, 0, g_mu_tilde.rows(), g_mu_tilde.cols()) = -1.0*g_mu_tilde;
  g_mu.block(g_mu_tilde.rows(), 0, dim1, dim1) = Eigen::MatrixXd::Identity(dim1, dim1);
  
  
  // truncate
  Eigen::VectorXd mu_trunc;
  Eigen::MatrixXd sigma_trunc, g_mu_trunc, g_sigma_trunc;
  if(dimdiff==1)
  {
    // not relevant for TVBS 
  }
  if(dimdiff==2)
  {double tol = 0.000001;
    struct TVBS_Matrix_two output = TVBS_truncate_bi_normal_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x, tol);
    mu_trunc = output.mat1;
    sigma_trunc = output.mat2;
    output = TVBS_truncate_bi_normal_gradient_cpp(mu.head<2>(), Sigma.topLeftCorner<2,2>(), x);
    g_mu_trunc = output.mat1;
    g_sigma_trunc = output.mat2;
  }
  
  Eigen::MatrixXd g_y_1 = Sigma_12*Sigma_11_inv*(mu_trunc-mu.head(dimdiff));
  Eigen::MatrixXd g_g_y(g_y_1.rows(), g_y_1.rows());
  g_g_y.setZero();
  g_g_y.diagonal() = g_y_1;
  
  
  Eigen::MatrixXd g_mu_tilde_new = g_mu_trunc*g_mu_tilde;
  g_mu.block(0, 0, dimdiff, g_mu.cols()) += g_mu_tilde_new.block(0, 0, dimdiff, g_mu_tilde_new.cols());
  
  Eigen::MatrixXd g_x_12 = kroneckerProduct(Id_mat, Sigma_11_inv*(mu_trunc-mu.head(dimdiff)));
  
  // reordering indices for matrix transposion
  {
    Eigen::MatrixXd indic_1_t = Eigen::VectorXd::LinSpaced(dim1*dimdiff,0,dim1*dimdiff-1);
    indic_1_t.resize(dimdiff, dim1);
    Eigen::MatrixXd indic_1 = indic_1_t.transpose();
    indic_1.resize(dim1*dimdiff,1);
    {
      Eigen::MatrixXd g_x_12_ordered(indic_1.rows(), g_x_12.cols());
      for(int i_rows = 0; i_rows<indic_1.rows(); i_rows++)
      {
        g_x_12_ordered.row(i_rows) = g_x_12.row(indic_1(i_rows));
      }
      g_x_12 = g_x_12_ordered;
    }
  }
  
  // set global variables
  omsymmetric_global = 1;
  omdiagonal_global = 0;
  
  if(condcovmeantrunc_global==1)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 0;
    
    Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcov,0,ddcov-1);
    Eigen::MatrixXd mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq,1);
    
    Eigen::MatrixXd indic_2_1_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcov, dimdiff*dim1+ddcov-1);
    indic_2_1_t.resize(dim1, dimdiff);
    Eigen::MatrixXd indic_2_1 = indic_2_1_t.transpose();
    
    Eigen::MatrixXd indic_2(mat_d_up.rows(), mat_d_up.cols()+indic_2_1.cols() );
    indic_2.block(0, 0, mat_d_up.rows(), mat_d_up.cols()) = mat_d_up;
    indic_2.block(0, mat_d_up.cols(), indic_2_1.rows(), indic_2_1.cols()) = indic_2_1;
    
    mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cov, ddcov+dimdiff*dim1, ddcov+dimdiff*dim1+dd1cov-1);
    mat_d_up = TVBS_vec_2_upper_diag_cpp(mat_n_seq,1);
    Eigen::MatrixXd indic_3(dim1, dimdiff+mat_d_up.cols());
    indic_3.setZero();
    indic_3.block(0, dimdiff, dim1, mat_d_up.cols()) = mat_d_up;
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_2.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(),1 );
    
    Eigen::MatrixXd g_x_11 = TVBS_g_inverse_cpp(Sigma_11)*TVBS_g_a_omega_b_cpp(Id_mat*Sigma_12, mu_trunc-mu.head(dimdiff));
    g_x_11 += g_mu_tilde_new.block(dimdiff, 0, (dimdiff+dimdiff*(dimdiff+1)/2)-(dimdiff+1)+1, g_mu_tilde_new.cols());
    
    Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows()+dd1cov, dim1);
    g_g_x_int.setZero();
    g_g_x_int.block(0, 0, g_x_11.rows(), dim1) = g_x_11;
    g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), dim1) = g_x_12;
    {
      Eigen::MatrixXd g_g_x_ordered(indic.size(), dim1);
      for(int i_rows=0; i_rows<g_g_x_int.rows(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
  }
  
  if(condcovmeantrunc_global==0)
  {
    xinvsymmetric_global = 1;
    xinvdiagonal_global = 0;
    xinvcorrelation_global = 1;
    
    Eigen::MatrixXd g_x_11 = TVBS_g_inverse_cpp(Sigma_11)*TVBS_g_a_omega_b_cpp(Id_mat*Sigma_12, mu_trunc-mu.head(dimdiff));
    
    if(dimdiff==2)
    {
      g_x_11 += g_mu_tilde_new.row(dimdiff+1);
    }
    
    // for reordering purposes:
    Eigen::MatrixXd indic_2, indic_3;
    if(ddcor==0)
    {
      Eigen::MatrixXd indic_2_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, ddcor+dim1*dimdiff-1);
      indic_2_t.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, dim1+1);
      indic_2_int.setZero();
      indic_2_int.block(0, 1, dimdiff, dim1) = indic_2_t.transpose();
      indic_2 = indic_2_int;
    }
    if(ddcor!=0)
    {
      Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(ddcor, 0, ddcor-1);
      Eigen::MatrixXd mat_nd_up_diag_zero = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      Eigen::MatrixXd indic_2_2_t = Eigen::VectorXd::LinSpaced(dimdiff*dim1, ddcor, dimdiff*dim1+ddcor-1);
      indic_2_2_t.resize(dim1, dimdiff);
      Eigen::MatrixXd indic_2_int(dimdiff, mat_nd_up_diag_zero.cols()+dim1);
      indic_2_int.block(0, 0, dimdiff, mat_nd_up_diag_zero.cols()) = mat_nd_up_diag_zero;
      indic_2_int.block(0, mat_nd_up_diag_zero.cols(), dimdiff, dim1) = indic_2_2_t.transpose();
      indic_2 = indic_2_int;
    }
    
    if(dd1cor==0)
    {
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+1);
      indic_3_int.setZero();
      indic_3 = indic_3_int;
    }
    if(dd1cor!=0)
    {
      Eigen::VectorXd mat_n_seq = Eigen::VectorXd::LinSpaced(dd1cor, ddcor+dimdiff*dim1, ddcor+dimdiff*dim1+ddcor-1);
      Eigen::MatrixXd mat_nd_up_diag_zero = TVBS_vec_2_upper_diag_cpp(mat_n_seq, 0);
      Eigen::MatrixXd indic_3_int(dim1, dimdiff+mat_nd_up_diag_zero.cols());
      indic_3_int.setZero();
      indic_3_int.block(0, dimdiff, mat_nd_up_diag_zero.rows(), mat_nd_up_diag_zero.cols()) = mat_nd_up_diag_zero;
      indic_3 = indic_3_int;
    }
    
    Eigen::MatrixXd indic_4(indic_2.rows()+indic_3.rows(), indic_3.cols());
    indic_4.block(0, 0, indic_2.rows(), indic_2.cols()) = indic_2;
    indic_4.block(indic_2.rows(), 0, indic_3.rows(), indic_3.cols()) = indic_3;
    
    Eigen::VectorXd indic = TVBS_lower_tri_entries_cpp(indic_4.transpose(),0);
    
    
    if((ddcor==0) && (dd1cor==0))
    {
      g_g_x = g_x_12;
    }
    if((ddcor==0) && (dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_12.rows()+dd1cor, dim1);
      g_g_x_int.setZero();
      g_g_x_int.block(0, 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      g_g_x = g_g_x_int;
    }
    if((ddcor!=0) && (dd1cor==0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows(), g_x_12.cols());
      g_g_x_int.block(0, 0, g_x_11.rows(), g_x_11.cols()) = g_x_11;
      g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      Eigen::MatrixXd g_g_x_ordered(indic.size(), g_x_12.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
    if((ddcor!=0) && (dd1cor!=0))
    {
      Eigen::MatrixXd g_g_x_int(g_x_11.rows()+g_x_12.rows()+dd1cor, g_x_12.cols());
      g_g_x_int.setZero();
      g_g_x_int.block(0, 0, g_x_11.rows(), g_x_11.cols()) = g_x_11;
      g_g_x_int.block(g_x_11.rows(), 0, g_x_12.rows(), g_x_12.cols()) = g_x_12;
      Eigen::MatrixXd g_g_x_ordered(indic.size(), g_x_12.cols());
      for(int i_rows=0; i_rows<indic.size(); i_rows++)
      {
        g_g_x_ordered.row(i_rows) = g_g_x_int.row(indic(i_rows));
      }
      g_g_x = g_g_x_ordered;
    }
    
  }
  
  if((cholesky_global==1) && (condcovmeantrunc_global==1))
  {
    g_g_x = TVBS_g_cholesky_cov_cpp(Sigma)*g_g_x;
  }
  if((cholesky_global==1) && (condcovmeantrunc_global==0))
  {
    g_g_x = TVBS_g_cholesky_cor_cpp(Sigma)*g_g_x;
  }
  g_c = g_mu_tilde_new.block(dimdiff+(dimdiff*(dimdiff+1))/2, 0, g_mu_tilde_new.rows()-(dimdiff+dimdiff*((dimdiff+1))/2), g_mu_tilde_new.cols());
  
  
  struct TVBS_Matrix_four output_final = {g_g_y, g_mu, g_g_x, g_c};
  return output_final;
}












// // // // // // // // // // // // // // // // // // // // // // // // // // // // //
//                                                                                  //
//    Functions to evaluate bi, tri- and quadro-variate NCDF and the gradients      //
//                                                                                  //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // //


// // // // // // // // //
//  bivariate           //
// // // // // // // // //

struct TVBS_double_three TVBS_grad_cdf_bvn_cpp(Eigen::VectorXd x_temp, Eigen::MatrixXd Rho)
{
  double rhotilde, tr1, tr2, pdf2, g_w1, g_w2;
  
  
  rhotilde = sqrt(1.0-pow(Rho(0,1), 2));
  
  //Rcout << "rhotilde" << rhotilde << "x_t" << x_temp << std::endl; 
  
  tr1 = (x_temp(1) - Rho(0,1)*x_temp(0))/rhotilde;
  tr2 = (x_temp(0) - Rho(0,1)*x_temp(1))/rhotilde;
  pdf2 = std_normal_pdf(x_temp(0))*std_normal_pdf(tr1)/rhotilde;
  g_w1 = std_normal_pdf(x_temp(0))*std_normal_cdf(tr1);
  g_w2 = std_normal_pdf(x_temp(1))*std_normal_cdf(tr2);
  
  struct TVBS_double_three output_final = {g_w1, g_w2, pdf2};
  
  return output_final;
}


struct TVBS_double_three TVBS_grad_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd w, Eigen::MatrixXd Rho)
{
  double bivar_cdf, univar_cdf, g_w1, g_w2, grho, gw1new;
  Eigen::VectorXd mu(2);
  mu.setZero();
  
  bivar_cdf = TVBS_pmvnorm_cpp(w, mu, Rho);
  univar_cdf = std_normal_cdf(w(0));
  
  if (univar_cdf<tol) { univar_cdf = tol; }
  
  struct TVBS_double_three output = TVBS_grad_cdf_bvn_cpp(w, Rho);
  g_w1 = output.d1;
  g_w2 = output.d2;
  grho = output.d3;
  
  gw1new = (univar_cdf*g_w1 - bivar_cdf*std_normal_pdf(w(0)))/pow(univar_cdf, 2);
  
  struct TVBS_double_three output_final = {gw1new, g_w2/univar_cdf, grho/univar_cdf};
  return output_final;
}




struct TVBS_Matrix_three TVBS_grad_non_cdf_bvn_by_cdfn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  Eigen::VectorXd Sigma_diag_sqrt, x_norm, g_w(2);
  Eigen::MatrixXd Cor_mat, g_b_corcov, g_omega_corcov, g_cov, g_x, g_mu;
  double rho, g_w1, g_w2, grho;
  
  Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  
  x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  
  Eigen::MatrixXd Sigma_diag_sqrt_inv(2,2);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_by_cdfn_cpp(x_norm.head(2), Cor_mat);
  g_w1 = output_3.d1;
  g_w2 = output_3.d2;
  grho = output_3.d3;
  g_w << g_w1, g_w2;
  
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  g_b_corcov = output_2.mat1;
  g_omega_corcov = output_2.mat2;
  
  g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  g_cov = g_b_corcov*g_w + g_omega_corcov*grho;
  g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  return output_final;
}






struct TVBS_Matrix_two TVBS_grad_cdf_tvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr)
{
  // initialize variables
  int enteredloop1, enteredloop2;
  double epst, pt, h1, h2, h3, r12, r13, r23, tvn;
  Eigen::VectorXd gg(6), d_tvn(6);
  
  // safeguard threshold against numerical errors 
  epst = pow(10, -7);
  
  pt = pi_global/2;
  h1 = x_norm(0);
  h2 = x_norm(1);
  h3 = x_norm(2);
  r12 = Corr(0,1);
  r13 = Corr(0,2);
  r23 = Corr(1,2);
  
  // initialize vector to store gradient output
  gg.setZero();
  
  // following are safeguards against numerical errors (implemented by a reordering)
  enteredloop1 = 0;     // control variables to check which safeguard was used
  enteredloop2 = 0;
  
  // check if one of the correlations is bigger than the other
  // (kind of optimal ordering)
  if( (abs(r12)-abs(r13)) > epst)
  {
    // 1
    h2 = h3;
    h3 = x_norm(1);
    r12 = r13;
    r13 = Corr(0,1);
    enteredloop1 = 1;
  }
  if( (abs(r13)-abs(r23)) > epst)
  {
    // 2
    h1 = h2;
    h2 = x_norm(0);
    r23 = r13;
    r13 = Corr(1,2);
    enteredloop2 = 1;
  }
  
  // take care of some cases with small values seperately
  if( (abs(h1) + abs(h2) + abs(h3)) < epst )
  {
    // 3
    gg(3) = ( 2* 1.0/sqrt(1.0-pow(r12,2))/pi_global )/8;
    gg(4) = ( 2* 1.0/sqrt(1.0-pow(r13,2))/pi_global )/8;
    gg(5) = ( 2* 1.0/sqrt(1.0-pow(r23,2))/pi_global )/8;
    d_tvn = gg;
  } 
  else if( (abs(r12) + abs(r13))<epst )
  {
    // 4
    double tempval1 = std_normal_cdf(h1);
    
    Eigen::VectorXd h23(2), mu_h23(2);
    mu_h23.setZero();
    h23 << h2, h3;
    Eigen::MatrixXd corr_h23(2,2);
    corr_h23.setOnes();
    corr_h23(0,1) = r23;
    corr_h23(1,0) = r23;
    double tempval2 = TVBS_pmvnorm_cpp(h23, mu_h23, corr_h23);
    
    gg(0) = std_normal_pdf(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h23, corr_h23);
    gg(1) = tempval1*temp_bvn_grad.d1;
    gg(2) = tempval1*temp_bvn_grad.d2;
    gg(5) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r13) + abs(r23))<epst )
  {
    // 5
    double tempval1 = std_normal_cdf(h3);
    
    Eigen::VectorXd h12(2), mu_h12(2);
    mu_h12.setZero();
    h12 << h1, h2;
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    double tempval2 = TVBS_pmvnorm_cpp(h12, mu_h12, corr_h12);
    
    gg(2) = std_normal_pdf(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h12, corr_h12);
    gg(0) = tempval1*temp_bvn_grad.d1;
    gg(1) = tempval1*temp_bvn_grad.d2;
    gg(3) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r12) + abs(r23))<epst )
  {
    // 6
    double tempval1 = std_normal_cdf(h2);
    
    Eigen::VectorXd h13(2), mu_h13(2);
    mu_h13.setZero();
    h13 << h1, h3;
    Eigen::MatrixXd corr_h13(2,2);
    corr_h13.setOnes();
    corr_h13(0,1) = r13;
    corr_h13(1,0) = r13;
    double tempval2 = TVBS_pmvnorm_cpp(h13, mu_h13, corr_h13);
    
    // tvn should have no influence on this function
    tvn = std_normal_cdf(h2)*tempval2;
    
    
    gg(1) = std_normal_pdf(h2)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h13, corr_h13);
    gg(0) = tempval1*temp_bvn_grad.d1;
    gg(2) = tempval1*temp_bvn_grad.d2;
    gg(4) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (1.0-r23)<epst )
  {
    // 7
    double h_23_min = std::min(h2,h3);
    
    Eigen::VectorXd h1_23_min(2), mu_h13(2);
    mu_h13.setZero();
    h1_23_min << h1, h_23_min;
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    
    // tvn should have no influence on this function
    tvn = TVBS_pmvnorm_cpp(h1_23_min, mu_h13, corr_h12);
    
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h1_23_min, corr_h12);
    
    gg(0) = temp_bvn_grad.d1;
    if(h2==h3)
    {
      gg(1) = temp_bvn_grad.d2;
      gg(2) = temp_bvn_grad.d2;
    }
    else if(h2<h3)
    {
      gg(1) = temp_bvn_grad.d2;
    }
    else
    {
      gg(2) = temp_bvn_grad.d2;
    }
    gg(3) = temp_bvn_grad.d3;
    d_tvn = gg;
    
  }
  else if( (r23+1.0<epst) && (h2>-1.0*h3) )
  {
    // 8
    Eigen::VectorXd h_12(2), h_1m3(2);
    h_12 << h1, h2;
    h_1m3 << h1, -1.0*h3;
    
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h_12, corr_h12);
    struct TVBS_double_three temp_bvn_grad2 = TVBS_grad_cdf_bvn_cpp(h_1m3, corr_h12);
    gg(0) = temp_bvn_grad.d1 - temp_bvn_grad2.d1;
    gg(1) = temp_bvn_grad.d2;
    gg(2) = temp_bvn_grad2.d2;
    gg(3) = temp_bvn_grad.d3 - temp_bvn_grad2.d3;
    d_tvn = gg;
    
  }
  else
  {
    // now the general case....
    Eigen::Vector2d h_23, mu_h23;
    mu_h23.setZero();
    h_23 << h2, h3;
    Eigen::Matrix2d corr_h23;
    corr_h23 << 1, r23, r23, 1;
    
    double cdfbvn_h2_h3 = TVBS_pmvnorm_cpp(h_23, mu_h23, corr_h23);
    
    double cdfn_h1 = std_normal_cdf(h1);
    tvn = cdfn_h1*cdfbvn_h2_h3;
    
    double d_a = std_normal_pdf(h2)*std_normal_cdf((h3-r23*h2)/sqrt(1.0-pow(r23,2)));
    double d_b = std_normal_pdf(h3)*std_normal_cdf((h2-r23*h3)/sqrt(1.0-pow(r23,2)));
    double d_corr = (exp(-0.5*(pow(h2,2) + pow(h3,2) - 2*r23 * h2 * h3   ) / (1.0-r23*r23) )) / sqrt(1.0-r23 * r23)/(2*pi_global);
    
    d_tvn << (std_normal_pdf(h1)*cdfbvn_h2_h3), cdfn_h1*d_a, (cdfn_h1*d_b), 0, 0, (cdfn_h1*d_corr);
    
    
    // bunch of auxillary variables
    Eigen::VectorXd wg(12), xg(12);
    wg << 0.1279381953467518, 0.1258374563468280, 0.1216704729278031, 0.1155056680537265, 0.1074442701159659, 0.09761865210411358, 0.08619016153195296, 0.07334648141108081, 0.05929858491543594, 0.04427743881742087, 0.02853138862893389, 0.01234122979998693;
    xg << 0.06405689286260559, 0.1911188674736164, 0.3150426796961635, 0.4337935076260450, 0.5454214713888396, 0.6480936519369754, 0.7401241915785546, 0.8200019859739028, 0.8864155270044012, 0.9382745520027329, 0.9747285559713096, 0.9951872199970214;
    double rua = asin(r12);
    double rub = asin(r13);
    double res = 0.0;
    double d_rua_r12 = 1.0/sqrt(1.0-pow(r12,2));
    double d_rub_r13 = 1.0/sqrt(1.0-pow(r13,2));
    double d_res_h1 = 0;
    double d_res_h2 = 0;
    double d_res_h3 = 0;
    double d_res_r12 = 0;
    double d_res_r13 = 0;
    double d_res_r23 = 0;
    
    // here comes the (numerical) magic
    for(int j_int = 0; j_int<12; j_int++)
    {
      double fc = 0;
      double d_fc_h1 = 0;
      double d_fc_h2 = 0;
      double d_fc_h3 = 0;
      double d_fc_r12 = 0;
      double d_fc_r13 = 0;
      double d_fc_r23 = 0;
      
      Eigen::VectorXd temp_sincs_grad = TVBS_sincs_grad_cpp(rua*(1.0-xg(j_int))/2);
      double r12t = temp_sincs_grad(2);
      double rr2 = temp_sincs_grad(3);
      double d_r12t_r12 = temp_sincs_grad(0)*( 1.0 - xg(j_int) )/2*d_rua_r12;
      double d_rr2_r12 = temp_sincs_grad(1)*( 1.0 - xg(j_int) )/2*d_rua_r12;
      
      temp_sincs_grad = TVBS_sincs_grad_cpp(rub*(1.0-xg(j_int))/2);
      double r13t = temp_sincs_grad(2);
      double rr3 = temp_sincs_grad(3);
      double d_r13t_r13 = temp_sincs_grad(0)*( 1.0 - xg(j_int) )/2*d_rub_r13;
      double d_rr3_r13 = temp_sincs_grad(1)*( 1.0 - xg(j_int) )/2*d_rub_r13;
      
      if( abs(rua)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h2, h3, r13t, r23, r12t, rr2 );
        fc += rua*temp_pntgnd_grad(7);
        d_fc_h1 += rua*temp_pntgnd_grad(0);
        d_fc_h2 += rua*temp_pntgnd_grad(1);
        d_fc_h3 += rua*temp_pntgnd_grad(2);
        d_fc_r12 += d_rua_r12*temp_pntgnd_grad(7) + rua*(temp_pntgnd_grad(5)*d_r12t_r12 + temp_pntgnd_grad(6)*d_rr2_r12 );
        d_fc_r13 += rua*(temp_pntgnd_grad(3)*d_r13t_r13);
        d_fc_r23 += rua*temp_pntgnd_grad(4);
        
      }
      if( abs(rub)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h3, h2, r12t, r23, r13t, rr3 );
        fc += rub*temp_pntgnd_grad(7);
        d_fc_h1 += rub*temp_pntgnd_grad(0);
        d_fc_h2 += rub*temp_pntgnd_grad(2);
        d_fc_h3 += rub*temp_pntgnd_grad(1);
        d_fc_r12 += rub*temp_pntgnd_grad(3)*d_r12t_r12;
        d_fc_r13 += d_rub_r13*temp_pntgnd_grad(7) + rub*(temp_pntgnd_grad(5)*d_r13t_r13 + temp_pntgnd_grad(6)*d_rr3_r13);
        d_fc_r23 += rub*temp_pntgnd_grad(4);
        
      }
      
      temp_sincs_grad = TVBS_sincs_grad_cpp( rua*( 1.0 + xg(j_int) )/2 );
      r12t = temp_sincs_grad(2);
      rr2 = temp_sincs_grad(3);
      d_r12t_r12 = temp_sincs_grad(0)*( 1.0 + xg(j_int) )/2*d_rua_r12;
      d_rr2_r12 = temp_sincs_grad(1)*( 1.0 + xg(j_int) )/2*d_rua_r12;
      
      temp_sincs_grad = TVBS_sincs_grad_cpp( rub*( 1.0 + xg(j_int) )/2 );
      r13t = temp_sincs_grad(2);
      rr3 = temp_sincs_grad(3);
      d_r13t_r13 = temp_sincs_grad(0)*( 1.0 + xg(j_int) )/2*d_rub_r13;
      d_rr3_r13 = temp_sincs_grad(1)*( 1.0 + xg(j_int) )/2*d_rub_r13;
      
      if( abs(rua)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h2, h3, r13t, r23, r12t, rr2 );
        fc += rua*temp_pntgnd_grad(7);
        d_fc_h1 += rua*temp_pntgnd_grad(0);
        d_fc_h2 += rua*temp_pntgnd_grad(1);
        d_fc_h3 += rua*temp_pntgnd_grad(2);
        d_fc_r12 += d_rua_r12*temp_pntgnd_grad(7) + rua*(temp_pntgnd_grad(5)*d_r12t_r12 + temp_pntgnd_grad(6)*d_rr2_r12 );
        d_fc_r13 += rua*(temp_pntgnd_grad(3)*d_r13t_r13);
        d_fc_r23 += rua*temp_pntgnd_grad(4);
      }
      if( abs(rub)>0 )
      {
        Eigen::VectorXd temp_pntgnd_grad = TVBS_pntgnd_grad_cpp( h1, h3, h2, r12t, r23, r13t, rr3 );
        fc += rub*temp_pntgnd_grad(7);
        d_fc_h1 += rub*temp_pntgnd_grad(0);
        d_fc_h2 += rub*temp_pntgnd_grad(2);
        d_fc_h3 += rub*temp_pntgnd_grad(1);
        d_fc_r12 += rub*temp_pntgnd_grad(3)*d_r12t_r12;
        d_fc_r13 += d_rub_r13*temp_pntgnd_grad(7) + rub*(temp_pntgnd_grad(5)*d_r13t_r13 + temp_pntgnd_grad(6)*d_rr3_r13);
        d_fc_r23 += rub*temp_pntgnd_grad(4);
      }
      
      res += wg(j_int)*fc;
      d_res_h1 += wg(j_int)*d_fc_h1;
      d_res_h2 += wg(j_int)*d_fc_h2;
      d_res_h3 += wg(j_int)*d_fc_h3;
      d_res_r12 += wg(j_int)*d_fc_r12;
      d_res_r13 += wg(j_int)*d_fc_r13;
      d_res_r23 += wg(j_int)*d_fc_r23;
    }
    
    tvn += res/(4*pi_global);
    Eigen::VectorXd d_tvn_add(6);
    d_tvn_add << d_res_h1, d_res_h2, d_res_h3, d_res_r12, d_res_r13, d_res_r23;
    d_tvn += d_tvn_add/(4*pi_global);
    
  }
  
  // reorder according to which case was present
  if(enteredloop2==1)
  {
    Eigen::VectorXd d_tvn_unsorted = d_tvn;
    d_tvn_unsorted(0) = d_tvn(1);
    d_tvn_unsorted(1) = d_tvn(0);
    d_tvn_unsorted(4) = d_tvn(5);
    d_tvn_unsorted(5) = d_tvn(4);
    d_tvn = d_tvn_unsorted;
  }
  if(enteredloop1==1)
  {
    Eigen::VectorXd d_tvn_unsorted = d_tvn;
    d_tvn_unsorted(1) = d_tvn(2);
    d_tvn_unsorted(2) = d_tvn(1);
    d_tvn_unsorted(3) = d_tvn(4);
    d_tvn_unsorted(4) = d_tvn(3);
    d_tvn = d_tvn_unsorted;
  }
  
  struct TVBS_Matrix_two output_final = {d_tvn.head(3), d_tvn.tail(3)};
  return output_final;
}












struct Hessian_trunc_struct
{
  Eigen::MatrixXd mat1, mat2, mat3, mat4, mat5;
};


//Hessian_trunc_struct 
struct Hessian_trunc_struct TVBS_Hessian_biv_norm_trunc(Eigen::VectorXd w, double rho)
{
  //Hessian_trunc Hess; // Hessian of lambda_{0,1}, Omega_{0,0},Omega_{0,1},Omega_{1,1} as a function of w_0, w_1, rho. (indexation according to CPP to match code) 
  // structure contains five matrices, one Hessian for each element. 
  Eigen::MatrixXd H_lambda0(3,3);
  H_lambda0.setZero();
  
  Eigen::MatrixXd H_lambda1(3,3);
  H_lambda1.setZero();
  
  Eigen::MatrixXd H_Omega00(3,3);
  H_Omega00.setZero();
  
  Eigen::MatrixXd H_Omega01(3,3);
  H_Omega01.setZero();
  
  Eigen::MatrixXd H_Omega11(3,3);
  H_Omega11.setZero();
  
  
  // calculate deltas 
  double cdf_0 = std_normal_cdf((w(1) - rho*w(0))/sqrt(1-rho*rho));
  double delta_0 = std_normal_pdf(w(0))*cdf_0;
  double cdf_1 = std_normal_cdf((w(0) - rho*w(1))/sqrt(1-rho*rho));
  double delta_1 = std_normal_pdf(w(1))*cdf_1;
  
  // now for derivatives 
  // gradients of delta_0, 
  Eigen::MatrixXd grad_delta_0(3,1);
  grad_delta_0.setZero();
  
  grad_delta_0(0,0) = -std_normal_pdf(w(0))*cdf_0*w(0) - std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_0(1,0) = std_normal_pdf(w(0))*std_normal_pdf((w(1) - rho*w(0))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_0(2,0) = std_normal_pdf(w(0))*normal_pdf((w(1) - rho*w(0)),(1-rho*rho))*(-w(0)+(w(1)-rho*w(0))*rho/(1-rho*rho));
  
  // gradients of delta_1. 
  Eigen::MatrixXd grad_delta_1(3,1);
  grad_delta_1.setZero();
  grad_delta_1(1,0) = -std_normal_pdf(w(1))*cdf_1*w(1) - std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))*rho/(sqrt(1-rho*rho));
  grad_delta_1(0,0) = std_normal_pdf(w(1))*std_normal_pdf((w(0) - rho*w(1))/sqrt(1-rho*rho))/(sqrt(1-rho*rho));
  grad_delta_1(2,0) = std_normal_pdf(w(1))*normal_pdf((w(0) - rho*w(1)),(1-rho*rho))*(-w(1)+(w(0)-rho*w(1))*rho/(1-rho*rho));  
  
  // Hessians 
  Eigen::MatrixXd H_delta0 = Hessian_delta(w(0),w(1),rho);
  Eigen::MatrixXd HH_delta1 = Hessian_delta(w(1),w(0),rho);
  Eigen::MatrixXd H_delta1 = HH_delta1;
  
  H_delta1(0,0)= HH_delta1(1,1);
  H_delta1(1,1)= HH_delta1(0,0);
  H_delta1(0,2) = HH_delta1(1,2);
  H_delta1(1,2) = HH_delta1(0,2);
  H_delta1(2,0) = HH_delta1(2,1);
  H_delta1(2,1) = HH_delta1(2,0);
  
  
  // now for the lambdas. 
  // calculate lambdas
  double Phi2 = biv_normal_cdf(w(0),w(1),rho);
  double lambda_0 = -(delta_0+rho*delta_1)/Phi2;
  double lambda_1 = -(delta_1+rho*delta_0)/Phi2;
  
  Eigen::VectorXd dPhi2 = grad_cdf(w(0),w(1),rho);
  Eigen::MatrixXd H_Phi2 = Hess_cdf(w(0),w(1),rho);
  
  Eigen::VectorXd drho(3);
  drho.setZero();
  drho(2)=1; 
  
  // Hessian for lambda_0 
  // three terms. 
  Eigen::VectorXd grad_lambda_0(3);
  grad_lambda_0.setZero();
  
  grad_lambda_0 = -(grad_delta_0 + rho*grad_delta_1)/Phi2 - drho*delta_1/Phi2 + (delta_0+rho*delta_1)*dPhi2/pow(Phi2,2);
  
  // first term
  H_lambda0 = -(H_delta0 + rho*H_delta1)/Phi2 + dPhi2*(grad_delta_0.transpose() + rho*grad_delta_1.transpose())/pow(Phi2,2); 
  
  H_lambda0 += -(drho*grad_delta_1.transpose())/Phi2;
  
  // second term.
  H_lambda0 += -(grad_delta_1* drho.transpose())/Phi2 + delta_1*(dPhi2 * drho.transpose())/pow(Phi2,2); 
  
  // third term 
  H_lambda0 +=  (grad_delta_0 + rho*grad_delta_1 + drho*delta_1) * dPhi2.transpose()/pow(Phi2,2);
  H_lambda0 +=  (delta_0+rho*delta_1) *(H_Phi2/pow(Phi2,2) - 2*dPhi2 * dPhi2.transpose()/pow(Phi2,3));
  
  
  
  // Hessian for lambda_1 
  // three terms. 
  Eigen::VectorXd grad_lambda_1(3);
  grad_lambda_1.setZero();
  
  grad_lambda_1 = -(grad_delta_1 + rho*grad_delta_0)/Phi2 - drho*delta_0/Phi2 + (delta_1+rho*delta_0)*dPhi2/pow(Phi2,2);
  
  
  // first term
  H_lambda1 = -(H_delta1 + rho*H_delta0)/Phi2 + dPhi2*(grad_delta_1.transpose() + rho*grad_delta_0.transpose())/pow(Phi2,2); 
  H_lambda1 += -(drho*grad_delta_0.transpose())/Phi2;
  
  // second term.
  H_lambda1 += -(grad_delta_0* drho.transpose())/Phi2 + delta_0*(dPhi2 * drho.transpose())/pow(Phi2,2); 
  
  // third term 
  H_lambda1 +=  (grad_delta_1 + rho*grad_delta_0 + drho*delta_0) * dPhi2.transpose()/pow(Phi2,2);
  H_lambda1 +=  (delta_1+rho*delta_0) *(H_Phi2/pow(Phi2,2) - 2*dPhi2 * dPhi2.transpose()/pow(Phi2,3));
  
  ///////////////////////////
  // Hessian for the Omegas 
  ///////////////////////////
  
  
  // turn to the Omegas. 
  double phi2 = biv_normal_pdf(w(0),w(1),rho);
  
  Eigen::VectorXd dphi2(3);
  dphi2.setZero();
  
  double detSigma = 1 - rho*rho;
  double bexp = w(0)*w(0) + w(1)*w(1) - 2*w(0)*w(1)*rho;
  
  Eigen::VectorXd ddetSigma(3);
  ddetSigma.setZero();
  ddetSigma(2)= -2*rho;
  
  dphi2(0) = (-phi2)*(w(0)-rho*w(1))/(detSigma);
  dphi2(1) = (-phi2)*(w(1)-rho*w(0))/(detSigma);
  dphi2(2) = phi2*(rho/(detSigma)+ (w(0)*w(1)*detSigma-(bexp*rho))/pow(detSigma,2) );
  
  
  // Omega_{00}
  double exp_0 = w(0)*delta_0 + rho*rho*w(1)*delta_1 - rho*(1-rho*rho)*phi2;
  Eigen::MatrixXd H_phi2 = Hess_pdf(w(0),w(1),rho);
  Eigen::MatrixXd H_exp_0(3,3);
  H_exp_0.setZero();
  
  Eigen::VectorXd dw0(3);
  dw0.setZero();
  dw0(0)=1; 
  
  Eigen::VectorXd dw1(3);
  dw1.setZero();
  dw1(1)=1; 
  
  // first term 
  H_exp_0 = w(0)*H_delta0 + dw0 *grad_delta_0.transpose() + grad_delta_0 * dw0.transpose(); 
  
  // second term 
  H_exp_0 += rho*rho*(w(1)*H_delta1 + dw1 * grad_delta_1.transpose() + grad_delta_1 * dw1.transpose()); // derive delta1 + w1; twice delta1
  H_exp_0 += w(1)*2*rho*(drho * grad_delta_1.transpose() + grad_delta_1 *drho.transpose()); // derive delta1 and rho 
  H_exp_0 += w(1)*delta_1*2* drho *drho.transpose(); // twice rho 
  H_exp_0 += delta_1*2*rho*(drho * dw1.transpose() + dw1 *drho.transpose()); // derive w1 and rho
  
  // third term 
  H_exp_0 += -rho*(1-rho*rho)*H_phi2 - (dphi2*drho.transpose() + drho*dphi2.transpose())*(1-3*rho*rho); // twice phi2; phi2 and rho. 
  H_exp_0(2,2) += 3*2*rho*phi2;
  
  // gradient of exp_0. 
  Eigen::VectorXd grad_exp_0(3);
  grad_exp_0.setZero();
  for (int k=0;k<3;k++){
    grad_exp_0(k) = w(0)*grad_delta_0(k) + rho*rho*w(1)*grad_delta_1(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_0(0) += delta_0; 
  grad_exp_0(1) += rho*rho*delta_1; 
  grad_exp_0(2) += 2*rho*w(1)*delta_1 - phi2*(1-3*rho*rho);
  
  
  //double omega_00 = 1 - exp_0/Phi2 - lambda_0*lambda_0; 
  
  H_Omega00 = -H_exp_0/Phi2 + (dPhi2 * grad_exp_0.transpose() +grad_exp_0 * dPhi2.transpose() )/pow(Phi2,2) + exp_0*H_Phi2/pow(Phi2,2) -2*exp_0*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega00 += - 2*(H_lambda0 *lambda_0 + grad_lambda_0* grad_lambda_0.transpose());
  
  // Omega_{11}
  double exp_1 = w(1)*delta_1 + rho*rho*w(0)*delta_0 - rho*(1-rho*rho)*phi2;
  Eigen::MatrixXd H_exp_1(3,3);
  H_exp_1.setZero();
  
  // first term 
  H_exp_1 = w(1)*H_delta1 + dw1 *grad_delta_1.transpose() + grad_delta_1 * dw1.transpose(); 
  
  // second term 
  H_exp_1 += rho*rho*(w(0)*H_delta0 + dw0 * grad_delta_0.transpose() + grad_delta_0 * dw0.transpose()); // derive delta1 + w1; twice delta1
  H_exp_1 += w(0)*2*rho*(drho * grad_delta_0.transpose() + grad_delta_0 *drho.transpose()); // derive delta1 and rho 
  H_exp_1 += w(0)*delta_0*2* drho *drho.transpose(); // twice rho 
  H_exp_1 += delta_0*2*rho*(drho * dw0.transpose() + dw0 *drho.transpose()); // derive w1 and rho
  
  // third term 
  H_exp_1 += -rho*(1-rho*rho)*H_phi2 - (dphi2*drho.transpose() + drho*dphi2.transpose())*(1-3*rho*rho); // twice phi2; phi2 and rho. 
  H_exp_1(2,2) += 3*2*rho*phi2;
  
  // gradient of exp_1. 
  Eigen::VectorXd grad_exp_1(3);
  grad_exp_1.setZero();
  for (int k=0;k<3;k++){
    grad_exp_1(k) = w(1)*grad_delta_1(k) + rho*rho*w(0)*grad_delta_0(k) - rho*(1-rho*rho)*dphi2(k);
  }
  grad_exp_1(0) += rho*rho*delta_0; 
  grad_exp_1(1) += delta_1; 
  grad_exp_1(2) += 2*rho*w(0)*delta_0 - phi2*(1-3*rho*rho);
  
  //double omega_11 = 1 - exp_1/Phi2 - lambda_1*lambda_1; 
  
  H_Omega11 = -H_exp_1/Phi2 + (dPhi2 * grad_exp_1.transpose() +grad_exp_1 * dPhi2.transpose() )/pow(Phi2,2) + exp_1*H_Phi2/pow(Phi2,2) -2*exp_1*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega11 += - 2*(H_lambda1 *lambda_1 + grad_lambda_1* grad_lambda_1.transpose());
  
  // Omega_{01}
  double exp_01 = rho*w(1)*delta_1 + rho*w(0)*delta_0  - (1-rho*rho)*phi2;
  Eigen::VectorXd grad_exp_01(3);
  grad_exp_01.setZero();
  for (int k=0;k<3;k++){
    grad_exp_01(k) = rho*w(1)*grad_delta_1(k) + rho*w(0)*grad_delta_0(k) - (1-rho*rho)*dphi2(k);
  }
  grad_exp_01(0) += rho*delta_0; 
  grad_exp_01(1) += rho*delta_1; 
  grad_exp_01(2) += w(0)*delta_0 + w(1)*delta_1 + (2*rho* phi2);
  
  
  Eigen::MatrixXd H_exp_01(3,3);
  H_exp_01.setZero();
  
  // first term 
  H_exp_01 = w(1)*(drho * grad_delta_1.transpose() + grad_delta_1 * drho.transpose()) +  delta_1*(drho*dw1.transpose() + dw1 * drho.transpose());
  H_exp_01 += rho*(dw1 * grad_delta_1.transpose() + grad_delta_1 * dw1.transpose());
  H_exp_01 += w(1)*rho*H_delta1; 
  
  // second term 
  H_exp_01 += w(0)*(drho * grad_delta_0.transpose() + grad_delta_0 * drho.transpose()) +  delta_0*(drho*dw0.transpose() + dw0 * drho.transpose());
  H_exp_01 += rho*(dw0 * grad_delta_0.transpose() + grad_delta_0 * dw0.transpose());
  H_exp_01 += w(0)*rho*H_delta0; 
  
  // third term 
  H_exp_01 += 2*rho*(dphi2*drho.transpose() + drho*dphi2.transpose()) + (rho*rho-1)*H_phi2; 
  H_exp_01(2,2) += 2*phi2;
  
  //double omega_01 = rho - exp_01/Phi2 - lambda_0*lambda_1; 
  
  // finally the last Hessian
  H_Omega01 = -H_exp_01/Phi2 + (dPhi2 * grad_exp_01.transpose() +grad_exp_01 * dPhi2.transpose() )/pow(Phi2,2) + exp_01*H_Phi2/pow(Phi2,2) -2*exp_01*dPhi2 * dPhi2.transpose()/pow(Phi2,3);
  H_Omega01 += - H_lambda0 *lambda_1 - grad_lambda_0* grad_lambda_1.transpose() - grad_lambda_1* grad_lambda_0.transpose() - lambda_0 * H_lambda1;
  
  
  // return matrices. 
  struct Hessian_trunc_struct output_final = {H_lambda0,H_lambda1, H_Omega00,H_Omega01,H_Omega11};
  return output_final;
}


//Hessian_trunc_struct 
Eigen::MatrixXd TVBS_Hessian_biv_norm_trunc_R(Eigen::VectorXd w, double rho, int out)
{
  Eigen::MatrixXd mat_out(3,3);
  
  struct Hessian_trunc_struct output_final = TVBS_Hessian_biv_norm_trunc(w, rho);
  
  mat_out = output_final.mat1;
  if (out==2){
    mat_out = output_final.mat2;
  }
  if (out==3){
    mat_out = output_final.mat3;
  }
  if (out==4){
    mat_out = output_final.mat4;
  }
  if (out==5){
    mat_out = output_final.mat5;
  }
  
  return mat_out;
}

Eigen::MatrixXd TVBS_Jacobian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma)
{
  Eigen::MatrixXd Jacob(5,7);
  Eigen::MatrixXd Jacob_norm(3,7); // Jacobian of normalizing transformation. 
  Jacob_norm.setZero(); 
  
  Eigen::MatrixXd D(2,2);
  D.setZero();
  
  D(0,0)= sqrt(Sigma(0,0));
  D(1,1)= sqrt(Sigma(1,1));
  double rho = Sigma(0,1)/(D(0,0)*D(1,1));
  
  Eigen::VectorXd tildex = D.inverse() *(x-alpha); 
  
  // w(0) = sigma(0,0)^(-1)*(x(0)-alpha(0))
  Jacob_norm(0,0) = 1/D(0,0); 
  Jacob_norm(0,2) = -Jacob_norm(0,0);
  Jacob_norm(0,4) =- tildex(0)/(2*Sigma(0,0)); 
  
  // w(1) = sigma(1,1)^(-1)*(x(1)-alpha(1))
  Jacob_norm(1,1) = 1/D(1,1); 
  Jacob_norm(1,3) = -Jacob_norm(1,1);
  Jacob_norm(1,6) = -tildex(1)/(2*Sigma(1,1)); 
  
  // rho = Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1))
  Jacob_norm(2,4) = -rho/(2*Sigma(0,0));
  Jacob_norm(2,5) = 1/(D(0,0)*D(1,1));
  Jacob_norm(2,6) = -rho/(2*Sigma(1,1));
  
  // calculate using the normalized values 
  Eigen::VectorXd lambda_Omega = TVBS_biv_norm_trunc(tildex,rho); 
  Eigen::MatrixXd Jacob_lambda_Omega = TVBS_Jacobian_biv_norm_trunc(tildex,rho);
  
  // rescale 
  Eigen::MatrixXd Scale(5,5);
  Scale.setZero(); 
  Scale(0,0) = D(0,0);
  Scale(1,1)= D(1,1);
  Scale(2,2)= Sigma(0,0);
  Scale(3,3) = D(0,0)*D(1,1);
  Scale(4,4) = Sigma(1,1); 
  
  Jacob = Scale * Jacob_lambda_Omega * Jacob_norm; 
  Jacob(0,4) += lambda_Omega(0)/(2*D(0,0));
  Jacob(1,6) += lambda_Omega(1)/(2*D(1,1));
  Jacob(2,4) += lambda_Omega(2);
  Jacob(3,4) += lambda_Omega(3)*(D(1,1)/(2*D(0,0)));
  Jacob(3,6) += lambda_Omega(3)*(D(0,0)/(2*D(1,1)));
  Jacob(4,6) += lambda_Omega(4);
  
  Jacob(0,2) += 1;
  Jacob(1,3) += 1; 
  
  return Jacob;
  
}






struct Hessian_trunc_struct TVBS_Hessian_biv_gen_trunc(Eigen::VectorXd x, Eigen::VectorXd alpha, Eigen::MatrixXd Sigma)
{
  Eigen::MatrixXd nH_lambda0(5,5), nH_lambda1(5,5), nH_Omega00(5,5), nH_Omega01(5,5), nH_Omega11(5,5); // normalized matrices 
  Eigen::MatrixXd H_lambda0(5,5), H_lambda1(5,5), H_Omega00(5,5), H_Omega01(5,5), H_Omega11(5,5); // general including scaling.
  
  Eigen::MatrixXd Jacob(5,7);
  Eigen::MatrixXd Jacob_norm(3,7); // Jacobian of normalizing transformation. 
  Jacob_norm.setZero(); 
  
  Eigen::MatrixXd D(2,2);
  D.setZero();
  
  D(0,0)= sqrt(Sigma(0,0));
  D(1,1)= sqrt(Sigma(1,1));
  double rho = Sigma(0,1)/(D(0,0)*D(1,1));
  
  Eigen::VectorXd tildex = D.inverse() *(x-alpha); 
  
  // w(0) = sigma(0,0)^(-1)*(x(0)-alpha(0))
  Jacob_norm(0,0) = 1/D(0,0); 
  Jacob_norm(0,2) = -Jacob_norm(0,0);
  Jacob_norm(0,4) =- tildex(0)/(2*Sigma(0,0)); 
  
  // w(1) = sigma(1,1)^(-1)*(x(1)-alpha(1))
  Jacob_norm(1,1) = 1/D(1,1); 
  Jacob_norm(1,3) = -Jacob_norm(1,1);
  Jacob_norm(1,6) = -tildex(1)/(2*Sigma(1,1)); 
  
  // rho = Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1))
  Jacob_norm(2,4) = -rho/(2*Sigma(0,0));
  Jacob_norm(2,5) = 1/(D(0,0)*D(1,1));
  Jacob_norm(2,6) = -rho/(2*Sigma(1,1));
  
  // calculate using the normalized values 
  Eigen::VectorXd lambda_Omega = TVBS_biv_norm_trunc(tildex,rho); 
  // Jacobian of lambdaj,Omegaij w.r.t normalized w, rho.: 5x3. 
  Eigen::MatrixXd Jacob_lambda_Omega = TVBS_Jacobian_biv_norm_trunc(tildex,rho);
  
  // Jacobian without rescaling for easier Hessian calculation: 5x 7.  
  Jacob = Jacob_lambda_Omega * Jacob_norm; 
  
  // H_fmat: Hessian for normalized quantities: 5 matrices of size 3x3 each. 
  struct Hessian_trunc_struct H_fmat = TVBS_Hessian_biv_norm_trunc(tildex, rho);
  
  // Hessian for normalisation of each entry separately. 
  Eigen::MatrixXd HT_1(7,7), HT_2(7,7), HT_3(7,7);
  HT_1.setZero();
  HT_2.setZero();
  HT_3.setZero();
  
  // HT_1 Hessian of w0 w.r.t. x,alpha, Sigma. 
  HT_1(0,4)= -1/(2*sqrt(Sigma(0,0))*Sigma(0,0));
  HT_1(2,4)= 1/(2*sqrt(Sigma(0,0))*Sigma(0,0));
  HT_1(4,4)= 3*tildex(0)/(4*Sigma(0,0)*Sigma(0,0));
  HT_1(4,0) = -0.5/(Sigma(0,0)*sqrt(Sigma(0,0)));
  HT_1(4,2) = 0.5/(Sigma(0,0)*sqrt(Sigma(0,0)));
  
  // HT_2 Hessian of w1 w.r.t. x,alpha, Sigma. 
  HT_2(1,6)= -1/(2*sqrt(Sigma(1,1))*Sigma(1,1));
  HT_2(3,6)= 1/(2*sqrt(Sigma(1,1))*Sigma(1,1));
  HT_2(6,6)= 3*tildex(1)/(4*Sigma(1,1)*Sigma(1,1));
  HT_2(6,1) = -0.5/(Sigma(1,1)*sqrt(Sigma(1,1)));
  HT_2(6,3) = 0.5/(Sigma(1,1)*sqrt(Sigma(1,1)));
  
  // HT_3 Hessian of rho w.r.t. x,alpha, Sigma. 
  HT_3(4,4)= 3*rho/(4*Sigma(0,0)*Sigma(0,0));
  HT_3(5,4) = -0.5/(Sigma(0,0)*sqrt(Sigma(0,0)*Sigma(1,1)));
  HT_3(6,4)= rho/(4*Sigma(0,0)*Sigma(1,1));
  HT_3(4,5) = HT_3(5,4);
  
  HT_3(6,6)= 3*rho/(4*Sigma(1,1)*Sigma(1,1));
  HT_3(5,6) = -0.5/(Sigma(1,1)*sqrt(Sigma(0,0)*Sigma(1,1)));
  HT_3(4,6)= rho/(4*Sigma(0,0)*Sigma(1,1));
  HT_3(6,5) = HT_3(5,6);
  
  // Hessian without rescaling:
  // Hessian for lambda_0
  nH_lambda0 = Jacob_norm.transpose() * H_fmat.mat1 * Jacob_norm;
  nH_lambda0 += HT_1* Jacob_lambda_Omega(0,0) +  HT_2* Jacob_lambda_Omega(0,1) +  HT_3* Jacob_lambda_Omega(0,2);
  
  // Hessian for lambda_1
  nH_lambda1 = Jacob_norm.transpose() * H_fmat.mat2 * Jacob_norm;
  nH_lambda1 += HT_1* Jacob_lambda_Omega(1,0) +  HT_2* Jacob_lambda_Omega(1,1) +  HT_3* Jacob_lambda_Omega(1,2);
  
  // Hessian for Omega_00
  nH_Omega00 = Jacob_norm.transpose() * H_fmat.mat3 * Jacob_norm;
  nH_Omega00 += HT_1* Jacob_lambda_Omega(2,0) +  HT_2* Jacob_lambda_Omega(2,1) +  HT_3* Jacob_lambda_Omega(2,2);
  
  // Hessian for Omega_01
  nH_Omega01 = Jacob_norm.transpose() * H_fmat.mat4 * Jacob_norm;
  nH_Omega01 += HT_1* Jacob_lambda_Omega(3,0) +  HT_2* Jacob_lambda_Omega(3,1) +  HT_3* Jacob_lambda_Omega(3,2);
  
  // Hessian for Omega_11
  nH_Omega11 = Jacob_norm.transpose() * H_fmat.mat5 * Jacob_norm;
  nH_Omega11 += HT_1* Jacob_lambda_Omega(4,0) +  HT_2* Jacob_lambda_Omega(4,1) +  HT_3* Jacob_lambda_Omega(4,2);
  
  // rescale 
  Eigen::MatrixXd DS_1(7,1),DS_2(7,1),DS_3(7,1),DS_4(7,1),DS_5(7,1);
  DS_1.setZero();
  DS_2.setZero();
  DS_3.setZero();
  DS_4.setZero();
  DS_5.setZero();
  
  DS_1(4,0) = 0.5/sqrt(Sigma(0,0));
  DS_2(6,0)= 0.5/sqrt(Sigma(1,1));
  DS_3(4,0) = 1;
  DS_4(4,0) = 0.5*sqrt(Sigma(1,1)/Sigma(0,0));
  DS_4(6,0) = 0.5*sqrt(Sigma(0,0)/Sigma(1,1));
  DS_5(6,0)= 1;
  
  Eigen::MatrixXd HS_1(7,7),HS_2(7,7),HS_3(7,7),HS_4(7,7),HS_5(7,7);
  HS_1.setZero();
  HS_2.setZero();
  HS_3.setZero();
  HS_4.setZero();
  HS_5.setZero();
  
  HS_1(4,4) = -DS_1(4,0)/(2*Sigma(0,0));
  HS_2(6,6) = -DS_2(6,0)/(2*Sigma(1,1));
  
  HS_4(4,4) = -DS_4(4,0)/(2*Sigma(0,0));
  HS_4(6,6) = -DS_4(6,0)/(2*Sigma(1,1));
  HS_4(4,6) = 0.25/sqrt(Sigma(0,0)*Sigma(1,1));
  HS_4(6,4) = 0.25/sqrt(Sigma(0,0)*Sigma(1,1));
  
  // HS_5 and HS_3 are zero. 
  H_lambda0 = lambda_Omega(0)*HS_1 + DS_1 * Jacob.block(0,0,1,7)  + Jacob.block(0,0,1,7).transpose() * DS_1.transpose() + nH_lambda0 *sqrt(Sigma(0,0));
  H_lambda1 = lambda_Omega(1)*HS_2 + DS_2 * Jacob.block(1,0,1,7)  + Jacob.block(1,0,1,7).transpose() * DS_2.transpose() + nH_lambda1 *sqrt(Sigma(1,1));
  
  H_Omega00 = lambda_Omega(2)*HS_3 + DS_3 * Jacob.block(2,0,1,7)  + Jacob.block(2,0,1,7).transpose() * DS_3.transpose() + nH_Omega00 * (Sigma(0,0));
  H_Omega01 = lambda_Omega(3)*HS_4 + DS_4 * Jacob.block(3,0,1,7)  + Jacob.block(3,0,1,7).transpose() * DS_4.transpose() + nH_Omega01 * sqrt(Sigma(0,0)*Sigma(1,1));
  H_Omega11 = lambda_Omega(4)*HS_5 + DS_5 * Jacob.block(4,0,1,7)  + Jacob.block(4,0,1,7).transpose() * DS_5.transpose() + nH_Omega11 * (Sigma(1,1));
  
  
  // return matrices. 
  struct Hessian_trunc_struct output_final = {H_lambda0,H_lambda1, H_Omega00,H_Omega01,H_Omega11};
  return output_final;
  
}


