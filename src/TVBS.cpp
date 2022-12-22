#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;
//TVBS Approximation
//INPUTS
//- x: a vector containing the upper integration limits
//- mu: a vector containing the means of the distribution
//- Sigma: a covariance matrix
//- ll: 0 if Probability should be the output, 1 if the logarithm should be returned
//OUTPUT
//- TVBS: (log)-probability


// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;                                               // 'maps' rather than copies


#include "toms462.h"                                          // allows bivariate normal cdf calculations
#include "TVBS_vdb.h"


#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <random>       // for random numbers



// // // // // // // // // // // // // // // // //
//                                              //
//    DECLARATION OF GLOBAL VARIABLES           //
//                                              //
// // // // // // // // // // // // // // // // //

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


// // // // // // // // // // // // // // // // //
//                                              //
//    STRUCTURES                                //
//                                              //
// // // // // // // // // // // // // // // // //

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




// // // // // // // // // // // // // // // // //
//                                              //
//    PDF AND CDF FUNCTIONS                     //
//                                              //
// // // // // // // // // // // // // // // // //

double TVBS_std_normal_cdf_cpp(double x)                                 // create a function for the standard normal cumulative distribution function
{
  if (x>6) { x = 6;}
  if (x<-6) { x = -6;}
  
  return 0.5 * (1.0 + std::erf(x / 1.4142135623730950));            // computes 1/2 * (1 + int_0^(x/sqrt(2)) exp(-t^2) dt = CDF of standard normal dist.
}



double TVBS_std_normal_pdf_cpp(double x)                                 // create a function for the standard normal probability distribution function
{  
  if (x>6) { x = 6;}
  if (x<-6) { x = -6;}  
  return 0.3989422804014327 * std::exp(-0.5 * x * x);                 // compute PDF 1/sqrt(2+pi+sigma^2)*exp(-(x-mu)/(2*sigma^2))
}

double TVBS_std_normal_cdf_inv_1(double y)
{
  // Beasley-Springer-Moro algorithm
  // From "Monte Carlo Methods in Financial Engineering" - Paul Glasserman (2003)
  // pp 67
  double r, x;
  
  y += - 0.5;
  
  if(abs(y)<0.42)
  {
    r = y*y;
    x = y*((((-25.44106049637)*r+41.39119773534)*r+(-18.61500062529))*r+(2.50662823884))/(((((3.13082909833)*r+(-21.06224101826))*r+(23.08336743743))*r+(-8.47351093090))*r+1);
  }
  else
  {
    r = y+0.5;
    if(y>0)
    {
      r = 0.5-y;
    }
    r = log(-log(r));
    x = 0.3374754822726147 + r*(0.9761690190917186+r*(0.1607979714918209+r*(0.0276438810333863+r*(0.0038405729373609+r*(0.0003951896511919+r*(0.0000321767881768+r*(0.0000002888167364+r*0.0000003960315187)))))));
    if(y<0)
    {
      x = -x;
    }
  }
       
  return x;
}


Eigen::VectorXd TVBS_std_normal_cdf_inv(Eigen::VectorXd y)
{
  // Beasley-Springer-Moro algorithm
  // From "Monte Carlo Methods in Financial Engineering" - Paul Glasserman (2003)
  // pp 67
  double r;
  
  
  Eigen::VectorXd x(y.size());
  y.array() += - 0.5;
  
  for(int j = 0; j < y.size(); j++ )
  {
    if(abs(y(j))<0.42)
    {
      r = y(j)*y(j);
      x(j) = y(j)*((((-25.44106049637)*r+41.39119773534)*r+(-18.61500062529))*r+(2.50662823884))/(((((3.13082909833)*r+(-21.06224101826))*r+(23.08336743743))*r+(-8.47351093090))*r+1);
    }
    else
    {
      r = y(j)+0.5;
      if(y(j)>0)
      {
        r = 1-y(j)-0.5;
      }
      r = log(-log(r));
      x(j) = 0.3374754822726147 + r*(0.9761690190917186+r*(0.1607979714918209+r*(0.0276438810333863+r*(0.0038405729373609+r*(0.0003951896511919+r*(0.0000321767881768+r*(0.0000002888167364+r*0.0000003960315187)))))));
      if(y(j)<0)
      {
        x(j) = -x(j);
      }
    }
  }
  
  return x;
}





// // // // // // // // // // // // // // // // //
//                                              //
//    Seperation Of Variables Method (SOV)      //
//                                              //
// // // // // // // // // // // // // // // // //




Eigen::VectorXi TVBS_v_Korobov_cpp(int h, int K, int N)
{
  Eigen::VectorXi output(K);
  output(0) = 1;
  for(int i=1; i<K; i++)
  {
    output(i) = (output(i-1)*h) % N;
  }
  return output;
}


Eigen::MatrixXd TVBS_L_N_shift_cpp(int h, int K, int N, int s = 0)
{
  Eigen::MatrixXd mat_n_seq(1,N);
  mat_n_seq = Eigen::VectorXd::LinSpaced(N,0,N-1);

  Eigen::VectorXd test_1 = TVBS_v_Korobov_cpp(h,K,N).cast <double> ();

  Eigen::MatrixXd test_2 = test_1*mat_n_seq.transpose()/N;
  
  std::mt19937_64 rng;
  uint64_t theSeed = s;
  std::seed_seq ss{uint32_t(theSeed & 0xffffffff), uint32_t(theSeed>>32)};
  rng.seed(ss);
  std::uniform_real_distribution<double> unif(0, 1);
  
  for(int i=0; i<K; i++)
  {
    double r = unif(rng);
    test_2.row(i) = test_2.row(i).array() + r;
  }
  
  return (test_2).unaryExpr([](const double x) { return std::fmod(x,1); });
}



// pre-define lattices
//Eigen::MatrixXd Lattice_N_shift_3_00(3,263), Lattice_N_shift_3_01(3,263), Lattice_N_shift_3_02(3,263), Lattice_N_shift_3_03(3,263), Lattice_N_shift_3_04(3,263), Lattice_N_shift_3_05(3,263), Lattice_N_shift_3_06(3,263), Lattice_N_shift_3_07(3,263), Lattice_N_shift_3_08(3,263), Lattice_N_shift_3_09(3,263), Lattice_N_shift_3_10(3,263), Lattice_N_shift_3_11(3,263);
Eigen::MatrixXd Lattice_N_shift_3_00_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 0*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_01_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 1*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_02_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 2*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_03_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 3*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_04_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 4*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_05_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 5*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_06_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 6*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_07_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 7*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_08_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 8*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_09_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+ 9*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_10_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+10*85).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_3_11_L = abs(2.0*TVBS_L_N_shift_cpp(85, 3, 263, 0+11*85).array()-1.0);

Eigen::MatrixXd Lattice_N_shift_3_00_R = 1.0-Lattice_N_shift_3_00_L.array();
Eigen::MatrixXd Lattice_N_shift_3_01_R = 1.0-Lattice_N_shift_3_01_L.array();
Eigen::MatrixXd Lattice_N_shift_3_02_R = 1.0-Lattice_N_shift_3_02_L.array();
Eigen::MatrixXd Lattice_N_shift_3_03_R = 1.0-Lattice_N_shift_3_03_L.array();
Eigen::MatrixXd Lattice_N_shift_3_04_R = 1.0-Lattice_N_shift_3_04_L.array();
Eigen::MatrixXd Lattice_N_shift_3_05_R = 1.0-Lattice_N_shift_3_05_L.array();
Eigen::MatrixXd Lattice_N_shift_3_06_R = 1.0-Lattice_N_shift_3_06_L.array();
Eigen::MatrixXd Lattice_N_shift_3_07_R = 1.0-Lattice_N_shift_3_07_L.array();
Eigen::MatrixXd Lattice_N_shift_3_08_R = 1.0-Lattice_N_shift_3_08_L.array();
Eigen::MatrixXd Lattice_N_shift_3_09_R = 1.0-Lattice_N_shift_3_09_L.array();
Eigen::MatrixXd Lattice_N_shift_3_10_R = 1.0-Lattice_N_shift_3_10_L.array();
Eigen::MatrixXd Lattice_N_shift_3_11_R = 1.0-Lattice_N_shift_3_11_L.array();


//Eigen::MatrixXd Lattice_N_shift_4_00(4,383), Lattice_N_shift_4_01(4,383), Lattice_N_shift_4_02(4,383), Lattice_N_shift_4_03(4,383), Lattice_N_shift_4_04(4,383), Lattice_N_shift_4_05(4,383), Lattice_N_shift_4_06(4,383), Lattice_N_shift_4_07(4,383), Lattice_N_shift_4_08(4,383), Lattice_N_shift_4_09(4,383), Lattice_N_shift_4_10(4,383), Lattice_N_shift_4_11(4,383);
Eigen::MatrixXd Lattice_N_shift_4_00_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 0*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_01_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 1*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_02_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 2*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_03_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 3*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_04_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 4*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_05_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 5*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_06_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 6*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_07_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 7*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_08_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 8*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_09_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+ 9*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_10_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+10*62).array()-1.0);
Eigen::MatrixXd Lattice_N_shift_4_11_L = abs(2.0*TVBS_L_N_shift_cpp(62, 4, 383, 0+11*62).array()-1.0);

Eigen::MatrixXd Lattice_N_shift_4_00_R = 1.0-Lattice_N_shift_4_00_L.array();
Eigen::MatrixXd Lattice_N_shift_4_01_R = 1.0-Lattice_N_shift_4_01_L.array();
Eigen::MatrixXd Lattice_N_shift_4_02_R = 1.0-Lattice_N_shift_4_02_L.array();
Eigen::MatrixXd Lattice_N_shift_4_03_R = 1.0-Lattice_N_shift_4_03_L.array();
Eigen::MatrixXd Lattice_N_shift_4_04_R = 1.0-Lattice_N_shift_4_04_L.array();
Eigen::MatrixXd Lattice_N_shift_4_05_R = 1.0-Lattice_N_shift_4_05_L.array();
Eigen::MatrixXd Lattice_N_shift_4_06_R = 1.0-Lattice_N_shift_4_06_L.array();
Eigen::MatrixXd Lattice_N_shift_4_07_R = 1.0-Lattice_N_shift_4_07_L.array();
Eigen::MatrixXd Lattice_N_shift_4_08_R = 1.0-Lattice_N_shift_4_08_L.array();
Eigen::MatrixXd Lattice_N_shift_4_09_R = 1.0-Lattice_N_shift_4_09_L.array();
Eigen::MatrixXd Lattice_N_shift_4_10_R = 1.0-Lattice_N_shift_4_10_L.array();
Eigen::MatrixXd Lattice_N_shift_4_11_R = 1.0-Lattice_N_shift_4_11_L.array();




double TVBS_pmvnorm_old_cpp(Eigen::VectorXd upper, Eigen::VectorXd mu, Eigen::MatrixXd Sigma, int N = 0, int h = 0, int s = 0)
{
  int k         = upper.rows();
  double output;
  
  if( k == 1 )
  {
    output       = TVBS_std_normal_cdf_cpp((upper(0) - mu(0)) / sqrt(Sigma(0,0)));
  }
  if( k == 2 )
  {
    output       = bivnor(-(upper(0)-mu(0))/sqrt(Sigma(0,0)), -(upper(1)-mu(1))/sqrt(Sigma(1,1)), Sigma(0,1)/sqrt(Sigma(0,0)*Sigma(1,1)));
  }
  else
  {
    
    // reorder for better numerical performance / accuracy
    // and compute Cholesky decomposition
    
    struct TVBS_Vec_Mat b_L = TVBS_reorder_chol_cpp(upper-mu, Sigma);
    Eigen::VectorXd b = b_L.vec1;
    Eigen::MatrixXd L = b_L.mat1;
    
    int M = 12;
    Eigen::VectorXd I_Ni(M); 
    
    // standard method
    if(N==0 && h==0 && s==0){
      
      // set parameters for QMC method
      if(k == 3)
      {
        N = 263;
        h = 85;
      }
      if(k > 3)
      {
        N = 383;
        h = 62;
      }
      
      Eigen::VectorXd Lattice_evaluations_00(N), Lattice_evaluations_01(N), Lattice_evaluations_02(N), Lattice_evaluations_03(N), Lattice_evaluations_04(N), Lattice_evaluations_05(N), Lattice_evaluations_06(N), Lattice_evaluations_07(N), Lattice_evaluations_08(N), Lattice_evaluations_09(N), Lattice_evaluations_10(N), Lattice_evaluations_11(N);
      
      if(k==3)
      {
        Eigen::VectorXd e_00_L(k), z_norm_inv_00_L(k-1), e_00_R(k), z_norm_inv_00_R(k-1);
        e_00_L.setZero(), z_norm_inv_00_L.setZero(), e_00_R.setZero(), z_norm_inv_00_R.setZero();
        e_00_L(0) = TVBS_std_normal_cdf_cpp(b(0)/L(0,0));
        e_00_R(0) = e_00_L(0);
        Eigen::VectorXd e_01_L = e_00_L, e_02_L = e_00_L, e_03_L = e_00_L, e_04_L = e_00_L, e_05_L = e_00_L, e_06_L = e_00_L, e_07_L = e_00_L, e_08_L = e_00_L, e_09_L = e_00_L, e_10_L = e_00_L, e_11_L = e_00_L;
        Eigen::VectorXd e_01_R = e_00_R, e_02_R = e_00_R, e_03_R = e_00_R, e_04_R = e_00_R, e_05_R = e_00_R, e_06_R = e_00_R, e_07_R = e_00_R, e_08_R = e_00_R, e_09_R = e_00_R, e_10_R = e_00_R, e_11_R = e_00_R;
        Eigen::VectorXd z_norm_inv_01_L(k-1), z_norm_inv_02_L(k-1), z_norm_inv_03_L(k-1), z_norm_inv_04_L(k-1), z_norm_inv_05_L(k-1), z_norm_inv_06_L(k-1), z_norm_inv_07_L(k-1), z_norm_inv_08_L(k-1), z_norm_inv_09_L(k-1), z_norm_inv_10_L(k-1), z_norm_inv_11_L(k-1);
        Eigen::VectorXd z_norm_inv_01_R(k-1), z_norm_inv_02_R(k-1), z_norm_inv_03_R(k-1), z_norm_inv_04_R(k-1), z_norm_inv_05_R(k-1), z_norm_inv_06_R(k-1), z_norm_inv_07_R(k-1), z_norm_inv_08_R(k-1), z_norm_inv_09_R(k-1), z_norm_inv_10_R(k-1), z_norm_inv_11_R(k-1);
        for(int i=0; i<N; i++)
        {
          for(int i2=1; i2<k; i2++)
          {
            z_norm_inv_00_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_00_L(i2-1)*Lattice_N_shift_3_00_L(i2-1,i), tol));
            e_00_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_00_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_00_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_00_R(i2-1)*Lattice_N_shift_3_00_R(i2-1,i), tol));
            e_00_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_00_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_01_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_01_L(i2-1)*Lattice_N_shift_3_01_L(i2-1,i), tol));
            e_01_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_01_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_01_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_01_R(i2-1)*Lattice_N_shift_3_01_R(i2-1,i), tol));
            e_01_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_01_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_02_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_02_L(i2-1)*Lattice_N_shift_3_02_L(i2-1,i), tol));
            e_02_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_02_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_02_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_02_R(i2-1)*Lattice_N_shift_3_02_R(i2-1,i), tol));
            e_02_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_02_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_03_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_03_L(i2-1)*Lattice_N_shift_3_03_L(i2-1,i), tol));
            e_03_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_03_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_03_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_03_R(i2-1)*Lattice_N_shift_3_03_R(i2-1,i), tol));
            e_03_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_03_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_04_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_04_L(i2-1)*Lattice_N_shift_3_04_L(i2-1,i), tol));
            e_04_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_04_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_04_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_04_R(i2-1)*Lattice_N_shift_3_04_R(i2-1,i), tol));
            e_04_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_04_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_05_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_05_L(i2-1)*Lattice_N_shift_3_05_L(i2-1,i), tol));
            e_05_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_05_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_05_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_05_R(i2-1)*Lattice_N_shift_3_05_R(i2-1,i), tol));
            e_05_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_05_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_06_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_06_L(i2-1)*Lattice_N_shift_3_06_L(i2-1,i), tol));
            e_06_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_06_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_06_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_06_R(i2-1)*Lattice_N_shift_3_06_R(i2-1,i), tol));
            e_06_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_06_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_07_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_07_L(i2-1)*Lattice_N_shift_3_07_L(i2-1,i), tol));
            e_07_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_07_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_07_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_07_R(i2-1)*Lattice_N_shift_3_07_R(i2-1,i), tol));
            e_07_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_07_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_08_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_08_L(i2-1)*Lattice_N_shift_3_08_L(i2-1,i), tol));
            e_08_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_08_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_08_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_08_R(i2-1)*Lattice_N_shift_3_08_R(i2-1,i), tol));
            e_08_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_08_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_09_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_09_L(i2-1)*Lattice_N_shift_3_09_L(i2-1,i), tol));
            e_09_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_09_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_09_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_09_R(i2-1)*Lattice_N_shift_3_09_R(i2-1,i), tol));
            e_09_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_09_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_10_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_10_L(i2-1)*Lattice_N_shift_3_10_L(i2-1,i), tol));
            e_10_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_10_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_10_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_10_R(i2-1)*Lattice_N_shift_3_10_R(i2-1,i), tol));
            e_10_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_10_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_11_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_11_L(i2-1)*Lattice_N_shift_3_11_L(i2-1,i), tol));
            e_11_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_11_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_11_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_11_R(i2-1)*Lattice_N_shift_3_11_R(i2-1,i), tol));
            e_11_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_11_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
          }
          Lattice_evaluations_00(i) = (e_00_L.prod()+e_00_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_01(i) = (e_01_L.prod()+e_01_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_02(i) = (e_02_L.prod()+e_02_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_03(i) = (e_03_L.prod()+e_03_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_04(i) = (e_04_L.prod()+e_04_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_05(i) = (e_05_L.prod()+e_05_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_06(i) = (e_06_L.prod()+e_06_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_07(i) = (e_07_L.prod()+e_07_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_08(i) = (e_08_L.prod()+e_08_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_09(i) = (e_09_L.prod()+e_09_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_10(i) = (e_10_L.prod()+e_10_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_11(i) = (e_11_L.prod()+e_11_R.prod())/sqrt(2.0*N);
          
        }
      } else {
        
        Eigen::VectorXd e_00_L(k), z_norm_inv_00_L(k-1), e_00_R(k), z_norm_inv_00_R(k-1);
        e_00_L.setZero(), z_norm_inv_00_L.setZero(), e_00_R.setZero(), z_norm_inv_00_R.setZero();
        e_00_L(0) = TVBS_std_normal_cdf_cpp(b(0)/L(0,0));
        e_00_R(0) = e_00_L(0);
        Eigen::VectorXd e_01_L = e_00_L, e_02_L = e_00_L, e_03_L = e_00_L, e_04_L = e_00_L, e_05_L = e_00_L, e_06_L = e_00_L, e_07_L = e_00_L, e_08_L = e_00_L, e_09_L = e_00_L, e_10_L = e_00_L, e_11_L = e_00_L;
        Eigen::VectorXd e_01_R = e_00_R, e_02_R = e_00_R, e_03_R = e_00_R, e_04_R = e_00_R, e_05_R = e_00_R, e_06_R = e_00_R, e_07_R = e_00_R, e_08_R = e_00_R, e_09_R = e_00_R, e_10_R = e_00_R, e_11_R = e_00_R;
        Eigen::VectorXd z_norm_inv_01_L(k-1), z_norm_inv_02_L(k-1), z_norm_inv_03_L(k-1), z_norm_inv_04_L(k-1), z_norm_inv_05_L(k-1), z_norm_inv_06_L(k-1), z_norm_inv_07_L(k-1), z_norm_inv_08_L(k-1), z_norm_inv_09_L(k-1), z_norm_inv_10_L(k-1), z_norm_inv_11_L(k-1);
        Eigen::VectorXd z_norm_inv_01_R(k-1), z_norm_inv_02_R(k-1), z_norm_inv_03_R(k-1), z_norm_inv_04_R(k-1), z_norm_inv_05_R(k-1), z_norm_inv_06_R(k-1), z_norm_inv_07_R(k-1), z_norm_inv_08_R(k-1), z_norm_inv_09_R(k-1), z_norm_inv_10_R(k-1), z_norm_inv_11_R(k-1);
        
        for(int i=0; i<N; i++)
        {
          for(int i2=1; i2<k; i2++)
          {
            z_norm_inv_00_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_00_L(i2-1)*Lattice_N_shift_4_00_L(i2-1,i), tol));
            e_00_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_00_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_00_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_00_R(i2-1)*Lattice_N_shift_4_00_R(i2-1,i), tol));
            e_00_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_00_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_01_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_01_L(i2-1)*Lattice_N_shift_4_01_L(i2-1,i), tol));
            e_01_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_01_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_01_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_01_R(i2-1)*Lattice_N_shift_4_01_R(i2-1,i), tol));
            e_01_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_01_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_02_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_02_L(i2-1)*Lattice_N_shift_4_02_L(i2-1,i), tol));
            e_02_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_02_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_02_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_02_R(i2-1)*Lattice_N_shift_4_02_R(i2-1,i), tol));
            e_02_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_02_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_03_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_03_L(i2-1)*Lattice_N_shift_4_03_L(i2-1,i), tol));
            e_03_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_03_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_03_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_03_R(i2-1)*Lattice_N_shift_4_03_R(i2-1,i), tol));
            e_03_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_03_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_04_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_04_L(i2-1)*Lattice_N_shift_4_04_L(i2-1,i), tol));
            e_04_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_04_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_04_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_04_R(i2-1)*Lattice_N_shift_4_04_R(i2-1,i), tol));
            e_04_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_04_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_05_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_05_L(i2-1)*Lattice_N_shift_4_05_L(i2-1,i), tol));
            e_05_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_05_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_05_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_05_R(i2-1)*Lattice_N_shift_4_05_R(i2-1,i), tol));
            e_05_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_05_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_06_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_06_L(i2-1)*Lattice_N_shift_4_06_L(i2-1,i), tol));
            e_06_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_06_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_06_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_06_R(i2-1)*Lattice_N_shift_4_06_R(i2-1,i), tol));
            e_06_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_06_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_07_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_07_L(i2-1)*Lattice_N_shift_4_07_L(i2-1,i), tol));
            e_07_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_07_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_07_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_07_R(i2-1)*Lattice_N_shift_4_07_R(i2-1,i), tol));
            e_07_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_07_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_08_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_08_L(i2-1)*Lattice_N_shift_4_08_L(i2-1,i), tol));
            e_08_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_08_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_08_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_08_R(i2-1)*Lattice_N_shift_4_08_R(i2-1,i), tol));
            e_08_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_08_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_09_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_09_L(i2-1)*Lattice_N_shift_4_09_L(i2-1,i), tol));
            e_09_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_09_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_09_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_09_R(i2-1)*Lattice_N_shift_4_09_R(i2-1,i), tol));
            e_09_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_09_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_10_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_10_L(i2-1)*Lattice_N_shift_4_10_L(i2-1,i), tol));
            e_10_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_10_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_10_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_10_R(i2-1)*Lattice_N_shift_4_10_R(i2-1,i), tol));
            e_10_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_10_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
            z_norm_inv_11_L(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_11_L(i2-1)*Lattice_N_shift_4_11_L(i2-1,i), tol));
            e_11_L(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_11_L.transpose().head(i2).array())).sum())/L(i2,i2)));
            z_norm_inv_11_R(i2-1) = TVBS_std_normal_cdf_inv_1(std::max(e_11_R(i2-1)*Lattice_N_shift_4_11_R(i2-1,i), tol));
            e_11_R(i2) = TVBS_std_normal_cdf_cpp(((b(i2)-((L.row(i2).head(i2).array())*(z_norm_inv_11_R.transpose().head(i2).array())).sum())/L(i2,i2)));
            
          }
          Lattice_evaluations_00(i) = (e_00_L.prod()+e_00_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_01(i) = (e_01_L.prod()+e_01_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_02(i) = (e_02_L.prod()+e_02_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_03(i) = (e_03_L.prod()+e_03_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_04(i) = (e_04_L.prod()+e_04_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_05(i) = (e_05_L.prod()+e_05_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_06(i) = (e_06_L.prod()+e_06_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_07(i) = (e_07_L.prod()+e_07_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_08(i) = (e_08_L.prod()+e_08_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_09(i) = (e_09_L.prod()+e_09_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_10(i) = (e_10_L.prod()+e_10_R.prod())/sqrt(2.0*N);
          Lattice_evaluations_11(i) = (e_11_L.prod()+e_11_R.prod())/sqrt(2.0*N);
          
        }
      }
      
      I_Ni(0)  = Lattice_evaluations_00.sum()/sqrt(2.0*N);
      I_Ni(1)  = Lattice_evaluations_01.sum()/sqrt(2.0*N);
      I_Ni(2)  = Lattice_evaluations_02.sum()/sqrt(2.0*N);
      I_Ni(3)  = Lattice_evaluations_03.sum()/sqrt(2.0*N);
      I_Ni(4)  = Lattice_evaluations_04.sum()/sqrt(2.0*N);
      I_Ni(5)  = Lattice_evaluations_05.sum()/sqrt(2.0*N);
      I_Ni(6)  = Lattice_evaluations_06.sum()/sqrt(2.0*N);
      I_Ni(7)  = Lattice_evaluations_07.sum()/sqrt(2.0*N);
      I_Ni(8)  = Lattice_evaluations_08.sum()/sqrt(2.0*N);
      I_Ni(9)  = Lattice_evaluations_09.sum()/sqrt(2.0*N);
      I_Ni(10) = Lattice_evaluations_10.sum()/sqrt(2.0*N);
      I_Ni(11) = Lattice_evaluations_11.sum()/sqrt(2.0*N);
      
      
    } else {
      // set parameters for QMC method
      if(N==0)
      {
        if(k == 3)
        {
          N = 263;
          h = 85;
        }
        if(k > 3)
        {
          N = 383;
          h = 62;
        }
      }
      
      if(h==0)
      {
        h = N/3;
      }
      
      Eigen::VectorXd Lattice_evaluations(N);
      Eigen::MatrixXd Lattice_N_shift(k,N);
      for(int j = 0; j < M; j++)
      {
        Lattice_N_shift = TVBS_L_N_shift_cpp(h, k, N, s+j*h);
        
        for(int i=0; i<N; i++)
        {
          Eigen::VectorXd Lattice_N_shift_i = abs(2.0*Lattice_N_shift.col(i).array()-1.0);
          Eigen::VectorXd Lattice_N_shift_i_2 = (1.0-abs(2.0*Lattice_N_shift.col(i).array()-1.0));
          Lattice_evaluations(i) = (TVBS_std_normal_cdf_inv(Lattice_N_shift_i , b, L) + TVBS_std_normal_cdf_inv(Lattice_N_shift_i_2, b, L))/sqrt(2.0*N);
        }

        I_Ni(j) = Lattice_evaluations.sum()/sqrt(2.0*N);
      }

    }
    
    output = I_Ni.sum()/(M*1.0);
  }
  
  
  return output;
}







// // // // // // // // //
//  trivariate          //
// // // // // // // // //

struct TVBS_Matrix_two TVBSo_grad_cdf_tvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr)
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
    double tempval1 = TVBS_std_normal_cdf_cpp(h1);
    
    Eigen::VectorXd h23(2), mu_h23(2);
    mu_h23.setZero();
    h23 << h2, h3;
    Eigen::MatrixXd corr_h23(2,2);
    corr_h23.setOnes();
    corr_h23(0,1) = r23;
    corr_h23(1,0) = r23;
    double tempval2 = TVBS_pmvnorm_old_cpp(h23, mu_h23, corr_h23);
    
    gg(0) = TVBS_std_normal_pdf_cpp(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h23, corr_h23);
    gg(1) = tempval1*temp_bvn_grad.d1;
    gg(2) = tempval1*temp_bvn_grad.d2;
    gg(5) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r13) + abs(r23))<epst )
  {
    // 5
    double tempval1 = TVBS_std_normal_cdf_cpp(h3);
    
    Eigen::VectorXd h12(2), mu_h12(2);
    mu_h12.setZero();
    h12 << h1, h2;
    Eigen::MatrixXd corr_h12(2,2);
    corr_h12.setOnes();
    corr_h12(0,1) = r12;
    corr_h12(1,0) = r12;
    double tempval2 = TVBS_pmvnorm_old_cpp(h12, mu_h12, corr_h12);
    
    gg(2) = TVBS_std_normal_pdf_cpp(h1)*tempval2;
    struct TVBS_double_three temp_bvn_grad = TVBS_grad_cdf_bvn_cpp(h12, corr_h12);
    gg(0) = tempval1*temp_bvn_grad.d1;
    gg(1) = tempval1*temp_bvn_grad.d2;
    gg(3) = tempval1*temp_bvn_grad.d3;
    d_tvn = gg;
  }
  else if( (abs(r12) + abs(r23))<epst )
  {
    // 6
    double tempval1 = TVBS_std_normal_cdf_cpp(h2);
    
    Eigen::VectorXd h13(2), mu_h13(2);
    mu_h13.setZero();
    h13 << h1, h3;
    Eigen::MatrixXd corr_h13(2,2);
    corr_h13.setOnes();
    corr_h13(0,1) = r13;
    corr_h13(1,0) = r13;
    double tempval2 = TVBS_pmvnorm_old_cpp(h13, mu_h13, corr_h13);
    
    // tvn should have no influence on this function
    tvn = TVBS_std_normal_cdf_cpp(h2)*tempval2;
    
    
    gg(1) = TVBS_std_normal_pdf_cpp(h2)*tempval2;
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
    tvn = TVBS_pmvnorm_old_cpp(h1_23_min, mu_h13, corr_h12);
    
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
    
    double cdfbvn_h2_h3 = TVBS_pmvnorm_old_cpp(h_23, mu_h23, corr_h23);
    
    double cdfn_h1 = TVBS_std_normal_cdf_cpp(h1);
    tvn = cdfn_h1*cdfbvn_h2_h3;
    
    double d_a = TVBS_std_normal_pdf_cpp(h2)*TVBS_std_normal_cdf_cpp((h3-r23*h2)/sqrt(1.0-pow(r23,2)));
    double d_b = TVBS_std_normal_pdf_cpp(h3)*TVBS_std_normal_cdf_cpp((h2-r23*h3)/sqrt(1.0-pow(r23,2)));
    double d_corr = (exp(-0.5*(pow(h2,2) + pow(h3,2) - 2*r23 * h2 * h3   ) / (1.0-r23*r23) )) / sqrt(1.0-r23 * r23)/(2*pi_global);
        
    d_tvn << (TVBS_std_normal_pdf_cpp(h1)*cdfbvn_h2_h3), cdfn_h1*d_a, (cdfn_h1*d_b), 0, 0, (cdfn_h1*d_corr);
    
    
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




// some differences to R due to trivariate CDF approximation
struct TVBS_Matrix_two TVBSo_grad_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor)
{
  Eigen::VectorXd mu_3_0(3), mu_2_0(2);
  mu_3_0.setZero();
  mu_2_0.setZero();
  double tri_var_cdf = TVBS_pmvnorm_old_cpp(w, mu_3_0, cor);
  double bi_var_cdf = TVBS_pmvnorm_old_cpp(w.head(2), mu_2_0, cor.block(0,0,2,2));
  
  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(w.head(2), cor.block(0,0,2,2));
  double g_w_1 = output_3.d1;
  double g_w_2 = output_3.d2;
  double g_rho = output_3.d3;
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_tvn_cpp(w, cor);
  Eigen::VectorXd g_w_tri = output_2.mat1;
  Eigen::VectorXd g_rho_tri = output_2.mat2;
  
  Eigen::VectorXd g_w_12(2);
  g_w_12 << g_w_1, g_w_2;
  Eigen::VectorXd g_w1_w2_new = (bi_var_cdf*g_w_tri.head(2).array() - tri_var_cdf*g_w_12.array())/pow(bi_var_cdf,2);
  
  
  double g_rho12_new = (bi_var_cdf*g_rho_tri(0) - tri_var_cdf*g_rho)/pow(bi_var_cdf,2);
  
  Eigen::VectorXd output_final_1(3), output_final_2(3);
  output_final_1 << g_w1_w2_new, g_w_tri(2)/bi_var_cdf;
  output_final_2 << g_rho12_new, g_rho_tri.segment(1,2).array()/bi_var_cdf;
  
  struct TVBS_Matrix_two output_final = {output_final_1, output_final_2};
  
  return output_final;
}




// numerical issues might arrise
struct TVBS_Matrix_three TVBSo_grad_non_cdf_tvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  // cap the x values over 6 and under -6

  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  Eigen::MatrixXd Sigma_diag_sqrt_inv(3,3);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_tvn_by_cdf_bvn_cpp(x_norm.head(3), Cor_mat.block(0,0,3,3));
  Eigen::VectorXd g_w = output_2.mat1;
  Eigen::VectorXd g_rho = output_2.mat2;
  
  struct TVBS_Matrix_two output_3 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  Eigen::MatrixXd g_b_cor_cov = output_3.mat1;
  Eigen::MatrixXd g_omega_cor_cov = output_3.mat2;
  
  
  Eigen::VectorXd g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  Eigen::VectorXd g_cov = g_b_cor_cov*g_w + g_omega_cor_cov*g_rho;
  Eigen::VectorXd g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  
  return output_final;
}





// // // // // // // // //
//  quadrovariate       //
// // // // // // // // //


// issues arrise with approximated trivariate cdf
struct TVBS_Matrix_two TVBSo_grad_cdf_qvn_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Corr)
{
  // initiate probability vector
  Eigen::VectorXd p(2);
  p.setZero();
  
  // include reordering here if wanted
  Eigen::VectorXd x_temp = x_norm;
  Eigen::MatrixXd Corr_temp = Corr;
  Eigen::VectorXd mu_temp(x_norm.size());
  mu_temp.setZero();
  
  // calculate first part of probability using trivariate normal distribution
  p(0) = TVBS_pmvnorm_old_cpp(x_temp.head(3), mu_temp.head(3), Corr_temp.block(0,0,3,3));
  
  // truncate on first two variables
  struct TVBS_Matrix_two output_2 = TVBS_mu_l_trunc_bivariate_cpp(mu_temp, Corr_temp, x_temp.head(2));
  Eigen::VectorXd mu_temp_1 = output_2.mat1;
  Eigen::MatrixXd sigma_temp_1 = output_2.mat2;
  
  
  // calculate second part by calculating P(4)=P(3,4)/P(3)
  double p_2_1 = TVBS_pmvnorm_old_cpp(x_temp.segment(2,2), mu_temp_1.segment(2,2), sigma_temp_1.block(2,2,2,2));
  double p_2_2 = TVBS_std_normal_cdf_cpp((x_temp(2)-mu_temp_1(2))/sqrt(sigma_temp_1(2,2)));
  p(1) = p_2_1/p_2_2;
  
  // update global condition numbers
  condcovsigtrunc_global = 0;
  condcovmeantrunc_global = 0;
  
  struct TVBS_Matrix_four output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(2,2), mu_temp, Corr_temp, x_temp.head(2));
  Eigen::MatrixXd g_y = output_4.mat1;
  Eigen::MatrixXd g_mu_mean = output_4.mat2;
  Eigen::MatrixXd g_x_mean = output_4.mat3;
  Eigen::MatrixXd g_c_mean = output_4.mat4;
  
  
  struct TVBS_Matrix_three output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Corr_temp, x_temp.head(2));
  Eigen::MatrixXd g_mu_cov = output_3.mat1;
  Eigen::MatrixXd g_x_cov = output_3.mat2;
  Eigen::MatrixXd g_c_cov = output_3.mat3;
  
  
  
  Eigen::MatrixXd g_cu_mu_l_mu_sig(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
  g_cu_mu_l_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
  g_cu_mu_l_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
  
  Eigen::MatrixXd g_cu_mu_l_c(g_c_mean.rows()+2,5);
  g_cu_mu_l_c.setZero();
  g_cu_mu_l_c.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
  g_cu_mu_l_c.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
  
  
  Eigen::VectorXd g_c(4);
  g_c.setZero();
  
  output_3 = TVBS_grad_non_cdf_bvn_by_cdfn_cpp(mu_temp_1.segment(2,2), sigma_temp_1.block(2,2,2,2), x_temp.segment(2,2));
  Eigen::VectorXd g_from_p_mu = output_3.mat1;
  Eigen::VectorXd g_from_p_cov = output_3.mat2;
  Eigen::VectorXd g_from_p_c = output_3.mat3;
  
  
  Eigen::VectorXd g_rho_rho = p(0)*(g_cu_mu_l_mu_sig.block(0,0,g_cu_mu_l_mu_sig.rows(),2)*g_from_p_mu + g_cu_mu_l_mu_sig.block(0,2,g_cu_mu_l_mu_sig.rows(),3)*g_from_p_cov);
  
  
  g_c = p(0)*(g_cu_mu_l_c.block(0,0,g_cu_mu_l_c.rows(),2)*g_from_p_mu + g_cu_mu_l_c.block(0,2,g_cu_mu_l_c.rows(),3)*g_from_p_cov);
  
  // getting gradient contribution directly from noncdfbvn for absiccae except the first two
  g_c.segment(2,2) = p(0)*g_from_p_c;
  
  output_2 = TVBSo_grad_cdf_tvn_cpp(x_temp.segment(0,3), Corr_temp.block(0,0,3,3));
  Eigen::VectorXd g_c_1 = output_2.mat1;
  Eigen::VectorXd g_rho_1 = output_2.mat2;
  
  Eigen::VectorXd g_rho_rho_add(6);
  g_rho_rho_add << g_rho_1.head(2), 0, g_rho_1(2), 0, 0;
  g_rho_rho += p(1)*g_rho_rho_add;
  
  g_c.head(g_c_1.size()) += p(1)*g_c_1;
  // restore order, if reordered earlier, according to temp1
  
  Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
  g_rho_rho_mat += g_rho_rho_mat.transpose();
  
  // reorder the g_rho_rho_mat according to temp1
  
  g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat, 0);
  
  struct TVBS_Matrix_two output_final = {g_c, g_rho_rho};
  return output_final;
}





struct TVBS_Matrix_two TVBSo_grad_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd w, Eigen::MatrixXd cor)
{
  Eigen::VectorXd mu_temp_4(4);
  mu_temp_4.setZero();

  double quad_var_cdf = TVBS_pmvnorm_old_cpp(w.head(4), mu_temp_4, cor.block(0,0,4,4));
  double bi_var_cdf = TVBS_pmvnorm_old_cpp(w.head(2), mu_temp_4.head(2), cor.block(0,0,2,2));
  
  
  struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(w.head(2),cor.block(0,0,2,2));
  double g_w_1 = output_3.d1;
  double g_w_2 = output_3.d2;
  double g_rho = output_3.d3;
  Eigen::VectorXd g_w(2);
  g_w << g_w_1, g_w_2;
  
  struct TVBS_Matrix_two output_2 = TVBS_grad_cdf_qvn_cpp(w, cor);
  Eigen::VectorXd g_w_quad = output_2.mat1;
  Eigen::VectorXd g_rho_quad = output_2.mat2;
  
  
  Eigen::VectorXd g_w1_w2_new = (bi_var_cdf*g_w_quad.head(2) - quad_var_cdf*g_w)/pow(bi_var_cdf,2);
  double g_rho12_new = (bi_var_cdf*g_rho_quad(0) - quad_var_cdf*g_rho)/pow(bi_var_cdf,2);
  
  Eigen::VectorXd output_final_1(4), output_final_2(6);
  output_final_1 << g_w1_w2_new, g_w_quad.segment(2,2)/bi_var_cdf;
  output_final_2 << g_rho12_new, g_rho_quad.segment(1,5)/bi_var_cdf;
  
  struct TVBS_Matrix_two output_final = {output_final_1, output_final_2};
  return output_final;
}





// strange behaviour with cholesky_global==1
struct TVBS_Matrix_three TVBSo_grad_non_cdf_qvn_by_cdf_bvn_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x)
{
  // normalize x and calculate correlation matrix from covariance matrix
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Sigma_diag_sqrt_inv(4,4);
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  
  // cap the x values over 6 and under -6
  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_qvn_by_cdf_bvn_cpp(x_norm.head(4), Cor_mat.block(0,0,4,4));
  Eigen::VectorXd g_w = output_2.mat1;
  Eigen::VectorXd g_rho = output_2.mat2;
  
  output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
  Eigen::MatrixXd g_b_cor_cov = output_2.mat1;
  Eigen::MatrixXd g_omega_cor_cov = output_2.mat2;
  
  Eigen::VectorXd g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  Eigen::VectorXd g_cov = g_b_cor_cov*g_w + g_omega_cor_cov*g_rho;
  Eigen::VectorXd g_x = -1.0*g_mu;
  
  if(cholesky_global==1)
  {
    g_cov = TVBS_g_cholesky_cov_cpp(cov)*g_cov;
  }
  
  
  struct TVBS_Matrix_three output_final = {g_mu, g_cov, g_x};
  
  return output_final;
  
}





struct TVBS_Vector_three TVBSo_pdf_mvna_tvbs_cpp(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{
  // determine dimension of the problem
  int m = x_norm.size();
  
  // calculate number of iterations
  int k_tilde = std::max((m-1)/2,1);
  
  // initialize variables for the output
  Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
  p.setZero();
  p_out.setZero();
  g_c.setZero();
  g_rho_rho.setZero();
  p_out.setZero();
  
  // initialize mu_zero variable
  Eigen::VectorXd mu_zero(m);
  mu_zero.setZero();
  
  
  // // // // // // // // // // //
  // reorder the dimensions     //
  // // // // // // // // // // //
  
  struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
  
  Eigen::MatrixXd Cor_mat_temp = output_reorder.mat1;        // we'll work with this one
  Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
  
  Eigen::VectorXd x_temp = output_reorder.vec1;
  
  Eigen::VectorXd temp1 = output_reorder.vec2;
  
  Eigen::VectorXd mu_temp(m);
  mu_temp.setZero();
  
  // // // // // // // // // // //
  // end of reordering          //
  // // // // // // // // // // //
  
  
  
  // calculate initial LDLT decomposition
  struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
  Eigen::MatrixXd L = output_2.mat1;
  Eigen::MatrixXd D = output_2.mat2;
  
  
  // deal with low dimensional cases directly
  if(m==2)
  {
    struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(x_temp, Cor_mat_temp.block(0,0,2,2));
    g_c << output_3.d1, output_3.d2;
    g_rho_rho(0) = output_3.d3;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==3)
  {
    struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_tvn_cpp(x_temp, Cor_mat_temp.block(0,0,3,3));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==4)
  {
    struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_qvn_cpp(x_temp, Cor_mat_temp.block(0,0,4,4));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
  // if the output should be for log
  if(log_out==1)
  {
    g_c /= p_out(0);
    g_rho_rho /= p_out(0);
    double p = log(p_out(0));
    p_out(0) = p;
  }
  }
  else 
  {
    // initiate some variables
    struct TVBS_Matrix_two output_2;
    struct TVBS_Matrix_three output_3;
    struct TVBS_Matrix_four output_4;
    Eigen::VectorXd mu_tilde;
    Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
    double p_k1_nomi, p_k1_denomi;
    Eigen::VectorXd g_from_p_mu, g_from_p_cov, g_from_p_c;
    
    
    // calculate P1 using the first four variables
    p(0) = TVBS_pmvnorm_old_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    Eigen::VectorXd g_cu_mul_from_pc(m);
    g_cu_mul_from_pc.setZero();
    
    
    // start here from k1=1 to make indexing easier
    for(int k1=1; k1<k_tilde; k1++)
    {
      // claculate the truncated mean and variance
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
      mu_tilde = output_2.mat1;
      omega = output_2.mat2;
      
      
      // set global variables according to the case (if they have already been truncated or not)
      if(k1==1)
      {
        condcovsigtrunc_global = 0;
        condcovmeantrunc_global = 0;
      }
      else
      {
        condcovsigtrunc_global = 1;
        condcovmeantrunc_global = 1;
      }
      
      
      output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(m-2*k1, m-2*k1), mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_y = output_4.mat1;
      g_mu_mean = output_4.mat2;
      g_x_mean = output_4.mat3;
      g_c_mean = output_4.mat4;
      
      
      output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_mu_cov = output_3.mat1;
      g_x_cov = output_3.mat2;
      g_c_cov = output_3.mat3;
      
      
      // update g_cu_mul_mu_sig and g_cu_mu_lc
      if(k1==1)
      {
        g_cu_mul_mu_sig.resize(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
        g_cu_mul_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        
        g_cu_mu_lc.resize(g_c_mean.rows()+m-2,g_c_mean.cols()+g_c_cov.cols());
        g_cu_mu_lc.setZero();
        g_cu_mu_lc.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
        
      }
      else if(k1>1)
      {
        g_cu_mul_mu_sig_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(), g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mul_mu_sig_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mul_mu_sig_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mul_mu_sig = g_cu_mul_mu_sig*g_cu_mul_mu_sig_rhs;
        
        // update g_cu_mu_lc
        g_cu_mu_lc_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(),g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mu_lc_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mu_lc_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mu_lc = g_cu_mu_lc*g_cu_mu_lc_rhs;
        
        g_cu_mu_lc.block((2*(k1-1)+1)-1,0,2,g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block((2*(k1-1)+1)-1,g_c_mean.cols(),2,g_c_cov.cols()) = g_c_cov;
      }
      
      
      // update pi_kp1
      mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
      
      
      // update the covariance matrix LDLT decomposition
      output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
      L = output_2.mat1;
      D = output_2.mat2;
      Cor_mat_temp = L*D*(L.transpose());
      
      
      if(m>=2*k1+4)
      {
        // regular case
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        
        // update gradients
        output_3 = TVBSo_grad_non_cdf_qvn_by_cdf_bvn_cpp(mu_temp.head(4), Cor_mat_temp.block(0,0,4,4), x_temp.segment((2*k1+1)-1,4));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 4+3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),4) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),4);
        g_rho_rho_add.block(0,4,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,4+3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(3*m-6*k1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,4+3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,(4*m-8*k1-2)-1,g_cu_mul_mu_sig.rows(),1);

        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),4)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,4) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),4+3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),4) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),4);
        g_c_add.block(0,4,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,4+3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(3*m-6*k1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,4+3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,(4*m-8*k1-2)-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
      }
      else
      {
        // last iteration, if uneven dimension
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        
        
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        
        // update gradients
        
        // // // 
        // Errors of e-6 magnitude occur here
        // // //
        output_3 = TVBSo_grad_non_cdf_tvn_by_cdf_bvn_cpp(mu_temp.head(3), Cor_mat_temp.block(0,0,3,3), x_temp.segment((2*k1+1)-1,3));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,3*m-6*k1-1,g_cu_mul_mu_sig.rows(),1);
        
        // // // 
        // Errors of e-4 magnitude occur here
        // // //
        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),3)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,3) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,3*m-6*k1-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
        
      }
      
    }
    
    
    
    output_2 = TVBSo_grad_cdf_qvn_cpp(x_temp.head(4), Cor_mat_ordered.block(0,0,4,4));
    Eigen::VectorXd g_w = output_2.mat1;
    Eigen::VectorXd g_rho = output_2.mat2;
    
    
    // adding contribution of rho12,rho13, and rho14 from initial cdfqvn function
    g_rho_rho.head(3) += g_rho.head(3)/p(0);
    
    
    // adding contribution of rho23 & rho24 from initial cdfqvn function
    g_rho_rho.segment(m-1,2) += g_rho.segment(3,2)/p(0);
    
    
    // adding contribution of rho34 from initial cdfqvn function
    g_rho_rho(2*m-2-1) += g_rho(5)/p(0);
    
    
    // adding contribution of all abscissa originating from the probability function
    g_c += g_cu_mul_from_pc;
    
    
    // inserting gradient contribution of first four abscissae directly from the cdfqvn function
    g_c.head(4) += g_w/p(0);
    
    
    if(log_out == 0)
    {
      p_out(0) = p.prod();
      g_c = (p_out(0)*g_c);
      g_rho_rho = (p_out(0)*g_rho_rho);
    }
    else
    {
      p_out(0) = p.array().log().sum();
    }
  }
  
  
  
  
  // // // // // // // // // // // // // // //
  // invert the reordering                  //
  // // // // // // // // // // // // // // //
  
  // g_rho_rho has to be transformed into matrix form, reordered and then back into vectorization
  
  Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
  Eigen::MatrixXd g_rho_rho_mat_new = g_rho_rho_mat + g_rho_rho_mat.transpose();
  
  //Eigen::MatrixXd g_rho_rho_mat_new(m,m);
  //g_rho_rho_mat_new.setZero();
  Eigen::VectorXd g_c_new(m);
  g_c_new.setZero();
  
  for(int i_row=0; i_row<m; i_row++)
  {
    g_c_new(temp1(i_row)) = g_c(i_row);
    for(int i_col=0; i_col<m; i_col++)
    {
      g_rho_rho_mat(temp1(i_row),temp1(i_col)) = g_rho_rho_mat_new(i_row,i_col);
    }
  }
  g_c = g_c_new;
  g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat,0);
  
  
  // // // // // // // // // // // // // // //
  // end of reorder invertion               //
  // // // // // // // // // // // // // // //
  
  struct TVBS_Vector_three output_final = {p_out, g_c, g_rho_rho};
  return output_final;
}






struct TVBS_Vector_four TVBSo_pdf_mvn_analytic_cpp(Eigen::VectorXd mu, Eigen::MatrixXd cov, Eigen::VectorXd x, int log_out = 0)
{
  // initialize output variables
  Eigen::VectorXd p, g_mu, g_cov, g_x; 
  
  // calculate normalized x vector
  Eigen::VectorXd Sigma_diag_sqrt = cov.diagonal().array().sqrt();
  Eigen::VectorXd x_norm = (x.array()-mu.array())/Sigma_diag_sqrt.array();
  
  
  // if x values are further away from mu than 6 sd => P is in {0,1}
  x_norm = x_norm.array().min(6);
  x_norm = x_norm.array().max(-6);
  
  
  // calculate correlation matrix from covariance matrix
  Eigen::MatrixXd Sigma_diag_sqrt_inv(Sigma_diag_sqrt.size(),Sigma_diag_sqrt.size());
  Sigma_diag_sqrt_inv.setZero();
  Sigma_diag_sqrt_inv.diagonal() = 1.0/Sigma_diag_sqrt.array();
  Eigen::MatrixXd Cor_mat = Sigma_diag_sqrt_inv*cov*Sigma_diag_sqrt_inv;
  Cor_mat.diagonal().setOnes();
  
  
  // adjust global variable
  if(optimal_global==1)
  {
    optimal_global = 0;
  }
  
  
  // get probability and gradient of the normalized problem
  struct TVBS_Vector_three output_3 = TVBSo_pdf_mvna_tvbs_cpp(x_norm, Cor_mat, log_out);
  p = output_3.v1;
  Eigen::VectorXd g_w = output_3.v2;
  Eigen::VectorXd g_rho = output_3.v3;
  
  
  
  // calculate the derivatives of the COVARIANCE matrix from the CORRELATION matrix
  if(covarr_global==1)
  {
    struct TVBS_Matrix_two output_2 = TVBS_grad_cor_cov_cpp(x_norm, Sigma_diag_sqrt, Cor_mat);
    Eigen::MatrixXd g_b_cor_cov = output_2.mat1;
    Eigen::MatrixXd g_omega_cor_cov = output_2.mat2;
    
    g_cov = g_b_cor_cov*g_w + g_omega_cor_cov*g_rho;
  }
  else
  {
    g_cov = g_rho;
  }
  
  
  // resize the gradient of the mean and x vector
  g_mu = -1.0*g_w.array()/Sigma_diag_sqrt.array();
  g_x = -1.0*g_mu;
  

  struct TVBS_Vector_four output_final = {p, g_mu, g_cov, g_x};
  return output_final;
}






// // // // // // // // //
//  TVBS + grad norm    //
// // // // // // // // //

//' TVBS approximation to multivariate Gaussian CDF. 
//' @description
//' The function computes the  CDF for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS) as well as the gradient. 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @param log_out
//' integer; probability or log of probability? 
//' @return 
//' double; (log of) probability. 
//' @keywords internal
//'
//[[Rcpp::export]]
Eigen::VectorXd TVBS(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{
  // determine dimension of the problem
  int m = x_norm.size();
  
  // calculate number of iterations
  int k_tilde = std::max((m-1)/2,1);
  
  
  // initialize variables for the output
  Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
  p.setZero();
  p_out.setZero();
  g_c.setZero();
  g_rho_rho.setZero();
  p_out.setZero();
  
  // initialize mu_zero variable
  Eigen::VectorXd mu_zero(m);
  mu_zero.setZero();
  
  
  
  
  // // // // // // // // // // //
  // reorder the dimensions     //
  // // // // // // // // // // //
  
  struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
  
  Eigen::MatrixXd Cor_mat_temp = output_reorder.mat1;        // we'll work with this one
  Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
  
  Eigen::VectorXd x_temp = output_reorder.vec1;
  
  Eigen::VectorXd temp1 = output_reorder.vec2;
  
  Eigen::VectorXd mu_temp(m);
  mu_temp.setZero();
  
  // // // // // // // // // // //
  // end of reordering          //
  // // // // // // // // // // //

  
  
  // calculate initial LDLT decomposition
  struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
  Eigen::MatrixXd L = output_2.mat1;
  Eigen::MatrixXd D = output_2.mat2;
  
  
  // deal with low dimensional cases directly
  if(m==2)
  {
    struct TVBS_double_three output_3 = TVBS_grad_cdf_bvn_cpp(x_temp, Cor_mat_temp.block(0,0,2,2));
    g_c << output_3.d1, output_3.d2;
    g_rho_rho(0) = output_3.d3;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==3)
  {
    struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_tvn_cpp(x_temp, Cor_mat_temp.block(0,0,3,3));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else if(m==4)
  {
    struct TVBS_Matrix_two output_2 = TVBSo_grad_cdf_qvn_cpp(x_temp, Cor_mat_temp.block(0,0,4,4));
    g_c = output_2.mat1;
    g_rho_rho = output_2.mat2;
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    // if the output should be for log
    if(log_out==1)
    {
      g_c /= p_out(0);
      g_rho_rho /= p_out(0);
      double p = log(p_out(0));
      p_out(0) = p;
    }
  }
  else 
  {

    // initiate some variables
    struct TVBS_Matrix_two output_2;
    struct TVBS_Matrix_three output_3;
    struct TVBS_Matrix_four output_4;
    Eigen::VectorXd mu_tilde;
    Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
    double p_k1_nomi, p_k1_denomi;
    Eigen::VectorXd g_from_p_mu, g_from_p_cov, g_from_p_c;
    
    
    // calculate P1 using the first four variables
    p(0) = TVBS_pmvnorm_old_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    
    Eigen::VectorXd g_cu_mul_from_pc(m);
    g_cu_mul_from_pc.setZero();
    
    
    // start here from k1=1 to make indexing easier
    for(int k1=1; k1<k_tilde; k1++)
    {
      // claculate the truncated mean and variance
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
      mu_tilde = output_2.mat1;
      omega = output_2.mat2;
      
      
      // set global variables according to the case (if they have already been truncated or not)
      if(k1==1)
      {
        condcovsigtrunc_global = 0;
        condcovmeantrunc_global = 0;
      }
      else
      {
        condcovsigtrunc_global = 1;
        condcovmeantrunc_global = 1;
      }
      
      
      output_4 = TVBS_g_cond_mean_trunc_cpp(Eigen::MatrixXd::Identity(m-2*k1, m-2*k1), mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_y = output_4.mat1;
      g_mu_mean = output_4.mat2;
      g_x_mean = output_4.mat3;
      g_c_mean = output_4.mat4;
      
      
      output_3 = TVBS_g_cond_cov_trunc_cpp(mu_temp, Cor_mat_temp, x_temp.segment((2*k1-1)-1,2));
      g_mu_cov = output_3.mat1;
      g_x_cov = output_3.mat2;
      g_c_cov = output_3.mat3;
      
      
      // update g_cu_mul_mu_sig and g_cu_mu_lc
      if(k1==1)
      {
        g_cu_mul_mu_sig.resize(g_x_mean.rows(), g_x_mean.cols()+g_x_cov.cols());
        g_cu_mul_mu_sig.block(0,0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig.block(0,g_x_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        
        g_cu_mu_lc.resize(g_c_mean.rows()+m-2,g_c_mean.cols()+g_c_cov.cols());
        g_cu_mu_lc.setZero();
        g_cu_mu_lc.block(0,0,g_c_mean.rows(),g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block(0,g_c_mean.cols(),g_c_cov.rows(),g_c_cov.cols()) = g_c_cov;
        
      }
      else if(k1>1)
      {
        g_cu_mul_mu_sig_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(), g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mul_mu_sig_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mul_mu_sig_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mul_mu_sig_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mul_mu_sig = g_cu_mul_mu_sig*g_cu_mul_mu_sig_rhs;
        
        // update g_cu_mu_lc
        g_cu_mu_lc_rhs.resize(g_mu_mean.rows()+g_x_mean.rows(),g_mu_mean.cols()+g_mu_cov.cols());
        g_cu_mu_lc_rhs.block(0,0,g_mu_mean.rows(),g_mu_mean.cols()) = g_mu_mean;
        g_cu_mu_lc_rhs.block(0,g_mu_mean.cols(),g_mu_cov.rows(),g_mu_cov.cols()) = g_mu_cov;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),0,g_x_mean.rows(),g_x_mean.cols()) = g_x_mean;
        g_cu_mu_lc_rhs.block(g_mu_mean.rows(),g_mu_mean.cols(),g_x_cov.rows(),g_x_cov.cols()) = g_x_cov;
        g_cu_mu_lc = g_cu_mu_lc*g_cu_mu_lc_rhs;
        
        g_cu_mu_lc.block((2*(k1-1)+1)-1,0,2,g_c_mean.cols()) = g_c_mean;
        g_cu_mu_lc.block((2*(k1-1)+1)-1,g_c_mean.cols(),2,g_c_cov.cols()) = g_c_cov;
      }
      
      
      // update pi_kp1
      mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
      
      
      // update the covariance matrix LDLT decomposition
      output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
      L = output_2.mat1;
      D = output_2.mat2;
      Cor_mat_temp = L*D*(L.transpose());
      
      
      if(m>=2*k1+4)
      {
        // regular case
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        
        // update gradients
        output_3 = TVBSo_grad_non_cdf_qvn_by_cdf_bvn_cpp(mu_temp.head(4), Cor_mat_temp.block(0,0,4,4), x_temp.segment((2*k1+1)-1,4));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 4+3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),4) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),4);
        g_rho_rho_add.block(0,4,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,4+3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(3*m-6*k1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,4+3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,(4*m-8*k1-2)-1,g_cu_mul_mu_sig.rows(),1);
        
        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),4)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,4) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),4+3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),4) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),4);
        g_c_add.block(0,4,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,4+3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(3*m-6*k1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,4+3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,(4*m-8*k1-2)-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
        
      }
      else
      {
        // last iteration, if uneven dimension
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        
        p(k1) = p_k1_nomi/p_k1_denomi;
        
        // update gradients
        
        // // // 
        // Errors of e-6 magnitude occur here
        // // //
        output_3 = TVBSo_grad_non_cdf_tvn_by_cdf_bvn_cpp(mu_temp.head(3), Cor_mat_temp.block(0,0,3,3), x_temp.segment((2*k1+1)-1,3));
        g_from_p_mu = output_3.mat1;
        g_from_p_cov = output_3.mat2;
        g_from_p_c = output_3.mat3;
        
        
        // update g_rho_rho
        g_rho_rho_add.resize(g_cu_mul_mu_sig.rows(), 3+2+1);
        g_rho_rho_add.block(0,0,g_rho_rho_add.rows(),3) = g_cu_mul_mu_sig.block(0,(m-2*k1+1)-1,g_cu_mul_mu_sig.rows(),3);
        g_rho_rho_add.block(0,3,g_rho_rho_add.rows(),2) = g_cu_mul_mu_sig.block(0,(2*m-4*k1+1)-1,g_cu_mul_mu_sig.rows(),2);
        g_rho_rho_add.block(0,3+2,g_rho_rho_add.rows(),1) = g_cu_mul_mu_sig.block(0,3*m-6*k1-1,g_cu_mul_mu_sig.rows(),1);
        
        // // // 
        // Errors of e-4 magnitude occur here
        // // //
        g_rho_rho += 1.0/p(k1)*( g_cu_mul_mu_sig.block(0,0,g_cu_mul_mu_sig.rows(),3)*g_from_p_mu +  g_rho_rho_add*g_from_p_cov);
        
        
        // update g_cu_mul_from_pc
        g_cu_mul_from_pc.segment((2*k1+1)-1,3) += g_from_p_c/p(k1);
        
        
        // update g_c
        g_c_add.resize(g_cu_mu_lc.rows(),3+2+1);
        g_c_add.block(0,0,g_cu_mu_lc.rows(),3) = g_cu_mu_lc.block(0,(m-2*k1+1)-1,g_cu_mu_lc.rows(),3);
        g_c_add.block(0,3,g_cu_mu_lc.rows(),2) = g_cu_mu_lc.block(0,(2*m-4*k1+1)-1,g_cu_mu_lc.rows(),2);
        g_c_add.block(0,3+2,g_cu_mu_lc.rows(),1) = g_cu_mu_lc.block(0,3*m-6*k1-1,g_cu_mu_lc.rows(),1);
        
        g_c += 1.0/p(k1)*( g_cu_mu_lc.block(0,0,g_cu_mu_lc.rows(),3)*g_from_p_mu + g_c_add*g_from_p_cov );
      }
      
    }
    
    
    
    output_2 = TVBSo_grad_cdf_qvn_cpp(x_temp.head(4), Cor_mat_ordered.block(0,0,4,4));
    Eigen::VectorXd g_w = output_2.mat1;
    Eigen::VectorXd g_rho = output_2.mat2;
    
    
    // adding contribution of rho12,rho13, and rho14 from initial cdfqvn function
    g_rho_rho.head(3) += g_rho.head(3)/p(0);
    
    
    // adding contribution of rho23 & rho24 from initial cdfqvn function
    g_rho_rho.segment(m-1,2) += g_rho.segment(3,2)/p(0);
    
    
    // adding contribution of rho34 from initial cdfqvn function
    g_rho_rho(2*m-2-1) += g_rho(5)/p(0);
    
    
    // adding contribution of all abscissa originating from the probability function
    g_c += g_cu_mul_from_pc;
    
    
    // inserting gradient contribution of first four abscissae directly from the cdfqvn function
    g_c.head(4) += g_w/p(0);
    
    
    if(log_out == 0)
    {
      p_out(0) = p.prod();
      g_c = (p_out(0)*g_c);
      g_rho_rho = (p_out(0)*g_rho_rho);
    }
    else
    {
      p_out(0) = p.array().log().sum();
    }
    
  }
  
  
  
  
  // // // // // // // // // // // // // // //
  // invert the reordering                  //
  // // // // // // // // // // // // // // //
  
  // g_rho_rho has to be transformed into matrix form, reordered and then back into vectorization
  
  Eigen::MatrixXd g_rho_rho_mat = TVBS_vec_2_upper_diag_cpp(g_rho_rho, 0);
  Eigen::MatrixXd g_rho_rho_mat_new = g_rho_rho_mat + g_rho_rho_mat.transpose();

  Eigen::VectorXd g_c_new(m);
  g_c_new.setZero();
  
  for(int i_row=0; i_row<m; i_row++)
  {
    g_c_new(temp1(i_row)) = g_c(i_row);
    for(int i_col=0; i_col<m; i_col++)
    {
      g_rho_rho_mat(temp1(i_row),temp1(i_col)) = g_rho_rho_mat_new(i_row,i_col);
    }
  }
  g_c = g_c_new;
  g_rho_rho = TVBS_lower_tri_entries_cpp(g_rho_rho_mat,0);
  
  // // // // // // // // // // // // // // //
  // end of reorder invertion               //
  // // // // // // // // // // // // // // //
  
  
  
  Eigen::VectorXd output_final(1+g_c.size()+g_rho_rho.size()); 
  output_final << g_c, g_rho_rho, p_out(0);
  return output_final;
}








// // // // // // // // //
//  TVBS - prob only    //
// // // // // // // // //
//' TVBS approximation to multivariate Gaussian CDF. 
//' @description
//' The function computes the  CDF for a multivariate Gaussian distribution according to the method of Chandra Bhat (TVBS). 
//' @param x_norm
//' nx1 vector of point to evaluate the CDF at. 
//' @param Cor_mat
//' nxn correlation matrix.
//' @param log_out
//' integer; probability or log of probability? 
//' @return 
//' double; (log of) probability. 
//' @keywords internal
//'
// [[Rcpp::export]]
double TVBS_p(Eigen::VectorXd x_norm, Eigen::MatrixXd Cor_mat, int log_out = 0)
{
  // determine dimension of the problem
  int m = x_norm.size();
  
  // calculate number of iterations
  int k_tilde = std::max((m-1)/2,1);
  
  
  // initialize variables for the output
  Eigen::VectorXd p(k_tilde), p_out(1), g_c(m), g_rho_rho((m*(m-1))/2);
  p.setZero();
  p_out.setZero();
  g_c.setZero();
  g_rho_rho.setZero();
  p_out.setZero();
  
  // initialize mu_zero variable
  Eigen::VectorXd mu_zero(m);
  mu_zero.setZero();
  
  
  
  
  // // // // // // // // // // //
  // reorder the dimensions     //
  // // // // // // // // // // //
  
  struct TVBS_Vec_Mat_Vec output_reorder = TVBS_reorder_sig_cpp(x_norm, Cor_mat);
  
  Eigen::MatrixXd Cor_mat_temp = output_reorder.mat1;        // we'll work with this one
  Eigen::MatrixXd Cor_mat_ordered = Cor_mat_temp;            // we'll keep this as the original one
  
  Eigen::VectorXd x_temp = output_reorder.vec1;
  
  Eigen::VectorXd temp1 = output_reorder.vec2;
  
  Eigen::VectorXd mu_temp(m);
  mu_temp.setZero();
  
  // // // // // // // // // // //
  // end of reordering          //
  // // // // // // // // // // //
  
  
  
  
  // calculate initial LDLT decomposition
  struct TVBS_Matrix_two output_2 = TVBS_LDLT_decomp_cpp(Cor_mat_temp);
  Eigen::MatrixXd L = output_2.mat1;
  Eigen::MatrixXd D = output_2.mat2;
  
  
  
  // deal with low dimensional cases directly
  if(m==2)
  {
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(2), Cor_mat_temp.block(0,0,2,2));
    
    // if the output should be for log
    if(log_out==1)
    {
      double p_temp = log(p_out(0));
      p_out(0) = p_temp;
    }
  }
  else if(m==3)
  {
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(3), Cor_mat_temp.block(0,0,3,3));
    
    // if the output should be for log
    if(log_out==1)
    {
      double p_temp = log(p_out(0));
      p_out(0) = p_temp;
    }
  }
  else if(m==4)
  {
    p_out(0) = TVBS_pmvnorm_old_cpp(x_temp, mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    // if the output should be for log
    if(log_out==1)
    {
      double p_temp = log(p_out(0));
      p_out(0) = p_temp;
    }
  }
  else 
  {
    // initiate some variables
    struct TVBS_Matrix_two output_2;
    struct TVBS_Matrix_three output_3;
    struct TVBS_Matrix_four output_4;
    Eigen::VectorXd mu_tilde;
    Eigen::MatrixXd omega, g_y, g_mu_mean, g_x_mean, g_c_mean, g_mu_cov, g_x_cov, g_c_cov, g_cu_mul_mu_sig, g_cu_mu_lc, g_cu_mul_mu_sig_rhs, g_cu_mu_lc_rhs, g_rho_rho_add, g_c_add;
    double p_k1_nomi, p_k1_denomi;

    
    // calculate P1 using the first four variables
    p(0) = TVBS_pmvnorm_old_cpp(x_temp.head(4), mu_zero.head(4), Cor_mat_temp.block(0,0,4,4));
    
    
    // start here from k1=1 to make indexing easier
    for(int k1=1; k1<k_tilde; k1++)
    {
      // claculate the truncated mean and variance
      output_2 = TVBS_truncate_bi_normal_cpp(mu_temp.head(2), D.block(0,0,2,2), x_temp.segment((2*k1-1)-1,2),tol);
      mu_tilde = output_2.mat1;
      omega = output_2.mat2;
      
      
      // set global variables according to the case (if they have already been truncated or not)
      if(k1==1)
      {
        condcovsigtrunc_global = 0;
        condcovmeantrunc_global = 0;
      }
      else
      {
        condcovsigtrunc_global = 1;
        condcovmeantrunc_global = 1;
      }
      
      
      
      // update pi_kp1
      mu_temp = mu_temp.segment(2,m-2*k1) + L.block(2,0,m-2*k1,2)*(mu_tilde-mu_temp.head(2));
      
      
      // update the covariance matrix LDLT decomposition
      output_2 = TVBS_LDLT_update_cpp(L,D,omega,2);
      L = output_2.mat1;
      D = output_2.mat2;
      Cor_mat_temp = L*D*(L.transpose());
      
      
      if(m>=2*k1+4)
      {
        // regular case
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,4), mu_temp.head(4), Cor_mat_temp.block(0,0,4,4));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,4));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        p(k1) = p_k1_nomi/p_k1_denomi;
        
      }
      else
      {
        // last iteration, if uneven dimension
        
        // calculate next probability
        p_k1_nomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,3), mu_temp.head(3), Cor_mat_temp.block(0,0,3,3));
        p_k1_nomi = std::max(p_k1_nomi, pow(tol,3));
        p_k1_denomi = TVBS_pmvnorm_old_cpp(x_temp.segment((2*k1+1)-1,2), mu_temp.head(2), Cor_mat_temp.block(0,0,2,2));
        p_k1_denomi = std::max(p_k1_denomi, pow(tol,2));
        
        
        p(k1) = p_k1_nomi/p_k1_denomi;
        
      }
      
    }
    
    
    if(log_out == 0)
    {
      p_out(0) = p.prod();
    }
    else
    {
      p_out(0) = p.array().log().sum();
    }
    
    
  }
  
  double p_final = p_out(0);
  
  return p_final;
}

