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

#include "toms462.h"                                          // allows bivariate normal cdf calculations


double tol = 0.000000001;
double bound_max = 6.0;


// [[Rcpp::plugins(cpp11)]]

double std_normal_cdf(double x)                          // create a function for the standard normal cumulative distribution function
{
  double m = 0;                                                     // define mean == 0
  double s = 1;    // define standard deviation sigma == 1
  
  if (x>bound_max) { x = bound_max;}
  if (x<-bound_max) { x = -bound_max;}
  return 0.5 * (1 + erf((x - m) / (s * sqrt(2.))));            // computes 1/2 * (1 + int_0^(x/sqrt(2)) exp(-t^2) dt = CDF of standard normal dist.
}


double normal_cdf(double x, double v)                          // create a function for the standard normal cumulative distribution function
{
  double s = sqrt(v);                                                     // define standard deviation sigma == 1
  double m = x/s;
  if (m>bound_max) { m = bound_max;}
  if (m<-bound_max) { m = -bound_max;}
  
  return 0.5 * (1 + erf(m / (sqrt(2.))));            // computes 1/2 * (1 + int_0^(x/sqrt(2)) exp(-t^2) dt = CDF of standard normal dist.
}


double std_normal_pdf(double x)                                 // create a function for the standard normal probability distribution function
{
  double m = 0;                                                     // define mean == 0
  double s = 1;                                                     // define standard deviation sigma == 1
  static const double inv_sqrt_2pi = 0.3989422804014327;            // save 1/(sqrt(2*pi)) as the constant inv_sqrt_2pi
  double a = (x - m) / s;                                           // save a = (x - mu)/sd <- substitution/normalization
  
  if (a>bound_max) { a = bound_max;}
  if (a<-bound_max) { a = -bound_max;}
  
  return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);                 // compute PDF 1/sqrt(2+pi+sigma^2)*exp(-(x-mu)/(2*sigma^2))
}


double normal_pdf(double x, double v)                                 // create a function for the standard normal probability distribution function
{
  double m = 0;                                                     // define mean == 0
  double s = sqrt(v);                                                     // define standard deviation sigma == 1
  static const double inv_sqrt_2pi = 0.3989422804014327;            // save 1/(sqrt(2*pi)) as the constant inv_sqrt_2pi
  double a = (x - m) / s;                                           // save a = (x - mu)/sd <- substitution/normalization
  
  if (a>bound_max) { a = bound_max;}
  if (a<-bound_max) { a = -bound_max;}
  
  return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);                 // compute PDF 1/sqrt(2+pi+sigma^2)*exp(-(x-mu)/(2*sigma^2))
}

double std_normal_cdf_inv_1(double y)
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

Eigen::VectorXd std_normal_cdf_inv(Eigen::VectorXd y)
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



double TVBS_biv_std_norm_pdf(Eigen::VectorXd x, double rho)                 // create a function for the bivariate standard normal probability distribution function
{
  double pi_global              = 3.14159265358979323846;
  
  if (x(0)>bound_max) { x(0) = bound_max;}
  if (x(0)<-bound_max) { x(0) = -bound_max;}
  
  if (x(1)>bound_max) { x(1) = bound_max;}
  if (x(1)<-bound_max) { x(1) = -bound_max;}
  
  if (rho>1-tol) {rho= 1-tol;}
  if (rho< tol-1) {rho = tol-1;}
  
  return exp(-(pow(x(0),2) - 2*rho*x(0)*x(1) + pow(x(1),2))/(2*(1.0-pow(rho,2))))/( 2*pi_global*sqrt(1.0-pow(rho,2)) );
}

//' biv_normal_cdf
//' @description
//' The function computes the bivariate Gaussian CDF. 
//' @param w0
//' double; x-coordinate 
//' @param w1
//' double; y-coordinate 
//' @param rho
//' double; correlation
//' @return 
//' double; cdf 
//' @keywords internal
//'
// [[Rcpp::export]]
double biv_normal_cdf(double x0, double x1, double r12)                                 // create a function for the standard normal probability distribution function
{
  
  double cd = 0;
  if (r12>1-tol){ r12 = 1-tol;}
  if (r12<-tol-1){ r12 =tol-1;}
  if (x0>bound_max){ x0 = bound_max;}
  if (x1>bound_max){ x1 = bound_max;}
  if (x0<-bound_max){ x0 = -bound_max;}
  if (x1<-bound_max){ x1 = -bound_max;}

  cd = bivnor(-x0, -x1, r12);     
  // 
  if (cd < tol){ cd = tol;}
  return cd;
  
}

//' biv_normal_pdf
//' @description
//' The function computes the bivariate Gaussian PDF. 
//' @param w0
//' double; x-coordinate 
//' @param w1
//' double; y-coordinate 
//' @param rho
//' double; correlation
//' @return 
//' double; pdf 
//' @keywords internal
//'
// [[Rcpp::export]]
double biv_normal_pdf(double x0, double x1, double r12)                                 // create a function for the standard normal probability distribution function
{
  static const double inv_2pi = 0.1591549;  /* 1/(2pi) */ 
    double pd = 0;
    if (r12>1-tol){ r12 = 1-tol;}
    if (r12<-tol-1){ r12 =tol-1;}
    if (x0>bound_max){ x0 = bound_max;}
    if (x1>bound_max){ x1 = bound_max;}
    if (x0<-bound_max){ x0 = -bound_max;}
    if (x1<-bound_max){ x1 = -bound_max;}

    pd =  inv_2pi/sqrt(1-r12*r12)*std::exp(-0.5*(x0*x0 + x1*x1 - 2*x0*x1*r12)/(1-r12*r12)); 
    // 
      return pd; 
}


Eigen::VectorXd grad_cdf(double x0, double x1, double r12)                                 // create a function for the standard normal probability distribution function
{
  Eigen::VectorXd grad(3);
  grad.setZero(); 
  
  // initialize 
  if (r12>1-tol){ r12 = 1-tol;}
  if (r12<-tol-1){ r12 =tol-1;}
  if (x0>bound_max){ x0 = bound_max;}
  if (x1>bound_max){ x1 = bound_max;}
  if (x0<-bound_max){ x0 = -bound_max;}
  if (x1<-bound_max){ x1 = -bound_max;}
  
  // derivatives 
  grad(0) = std_normal_pdf(x0)*normal_cdf(x1-r12*x0,1-r12*r12);
  grad(1) = std_normal_pdf(x1)*normal_cdf(x0-r12*x1,1-r12*r12);
  grad(2) = biv_normal_pdf(x0,x1,r12);
  // 
    return grad;   
}


Eigen::MatrixXd Hess_cdf(double x0, double x1, double r12)                                 // create a function for the standard normal probability distribution function
{
  //static const double inv_2pi = 0.1591549;  /* 1/(2pi) */ 
    Eigen::MatrixXd Hess(3,3);
  Hess.setZero(); 
  
  // initialize 
  if (r12>1-tol){ r12 = 1-tol;}
  if (r12<-tol-1){ r12 =tol-1;}
  if (x0>bound_max){ x0 = bound_max;}
  if (x1>bound_max){ x1 = bound_max;}
  if (x0<-bound_max){ x0 = -bound_max;}
  if (x1<-bound_max){ x1 = -bound_max;}
  
  // gradient needed?
    Eigen::VectorXd gr(3);
  gr.setZero();
  gr = grad_cdf(x0,x1,r12);
  
  // Hessian 
  Hess(0,0) = gr(0)*(-x0) + biv_normal_pdf(x0,x1,r12)*(-r12);
  Hess(0,1) = std_normal_pdf(x0)*normal_pdf(x1-r12*x0,1-r12*r12);
  Hess(0,2) = gr(2)*(-x0+x1*r12)/(1-r12*r12);
  Hess(1,0)=Hess(0,1);
  Hess(1,1) = gr(1)*(-x1) + biv_normal_pdf(x0,x1,r12)*(-r12);
  Hess(1,2) = gr(2)*(-x1+x0*r12)/(1-r12*r12);
  Hess(2,0) = Hess(0,2);
  Hess(2,1)= Hess(1,2);
  // last with respect to two times the corr.
  double onemr2 = 1-r12*r12;
  double he = x0*x1*onemr2 - (x0*x0+x1*x1 - 2*r12*x0*x1)*r12;
  Hess(2,2) = biv_normal_pdf(x0,x1,r12)*(he/(onemr2*onemr2)+r12/onemr2);
  
  // 
    return Hess;   
}



double biv_gen_pdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)                 // create a function for the bivariate standard normal probability distribution function
{
  double pi     = 3.14159265358979323846;
  
  double expo = b.transpose() * Sigma.inverse() * b;
  double deter = Sigma(0,0)*Sigma(1,1) - Sigma(0,1)*Sigma(1,0);
  return exp(-expo/2)/( 2*pi*sqrt(deter));
}

//' Hess_pdf
//' @description
//' The function computes the Hessian of the bivariate Gaussian CDF. 
//' @param w0
//' double; x-coordinate 
//' @param w1
//' double; y-coordinate 
//' @param rho
//' double; correlation
//' @return 
//' matrix; Hessian of pdf 
//' @keywords internal
//'
// [[Rcpp::export]]
Eigen::MatrixXd Hess_pdf(double w0, double w1, double rho)                                 // create a function for the standard normal probability distribution function
{
  //static const double inv_2pi = 0.1591549;  /* 1/(2pi) */ 
    Eigen::MatrixXd Hess(3,3);
  Hess.setZero(); 
  
  // initialize 
  if (rho>1-tol){ rho = 1-tol;}
  if (rho<-tol-1){ rho =tol-1;}
  if (w0>bound_max){ w0 = bound_max;}
  if (w1>bound_max){ w1 = bound_max;}
  if (w0<-bound_max){ w0 = -bound_max;}
  if (w1<-bound_max){ w1 = -bound_max;}
  
  double phi2 = biv_normal_pdf(w0,w1,rho);
  
  // gradient
  Eigen::VectorXd grad(3);
  grad.setZero();
  
  // deriv 
  double exppo = w0*w0 + w1*w1 - rho*w0*w1*2;
  grad(0) = -phi2*(w0-rho*w1)/(1-rho*rho);
  grad(1) = -phi2*(w1-rho*w0)/(1-rho*rho);
  grad(2) = phi2*(rho/(1-rho*rho) - (-w0*w1*(1-rho*rho) + exppo*rho)/pow(1-rho*rho,2));
  
  // Hessian 
  Hess(0,0) = -grad(0)*(w0-rho*w1)/(1-rho*rho) - phi2/(1-rho*rho);
  Hess(0,1) = -grad(1)*(w0-rho*w1)/(1-rho*rho) + phi2*rho/(1-rho*rho);
  Hess(1,0) = Hess(0,1);
  Hess(1,1) = -grad(1)*(w1-rho*w0)/(1-rho*rho) - phi2/(1-rho*rho);
  Hess(0,2) = -grad(2)*(w0-rho*w1)/(1-rho*rho) - phi2*(-w1*(1-rho*rho)+(w0-rho*w1)*2*rho)/pow(1-rho*rho,2);
  Hess(2,0)= Hess(0,2); 
  Hess(1,2) = -grad(2)*(w1-rho*w0)/(1-rho*rho) - phi2*(-w0*(1-rho*rho)+(w1-rho*w0)*2*rho)/pow(1-rho*rho,2);
  Hess(2,1) = Hess(1,2); 
  Hess(2,2) = grad(2)*(rho/(1-rho*rho) - (-w0*w1*(1-rho*rho) + exppo*rho)/pow(1-rho*rho,2));
  double f = -w0*w1*(1-rho*rho) + exppo*rho;
  double fstrich = w0*w1*2*rho -2*w0*w1*rho + exppo;
  Hess(2,2) += phi2*((1+rho*rho)/pow(1-rho*rho,2)-(fstrich*(1-rho*rho)+4*f*rho)/pow(1-rho*rho,3));
  
  return Hess;   
}

////////////////////
//// bivariate   //
//////////////////
  
// Jacobian of normalized (w1,w2,rho) with respect to (b1,b2,sigmax2,sigmay2,sigmaxy) 
Eigen::MatrixXd J_wR_bSig(Eigen::VectorXd b, Eigen::MatrixXd Sigma)    
{
  double sigmax2 = Sigma(0,0);
  double sigmay2 = Sigma(1,1);
  
  double sigmax = sqrt(sigmax2);
  double sigmay = sqrt(sigmay2);
  
  double rho = Sigma(0,1)/(sigmax*sigmay);
  
  Eigen::MatrixXd J_r(3,5);
  J_r.setZero();
  
  J_r(0,0)=1/sigmax;
  J_r(1,1)=1/sigmay;
  J_r(0,2)= -b(0)/(2*sigmax*sigmax2);
  J_r(1,3)= -b(1)/(2*sigmay*sigmay2);
  J_r(2,2)= -rho/(2*sigmax2);
  J_r(2,3)= -rho/(2*sigmay2);
  J_r(2,4) = 1/(sigmax*sigmay);
  
  // return Jacobian 
  return J_r;
}


//' general bivariate normal cdf calculation 
//' @description
//' The function computes the bivariate Gaussian CDF. 
//' @param b
//' vector; point to evaluate CDF 
//' @param Sigma
//' 2x2 matrix; variance matrix.
//' @return 
//' double; CDF
//' @keywords internal 
// [[Rcpp::export]]
double biv_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  double sigmax2 = Sigma(0,0);
  double sigmay2 = Sigma(1,1);
  
  double sigmax = sqrt(sigmax2);
  double sigmay = sqrt(sigmay2);
  
  double rho = Sigma(0,1)/(sigmax*sigmay);
  
  double x0 = b(0)/sigmax;
  double x1 = b(1)/sigmay; 
  
  double p = biv_normal_cdf(x0,x1,rho);
  return p;
}

//' general bivariate normal gradient cdf calculation 
//' @description
//' The function computes the gradient of the bivariate Gaussian CDF. 
//' @param b
//' vector; point to evaluate CDF 
//' @param Sigma
//' 2x2 matrix; variance matrix.
//' @return 
//' vector; gradient of CDF
//' @keywords internal 
// [[Rcpp::export]]
Eigen::VectorXd grad_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  double sigmax2 = Sigma(0,0);
  double sigmay2 = Sigma(1,1);
  
  double sigmax = sqrt(sigmax2);
  double sigmay = sqrt(sigmay2);
  
  double rho = Sigma(0,1)/(sigmax*sigmay);
  
  double x0 = b(0)/sigmax;
  double x1 = b(1)/sigmay; 
  
  Eigen::VectorXd gr = grad_cdf(x0,x1,rho); // calculate gradient 
  Eigen::MatrixXd J = J_wR_bSig(b, Sigma); // calculate Jacobian
  
  Eigen::VectorXd grad = gr.transpose() * J; 
  
  return grad;
}

//' general bivariate normal Hessian of cdf calculation 
//' @description
//' The function computes the Hessian of the bivariate Gaussian CDF. 
//' @param b
//' vector; point to evaluate CDF 
//' @param Sigma
//' 2x2 matrix; variance matrix.
//' @return 
//' 5x5 matrix; Hessian of CDF 
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd Hess_gen_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  double sigmax2 = Sigma(0,0);
  double sigmay2 = Sigma(1,1);
  
  double sigmax = sqrt(sigmax2);
  double sigmay = sqrt(sigmay2);
  
  double rho = Sigma(0,1)/(sigmax*sigmay);
  
  double x0 = b(0)/sigmax;
  double x1 = b(1)/sigmay; 
  
  Eigen::VectorXd gr = grad_cdf(x0,x1,rho); // calculate gradient for normalized version
  Eigen::MatrixXd J = J_wR_bSig(b, Sigma); // calculate Jacobian
  
  Eigen::MatrixXd Hess(5,5); 
  Eigen::MatrixXd H = Hess_cdf(x0,x1,rho); // calculate Hessian for normalized version
  
  // twice derivatives in normalized cdf
  Hess = J.transpose() * H * J; 
  
  // add gradient times derviative of Jacobian.   
  // b_1
  Hess(0,2) -= gr(0) / (2 * sigmax2 * sigmax);
  
  // b_2
  Hess(1,3) -= gr(1) /(2*  sigmay2 * sigmay);
  Hess(2,0) = Hess(0,2);
  Hess(3,1) = Hess(1,3);
  
  // sigmaxy
  Hess(4,2) -= gr(2) /(2*  sigmay * sigmax2 * sigmax);
  Hess(4,3) -= gr(2) /(2*  sigmay * sigmay2 * sigmax);
  
  Hess(2,4)=Hess(4,2);
  Hess(3,4)=Hess(4,3);
  
  // sigmax2 
  Hess(2,2) += ((gr(0) * b(0)) + (gr(2)*rho * sigmax)) * 3/(4*pow(sigmax,5));
  
  Hess(2,3) += gr(2) * rho/(4*sigmax2 * sigmay2);
  Hess(3,2)=Hess(2,3);
  
  // sigmay2 
  Hess(3,3) += ((gr(1) * b(1)) + (gr(2)*rho * sigmay)) * 3/(4*pow(sigmay,5));
  
  
  return Hess;
}


/////////////////////////////////////7
/// trivariate    ///
//////////////////////
  
  double tvn_gen_pdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)                 // create a function for the bivariate standard normal probability distribution function
{
  double pi     = 3.14159265358979323846;
  
  double expo = b.transpose() * Sigma.inverse() * b;
  double deter = Sigma.determinant();
  return exp(-expo/2)/( sqrt(pow(2*pi,3)*deter));
}

//' general trivariate normal cdf calculation 
//' @description
//' The function computes the gradient of the trivariate Gaussian CDF. 
//' @param b
//' vector; point to evaluate CDF 
//' @param Sigma
//' 3x3 matrix; variance matrix.
//' @return 
//' 9x1; gradient of the CDF 
//' @keywords internal 
// [[Rcpp::export]]
Eigen::VectorXd grad_gen_tvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  Eigen::VectorXd grad(9);
  grad.setZero();
  
  // derivative with respect to b_i. 
  
  Eigen::VectorXd bh=b;
  Eigen::MatrixXd Sigmah = Sigma;
  
  for (int i=0;i<3;i++){
    double f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
    Eigen::MatrixXd R21 = Sigmah.block(1,0,2,1) / Sigmah(0,0); 
    Eigen::VectorXd mub = bh.block(1,0,2,1) - R21 * bh(0);
    Eigen::MatrixXd Sigb = Sigmah.block(1,1,2,2) - R21 * Sigmah(0,0) * R21.transpose();  
    double f2 = biv_gen_cdf( mub, Sigb); 
    
    grad(i) = f1* f2 /sqrt(Sigmah(0,0));
    
    // rotate 
    Eigen::MatrixXd ROT(3,3);
    ROT.setZero();
    
    ROT(0,1) =1;
    ROT(1,2) =1;
    ROT(2,0)=1;
    
    bh = ROT * bh;
    Sigmah = ROT*Sigmah * ROT.transpose();
    
  }
  
  // derivative with respect to rho_{kl}
  Eigen::MatrixXd Gr(3,3);
  Gr.setZero();
  
  Eigen::MatrixXd D(3,3);
  D.setZero();
  for (int d=0;d<3;d++){ 
    D(d,d)= sqrt(Sigma(d,d));
  }
  bh=D.inverse()* b;
  Sigmah = D.inverse()* Sigma * D.inverse();
  
  
  for (int a=0;a<3;a++){
    for (int b=a+1;b<3;b++){
      
      
      Eigen::MatrixXd Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,1) * Sigmah.block(2,0,1,2)/ Sigmah(2,2);
      Eigen::VectorXd bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,1) / Sigmah(2,2) *bh(2);
      double f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
      
      Eigen::MatrixXd Sh = Sigmah.block(0,0,2,2);
      Eigen::MatrixXd R12 = Sigmah.block(2,0,1,2) * Sh.inverse(); 
      double mub = bh(2) - (R12(0) * bh(0)) - (R12(1)* bh(1));
      Eigen::MatrixXd Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
      double f2 = std_normal_cdf(mub/sqrt(Sigmah(2,2) - Sigb(0,0)));
      
      
      
      Gr(a,b) = f1* f2;
      
      
      // rotate second entry to the back. 
      Eigen::MatrixXd ROT(3,3);
      ROT.setZero();
      
      ROT(0,0) =1;
      ROT(1,2) =1;
      ROT(2,1)=1;
      
      bh = ROT * bh;
      Sigmah = ROT*Sigmah * ROT.transpose();
      
      
    }
    // rotate first entry to last entry. 
    Eigen::MatrixXd ROT(3,3);
    ROT.setZero();
    
    ROT(0,1) =1;
    ROT(1,2) =1;
    ROT(2,0)=1;
    bh = ROT * bh;
    Sigmah = ROT*Sigmah * ROT.transpose();
  }
  
  
  
  // fill in results into grad    
  grad(4) = Gr(0,1) / (D(0,0)* D(1,1)); // entry sigma(1,2)
  grad(5) = Gr(0,2) / (D(0,0)* D(2,2)); // entry sigma(1,3)
  grad(7) = Gr(1,2) / (D(1,1)* D(2,2)); // entry sigma(2,3)
  
  // finally diagonal entries
  grad(3) = (-grad(0)*b(0)-grad(4)*Sigma(0,1) - grad(5)*Sigma(0,2))/(2*Sigma(0,0)); 
  grad(6) = (-grad(1)*b(1)-grad(4)*Sigma(0,1) - grad(7)*Sigma(1,2))/(2*Sigma(1,1)); 
  grad(8) = (-grad(2)*b(2)-grad(7)*Sigma(1,2) - grad(5)*Sigma(0,2))/(2*Sigma(2,2)); 
  
  return grad;
}  

//' general trivariate normal cdf calculation 
//' @description
//' The function computes the Hessian of the trivariate Gaussian CDF. 
//' @param b
//' vector; point to evaluate CDF 
//' @param Sigma
//' 3x3 matrix; variance matrix.
//' @return 
//' 9x9; Hessian of the CDF 
//' @keywords internal 
// [[Rcpp::export]]
Eigen::MatrixXd Hessian_gen_tvn_cdf(Eigen::VectorXd b, Eigen::MatrixXd Sigma)
{
  // initialize 
  Eigen::MatrixXd Hess(9,9);
  Hess.setZero();
  Eigen::VectorXd bh = b;
  Eigen::MatrixXd Sigmah = Sigma;
  
  
  // cycle over rows. 
  // first with respect to b1
  
  double f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  double f1sqrt = f1/sqrt(Sigmah(0,0));
  
  Eigen::MatrixXd R21 = Sigmah.block(1,0,2,1) / Sigmah(0,0); 
  Eigen::VectorXd mub = bh.block(1,0,2,1) - R21 * bh(0);
  Eigen::MatrixXd Sigb = Sigmah.block(1,1,2,2) - R21 * Sigmah(0,0) * R21.transpose();  
  Eigen::VectorXd gr = grad_gen_cdf(mub, Sigb);
  double f2 = biv_gen_cdf( mub, Sigb); 
  
  // second deriv w.r.t. bi 
  Hess(0,0) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2))/Sigmah(0,0);
  
  Hess(0,1) = f1sqrt*(gr(0));
  Hess(0,2) = f1sqrt*(gr(1));
  
  // derivative with respect to sigma_{11}. 
  Eigen::VectorXd dSigb_dSig(5);
  dSigb_dSig.setZero();
  
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(3) = Sigmah(0,2)*Sigmah(0,2);
  
  
  double pr = (gr.transpose() * dSigb_dSig);
  Hess(0,3) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  
  // derivative w.r.t sigma_{1,2} 
  Hess(0,4) = -f1sqrt*( gr(0)*bh(0)+ gr(2)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,3}
  Hess(0,5) =  -f1sqrt*( gr(1)*bh(0) + gr(3)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{2,2}
  Hess(0,6)= f1sqrt*gr(2);
  
  // derivative w.r.t sigma_{2,3} 
  Hess(0,7) = f1sqrt*gr(4);
  
  // derivative w.r.t sigma_{3,3}
  Hess(0,8) = f1sqrt*gr(3);
  
  
  // next up: derivativves with respect to b2. 
  Eigen::MatrixXd ROT(3,3);
  ROT.setZero();
  
  ROT(0,1) =1;
  ROT(1,2) =1;
  ROT(2,0)=1;
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new ordering: (1,2,0).
  f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  f1sqrt = f1/sqrt(Sigmah(0,0));
  
  R21 = Sigmah.block(1,0,2,1) / Sigmah(0,0); 
  mub = bh.block(1,0,2,1) - R21 * bh(0);
  Sigb = Sigmah.block(1,1,2,2) - R21 * Sigmah(0,0) * R21.transpose();  
  gr = grad_gen_cdf(mub, Sigb);
  f2 = biv_gen_cdf( mub, Sigb); 
  
  // second deriv w.r.t. bi 
  Hess(1,1) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2))/Sigmah(0,0);
  
  Hess(1,2) = f1sqrt*(gr(0));
  Hess(1,0) = f1sqrt*(gr(1));
  
  // derivative with respect to sigma_{22}. 
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(3) = Sigmah(0,2)*Sigmah(0,2);
  
  
  pr = (gr.transpose() * dSigb_dSig);
  Hess(1,6) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  
  // derivative w.r.t sigma_{2,3} 
  Hess(1,7) = -f1sqrt*( gr(0)*bh(0)+ gr(2)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{2,1}
  Hess(1,4) =  -f1sqrt*( gr(1)*bh(0) + gr(3)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{3,3}
  Hess(1,8)= f1sqrt*gr(2);
  
  // derivative w.r.t sigma_{3,1} 
  Hess(1,5) = f1sqrt*gr(4);
  
  // derivative w.r.t sigma_{1,1}
  Hess(1,3) = f1sqrt*gr(3);
  
  
  // derivative with respect to b3. 
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new ordering: (1,2,0).
  f1 = std_normal_pdf(bh(0)/sqrt(Sigmah(0,0)));
  f1sqrt = f1/sqrt(Sigmah(0,0));
  
  R21 = Sigmah.block(1,0,2,1) / Sigmah(0,0); 
  mub = bh.block(1,0,2,1) - R21 * bh(0);
  Sigb = Sigmah.block(1,1,2,2) - R21 * Sigmah(0,0) * R21.transpose();  
  gr = grad_gen_cdf(mub, Sigb);
  f2 = biv_gen_cdf( mub, Sigb); 
  
  // second deriv w.r.t. b3 
  Hess(2,2) = -f1sqrt*f2*bh(0)/(Sigmah(0,0)) - f1sqrt*(gr(0)*Sigmah(0,1)+gr(1)*Sigmah(0,2))/Sigmah(0,0);
  
  Hess(2,0) = f1sqrt*(gr(0));
  Hess(2,1) = f1sqrt*(gr(1));
  
  // derivative with respect to sigma_{33}. 
  dSigb_dSig(0)= bh(0)*Sigmah(0,1);
  dSigb_dSig(1)= bh(0)*Sigmah(0,2);
  dSigb_dSig(2) = Sigmah(0,1)*Sigmah(0,1);
  dSigb_dSig(4) = Sigmah(0,2)*Sigmah(0,1);
  dSigb_dSig(3) = Sigmah(0,2)*Sigmah(0,2);
  
  
  pr = (gr.transpose() * dSigb_dSig);
  Hess(2,8) = f1sqrt*f2*bh(0)*bh(0)/(2*Sigmah(0,0)*Sigmah(0,0)) + f1sqrt/pow(Sigmah(0,0),2)* pr - f1sqrt*f2/(2*Sigmah(0,0));
  
  // derivative w.r.t sigma_{3,1} 
  Hess(2,5) = -f1sqrt*( gr(0)*bh(0)+ gr(2)*2*Sigmah(0,1)+gr(4)*Sigmah(0,2))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{3,2}
  Hess(2,7) =  -f1sqrt*( gr(1)*bh(0) + gr(3)*2*Sigmah(0,2)+gr(4)*Sigmah(0,1))/Sigmah(0,0);
  
  // derivative w.r.t. sigma_{1,1}
  Hess(2,3)= f1sqrt*gr(2);
  
  // derivative w.r.t sigma_{1,2} 
  Hess(2,4) = f1sqrt*gr(4);
  
  // derivative w.r.t sigma_{2,2}
  Hess(2,6) = f1sqrt*gr(3);
  
  // use symmetry of Hessian
  
  Eigen::MatrixXd He = Hess + Hess.transpose();
  
  He.block(0,0,3,3) = He.block(0,0,3,3)/2; 
  
  ////////////////////////////////////////
  /// done with derivatives involving b //
  ////////////////////////////////////////
    
  // derivative with respect to rho_{kl}
  Eigen::MatrixXd Gr(3,3);
  Gr.setZero();
  
  Eigen::MatrixXd D(3,3);
  D.setZero();
  for (int d=0;d<3;d++){ 
    D(d,d)= sqrt(Sigma(d,d));
  }
  bh=D.inverse()* b;
  Sigmah = D.inverse()* Sigma * D.inverse();
  
  // order (1,2,3)
  Eigen::MatrixXd Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,1) * Sigmah.block(2,0,1,2)/ Sigmah(2,2);
  Eigen::VectorXd bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,1) / Sigmah(2,2) *bh(2);
  f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
  
  Eigen::MatrixXd Sh = Sigmah.block(0,0,2,2);
  Eigen::MatrixXd R12 = Sigmah.block(2,0,1,2) * Sh.inverse(); 
  double mu_b = bh(2) - (R12(0) * bh(0)) - (R12(1)* bh(1));
  Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
  f2 = std_normal_cdf(mu_b/sqrt(Sigmah(2,2) - Sigb(0,0)));
  
  Gr(0,1) = f1* f2;
  
  // now start Hessian calculations. 
  // f1 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  Eigen::MatrixXd df1(1,6);
  df1.setZero();
  
  double detSigma = Sigmah(0,0)*Sigmah(1,1) - Sigmah(0,1)*Sigmah(1,0);
  double bexp = bh(0)*bh(0)*Sigmah(1,1) + bh(1)*bh(1)*Sigmah(0,0) - 2*bh(0)*bh(1)*Sigmah(0,1);
  
  Eigen::MatrixXd ddetSigma(1,6);
  ddetSigma.setZero();
  
  ddetSigma(0,0)= Sigmah(1,1);
  ddetSigma(0,1)= -2*Sigmah(0,1);
  ddetSigma(0,3) = Sigmah(0,0);
  
  df1 = f1*(-1/(2*detSigma)+bexp/(2*detSigma*detSigma)) * ddetSigma;
  df1(0,0) += (-f1)*(bh(1)*bh(1)/(2*(detSigma)));
  df1(0,1) += (f1)*(bh(0)*bh(1)/((detSigma)));
  df1(0,3) += (-f1)*(bh(0)*bh(0)/(2*(detSigma)));
  
  
  // f2 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  double sig_bar = Sigmah(2,2) - Sigb(0,0); 
  Eigen::MatrixXd df2(1,6);
  df2.setZero();
  
  Eigen::MatrixXd db_bar(1,6);
  db_bar.setZero();
  
  Eigen::MatrixXd dsig_bar(1,6);
  dsig_bar.setZero();
  
  Eigen::MatrixXd dR_12(2,6);
  dR_12.setZero();
  
  dR_12(0,0) = -(Sigmah(0,2)*Sigmah(1,1)-Sigmah(1,2)*Sigmah(0,1))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(0,1) = -2*dR_12(0,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(0,3) = dR_12(0,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(0,1) += (-1/detSigma)*Sigmah(1,2);
  dR_12(0,2) += Sigmah(1,1)/detSigma;
  dR_12(0,3) += Sigmah(0,2)/detSigma;
  dR_12(0,4) += (-1/detSigma)*Sigmah(0,1);
  
  
  dR_12(1,0)= -(-Sigmah(0,2)*Sigmah(0,1)+Sigmah(1,2)*Sigmah(0,0))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(1,1) = -2*dR_12(1,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(1,3) = dR_12(1,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(1,0) += Sigmah(1,2)/detSigma;
  dR_12(1,1) += -Sigmah(0,2)/detSigma;
  dR_12(1,2) += -Sigmah(0,1)/detSigma;
  dR_12(1,4) += Sigmah(0,0)/detSigma;
  
  
  
  db_bar = -bh(0)*dR_12.block(0,0,1,6) - bh(1)* dR_12.block(1,0,1,6);
  dsig_bar.block(0,0,1,6) = - Sigmah(0,2)*dR_12.block(0,0,1,6) - Sigmah(1,2)*dR_12.block(1,0,1,6);
  dsig_bar(0,2) += -R12(0,0);
  dsig_bar(0,4) += -R12(1,0);
  
  dsig_bar(0,5) = 1; 
  
  df2 = normal_pdf(mu_b,sig_bar)*(db_bar - mu_b*dsig_bar/(2*sig_bar)); 
  
  
  
  // combine the pieces 
  // grad = f1*f2 / sqrt(sigma_11 *sigma_22)
  He.block(4,3,1,6) = (df1*f2 + f1*df2)/(D(0,0)*D(1,1));
  
  He(4,4) = He(4,4)/(D(0,0)*D(1,1));
  He(4,5) = He(4,5)/(D(0,0)*D(2,2));
  He(4,7) = He(4,7)/(D(2,2)*D(1,1));
  // derivative w.r.t. sigma_11 reconstructed from others. 
  He(4,3) = -(He(4,0)*b(0) + He(4,4)*Sigma(0,1)+ He(4,5)*Sigma(0,2)+ (Gr(0,1) / (D(0,0)* D(1,1))))/(2*Sigma(0,0));
  He(4,6) = -(He(4,1)*b(1) + He(4,4)*Sigma(0,1)+ He(4,7)*Sigma(1,2)+ (Gr(0,1) / (D(0,0)* D(1,1))))/(2*Sigma(1,1));
  
  He(4,8) = He(4,8)/Sigma(2,2);
  
  ///////////////////////////////////////////////
    //// derivative w.r.t. sigma_13             ///
    ///////////////////////////////////////////////
    
    // rotate second entry to the back. 
  ROT.setZero();
  
  ROT(0,0) =1;
  ROT(1,2) =1;
  ROT(2,1)=1;
  
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  
  // new order: (1,3,2) 
  
  Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,1) * Sigmah.block(2,0,1,2)/ Sigmah(2,2);
  bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,1) / Sigmah(2,2) *bh(2);
  f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
  
  Sh = Sigmah.block(0,0,2,2);
  R12 = Sigmah.block(2,0,1,2) * Sh.inverse(); 
  mu_b = bh(2) - (R12(0) * bh(0)) - (R12(1)* bh(1));
  Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
  f2 = std_normal_cdf(mu_b/sqrt(Sigmah(2,2) - Sigb(0,0)));
  
  Gr(0,2) = f1* f2;
  
  // now start Hessian calculations. 
  // f1 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  df1.setZero();
  
  detSigma = Sigmah(0,0)*Sigmah(1,1) - Sigmah(0,1)*Sigmah(1,0);
  bexp = bh(0)*bh(0)*Sigmah(1,1) + bh(1)*bh(1)*Sigmah(0,0) - 2*bh(0)*bh(1)*Sigmah(0,1);
  
  ddetSigma.setZero();
  
  ddetSigma(0,0)= Sigmah(1,1);
  ddetSigma(0,1)= -2*Sigmah(0,1);
  ddetSigma(0,3) = Sigmah(0,0);
  
  df1 = f1*(-1/(2*detSigma)+bexp/(2*detSigma*detSigma)) * ddetSigma;
  df1(0,0) += (-f1)*(bh(1)*bh(1)/(2*(detSigma)));
  df1(0,1) += (f1)*(bh(0)*bh(1)/((detSigma)));
  df1(0,3) += (-f1)*(bh(0)*bh(0)/(2*(detSigma)));
  
  // f2 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  sig_bar = Sigmah(2,2) - Sigb(0,0); 
  df2.setZero();
  db_bar.setZero();
  dsig_bar.setZero();
  dR_12.setZero();
  
  dR_12(0,0) = -(Sigmah(0,2)*Sigmah(1,1)-Sigmah(1,2)*Sigmah(0,1))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(0,1) = -2*dR_12(0,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(0,3) = dR_12(0,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(0,1) += (-1/detSigma)*Sigmah(1,2);
  dR_12(0,2) += Sigmah(1,1)/detSigma;
  dR_12(0,3) += Sigmah(0,2)/detSigma;
  dR_12(0,4) += (-1/detSigma)*Sigmah(0,1);
  
  dR_12(1,0)= -(-Sigmah(0,2)*Sigmah(0,1)+Sigmah(1,2)*Sigmah(0,0))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(1,1) = -2*dR_12(1,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(1,3) = dR_12(1,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(1,0) += Sigmah(1,2)/detSigma;
  dR_12(1,1) += -Sigmah(0,2)/detSigma;
  dR_12(1,2) += -Sigmah(0,1)/detSigma;
  dR_12(1,4) += Sigmah(0,0)/detSigma;
  
  
  db_bar = -bh(0)*dR_12.block(0,0,1,6) - bh(1)* dR_12.block(1,0,1,6);
  dsig_bar.block(0,0,1,6) = - Sigmah(0,2)*dR_12.block(0,0,1,6) - Sigmah(1,2)*dR_12.block(1,0,1,6);
  dsig_bar(0,2) += -R12(0,0);
  dsig_bar(0,4) += -R12(1,0);
  
  dsig_bar(0,5) = 1; 
  
  df2 = normal_pdf(mu_b,sig_bar)*(db_bar - mu_b*dsig_bar/(2*sig_bar)); 
  
  
  
  // combine the pieces 
  // grad = f1*f2 / sqrt(sigma_11 *sigma_33)
  Eigen::VectorXd gr_df(6);
  gr_df.setZero();
  gr_df = (df1*f2 + f1*df2)/(D(0,0)*D(2,2));
  
  He(5,3)= df1(0)*f2 + f1*df2(0);
  He(5,4)= df1(2)*f2 + f1*df2(2);
  He(5,5)= df1(1)*f2 + f1*df2(1);
  He(5,6) = df1(5)*f2 + f1*df2(5);
  He(5,7)= df1(4)*f2 + f1*df2(4);
  He(5,8) = df1(3)*f2 + f1*df2(3);
  
  He.block(5,3,1,6)=He.block(5,3,1,6)/(D(0,0)*D(2,2));
  
  He(5,4) = He(5,4)/(D(0,0)*D(1,1));
  He(5,5) = He(5,5)/(D(0,0)*D(2,2));
  He(5,7) = He(5,7)/(D(2,2)*D(1,1));
  // derivative w.r.t. sigma_11 reconstructed from others. 
  He(5,3) = -(He(5,0)*b(0) + He(5,4)*Sigma(0,1)+ He(5,5)*Sigma(0,2)+ (Gr(0,2) / (D(0,0)* D(2,2))))/(2*Sigma(0,0));
  He(5,8) = -(He(5,2)*b(2) + He(5,5)*Sigma(0,2)+ He(5,7)*Sigma(1,2)+ (Gr(0,2) / (D(0,0)* D(2,2))))/(2*Sigma(2,2));
  
  He(5,6) = He(5,6)/Sigma(1,1);
  
  
  ////////////////////////////////////
    /// derivative w.r.t. sigma_23.
  ////////////////////////////////////
    
    // rotate for last ordering 
  ROT.setZero();
  
  ROT(2,0) =1;
  ROT(1,1) =1;
  ROT(0,2)=1;
  
  bh = ROT * bh;
  Sigmah = ROT*Sigmah * ROT.transpose();
  // new ordering now: (2,3,1).   
  
  Sigbb = Sigmah.block(0,0,2,2) - Sigmah.block(0,2,2,1) * Sigmah.block(2,0,1,2)/ Sigmah(2,2);
  bbb = bh.block(0,0,2,1) - Sigmah.block(0,2,2,1) / Sigmah(2,2) *bh(2);
  f1 = biv_gen_pdf(bh.block(0,0,2,1),Sigmah.block(0,0,2,2));
  
  Sh = Sigmah.block(0,0,2,2);
  R12 = Sigmah.block(2,0,1,2) * Sh.inverse(); 
  mu_b = bh(2) - (R12(0) * bh(0)) - (R12(1)* bh(1));
  Sigb = R12 * Sigmah.block(0,0,2,2) * R12.transpose();  
  f2 = std_normal_cdf(mu_b/sqrt(Sigmah(2,2) - Sigb(0,0)));
  
  Gr(1,2) = f1* f2;
  
  // now start Hessian calculations. 
  // f1 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  df1.setZero();
  
  detSigma = Sigmah(0,0)*Sigmah(1,1) - Sigmah(0,1)*Sigmah(1,0);
  bexp = bh(0)*bh(0)*Sigmah(1,1) + bh(1)*bh(1)*Sigmah(0,0) - 2*bh(0)*bh(1)*Sigmah(0,1);
  
  ddetSigma.setZero();
  
  ddetSigma(0,0)= Sigmah(1,1);
  ddetSigma(0,1)= -2*Sigmah(0,1);
  ddetSigma(0,3) = Sigmah(0,0);
  
  df1 = f1*(-1/(2*detSigma)+bexp/(2*detSigma*detSigma)) * ddetSigma;
  df1(0,0) += (-f1)*(bh(1)*bh(1)/(2*(detSigma)));
  df1(0,1) += (f1)*(bh(0)*bh(1)/((detSigma)));
  df1(0,3) += (-f1)*(bh(0)*bh(0)/(2*(detSigma)));
  
  // f2 derived with respect to (s_11,s_12,s_13,s_22,s_23,s_33)
  sig_bar = Sigmah(2,2) - Sigb(0,0); 
  df2.setZero();
  db_bar.setZero();
  dsig_bar.setZero();
  dR_12.setZero();
  
  dR_12(0,0) = -(Sigmah(0,2)*Sigmah(1,1)-Sigmah(1,2)*Sigmah(0,1))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(0,1) = -2*dR_12(0,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(0,3) = dR_12(0,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(0,1) += (-1/detSigma)*Sigmah(1,2);
  dR_12(0,2) += Sigmah(1,1)/detSigma;
  dR_12(0,3) += Sigmah(0,2)/detSigma;
  dR_12(0,4) += (-1/detSigma)*Sigmah(0,1);
  
  
  dR_12(1,0)= -(-Sigmah(0,2)*Sigmah(0,1)+Sigmah(1,2)*Sigmah(0,0))/(detSigma*detSigma)*Sigmah(1,1);
  dR_12(1,1) = -2*dR_12(1,0)/Sigmah(1,1)*Sigmah(0,1);
  dR_12(1,3) = dR_12(1,0)/Sigmah(1,1)*Sigmah(0,0);
  dR_12(1,0) += Sigmah(1,2)/detSigma;
  dR_12(1,1) += -Sigmah(0,2)/detSigma;
  dR_12(1,2) += -Sigmah(0,1)/detSigma;
  dR_12(1,4) += Sigmah(0,0)/detSigma;
  
  
  
  db_bar = -bh(0)*dR_12.block(0,0,1,6) - bh(1)* dR_12.block(1,0,1,6);
  dsig_bar.block(0,0,1,6) = - Sigmah(0,2)*dR_12.block(0,0,1,6) - Sigmah(1,2)*dR_12.block(1,0,1,6);
  dsig_bar(0,2) += -R12(0,0);
  dsig_bar(0,4) += -R12(1,0);
  
  dsig_bar(0,5) = 1; 
  
  
  df2 = normal_pdf(mu_b,sig_bar)*(db_bar - mu_b*dsig_bar/(2*sig_bar)); 
  
  
  
  // combine the pieces 
  // grad = f1*f2 / sqrt(sigma_11 *sigma_33)
  
  He(7,3)= df1(5)*f2 + f1*df2(5);
  He(7,4)= df1(2)*f2 + f1*df2(2);
  He(7,5)= df1(4)*f2 + f1*df2(4);
  He(7,6) = df1(0)*f2 + f1*df2(0);
  He(7,7)= df1(1)*f2 + f1*df2(1);
  He(7,8) = df1(3)*f2 + f1*df2(3);
  
  He.block(7,3,1,6)=He.block(7,3,1,6)/(D(1,1)*D(2,2));
  
  He(7,4) = He(7,4)/(D(0,0)*D(1,1));
  He(7,5) = He(7,5)/(D(0,0)*D(2,2));
  He(7,7) = He(7,7)/(D(2,2)*D(1,1));
  // derivative w.r.t. sigma_11 reconstructed from others. 
  He(7,6) = -(He(7,1)*b(1) + He(7,4)*Sigma(0,1)+ He(7,7)*Sigma(1,2)+ (Gr(1,2) / (D(1,1)* D(2,2))))/(2*Sigma(1,1));
  He(7,8) = -(He(7,2)*b(2) + He(7,5)*Sigma(0,2)+ He(7,7)*Sigma(1,2)+ (Gr(1,2) / (D(1,1)* D(2,2))))/(2*Sigma(2,2));
  
  He(7,3) = He(7,3)/Sigma(0,0);
  
  // use symmetry to obtain entries in rows 3,6 and 8 and columns 4,5 and 7. 
  He(3,4)=He(4,3);
  He(3,5)= He(5,3);
  He(3,7)= He(7,3);
  He(6,4)= He(4,6);
  He(6,5)= He(5,6);
  He(6,7)= He(7,6);
  He(8,4)= He(4,8);
  He(8,7)= He(7,8);
  He(8,5)= He(5,8);
  
  ////////////////////////////////////////////////
    ///////////////////////////////////////////////
    // fill in results into grad   
  Eigen::VectorXd grad(9);
  grad.setZero();
  grad = grad_gen_tvn_cdf(b, Sigma);
  
  // use the gradients to fill in the gaps for sigma_11, sigma_22  and sigma_33.
  
  He(3,3) = (-He(0,3)*b(0) - He(4,3)*Sigma(0,1) - He(5,3)*Sigma(0,2) - 2*grad(3))/(2*Sigma(0,0));
  He(3,6) = (-He(0,6)*b(0) - He(4,6)*Sigma(0,1) - He(5,6)*Sigma(0,2))/(2*Sigma(0,0));
  He(3,8) = (-He(0,8)*b(0) - He(4,8)*Sigma(0,1) - He(5,8)*Sigma(0,2))/(2*Sigma(0,0));
  
  He(6,3) = He(3,6);
  He(8,3) = He(3,8); 
  
  He(6,6) = (-He(1,6)*b(1) - He(4,6)*Sigma(0,1) - He(7,6)*Sigma(1,2) - 2*grad(6))/(2*Sigma(1,1));
  He(6,8) = (-He(1,8)*b(1) - He(4,8)*Sigma(0,1) - He(7,8)*Sigma(1,2))/(2*Sigma(1,1));
  He(8,6)= He(6,8); 
  
  He(8,8) = (-He(2,8)*b(2)-He(7,8)*Sigma(1,2) - He(5,8)*Sigma(0,2)- 2*grad(8))/(2*Sigma(2,2)); 
  
  
  // return statement
  return He;
  
}
