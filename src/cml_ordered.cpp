// // // // // // // // // // // // // // // // //
//                                              //
//    Oredered PROBIT CML WITH APPROXIMATED MVNCDF        //
//                                              //
// // // // // // // // // // // // // // // // //
// Does not include the prediction part and hence uses
// data that has been treated before to subtract
// the X-part for the actual choice.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <omp.h>

using namespace Rcpp;



// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

#include "SJ_vdb.h"
#include "TVBS.h"
#include "ME.h"
#include "macml_grad_Hessian.h"

#include "toms462.h"                                         // allows bivariate normal cdf calculations


//////////////////////////////////////////////////////////////////////
///// ordered probit calculations using CML criterion function  /////
////////////////////////////////////////////////////////////////////

//' MACML function for ordered probit models
//' @description
//' Computes the approximate composite likelihood for the ordered probit case.
//' @param theta
//' parameter vector
//' @param data_obj
//' \link{data_cl} organised as a data_cl object containing (a) list with elements contains the nxk matrix X and the nx1 vector y (b) N number of deciders (c) number of choice situations per decider
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the ordered probit estimation.
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//' @export
//'
// [[Rcpp::export]]
NumericVector ll_macml_o(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control){
  //............................................................//
  // UNPACK ESTIMATION CONTROLS
  //............................................................//

  std::string approx_method = "SJ";
  int cml_pair_type = 1;
  int hess = 0;
  int el = 0;

  if(control.containsElementNamed("approx_method")){
    CharacterVector approx_method_temp = control["approx_method"];
    approx_method = Rcpp::as<std::string>(approx_method_temp);
  }

  if(control.containsElementNamed("control_weights")){
    Rcpp::List control_weights = control["control_weights"];
    if(control_weights.containsElementNamed("cml_pair_type")){
      Rcpp::List cml_pair_type_list = control_weights["cml_pair_type"];
      if(cml_pair_type_list.containsElementNamed("pair_type")){
        cml_pair_type = cml_pair_type_list["pair_type"];
      }
    }
  }


  if(control.containsElementNamed("hess")){
    hess = control["hess"];
  }

  if(control.containsElementNamed("el")){
    el = control["el"];
  }


  //............................................................//
  // Safeguard for approx_method input
  //............................................................//

  if(approx_method == "tvbs" || approx_method == "Tvbs"){
    approx_method = "TVBS";
  }
  if(approx_method == "me" || approx_method == "ME"){
    approx_method = "ME";
  }
  if(approx_method == "sj" || approx_method == "SJ"){
    approx_method = "SJ";
  }

  if(approx_method != "SJ" && approx_method != "ME" && approx_method != "TVBS"){
    approx_method = "TVBSv2";
  }


  //...........................................................//
  // set up output numeric vector
  //...........................................................//

  NumericVector out;
  double ll=0;

  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//


  int n=0;
  int N= data_obj["N"];
  int lthb = mod["lthb"];
  int lthO = mod["lthO"];
  int lthL = mod["lthL"];
  int lRE = mod["lRE"];
  int alt = mod["alt"];
  Eigen::VectorXd Tp = data_obj["Tp"];
  n = Tp(0);
  Eigen::MatrixXd Hb = mod["Hb"];
  Eigen::VectorXd fb = mod["fb"];
  Eigen::VectorXd fL = mod["fL"];
  Eigen::MatrixXd HL = mod["HL"];
  Eigen::VectorXd fO = mod["fO"];
  Eigen::MatrixXd HO = mod["HO"];


  //............................................................//
  // RECOVER PARAMETERS
  //............................................................//

  // Recover the vector b
  Eigen::VectorXd b = ((Hb*theta.head(lthb)) + fb);
  int lthbb = b.rows();


  // Recover the Sigma matrix
  Eigen::MatrixXd ELt =  elimmat(n).transpose();
  Eigen::MatrixXd Sigma_chol_vech;

  // Without free parameters only fL is used
  if(lthL==0){
    Sigma_chol_vech = fL;
  }else{
    Sigma_chol_vech = ((HL*theta.segment((lthO+lthb), lthL)) + fL);
  }

  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();
  Eigen::MatrixXd Sigma_chol = ELt*Sigma_chol_vech;
  Sigma_chol.resize(n,n);
  Eigen::MatrixXd Sigma_chol_t = Sigma_chol.transpose();
  Eigen::MatrixXd Sigma = Sigma_chol*Sigma_chol_t;


  // Recover the matrix Omega
  Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
  Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t, Omega;
  int lth_Omega_chol_vech;

  // check if we have mixed parameters
  if(lRE!=0){
    // Without free parameters only fO is used
    if(lthO==0){
      Omega_chol_vech = fO;
    }else{
      //Omega_chol_vech = ((HO*o) + fO);
      Omega_chol_vech = ((HO*theta.segment((lthb), lthO)) + fO);
    }
    lth_Omega_chol_vech = Omega_chol_vech.rows();
    Omega_chol = EOt*Omega_chol_vech;
    Omega_chol.resize(lRE,lRE);
    Omega_chol_t = Omega_chol.transpose();
    Omega = Omega_chol*Omega_chol_t;
  }else{
    lth_Omega_chol_vech = 0;
    Omega_chol_vech = fO;
    Omega_chol = fO;
    Omega_chol_t = fO;
    Omega = fO;
  }

  // get the cutoff points

  Eigen::VectorXd dtauk = theta.tail(alt-1);

  //............................................................//
  // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT
  //............................................................//

  // deriv of xb,vech(Omega),vech(Sigma) with respect to theta.
  int m_l =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech); //compute the dimension of the matrix from the length of the vector

  Eigen::MatrixXd L_JcholL = elimmat(m_l);
  Eigen::MatrixXd D_JcholL  = duplmat(m_l);
  Eigen::MatrixXd K_JcholL  = commmat(m_l,m_l);

  Eigen::VectorXd HlL = Sigma_chol_vech.array();
  Eigen::VectorXd HlL_n = Sigma_chol_vech.array();

  int m_o =  -0.5 + sqrt(0.25 + 2*lth_Omega_chol_vech); //compute the dimension of the matrix from the length of the vector
  Eigen::MatrixXd L_JcholO = elimmat(m_o);
  Eigen::MatrixXd D_JcholO  = duplmat(m_o);
  Eigen::MatrixXd K_JcholO  = commmat(m_o,m_o);
  //This is used later in the computation to compute the Hessian of the Cholesky
  Eigen::VectorXd HoO = Omega_chol_vech.array();
  Eigen::VectorXd HoO_n = Omega_chol_vech.array();

  Eigen::MatrixXd JJ1(lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech,lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech);
  JJ1.setZero();

  Eigen::MatrixXd JJ1_n(lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech,lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech);
  JJ1_n.setZero();

  if(lRE!=0){
    JJ1 = J_bvOvSig_bLOLSig(HoO, HlL, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);
  }

  if(lRE == 0){
    JJ1 = J_bvSig_bLSig(HlL, lthbb, L_JcholL, K_JcholL, D_JcholL);
  }


  Eigen::MatrixXd JJ2(lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech,lthb+lthO+lthL);
  JJ2.setZero();
  JJ2.block(0,0,lthbb,lthb) = Hb;
  if(lRE!=0){
    JJ2.block(lthbb,lthb,lth_Omega_chol_vech,lthO) = HO;
  }
  JJ2.block(lthbb+lth_Omega_chol_vech,lthb+lthO,lth_Sigma_chol_vech,lthL) = HL;

  Eigen::MatrixXd IIl(lth_Sigma_chol_vech,lth_Sigma_chol_vech);
  IIl.setIdentity();

  //............................................................//
  // INITIALIZE GRADIENT OUTPUT
  //............................................................//

  // freeth denotes the number of free parameters for b, LO and LSig.
  int freeth = Hb.cols() +HL.cols() + HO.cols() + alt-1; // bsum, Osum, Lsum;


  Eigen::VectorXd grad(freeth);
  grad.setZero();
  Eigen::VectorXd grad_n(freeth);
  grad_n.setZero();
  Eigen::MatrixXd Hess(freeth,freeth), Hess_n(freeth,freeth);
  Hess.setZero();

  Eigen::MatrixXd Hess_approx(freeth,freeth);
  Hess_approx.setZero();

  Eigen::MatrixXd Hess_approx_n(freeth,freeth);
  Hess_approx_n.setZero();

  Eigen::MatrixXd Hess_approx2(freeth,freeth);
  Hess_approx2.setZero();

  Eigen::MatrixXd Hess_approx_n2(freeth,freeth);
  Hess_approx_n2.setZero();

  Eigen::MatrixXd I_alt(alt-1,alt-1);
  I_alt.setIdentity();

  double ll_n = 0;

  for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
    // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker

    ll_n= 0;
    grad_n.setZero();
    Hess_n.setZero();
    Hess_approx_n.setZero();
    Hess_approx_n2.setZero();

    Rcpp::List data = data_obj["data"];
    Rcpp::List data_n = data[i_n]; //unpack the data
    Eigen::MatrixXd X_n = data_n["X"];
    n = X_n.rows();
    Eigen::VectorXd y_n = data_n["y"];

    // compute X*b.
    Eigen::VectorXd xb(n);
    xb.setZero();

    xb = X_n * b;

    // compute Lambda.
    Eigen::MatrixXd Lambda(n,n);
    Lambda.setZero();
    Lambda = Sigma;

    Eigen::MatrixXd xRE(n,lRE);

    if(lRE!=0){
      xRE = X_n.block(0,0,n,lRE);
      Lambda+= xRE * Omega * xRE.transpose();
    }


    ll_n = prob_ordered_CML(xb, y_n, Lambda, alt, dtauk, cml_pair_type);


    Eigen::MatrixXd JJ0(n+lth_Sigma_chol_vech,lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech);
    JJ0.setZero();

    JJ0.block(0,0,n,lthbb)= X_n;

    if(lRE!=0){
      JJ0.block(n,lthbb,lth_Sigma_chol_vech,lth_Omega_chol_vech) = L_JcholL * kroneckerProduct(xRE,xRE) *  D_JcholO;
    }

    JJ0.block(n,lthbb+lth_Omega_chol_vech,lth_Sigma_chol_vech,lth_Sigma_chol_vech) = IIl;


    Eigen::MatrixXd J = JJ0*JJ1*JJ2; // d(b,vech(Omega),vech(Sigma)/dtheta = d(b,vech(Omega),vech(Sigma))/d(thb,vech(LO),vech(LS))* d(thb...)/dtheta.

    Eigen::VectorXd gradCML = grad_ordered_CML(xb, y_n, Lambda, alt, dtauk, cml_pair_type);


    Eigen::MatrixXd J_full(J.cols()+alt-1,J.rows()+alt-1);
    J_full.setZero();

    J_full.block(0,0,J.cols(),J.rows()) = J.transpose();
    J_full.block(J.cols(),J.rows(),alt-1,alt-1) = I_alt;

    grad_n = J_full * gradCML;


    /////////////////////////////////////////////////////////
    //// Hessian     ///
    //////////////////////////////////////////////////////////

    if (hess == 1){
      Eigen::MatrixXd hessCML = hess_ordered_CML(xb, y_n, Lambda, alt, dtauk, cml_pair_type);

      Hess_n = J_full * hessCML* J_full.transpose();

      HlL_n.setZero();
      Eigen::MatrixXd Zbb(lthbb,lthbb);
      Zbb.setZero();

      J_full.setZero();

      // parameters corresponding to Omega.
      for (int j=0;j<lth_Omega_chol_vech;j++){

        JJ1_n.setZero();
        HoO_n.setZero();
        HoO_n(j)=1;

        if(lRE!=0){
          JJ1_n = J_bvOvSig_bLOLSig(HoO_n, HlL_n, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);
        }
        JJ1_n.block(0,0,lthbb,lthbb)=Zbb;

        J = JJ0*JJ1_n*JJ2;
        J_full.block(0,0,J.cols(),J.rows()) = J.transpose();
        Hess_n.block(0,lthb,Hess.rows(),lthO) += J_full * gradCML * HO.block(j,0,1,lthO);
      }

      HoO_n.setZero();
      HlL_n.setZero();
      J_full.setZero();

      // parameters corresponding to Sigma.
      for (int j=0;j<lth_Sigma_chol_vech;j++){

        JJ1_n.setZero();
        HlL_n.setZero();
        HlL_n(j)=1;

        if(lRE!=0){
          JJ1_n = J_bvOvSig_bLOLSig(HoO_n, HlL_n, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);
        }
        if(lRE == 0){
          JJ1_n = J_bvSig_bLSig(HlL_n, lthbb, L_JcholL, K_JcholL, D_JcholL);
        }

        JJ1_n.block(0,0,lthbb,lthbb)=Zbb;

        J = JJ0*JJ1_n*JJ2;
        J_full.block(0,0,J.cols(),J.rows()) = J.transpose();
        Hess_n.block(0,lthb+lthO,Hess_n.rows(),lthL) += J_full * gradCML * HL.block(j,0,1,lthL);
      }
      Hess_n = 0.5*(Hess_n+Hess_n.transpose());
    }

    // return value.
    ll += ll_n; //it is already negative in the pairs loop

    grad += grad_n;
    Hess_approx2 += Hess_approx_n2;

    if (hess == 1){
      Hess += Hess_n;//it is already negative in the pairs loop
    }
    Hess_approx += grad_n * grad_n.transpose();//is positive throughout because it is only used this way

  } //end N-loop


  //............................................................//
  // Format the output
  //............................................................//
  out = -ll;
  out.attr("gradient") = -(grad.array());

  if(hess == 1){
    out.attr("hessian") =-(Hess.array());
  }

  out.attr("hessian1") =(Hess_approx.array());
  out.attr("hessian2") =(Hess_approx.array());
  return out;

}

Eigen::VectorXd int2bin(int num,int noBits){

  // set up return list
  Eigen::VectorXd out(noBits);
  for (int j=0;j<noBits;j++){
    out[j]=0;
    if (num % 2 == 1){
      out[j]=1;
    }
    num = (num - (num % 2))/2;
  }

  // return value

  return out;
}



//////////////////////////////////////////////////////////////////////
///// ordered probit calculations using likelihood  /////
////////////////////////////////////////////////////////////////////
//' Probability to multivariate ordered probit model.
//' @description
//' The function computes the probability of a number of choices in repeated ordered probit choices.
//' @param xb
//' nx1 vector of systematic part
//' @param y_n
//' nx1 integer vector of choices.
//' @param Lambda
//' Covariance matrix for random utilities including the random coefficients.
//' @param alt
//' integer; number of choice alternatives.
//' @param dtauk
//' (alt-1) real vector of increments of category limits obtained as the cumulative sum of the vector.
//' @param approx_method
//' string: approximation method, either "SJ" or "ME"
//' @return
//' double; log of probability.
//' @keywords internal
//' 
 // [[Rcpp::export]]
double prob_ordered(Eigen::VectorXd xb, Eigen::VectorXd y_n, Eigen::MatrixXd Lambda, int alt, Eigen::VectorXd dtauk, std::string approx_method){
  double ll =0;
  double pr = 0;


  // tauk's.
  Eigen::VectorXd tauk(alt-1);
  tauk.setZero();

  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    tauk(j)=tauk(j-1)+std::exp(dtauk(j));
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


  // find upper an lower bounds
  int n=xb.size();
  Eigen::MatrixXd upper_lower(n,2);
  upper_lower.setZero();

  for (int i_n=0;i_n<n;i_n++){

    upper_lower(i_n,0) = tauk(0)-10.0;
    upper_lower(i_n,1)=  tauk(alt-2)+10.0;

    if (y_n(i_n)<1){ y_n(i_n)=1;}
    if (y_n(i_n)>alt){ y_n(i_n)=alt;}

    if (y_n(i_n)==1){
      upper_lower(i_n,1) =tauk(0);
    }
    if (y_n(i_n)==alt){
      upper_lower(i_n,0) = tauk(alt-2);
    }

    if (y_n(i_n)>1){
      if(y_n(i_n)<alt){
        upper_lower(i_n,0) = tauk(y_n(i_n)-2);
        upper_lower(i_n,1) = tauk(y_n(i_n)-1);
      }
    }

  }


  // choice probability is the sum over all upper,lower bounds
  // cycle over all
  Eigen::VectorXd bound(n);

  for (int i_j=0;i_j<pow(2,n);i_j++){
    Eigen::VectorXd ind_upper_lower = int2bin(i_j,n);
    // bound sets the
    bound.setZero();

    double sign_ind=1;

    for (int i_n=0;i_n<n;i_n++){
      if (ind_upper_lower(i_n) == 1){
        bound(i_n) = upper_lower(i_n,1)/lambda(i_n);
      } else {
        sign_ind = sign_ind * (-1.0);
        bound(i_n) = upper_lower(i_n,0)/lambda(i_n);
      }
    }

    double pr_n=1;

    if (n>1){ // SJ and ME only used for n>1!!
      if (approx_method == "SJ"){
        pr_n = exp(SJ(bound-x_norm,Lambda_cor));
      } else {
        pr_n = exp(ME(bound-x_norm,Lambda_cor));
      }

    } else {
      pr_n = std_normal_cdf(bound(0)-x_norm(0));
    }

    if (sign_ind<0 ){
      pr -=pr_n;
    } else {
      pr += pr_n;
    }

  }
  ll = log(pr);
  return ll;
}



// // // // // // // // // // // // // // // // //
//                                              //
//    PREDICTION CALCULATION                    //
//                                              //
// // // // // // // // // // // // // // // // //

// input:
//  - theta           VectorXd    parameter vector (to optimize)
//  - data            List        observed data for one decider.
//  - approx_method   string      approximation method for the MVNCDF ("SJ", "ME" or "TVBS")

//' Choice probabilities for ordered probit models
//' @description
//' Computes the approximate choice probabilities for the ordered probit case.
//' @param theta
//' parameter vector
//' @param Xn
//' matrix of regressors
//' @param yn
//' vector of responses
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @return
//' A matrix, containing the predicted choice probabilities for each choice marginally (ignoring correlations)
//' @export
//'
// [[Rcpp::export]]
Eigen::MatrixXd pred_probit_ordered_approx(Eigen::VectorXd theta, Eigen::MatrixXd Xn, Eigen::VectorXd yn, Rcpp::List mod)
{


  std::string approx_method = "SJ";


  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//


  int lthb = mod["lthb"];
  int lthO = mod["lthO"];
  int lthL = mod["lthL"];
  int lRE = mod["lRE"];
  int alt = mod["alt"];
  Eigen::MatrixXd Hb = mod["Hb"];
  Eigen::VectorXd fb = mod["fb"];
  Eigen::VectorXd fL = mod["fL"];
  Eigen::MatrixXd HL = mod["HL"];
  Eigen::VectorXd fO = mod["fO"];
  Eigen::MatrixXd HO = mod["HO"];



  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//
  int n = yn.rows();

  Eigen::MatrixXd prob_predict(n,alt);
  prob_predict.setZero();


  //............................................................//
  // RECOVER PARAMETERS
  //............................................................//

  // Recover the vector b
  Eigen::VectorXd b = ((Hb*theta.head(lthb)) + fb);
  int lthbb = b.rows();

  // Recover the Sigma matrix
  Eigen::MatrixXd ELt =  elimmat(n).transpose();
  Eigen::MatrixXd Sigma_chol_vech;

  // Without free parameters only fL is used
  if(lthL==0){
    Sigma_chol_vech = fL;
  }else{
    Sigma_chol_vech = ((HL*theta.segment((lthO+lthb), lthL)) + fL);
  }

  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();
  Eigen::MatrixXd Sigma_chol = ELt*Sigma_chol_vech;
  Sigma_chol.resize(n,n);
  Eigen::MatrixXd Sigma_chol_t = Sigma_chol.transpose();
  Eigen::MatrixXd Sigma = Sigma_chol*Sigma_chol_t;


  // Recover the matrix Omega
  Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
  Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t, Omega;
  int lth_Omega_chol_vech;

  // check if we have mixed parameters
  if(lRE!=0){
    // Without free parameters only fO is used
    if(lthO==0){
      Omega_chol_vech = fO;
    }else{
      //Omega_chol_vech = ((HO*o) + fO);
      Omega_chol_vech = ((HO*theta.segment((lthb), lthO)) + fO);
    }
    lth_Omega_chol_vech = Omega_chol_vech.rows();
    Omega_chol = EOt*Omega_chol_vech;
    Omega_chol.resize(lRE,lRE);
    Omega_chol_t = Omega_chol.transpose();
    Omega = Omega_chol*Omega_chol_t;
  }else{
    lth_Omega_chol_vech = 0;
    Omega_chol_vech = fO;
    Omega_chol = fO;
    Omega_chol_t = fO;
    Omega = fO;
  }

  // get the cutoff points

  Eigen::VectorXd dtauk = theta.tail(alt-1);
  Eigen::VectorXd tauk(alt-1);
  tauk.setZero();

  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    tauk(j)=tauk(j-1)+std::exp(dtauk(j));
  }

  //............................................................//
  // END of setup
  //............................................................//

  //............................................................//
  // CYCLE over choices to calculate the  probs.
  //............................................................//

  Eigen::VectorXd xb(n);
  xb.setZero();

  xb = Xn * b;

  // compute Lambda.
  Eigen::MatrixXd Lambda(n,n);
  Lambda.setZero();
  Lambda = Sigma;

  Eigen::MatrixXd xRE(n,lRE);

  if(lRE!=0){
    xRE = Xn.block(0,0,n,lRE);
    Lambda+= xRE * Omega * xRE.transpose();
  }

  double xbn = 1.0;

  for(int i_choice_1 = 0; i_choice_1 < n; i_choice_1++){
    xbn = xb(i_choice_1) / Lambda(i_choice_1,i_choice_1);
    prob_predict(i_choice_1,0) = std_normal_cdf(tauk(0)-xbn);
    for (int y=1; y < alt-1; y++){
      prob_predict(i_choice_1,y) = std_normal_cdf(tauk(y)-xbn)-std_normal_cdf(tauk(y-1)-xbn);
    }

    // last probability
    prob_predict(i_choice_1,alt-1) = 1.0- std_normal_cdf(tauk(alt-2)-xbn);

  } //end i_choice_1 loop over choice decisions.


  //............................................................//
  // Format the output
  //............................................................//

  return  prob_predict;

}


