// // // // // // // // // // // // // // // // //
//                                              //
//    PROBIT CML WITH APPROXIMATED MVNCDF        //
//    functions for comparing the asymptotic variance //
//                                              //
// // // // // // // // // // // // // // // // //

#include <RcppEigen.h>
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



//////////////////////
// define structure //
//////////////////////

struct Contrib
{
  double ll;
  Eigen::VectorXd grad;
  Eigen::MatrixXd Hess;
};

struct Contrib_Var
{
  double ll;
  Eigen::MatrixXd Hess2;
  Eigen::MatrixXd Hess1;
  Eigen::MatrixXd Hess;
};



// // // // // // // // // // // // // // // // //
//                                              //
//    PROBIT LIKELIHOOD: contribution of one person  //
//                                              //
// // // // // // // // // // // // // // // // //

Contrib ll_probit_contrib(Rcpp::List X_n, Eigen::VectorXd y_n, Eigen::VectorXd theta, Rcpp::List mod, std::string approx_method)
{

  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//

  //int N = mod["N"];
  int lthb = mod["lthb"];
  int lthO = mod["lthO"];
  int lthL = mod["lthL"];
  int lRE = mod["lRE"];
  int alt = mod["alt"];
  //Eigen::VectorXd Tp = mod["Tp"];
  Eigen::MatrixXd Hb = mod["Hb"];
  Eigen::VectorXd fb = mod["fb"];
  Eigen::VectorXd fL = mod["fL"];
  Eigen::MatrixXd HL = mod["HL"];
  Eigen::VectorXd fO = mod["fO"];
  Eigen::MatrixXd HO = mod["HO"];


  //............................................................//
  // Safeguard for approx_method input
  //............................................................//

  if(approx_method == "me" || approx_method == "ME"){
    approx_method = "ME";
  }
  if(approx_method == "sj" || approx_method == "SJ"){
    approx_method = "SJ";
  }

  if(approx_method != "SJ" && approx_method != "ME"){
    approx_method = "SJ";
  }


  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//

  Contrib out;

  double ll = 0;

  //............................................................//
  // RECOVER PARAMETERS
  //............................................................//

  // Recover the vector b
  Eigen::VectorXd b = ((Hb*theta.head(lthb)) + fb);
  int lthbb = b.rows();

  // Recover the Sigma matrix
  Eigen::MatrixXd ELt =  elimmat(alt).transpose();
  Eigen::MatrixXd Sigma_chol_vech;

  // Without free parameters only fL is used
  if(lthL==0){
    Sigma_chol_vech = fL;
  }else{
    Sigma_chol_vech = ((HL*theta.segment((lthO+lthb), lthL)) + fL);
  }

  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();
  Eigen::MatrixXd Sigma_chol = ELt*Sigma_chol_vech;
  Sigma_chol.resize(alt,alt);
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

  //............................................................//
  // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT AND HESSIAN
  //............................................................//

  Eigen::MatrixXd L_JOm = elimmat(lRE); //only for pairwise with equal number of alts for all Decision Makers
  Eigen::MatrixXd D_JOm  = duplmat(lRE);

  int m_l =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech); //compute the dimension of the matrix from the length of the vector
  int m_o =  -0.5 + sqrt(0.25 + 2*lth_Omega_chol_vech); //compute the dimension of the matrix from the length of the vector

  Eigen::MatrixXd L_JcholL = elimmat(m_l);
  Eigen::MatrixXd D_JcholL  = duplmat(m_l);
  Eigen::MatrixXd K_JcholL  = commmat(m_l,m_l);

  Eigen::MatrixXd L_JcholO = elimmat(m_o);
  Eigen::MatrixXd D_JcholO  = duplmat(m_o);
  Eigen::MatrixXd K_JcholO  = commmat(m_o,m_o);

  //This is used later in the computation to compute the Hessian of the Cholesky which does not depend on the specific value of the parameters
  Eigen::VectorXd HoO = Omega_chol_vech.array() + 1;
  Eigen::VectorXd HlL = Sigma_chol_vech.array() + 1;

  Eigen::MatrixXd HJ2 = J_bvOvSig_bLOLSig(HoO, HlL, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);

  //............................................................//
  // END of setup
  //............................................................//

  // freeth counts the free parameters for b, LOm and LSig.
  int freeth = Hb.rows() +HL.rows() + HO.rows(); // bsum, Osum, Lsum;

  //............................................................//
  // INITIALIZE GRADIENT OUTPUT
  //............................................................//

  Eigen::MatrixXd Hess(lthb+lthO+lthL,lthb+lthO+lthL);
  Hess.setZero();
  Eigen::MatrixXd Hess_approx(lthb+lthO+lthL,lthb+lthO+lthL);
  Hess_approx.setZero();
  Eigen::MatrixXd Hess_approx2(lthb+lthO+lthL,lthb+lthO+lthL);
  Hess_approx2.setZero();
  Eigen::VectorXd grad(lthb+lthO+lthL);
  grad.setZero();

  int Tpn = y_n.rows();        //    .length();

  //save the data of the current decision maker
  int alt2 = Tpn*(alt-1);

  // setup X, M, M_Sig.
  Eigen::MatrixXd X(alt2,lthbb);
  X.setZero();

  Eigen::MatrixXd M_Sig(alt2,alt2);
  M_Sig.setZero();


  // NOW LOOP OVER THE observations MATRIX
  for(int jj = 0; jj < Tpn; jj++){ // loop jj over observations.

    //............................................................//
    // Data Prep and Lin Utilities
    //............................................................//
    //
    //  extract the chosen alternatives for the current pair of choice occasions
    int y1;
    y1 = y_n[jj];

    // extract the covariance matrix for the current pair of choice occasions
    Eigen::MatrixXd X1(alt,lthbb);
    X1.setZero();
    X1 = X_n[jj];

    Eigen::MatrixXd M1 = makeM((alt),(y1-1));
    X.block((jj)*(alt-1),0,alt-1,lthbb) = M1*X1;

    Eigen::MatrixXd M1t = M1.transpose();
    Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;

    //............................................................//
    // Covariance Matrix
    //............................................................//

    M_Sig.block((jj)*(alt-1),(jj)*(alt-1),alt-1,alt-1) = Sigma_norm_1;
  } // end jj loop over observations.

  // find the NA elements

  //NumericVector x1c = X.block(0,0,X.rows(),1);
  Eigen::MatrixXd X_red = X;
  Eigen::MatrixXd M_Sig_red = M_Sig;

  IntegerVector ind = NA_ind(X.block(0,0,X.rows(),1));

  int cur = ind.size();

  for (int ij=0;ij<cur;ij++){
    X_red.row(ij) = X.row(ind[ij]);
    M_Sig_red.row(ij) = M_Sig.row(ind[ij]);
  }

  for (int ij=0;ij<cur;ij++){
    M_Sig_red.col(ij) = M_Sig_red.col(ind[ij]);
  }

  // // unstandardized lin-utilties (upper limits of the MVNCDF)

  Eigen::VectorXd x_upper(cur);
  x_upper = (-1)*(X_red.block(0,0,cur,X.cols())*b);

  // // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
  Eigen::MatrixXd X_RE(cur,lRE);
  X_RE = X_red.block(0,0,cur, lRE);
  Eigen::MatrixXd X_RE_t = X_RE.transpose();

  // Compute combined covariance matrix
  // Changed Sig to Lambda for combined Covariance
  Eigen::MatrixXd Lambda(cur,cur);
  Lambda = M_Sig_red.block(0,0,cur,cur);
  if(lRE!=0){
    Lambda = Lambda + X_RE*Omega*X_RE_t;
  }

  Eigen::MatrixXd lambda(cur,1);
  lambda = Lambda.diagonal().array().sqrt();
  Eigen::MatrixXd lambda_t = lambda.transpose();

  // ............................................................//
  //  Normalize upper limits and covariance matrix (to correlation matrix)
  // ............................................................//
  //
  Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
  Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
  Eigen::MatrixXd Lambda_cor(cur,cur);
  Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();

  //............................................................//
  // Likelihood and Gradient
  //............................................................//
  //

  if (approx_method == "SJ"){
    ll= -SJ(x_norm,Lambda_cor); // added the weight
  } else {
    ll= -ME(x_norm,Lambda_cor); // added the weight
  }

  Eigen::VectorXd grad_output;
  grad_output.setZero();
  int lth_output;

  if(approx_method == "SJ"){
    grad_output = dlcond(x_norm,Lambda_cor);
  }

  if(approx_method == "ME"){
    grad_output = dlcond_ME(x_norm,Lambda_cor);
  }

  if(approx_method == "TVBSv2"){
    grad_output = TVBS_grad(x_norm,Lambda_cor,1);
  }

  if(approx_method == "TVBS"){
    grad_output = TVBS(x_norm,Lambda_cor,1);
  }

  lth_output = grad_output.size();
  ll = -grad_output[lth_output-1];
  //............................................................//
  // Gradient
  //............................................................//
  //

  Eigen::VectorXd grad_MVNCDF_approx;
  grad_MVNCDF_approx = grad_output.head((lth_output-1));

  //............................................................//
  // START: Manipulations to go from MVNCDF approximation gradient
  //        to gradient of the original problem
  //............................................................//
  //
  // calculate the Jacobian matrix of theta with respect to b and vech(Lambda)
  // // calculate the Jacobian matrix of w,vech(Lambda_cor) with respect to b and vech(Lambda)
  Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta(X_red.block(0,0,cur,X.cols()), b, Lambda);

  // // partial log prob / partial (b, vech(Lambda))
  Eigen::VectorXd g_bvechLambda = grad_MVNCDF_approx.transpose()*J_wLambda_cor_bvechLambda; //column to row

  // partial (b, vech(Lambda), vech(X Omega X')) / partial (b, vech(Omega), vech(Sigma))
  Eigen::MatrixXd L = elimmat(cur);

  Eigen::MatrixXd J_dM_dbOmSig =  J_M_bOmSig_red(alt,ind, y_n, lthbb, X_RE.block(0,0,cur,X_RE.cols()),L ,D_JOm);

  // // partial (b, vech(Omega), vech(Sigma)) / partial (b, vech(L_Omega),vech(L_Sigma))
  Eigen::MatrixXd J_chol = J_bvOvSig_bLOLSig(Omega_chol_vech, Sigma_chol_vech, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);

  Eigen::VectorXd g_bOmSig  = g_bvechLambda.transpose()*J_dM_dbOmSig;

  Eigen::VectorXd g_bLOLSig = g_bOmSig.transpose()*J_chol;


  Eigen::MatrixXd gradbb = g_bLOLSig.head(lthbb);
  // Hb ... partial b / partial th_b.
  Eigen::MatrixXd gradb =  Hb.transpose()*gradbb;

  // HL ... partial vech(L) / partial th_Sig.
  Eigen::VectorXd gradL1 = g_bLOLSig.segment((lth_Omega_chol_vech+lthbb), lth_Sigma_chol_vech);
  Eigen::MatrixXd gradL = HL.transpose()*gradL1; //pick only the relevant parameters

  // HO ... partial vech(L_Omega) / partial th_Om.
  Eigen::VectorXd gradO1 = g_bLOLSig.segment((lthbb), lth_Omega_chol_vech);
  Eigen::MatrixXd gradO = HO.transpose()*gradO1; //pick only the relevant parameters

  // fill in the gradients, depending on which is present.
  if(lthL == 0){
    if(lthO == 0){
      grad << gradb;
    } else {
      grad << gradb,gradO;
    }
  } else {
    if(lthO == 0){
      grad << gradb, gradL;
    } else {
      grad << gradb,gradO, gradL;
    }
  }

  //............................................................//
  // END of gradient manipulations
  //............................................................//

  // add gradient of the current pair to the gradient of the decicion maker
  grad = (-1)*grad;

  //............................................................//
  //Hessian
  //............................................................//


  //First step
  Eigen::MatrixXd hes;
  if (approx_method == "SJ"){
    hes = SJ_hess_new(x_norm, Lambda_cor);
  } else {
    hes = ME_hess_new(x_norm, Lambda_cor);
  }
  Eigen::MatrixXd H_th = H_theta(X_red.block(0,0,cur,X.cols()), b, Lambda);

  Eigen::MatrixXd Ide(H_th.cols(), H_th.cols());
  Ide.setIdentity();
  Eigen::MatrixXd H_bvechLambda = J_wLambda_cor_bvechLambda.transpose()*hes*J_wLambda_cor_bvechLambda + Eigen::kroneckerProduct(grad_MVNCDF_approx.transpose(), Ide)*H_th;
  //Second step
  Eigen::MatrixXd H_bOmSig = J_dM_dbOmSig.transpose()*H_bvechLambda*J_dM_dbOmSig;

  //Third Step
  //b identity
  //Sigma_chol_vech & oL Cholesky
  Eigen::MatrixXd H_bLOLSig = H_bvOvSig_bLOLSig(HJ2, lthbb, lth_Omega_chol_vech, lth_Sigma_chol_vech);

  Eigen::SparseMatrix<double> sH_bLOLSig = H_bLOLSig.sparseView();
  Eigen::SparseMatrix<double> Ide1(H_bLOLSig.cols(), H_bLOLSig.cols());
  Ide1.setIdentity();

  Eigen::SparseMatrix<double> temp = Eigen::kroneckerProduct(g_bOmSig.transpose(),Ide1)*sH_bLOLSig;

  Eigen::MatrixXd H_t = Eigen::MatrixXd(temp);

  Eigen::MatrixXd Hess_n = J_chol.transpose()*H_bOmSig*J_chol+ H_t;


  // add multiplication with HO, HL und Hb.

  Eigen::MatrixXd J_H = Eigen::MatrixXd(lthb+lthO+lthL,freeth);
  J_H.setZero();

  J_H.block(0,0,lthb,lthbb) = Hb.transpose();
  J_H.block(lthb,lthbb,lthO,lth_Omega_chol_vech) = HO.transpose();
  J_H.block(lthb+lthO,lthbb+lth_Omega_chol_vech,lthL,lth_Sigma_chol_vech) = HL.transpose();

  Hess = (-1) * J_H * Hess_n * J_H.transpose();//it is already negative in the pairs loop

  //............................................................//
  // Format the output
  //............................................................//
  out.ll = ll;
  out.grad = grad;
  out.Hess = Hess;

  return out;

}


/////////////////////////////////////////
// Function to convert a number into a //
// vector of choice values             //
/////////////////////////////////////////

//' providing categories for integers for cycling through upper and lower bounds
//' @description
//' converts integers into vector of category indices
//' @param num
//' integer index
//' @param cat
//' number of alternatives
//' @param noBits
//' integer number of bits
//' @return
//' vector of categories
//' @keywords internal
//'
// [[Rcpp::export]]
Eigen::VectorXd int2cats(int num, int cat,int noBits){

  // set up return list
  Eigen::VectorXd out(noBits);
  for (int j=0;j<noBits;j++){
    out[j]=num % cat + 1;
    num = int(num /cat);
  }
  // return value

  return out;
}



// evaluating the variance components.

Contrib ll_probit_exp_contrib(Rcpp::List X_n, Eigen::VectorXd theta, Rcpp::List mod, std::string approx_method){

  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//

  Contrib out;

  double ll = 0;

  int lth = theta.size();
  Eigen::VectorXd grad(lth);
  grad.setZero();

  Eigen::MatrixXd Hess(lth,lth);
  Hess.setZero();

  Contrib lln;

  //.......................................//
  //  loop over y-values
  //.......................................//

  Eigen::VectorXd yn(X_n.size());
  yn.setZero();
  int alt = mod["alt"];

  int Tp = pow(alt,X_n.size());

  for (int jj=0; jj<Tp;jj++){
    yn = int2cats(jj, alt,yn.size());

    Contrib lln = ll_probit_contrib(X_n,yn,theta,mod,approx_method);
    // update
    double prob = exp((-1)*lln.ll);

    ll += prob * lln.ll;
    grad += prob * lln.grad;
    Hess += prob * lln.Hess;
  }
  //////////////////////
  // return results
  //////////////////////
  out.ll = ll;
  out.grad = grad;
  out.Hess = Hess;

  return out;
}

//' Calculation of contribution of one person to the approximated probit log likelihood calculation.
//' @description
//' Computes the contribution of one person to the log likelihood.
//' @param X_n
//' X list of regressor matrices
//' @param theta
//' parameter vector
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param pairwise
//' integer specifying the pairwise option
//' @param approx_method
//' string specifying the approx_method
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector ll_probit_person(Rcpp::List X_n, Eigen::VectorXd theta, Rcpp::List mod, int pairwise, std::string approx_method){

  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//

  NumericVector out;

  double ll;

  int lth = theta.size();
  Eigen::VectorXd grad(lth), gradf(lth);
  gradf.setZero();
  grad.setZero();
  Eigen::MatrixXd Hess(lth,lth);
  Hess.setZero();
  Eigen::MatrixXd Hess1(lth,lth);
  Hess1.setZero();
  Eigen::MatrixXd Hess2(lth,lth);
  Hess2.setZero();
  Contrib lln,llm;

  //.......................................//
  //  loop over y-values
  //.......................................//

  Eigen::VectorXd yn(X_n.size());
  yn.setZero();
  int alt = mod["alt"];

  int Tp = pow(alt,X_n.size());

  Eigen::VectorXd pr(Tp);
  pr.setZero();

    // calculate probability of happening
  double prob, prob_full;

  prob_full = 0;
  for (int jj=0; jj<Tp;jj++){
    yn = int2cats(jj, alt,yn.size());

    lln = ll_probit_contrib(X_n, yn, theta, mod, approx_method);
    prob = exp((-1)*lln.ll);
    pr[jj]= prob;
    prob_full += prob;


    if (pairwise == 0){
      ///////////////////////
      // probit case
      ///////////////////////
      grad.setZero();
      ll += lln.ll;
      grad = lln.grad;
      gradf += (prob * grad);
      Hess += (prob* lln.Hess);
      Hess1 += (prob* grad * grad.transpose());
    }

    if (pairwise == 1){
      ////////////////////////
      // full pairwise
      ////////////////////////
      int lX = X_n.size();


      Rcpp::List Xn(2); // list of pairs of decisions.
      Eigen::VectorXd yy(2);
      yy.setZero();
      grad.setZero();

      for (int jm=0;jm<lX-1;jm++){
        for (int im=(jm+1);im<lX;im++){
          Xn[0] = X_n[jm];
          Xn[1] = X_n[im];
          yy[0] = yn[jm];
          yy[1] = yn[im];

          llm = ll_probit_contrib(Xn, yy,theta, mod, approx_method);
          ll += llm.ll;
          gradf += (prob* llm.grad);
          grad += llm.grad;
          Hess += (prob* llm.Hess);
          Hess2 += (prob* llm.grad * llm.grad.transpose());
        }
      }
      Hess1 += (prob* grad * grad.transpose());
    }

    if (pairwise == 2){
      //////////////////////////////
      // adjacent pairwise
      //////////////////////////////

      Rcpp::List Xn(2); // list of adjacent pairs of decisions.
      Eigen::VectorXd yy(2);

      grad.setZero();
      int lX = X_n.size();
      for (int jm=0;jm<lX-1;jm++){
        Xn[0] = X_n[jm];
        Xn[1] = X_n[jm+1];
        yy(0)= yn(jm);
        yy(1) = yn(jm+1);

        lln = ll_probit_contrib(Xn,yy, theta, mod, approx_method);
        ll += lln.ll;
        grad += lln.grad;
        Hess += (prob* lln.Hess);
        Hess2 += (prob* lln.grad * lln.grad.transpose());
      }
      gradf += (prob * grad);
      Hess1 += (prob* grad * grad.transpose());
    }

    if (pairwise == 3){
      //////////////////////////////
      //  marginal + full pairwise
      //////////////////////////////

      Rcpp::List Xn(2); // list of pairs of decisions.
      Eigen::VectorXd yy(1);

      grad.setZero();

      int lX = X_n.size();
      for (int jm=0;jm<lX;jm++){


        Xn[0] = X_n[jm];
        Xn[1] = X_n[jm];
        yy[0]= yn(jm);

        llm = ll_probit_contrib(Xn, yy,theta, mod, approx_method);
        ll += llm.ll;
        grad += llm.grad;
        Hess2 += prob*llm.grad * llm.grad.transpose();
        Hess += (prob* llm.Hess);
      }


      Rcpp::List Xnf(2); // list of pairs of decisions.
      Eigen::VectorXd yyf(2);

      for (int jm=0;jm<lX-1;jm++){
        Xnf[0] = X_n[jm];
        Xnf[1] = X_n[jm+1];
        yyf(0)= yn(jm);
        yyf(1) = yn(jm+1);

        llm = ll_probit_contrib(Xnf,yyf, theta, mod, approx_method);
        ll += llm.ll;
        grad += llm.grad;
        Hess2 += prob*llm.grad * llm.grad.transpose();
        Hess += (prob* llm.Hess);
      }


      gradf += (prob * grad);
      Hess1 += (prob* grad * grad.transpose());
    }

  }

  gradf= gradf/prob_full;
  Hess = Hess/prob_full;
  Hess1 = Hess1/prob_full;
  Hess2 = Hess2/prob_full;

  //////////////////////
  // return results
  //////////////////////
  out = ll;
  out.attr("prob") = (pr.array());
  out.attr("gradient") = (gradf.array());
  //out.attr("grad_m") = (grad_m.array());
  out.attr("hessian1") = (Hess1.array());
  out.attr("hessian") =(Hess.array());
  out.attr("hessian2") =(Hess2.array());
  return out;
}


//' Calculation of contribution of one person to the approximated probit log likelihood calculation.
//' @description
//' Computes the contribution of one person to the log likelihood.
//' @param X_n
//' X_n list of regressor matrices
//' @param y_n
//' vector of choices
//' @param theta
//' parameter vector
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param approx_method
//' string specifying the approximation method
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector ll_probit_contrib_R(Rcpp::List X_n, Eigen::VectorXd y_n, Eigen::VectorXd theta, Rcpp::List mod, std::string approx_method){

  Contrib llm;
  NumericVector out;

  llm = ll_probit_contrib(X_n, y_n,theta, mod, approx_method);


  out = llm.ll;
  out.attr("gradient") = (llm.grad.array());
  //out.attr("grad_m") = (grad_m.array());
  out.attr("hessian") =(llm.Hess.array());

return out;
}
