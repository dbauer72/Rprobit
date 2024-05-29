// // // // // // // // // // // // // // // // //
//                                              //
//    PROBIT CML WITH APPROXIMATED MVNCDF        //
//                                              //
// // // // // // // // // // // // // // // // //
// Does not include the prediction part and hence uses
// data that has been treated before to subtract
// the X-part for the actual choice.

#include <RcppEigen.h>
//#include <omp.h>

using namespace Rcpp;



// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision


#include "distrib.h"
#include "SJ_vdb.h"
#include "TVBS_vdb.h"
#include "TVBS.h"
#include "ME.h"

#include "macml_grad_Hessian.h"
#include "toms462.h"                                         // allows bivariate normal cdf calculations


// provide vectorization of upper diagonal of a correlation matrix.
Eigen::VectorXd vectorize_lowerdiag(Eigen::MatrixXd R, int n){

  int cur=0;
  int N = n*(n+1)/2; // full number of parameters
  Eigen::VectorXd vecR(N);
  vecR.setZero();

  for (int jj=0;jj<n;jj++){
    vecR.segment(cur,jj+1) = R.block(jj,0,1,jj+1).transpose();
    cur += jj+1;
  }

  return vecR;
}


//' generate identity matrix of size 'alt' that adds a zero y'th row.
//' @description
//' generate identity matrix of size 'alt' that adds a zero y'th row.
//' @param y
//' integer; row to zero out
//' @param alt
//' integer; numbre of alternatives
//' @return
//' matrix; altxalt
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd makeAM(int y,int alt){

  Eigen::MatrixXd AM(alt+1,alt);
  AM.setZero();

  for (int jj=0; jj<alt;jj++){
    if (y>jj+1){
      AM(jj,jj) = 1;
    }
    if (y<=jj+1){
      AM(jj+1,jj)=1;
    }
  }

  return AM;
}

// // // // // // // // // // // // // // // // //
//                                              //
//    LIKELIHOOD CONTRIBUTION                   //
//                                              //
// // // // // // // // // // // // // // // // //

// input:
//  - theta           VectorXd    parameter vector (to optimize)
//  - data            List        observed data (one list entry for each decision maker)
//  - approx_method   string      approximation method for the MVNCDF ("SJ", "ME" or "TVBS")
//  - cml_pair_type   int         0=full pairwise, 1=adjacent pairs looped, 2=adjacent pairs chained (not yet correctly implemented -> needs weights)

//' Calculation of of approximated composite likelihood calculation.
//' @description
//' Computes the CML.
//' @param theta
//' parameter vector
//' @param data_obj
//' \link{data_cl} object
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the probit estimation.
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector ll_macml_norm_est(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control)
{

  //............................................................//
  // UNPACK ESTIMATION CONTROLS
  //............................................................//

  std::string approx_method = "SJ";
  int cml_pair_type = 1;
  int hess = 0;
  int el = 0;
  Rcpp::List pairs_list;

  if(control.containsElementNamed("approx_method")){
    CharacterVector approx_method_temp = control["approx_method"];
    approx_method = Rcpp::as<std::string>(approx_method_temp);
  }

  if(control.containsElementNamed("control_weights")){
    if(control["control_weights"] != R_NilValue){
      Rcpp::List control_weights = control["control_weights"];
      if(control_weights.containsElementNamed("cml_pair_type")){
        Rcpp::List cml_pair_type_list = control_weights["cml_pair_type"];
        if(cml_pair_type_list.containsElementNamed("pair_type")){
          cml_pair_type = cml_pair_type_list["pair_type"];
        }
      }
    }
  }

  if(control.containsElementNamed("hess")){
    hess = control["hess"];
  }

  if(control.containsElementNamed("el")){
    el = control["el"];
  }

  if(control.containsElementNamed("pairs_list")){
    if(control["pairs_list"] != R_NilValue){
      pairs_list = control["pairs_list"];
    } else {
      pairs_list = Rcpp::List::create(Named("NP_LIST") = 0);
    }
  } else {
    pairs_list = Rcpp::List::create(Named("NP_LIST") = 0);
  }

  if (pairs_list.size() == 0){
    pairs_list = Rcpp::List::create(Named("NP_LIST") = 0);
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

  if(approx_method != "SJ" && approx_method != "ME" && approx_method != "TVBS" && approx_method != "TVBSv2"){
    approx_method = "SJ";
  }

  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//

  int N = data_obj["N"];
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
  
  Eigen::VectorXd Tp = data_obj["Tp"];
  int Tp1 = Tp(0);
  
  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//

  NumericVector out;

  double ll = 0;

  int alt2m1 = 2*(alt-1);
  int Nrows = 0;
  for (int jj =0;jj<N;jj++){
    Nrows += Tp(jj)*(Tp(jj)-1)/2;
  }

  int pars = alt2m1 + (alt2m1+2)*(alt2m1+3)/2 + 3;
  Eigen::MatrixXd prob_calcs(Nrows,pars);
  prob_calcs.setZero();


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
    //Sigma_chol_vech = ((HL*l) + fL);
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

  int alt2 = alt+alt-2;
  Eigen::MatrixXd L_JOm = elimmat(alt2); //only for pairwise with equal number of alts for all Decision Makers
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


  // freeth denotes the number of free parameters for b, LO and LSig.
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

  int cur_lp =0;
  int cur_norm_prob= 0;

  for(int i_n = 0; i_n < N; i_n++){ //loop over the sample

    // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker
    double ll_n = 0;

    Eigen::VectorXd grad_n(lthb+lthO+lthL);
    grad_n.setZero();

    Eigen::MatrixXd Hess_n(freeth,freeth);
    Hess_n.setZero();

    Eigen::MatrixXd Hess_approx_n(lthb+lthO+lthL,lthb+lthO+lthL);
    Hess_approx_n.setZero();


    // save the data of the current decision maker
    Rcpp::List data = data_obj["data"];
    Rcpp::List data_n = data[i_n]; //unpack the data
    Rcpp::List X_n = data_n["X"];
    Eigen::VectorXd y_n = data_n["y"];

    //............................................................//
    // CML SETUP
    //............................................................//



    //............................................................//
    // START: Build the pairs structure
    //............................................................//

    // determine the panel length
    int Tp_n = y_n.rows();
    int lth_pairs_n = (Tp_n-1)*Tp_n/2;

    // Define a matrix that allows to specify the pairs (possible extension to higher order CMLs)
    Eigen::MatrixXd pairs_n;
    int pairs_check = 0;
    
    // Check if pairs list was provided
    if(pairs_list.containsElementNamed("NP_LIST")){
      pairs_check = 0;
    } else {
      Eigen::MatrixXd pairs_n_from_list = pairs_list[i_n];
      if(pairs_n_from_list.cols()>2){
        pairs_n = pairs_n_from_list;
        lth_pairs_n = pairs_n.rows();
        pairs_check = 1;
        
      }
    }
    
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


    //............................................................//
    // END: Build the pairs structure
    //............................................................//

    //............................................................//
    // deal with observations, where only one choice is observed  //
    //............................................................//

    if (lth_pairs_n == 0){

      //save the data of the current decision maker
      int alt2 = (alt-1);

      // setup X, M, M_Sig.
      Eigen::MatrixXd X(alt2,lthbb);
      X.setZero();

      Eigen::MatrixXd M_Sig(alt2,alt2);
      M_Sig.setZero();

      //............................................................//
      // Data Prep and Lin Utilities
      //............................................................//
      //


      //  extract the chosen alternatives for the current pair of choice occasions
      int y1;
      y1 = y_n[0];

      // extract the covariance matrix for the current pair of choice occasions
      X = X_n[0];


      Eigen::MatrixXd M1 = makeM((alt),(y1-1));
      Eigen::MatrixXd M1t = M1.transpose();
      Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;

      //............................................................//
      // Covariance Matrix
      //............................................................//

      M_Sig = Sigma_norm_1;

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


      // unstandardized lin-utilties (upper limits of the MVNCDF)
      Eigen::VectorXd x_upper = (-1)*(X_red.block(0,0,cur,X_red.cols())*b);

            // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
      Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
      Eigen::MatrixXd X_RE_t = X_RE.transpose();

      // // Compute combined covariance matrix
      // // Changed Sig to Lambda for combined Covariance
      Eigen::MatrixXd Lambda = M_Sig_red.block(0,0,cur,cur);
      if(lRE!=0){
        Lambda = Lambda + X_RE*Omega*X_RE_t;
      }

      Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
      Eigen::MatrixXd lambda_t = lambda.transpose();



      // ............................................................//
      //  Normalize upper limits and covariance matrix (to correlation matrix)
      // ............................................................//
      //
      Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
      Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
      Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();

      //............................................................//
      // Likelihood and Gradient
      //............................................................//
      //
      //


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
      ll_n -= grad_output[lth_output-1];
      //............................................................//
      // Gradient
      //............................................................//
      //

      Eigen::VectorXd grad_MVNCDF_approx;
      grad_MVNCDF_approx = grad_output.head((lth_output-1));

      //
      // initialize gradient contribution of the current pair
      Eigen::VectorXd grad_n_pair(lthb+lthO+lthL);
      grad_n_pair.setZero();
      //
      //............................................................//
      // START: Manipulations to go from MVNCDF approximation gradient
      //        to gradient of the original problem
      //............................................................//
      //
      // // calculate the Jacobian matrix of w,vech(Lambda_cor) with respect to b and vech(Lambda)
      Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta(X_red.block(0,0,cur,X_red.cols()), b, Lambda);

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
          grad_n_pair << gradb;
        } else {
          grad_n_pair << gradb,gradO;
        }
      } else {
        if(lthO == 0){
          grad_n_pair << gradb, gradL;
        } else {
          grad_n_pair << gradb,gradO, gradL;
        }
      }


      //............................................................//
      // END of gradient manipulations
      //............................................................//

      // add gradient of the current pair to the gradient of the decicion maker
      grad_n -= grad_n_pair;

      //............................................................//
      //Hessian
      //............................................................//
      // Hessian estimation based on gradients 
      Hess_approx_n  = grad_n_pair * grad_n_pair.transpose();
      Hess_approx += Hess_approx_n;//is positive throughout because it is only used this way 
      Hess_approx2 += Hess_approx_n;//is positive throughout because it is only used this way 

      if(hess == 1){
           //............................................................//
          // exact hessian
          //............................................................//

          //First step
          Eigen::MatrixXd hes;
          if (approx_method == "SJ"){
            hes  = SJ_hess_new(x_norm, Lambda_cor);
          }
          if (approx_method == "ME") {
            hes  = ME_hess_new(x_norm, Lambda_cor);
          }
          if (approx_method == "TVBS") {
            hes  = TVBS_hess_new(x_norm, Lambda_cor);
          }
          if (approx_method == "TVBSv2") {
            hes  = TVBS_hess_new(x_norm, Lambda_cor);
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

          Eigen::MatrixXd Hess_n_pair = J_chol.transpose()*H_bOmSig*J_chol+ H_t;

          // update Hessian
          //

          Hess_n -= Hess_n_pair;

          //............................................................//
          // approximated hessian
          //............................................................//


      } //end if hessian

    }


    // NOW LOOP OVER THE PAIRS MATRIX to compute the CML
    for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){


      // only calculate the pair if the weight is !=0
      if(pairs_n(i_pair,2)!=0){
        //............................................................//
        // Data Prep and Lin Utilities
        //............................................................//

        // extract the chosen alternatives for the current pair of choice occasions
        int y1, y2;
        y1 = y_n[pairs_n(i_pair,0)];
        y2 = y_n[pairs_n(i_pair,1)];

        Eigen::VectorXd y12(2);
        y12(0)= y1;
        y12(1)= y2;


        // extract the covariance matrix for the current pair of choice occasions
        Eigen::MatrixXd X1, X2;
        X1 = X_n[pairs_n(i_pair,0)];
        X2 = X_n[pairs_n(i_pair,1)];


        // Joint matrix for pairs
        int alt1 = X1.rows(); // Not needed yet but here we are flexible w.r.t to the number of choice alternatives in each time period
        // do not expect this flexiblity beyond this point

        int alt2 = X2.rows();
        Eigen::MatrixXd X_pair_norm((alt1+alt2),lthbb);
        X_pair_norm.setZero();
        X_pair_norm << X1, X2;

        //............................................................//
        // Covariance Matrix
        //............................................................//

        // use the matrices M1 and M2 to obtain the covariance matrix Sigma after normalizing for the chosen alternative
        Eigen::MatrixXd M1 = makeM((alt1+1),(y1-1));
        Eigen::MatrixXd M2 = makeM((alt2+1), (y2-1));
        Eigen::MatrixXd M1t = M1.transpose();
        Eigen::MatrixXd M2t = M2.transpose();

        Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;
        Eigen::MatrixXd Sigma_norm_2 = M2*Sigma*M2t;


        // compute joint Covariance matrix for the current pair
        Eigen::MatrixXd Sigma_norm_pair((alt1+alt2), (alt1+alt2));
        Sigma_norm_pair.setZero();
        Sigma_norm_pair.block(0,0,alt1,alt1) = Sigma_norm_1;
        Sigma_norm_pair.block(alt1,alt1,alt2,alt2) = Sigma_norm_2;

        Eigen::MatrixXd Sigma_norm_full((alt1+alt2+2), (alt1+alt2+2));
        Sigma_norm_full.setZero();
        Sigma_norm_full.block(0,0,alt1+1,alt1+1) = Sigma;
        Sigma_norm_full.block(alt1+1,alt1+1,alt2+1,alt2+1) = Sigma;

        // find the NA elements

        Eigen::MatrixXd X_red =  X_pair_norm;
        Eigen::MatrixXd M_Sig_red =  Sigma_norm_pair;

        IntegerVector ind = NA_ind(X_red.block(0,0,X_red.rows(),1));


        int cur = ind.size();

        for (int ij=0;ij<cur;ij++){
          X_red.row(ij) = X_pair_norm.row(ind[ij]);
          M_Sig_red.row(ij) = Sigma_norm_pair.row(ind[ij]);
        }

        for (int ij=0;ij<cur;ij++){
          M_Sig_red.col(ij) = M_Sig_red.col(ind[ij]);
        }

        // // unstandardized lin-utilties (upper limits of the MVNCDF)

        Eigen::VectorXd x_upper = (-1)*(X_red.block(0,0,cur,X_pair_norm.cols())*b);


        // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
        Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
        Eigen::MatrixXd X_RE_t = X_RE.transpose();

        // Compute combined covariance matrix
        // Changed Sig to Lambda for combined Covariance
        Eigen::MatrixXd Lambda = M_Sig_red.block(0,0,cur,cur);
        if(lRE!=0){
          Lambda = Lambda + X_RE*Omega*X_RE_t;
        }

        // fill in full variance formula
        Eigen::MatrixXd AM1(alt1+1,alt1), AM2(alt2+1,alt2);
        AM1.setZero();
        AM2.setZero();

        AM1 = makeAM(y1,alt1);
        AM2 = makeAM(y2,alt2);

        Eigen::MatrixXd AM(alt1+alt2+2,alt1+alt2);
        AM.setZero();
        AM.block(0,0,alt1+1,alt1) = AM1;
        AM.block(alt1+1,alt1,alt2+1,alt2) = AM2;


        Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
        Eigen::MatrixXd lambda_t = lambda.transpose();

        //............................................................//
        // Normalize upper limits and covariance matrix (to correlation matrix)
        //............................................................//

        Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
        Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
        Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();

        //............................................................//
        // write out the current bn and R and choice pair             //
        //............................................................//

        int sizealt = alt2m1+2;
        int params = alt2m1+sizealt*(sizealt+1)/2;
        prob_calcs.block(cur_norm_prob,0,1,alt2m1) = x_upper.transpose();
        prob_calcs.block(cur_norm_prob,alt2m1,1,sizealt*(sizealt+1)/2) = vectorize_lowerdiag(Sigma_norm_full,sizealt).transpose();
        prob_calcs(cur_norm_prob,params)=y1;
        prob_calcs(cur_norm_prob,params+1)=y2;



        //............................................................//
        // Likelihood and Gradient
        //............................................................//

        double lp=0;
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
        lp =  grad_output[lth_output-1];

        ll_n -= lp*pairs_n(i_pair,2);
        prob_calcs(cur_norm_prob,params+2)=lp;

        cur_norm_prob += 1;

        //............................................................//
        // Gradient
        //............................................................//
        //

        Eigen::VectorXd grad_MVNCDF_approx;
        grad_MVNCDF_approx = grad_output.head((lth_output-1));

        // initialize gradient contribution of the current pair
        Eigen::VectorXd grad_n_pair(lthb+lthO+lthL);
        grad_n_pair.setZero();



        //............................................................//
        // START: Manipulations to go from MVNCDF approximation gradient
        //        to gradient of the original problem
        //............................................................//


        // // calculate the Jacobian matrix of w,vech(Lambda_cor) with respect to b and vech(Lambda)
        Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta(X_red.block(0,0,cur,X_red.cols()), b, Lambda);

        // // partial log prob / partial (b, vech(Lambda))
        Eigen::VectorXd g_bvechLambda = grad_MVNCDF_approx.transpose()*J_wLambda_cor_bvechLambda; //column to row

        // partial (b, vech(Lambda), vech(X Omega X')) / partial (b, vech(Omega), vech(Sigma))
        Eigen::MatrixXd L = elimmat(cur);

        Eigen::MatrixXd J_dM_dbOmSig =  J_M_bOmSig_red(alt1+1,ind, y12, lthbb, X_RE.block(0,0,cur,X_RE.cols()),L ,D_JOm);

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
            grad_n_pair << gradb;
          } else {
            grad_n_pair << gradb,gradO;
          }
        } else {
          if(lthO == 0){
            grad_n_pair << gradb, gradL;
          } else {
            grad_n_pair << gradb,gradO, gradL;
          }
        }

        //............................................................//
        // END of gradient manipulations
        //............................................................//



        // add gradient of the current pair to the gradient of the decicion maker
        grad_n -= grad_n_pair*pairs_n(i_pair,2);  // added the weight

        
        // Hessian estimation based on gradients 
        Hess_approx_n  = grad_n_pair * grad_n_pair.transpose();
        Hess_approx2 += pairs_n(i_pair,2)*Hess_approx_n;//is positive throughout because it is only used this way 


        //............................................................//
        //Hessian
        //............................................................//

        if(hess == 1){

            //............................................................//
            // exact hessian
            //............................................................//

            //First step
            Eigen::MatrixXd hes;
            if (approx_method == "SJ"){
              hes  = SJ_hess_new(x_norm, Lambda_cor);
            }
            if (approx_method == "ME") {
              hes  = ME_hess_new(x_norm, Lambda_cor);
            }
            if (approx_method == "TVBS") {
              hes  = TVBS_hess_new(x_norm, Lambda_cor);
            }
            if (approx_method == "TVBSv2") {
              hes  = TVBS_hess_new(x_norm, Lambda_cor);
            }

            Eigen::MatrixXd H_th = H_theta(X_red.block(0,0,cur,X_red.cols()), b, Lambda);

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

            Eigen::MatrixXd Hess_n_pair = J_chol.transpose()*H_bOmSig*J_chol+ H_t;

            // update Hessian
            Hess_n -= Hess_n_pair*pairs_n(i_pair,2);  // added the weight



        } //end if hessian


      } // end if-weight!=0

    } //end pairs-loop


    grad += grad_n;//it is already negative in the pairs loop
    ll += ll_n; //it is already negative in the pairs loop
   
   Hess_approx_n  = grad_n * grad_n.transpose();
   Hess_approx += Hess_approx_n;//is positive throughout because it is only used this way 
   
   if (hess == 1){



      //Hess_approx_n  = grad_n * grad_n.transpose();

      Eigen::MatrixXd J_H = Eigen::MatrixXd(lthb+lthO+lthL,freeth);
      J_H.setZero();

      J_H.block(0,0,lthb,lthbb) = Hb.transpose();
      J_H.block(lthb,lthbb,lthO,lth_Omega_chol_vech) = HO.transpose();
      J_H.block(lthb+lthO,lthbb+lth_Omega_chol_vech,lthL,lth_Sigma_chol_vech) = HL.transpose();

      Hess += J_H * Hess_n * J_H.transpose();//it is already negative in the pairs loop
      //Hess_approx += Hess_approx_n;//is positive throughout because it is only used this way 
    }


  } //end N-loop


  //............................................................//
  // Format the output
  //............................................................//
  out = ll;
  out.attr("gradient") = (grad.array());
  out.attr("prob_calcs") = prob_calcs.array();
  out.attr("hessian2") =(Hess_approx2.array());
  out.attr("hessian1") =(Hess_approx.array());
 

  if(hess == 1){
    out.attr("hessian") =(Hess.array());
  }

  return out;

}







