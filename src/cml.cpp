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

#include "SJ_vdb.h"
#include "TVBS.h"
#include "TVBS_vdb.h"
#include "ME.h"
#include "macml_grad_Hessian.h"
#include "distrib.h"

#include "toms462.h"                                         // allows bivariate normal cdf calculations

// // // // // // // // // // // // // // // // // // // // // //
//                                                             //
//    marginal CML                                             //
//                                                             //
// // // // // // // // // // // // // // // // // // // // // //

// input:
//  - theta           VectorXd    parameter vector (to optimize)
//  - data            List        observed data (one list entry for each decision maker)
//  - approx_method   string      approximation method for the MVNCDF ("SJ", "ME" or "TVBS")
//  - cml_pair_type   int         0=full pairwise, 1=adjacent pairs looped, 2=adjacent pairs chained (not yet correctly implemented -> needs weights)

//' Calculation of Approximated Composite Marginal Likelihood
//'
//' @description
//' This function computes the marginal likelihood.
//'
//' @param theta
//' parameter vector
//' @param data_obj
//' data_cl object
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the probit estimation.
//'
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector ll_macml_marginal(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   int cml_pair_type = 0;
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
   
   // if a different cml_pair_type is used: For including the marginals we added 10 -> substract here. 
   if (cml_pair_type > 9){
     cml_pair_type -= 10;
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
   
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   NumericVector out;
   Eigen::VectorXd llEL(N);
   Rcpp::List gradEL(N);
   Rcpp::List HessEL(N);
   
   double ll = 0;
   //Eigen::MatrixXd prob_predict_pairwise(alt,alt);
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
   
   //int alt2 = alt+alt-2;
   //Eigen::MatrixXd L_JOm = elimmat(alt2); //only for pairwise with equal number of alts for all Decision Makers
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
   
   int cur_lp =1;
   double ll_n = 0;
   
   for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
     
     // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker
     ll_n = 0.0;
     
     Eigen::VectorXd grad_n(lthb+lthO+lthL);
     grad_n.setZero();
     
     Eigen::MatrixXd Hess_n(freeth,freeth);
     Hess_n.setZero();
     
     Eigen::MatrixXd Hess_approx_n(lthb+lthO+lthL,lthb+lthO+lthL);
     Hess_approx_n.setZero();
     
     
     // save the data of the current decision maker
     Rcpp::List data = data_obj["data"]; // extract the data from the data object
     Rcpp::List data_n = data[i_n]; //unpack the data
     Rcpp::List X_n = data_n["X"];
     Eigen::VectorXd y_n = data_n["y"];
     
     
     
     //............................................................//
     // CML SETUP
     //............................................................//
     
     // determine the panel length
     int Tp_n = y_n.rows();
     
     //............................................................//
     // END: Build the pairs structure
     //............................................................//
     
     //............................................................//
     // marginals                                              //
     //............................................................//
     
     for (int mj = 0; mj < Tp_n;mj++){
       
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
       y1 = y_n[mj];
       
       // extract the covariance matrix for the current pair of choice occasions
       X = X_n[mj];
       
       
       Eigen::MatrixXd M1 = makeM((alt),(y1-1));
       Eigen::MatrixXd M1t = M1.transpose();
       Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;
       
       //............................................................//
       // Covariance Matrix
       //............................................................//
       
       M_Sig = Sigma_norm_1;
       
       // find the NA elements
       
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
       //
       // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
       Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
       Eigen::MatrixXd X_RE_t = X_RE.transpose();
       //
       
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
       
       //
       //............................................................//
       // Gradient
       //............................................................//
       //
       
       Eigen::VectorXd grad_MVNCDF_approx;
       grad_MVNCDF_approx = grad_output.head(lth_output-1);
       
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
       
       // add gradient of the current pair to the gradient of the decision maker
       grad_n -= grad_n_pair;
       
       //............................................................//
       //Hessian
       //............................................................//
       Hess_approx2 += grad_n_pair * grad_n_pair.transpose();
       if(hess == 1){
         
         //............................................................//
         // exact hessian
         //............................................................//
         
         //First step
         Eigen::MatrixXd hes;
         
         if (approx_method == "SJ"){
           hes  = SJ_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "ME"){
           hes  = ME_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "TVBSv2"){
           hes  = TVBS_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "TVBS"){
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
         Hess_n -= Hess_n_pair;
         
       } //end if hessian
       
     }
     
     grad += grad_n;//it is already negative in the pairs loop
     ll += ll_n; //it is already negative in the pairs loop
     
     if(el == 1){
       llEL[i_n] = ll_n;
       gradEL[i_n] = grad_n;
     }
     
     Hess_approx_n  = grad_n * grad_n.transpose();
     Hess_approx += Hess_approx_n;//is positive throughout because it is only used this way
     
     if(el == 1){
       HessEL[i_n] = Hess_approx_n;
     }
     
     if (hess == 1){
       
       Eigen::MatrixXd J_H = Eigen::MatrixXd(lthb+lthO+lthL,freeth);
       J_H.setZero();
       
       J_H.block(0,0,lthb,lthbb) = Hb.transpose();
       J_H.block(lthb,lthbb,lthO,lth_Omega_chol_vech) = HO.transpose();
       J_H.block(lthb+lthO,lthbb+lth_Omega_chol_vech,lthL,lth_Sigma_chol_vech) = HL.transpose();
       
       Hess += J_H * Hess_n * J_H.transpose();//it is already negative in the pairs loop
       
       if(el == 1){
         HessEL[i_n] = J_H * Hess_n * J_H.transpose();
       }
     }
     
   } //end N-loop
   
   
   //............................................................//
   // Format the output
   //............................................................//
   out = ll;
   out.attr("gradient") = (grad.array());
   out.attr("hessian2") =(Hess_approx2.array());
   out.attr("hessian1") =(Hess_approx.array());
   
   if(hess == 1){
     out.attr("hessian") =(Hess.array());
   }
   if(el == 1){
     out.attr("llEL")   = llEL;
     out.attr("gradEL") = gradEL;
     out.attr("HessEL") = HessEL;
   }
   
   return out;
 }
 
 
 
 
 // // // // // // // // // // // // // // // // //
 //                                              //
 //    LIKELIHOOD CONTRIBUTION                   //
 //                                              //
 // // // // // // // // // // // // // // // // //
 
 // input:
 //  - theta           VectorXd        parameter vector (to optimize)
 //  - data_obj        data_obj        observed data 
 //  - mod             mod_cl object   describes the model 
 //  - control         control_cl object   contains info on various parameter, e.g. 0=full pairwise, 1=adjacent pairs looped, 2=adjacent pairs chained (not yet correctly implemented -> needs weights)
 

//' Calculation of Approximated Composite Marginal Likelihood
//'
//' @description
//' This function computes the composite marginal likelihood.
//'
//' @param theta
//' parameter vector
//' @param data_obj
//' data_cl object
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the probit estimation.
//'
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and 
//' (if specified in the controls) the Hessian.
//' For el==1 the log CML contribution as well as gradients for each decider is contained in the output.
//' @export
//'
// [[Rcpp::export]]
NumericVector ll_macml(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   int cml_pair_type = 0;
   int hess = 0;
   int el = 0;
   Rcpp::List pairs_list;
   
   //Rcout << "start " <<  std::endl; 
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
   
   // initialize output 
   NumericVector out;
   
   //Rcout << cml_pair_type << std::endl; 
   
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
   
   //Rcout << "N:" << N << std::endl; 
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   Eigen::VectorXd llEL(N);
   llEL.setZero();
   Rcpp::List gradEL(N);
   Rcpp::List HessEL(N);
   
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
   
   int cur_lp =1;
   double ll_n = 0;
   
   for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
     
     
     // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker
     ll_n = 0.0;
     
     Eigen::VectorXd grad_n(lthb+lthO+lthL);
     grad_n.setZero();
     
     Eigen::MatrixXd Hess_n(freeth,freeth);
     Hess_n.setZero();
     
     Eigen::MatrixXd Hess_approx_n(lthb+lthO+lthL,lthb+lthO+lthL);
     Hess_approx_n.setZero();
     
     
     // save the data of the current decision maker
     Rcpp::List data = data_obj["data"]; // extract the data from the data object
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
     
     //// normalize to obtain similar likelihood as in full probit likelihood calculations. //
     //// reweight by Tp/sum of weights*2. 
     //
     //// sum of weights 
     //double sum_of_weights = 0;
     //
     //for (int js=0;js<pairs_n.rows();js++){
     //  sum_of_weights += pairs_n(js,2);
     //}
     //
     //sum_of_weights = sum_of_weights * 2; 
     //if (sum_of_weights == 0){ 
     //  sum_of_weights = 1.0;
     //}
     //
     //for (int js=0;js<pairs_n.rows();js++){
     //  pairs_n(js,2) = pairs_n(js,2)*Tp_n / sum_of_weights;
     //}
     
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
       //
       // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
       Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
       Eigen::MatrixXd X_RE_t = X_RE.transpose();
       //
       
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
       
       Eigen::VectorXd grad_MVNCDF_approx;
       grad_MVNCDF_approx = grad_output.head(lth_output-1);
       
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
       
       // add gradient of the current pair to the gradient of the decision maker
       grad_n -= grad_n_pair;
       
       //............................................................//
       //Hessian
       //............................................................//
       Hess_approx2 += grad_n_pair * grad_n_pair.transpose();
       if(hess == 1){
         
         //............................................................//
         // exact hessian
         //............................................................//
         
         //First step
         Eigen::MatrixXd hes;
         
         if (approx_method == "SJ"){
           hes  = SJ_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "ME"){
           hes  = ME_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "TVBSv2"){
           hes  = TVBS_hess_new(x_norm, Lambda_cor);
         }
         if (approx_method == "TVBS"){
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
         Hess_n -= Hess_n_pair;
         
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
         // TODO: if we want to allow different numbers of choice alternatives for different time points, this would be rather difficult, because we would also have to know which
         //        choice alternatives have to be excluded to cut the Lambda matrix and the theta vector accordingly
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
         
         Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
         Eigen::MatrixXd lambda_t = lambda.transpose();
         
         //............................................................//
         // Normalize upper limits and covariance matrix (to correlation matrix)
         //............................................................//
         
         Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
         Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
         Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
         
         //............................................................//
         // Likelihood and Gradient
         //............................................................//
         
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
         ll_n -= grad_output[lth_output-1]*pairs_n(i_pair,2);
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
         
         
         //............................................................//
         // Hessian
         //............................................................//
         
         Hess_approx2 += grad_n_pair * grad_n_pair.transpose()*pairs_n(i_pair,2);
         
         if(hess == 1){
           
           //............................................................//
           // exact hessian
           //............................................................//
           
           //First step
           Eigen::MatrixXd hes;
           if (approx_method == "SJ"){
             hes  = SJ_hess_new(x_norm, Lambda_cor);
           } else if (approx_method == "ME"){
             hes  = ME_hess_new(x_norm, Lambda_cor);
           } else { // same calculation for TVBS and TVBSv2
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
     
     if(el == 1){
       llEL[i_n] = ll_n;
       gradEL[i_n] = grad_n;
     }
     
     Hess_approx_n  = grad_n * grad_n.transpose();
     Hess_approx += Hess_approx_n;//is positive throughout because it is only used this way
     
     if(el == 1){
       HessEL[i_n] = Hess_approx_n;
     }
     
     if (hess == 1){
       
       Eigen::MatrixXd J_H = Eigen::MatrixXd(lthb+lthO+lthL,freeth);
       J_H.setZero();
       
       J_H.block(0,0,lthb,lthbb) = Hb.transpose();
       J_H.block(lthb,lthbb,lthO,lth_Omega_chol_vech) = HO.transpose();
       J_H.block(lthb+lthO,lthbb+lth_Omega_chol_vech,lthL,lth_Sigma_chol_vech) = HL.transpose();
       
       Hess += J_H * Hess_n * J_H.transpose();//it is already negative in the pairs loop
       
       if(el == 1){
         HessEL[i_n] = J_H * Hess_n * J_H.transpose();
       }
     }
     
   } //end N-loop
   
   
   //............................................................//
   // Format the output
   //............................................................//
   out = ll;
   out.attr("gradient") = (grad.array());
   out.attr("hessian2") =(Hess_approx2.array());
   out.attr("hessian1") =(Hess_approx.array());
   
   if(hess == 1){
     out.attr("hessian") =(Hess.array());
   }
   if(el == 1){
     out.attr("llEL") = llEL;
     out.attr("gradEL") = gradEL;
     out.attr("HessEL") = HessEL;
   }
   
   
   
   return out;
   
 }
 
 
 // // // // // // // // // // // // // // // // //
 //                                              //
 //    PROBIT LIKELIHOOD                         //
 //                                              //
 // // // // // // // // // // // // // // // // //
 
 // input:
 //  - theta           VectorXd    parameter vector (to optimize)
 //  - data            List        observed data (one list entry for each decision maker)
 //  - approx_method   string      approximation method for the MVNCDF ("SJ", "ME" or "TVBS")
 //  - hess           int          calculate the Hessian?
 

//' Calculation of of approximated probit log likelihood calculation.
//' @description
//' Computes the log likelihood.
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
//' @export
//'
// [[Rcpp::export]]
NumericVector ll_probit(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   int hess = 0;
   int el = 0;
   
   if(control.containsElementNamed("approx_method")){
     CharacterVector approx_method_temp = control["approx_method"];
     approx_method = Rcpp::as<std::string>(approx_method_temp);
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
   //Eigen::VectorXd Tp = data_obj["Tp"];
   Eigen::MatrixXd Hb = mod["Hb"];
   Eigen::VectorXd fb = mod["fb"];
   Eigen::VectorXd fL = mod["fL"];
   Eigen::MatrixXd HL = mod["HL"];
   Eigen::VectorXd fO = mod["fO"];
   Eigen::MatrixXd HO = mod["HO"];
   
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   NumericVector out;
   Rcpp::List out_prob_predict(1);
   
   Eigen::VectorXd llEL(N);
   llEL.setZero();
   Rcpp::List gradEL(N);
   Rcpp::List HessEL(N);
   
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
   Eigen::MatrixXd Hess_approx(lthb+lthO+lthL,lthb+lthO+lthL), Hess_approx2(lthb+lthO+lthL,lthb+lthO+lthL);
   Hess_approx.setZero();
   Hess_approx2.setZero();
   Eigen::VectorXd grad(lthb+lthO+lthL);
   grad.setZero();
   double ll_n=0.0;
   for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
     
     // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker
     ll_n = 0.0;
     
     Eigen::VectorXd grad_n(lthb+lthO+lthL);
     grad_n.setZero();
     
     Eigen::MatrixXd Hess_n(freeth,freeth);
     Hess_n.setZero();
     
     Eigen::MatrixXd Hess_approx_n(freeth,freeth);
     Hess_approx_n.setZero();
     
     
     //save the data of the current decision maker
     
     Rcpp::List data = data_obj["data"]; // extract the data from the data object
     Rcpp::List data_n = data[i_n]; //unpack the data
     Rcpp::List X_n = data_n["X"];
     Eigen::VectorXd y_n = data_n["y"];
     
     int Tp_n = y_n.rows();
     int alt2 = Tp_n*(alt-1);
     
     // setup X, M, M_Sig.
     Eigen::MatrixXd X(alt2,lthbb);
     X.setZero();
     
     Eigen::MatrixXd M_Sig(alt2,alt2);
     M_Sig.setZero();
     
     
     // NOW LOOP OVER THE observations MATRIX
     for(int jj = 0; jj < Tp_n; jj++){ // loop jj over observations.
       
       //............................................................//
       // Data Prep and Lin Utilities
       //............................................................//
       //
       //  extract the chosen alternatives for the current pair of choice occasions
       int y1;
       y1 = y_n[jj];
       
       // extract the covariance matrix for the current pair of choice occasions
       Eigen::MatrixXd X1;
       X1 = X_n[jj];
       
       X.block((jj)*(alt-1),0,alt-1,lthbb) = X1;
       
       Eigen::MatrixXd M1 = makeM((alt),(y1-1));
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
     
     Eigen::VectorXd x_upper = (-1)*(X_red.block(0,0,cur,X.cols())*b);
     //
     //Rcout << "x_upper" << x_upper << std::endl;
     // // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
     Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
     Eigen::MatrixXd X_RE_t = X_RE.transpose();
     
     // Compute combined covariance matrix
     // Changed Sig to Lambda for combined Covariance
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
     ll_n = -grad_output[lth_output-1];
     
     //............................................................//
     // Gradient
     //............................................................//
     //
     
     Eigen::VectorXd grad_MVNCDF_approx;
     grad_MVNCDF_approx = grad_output.head(lth_output-1);
     
     // initialize gradient contribution of the current pair
     Eigen::VectorXd grad_n_pair(lthb+lthO+lthL);
     grad_n_pair.setZero();
     //
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
     
     // add gradient of the current pair to the gradient of the decision maker
     grad_n = (-1)* grad_n_pair;
     
     
     //............................................................//
     //Hessian approximation
     //............................................................//
     
     Hess_approx_n  = grad_n * grad_n.transpose();
     Hess_approx2 += Hess_approx_n;//is positive throughout because it is only used this way
     
     
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
       } else if (approx_method == "ME"){
         hes  = ME_hess_new(x_norm, Lambda_cor);
       } else { // TVBSv2 or TVBS
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
       Hess_n -= Hess_n_pair;  // added the weight
       
     } //end if hessian
     
     
     grad += grad_n;//it is already negative in the pairs loop
     
     ll += ll_n; //it is already negative in the pairs loop
     
     if(el == 1){
       llEL[i_n]= ll_n; 
       gradEL[i_n] = grad_n;
       HessEL[i_n] = Hess_approx_n;
     }
     
     // add multiplication with HO, HL und Hb.
     if (hess == 1){
       Eigen::MatrixXd J_H = Eigen::MatrixXd(lthb+lthO+lthL,freeth);
       J_H.setZero();
       
       J_H.block(0,0,lthb,lthbb) = Hb.transpose();
       J_H.block(lthb,lthbb,lthO,lth_Omega_chol_vech) = HO.transpose();
       J_H.block(lthb+lthO,lthbb+lth_Omega_chol_vech,lthL,lth_Sigma_chol_vech) = HL.transpose();
       
       Hess += J_H * Hess_n * J_H.transpose();//it is already negative in the pairs loop
       
       
       if(el == 1){
         HessEL[i_n] = J_H* Hess_n *J_H.transpose();
       }
       
     }
     
   } //end N-loop
   
   
   
   //............................................................//
   // Format the output
   //............................................................//
   out = ll;
   out.attr("gradient") = (grad.array());
   out.attr("hessian2") = (Hess_approx2.array());
   
   if(el == 1){
     out.attr("llEL") = llEL;
     out.attr("gradEL") = gradEL;
     out.attr("HessEL") = HessEL;
   }
   
   if(hess == 1){
     out.attr("hessian") = (Hess.array());
     out.attr("hessian1") = (Hess_approx2.array());
   }
   
   
   return out;
   
 }
 
 
 
 
 // // // // // // // // // // // // // // // // // // // // // // 
 //    calculating probabilities                                //
 // // // // // // // // // // // // // // // // // // // // // // 
 
 double cal_choice_probs_1_nograd(Eigen::MatrixXd X, int alt, int y1, Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::MatrixXd Omega, std::string approx_method)
 {
   double out = 0;
   
   int alt2 = (alt-1);
   Eigen::MatrixXd M_Sig(alt2,alt2);
   M_Sig.setZero();
   
   Eigen::MatrixXd M1 = makeM((alt),(y1-1));
   Eigen::MatrixXd M1t = M1.transpose();
   Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;
   
   int lthbb = X.cols(); 
   
   //............................................................//
   // Covariance Matrix
   //............................................................//
   
   M_Sig = Sigma_norm_1;
   
   // find the NA elements
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
   
   //
   // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
   int lRE = Omega.rows(); 
   
   Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
   Eigen::MatrixXd X_RE_t = X_RE.transpose();
   //
   
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
   
   
   if(approx_method == "SJ"){
     out = SJ(x_norm,Lambda_cor);
   }
   
   if(approx_method == "ME"){
     out = ME(x_norm,Lambda_cor);
   }
   
   if(approx_method == "TVBSv2"){
     out = TVBS_pv2(x_norm,Lambda_cor,1);
   }
   
   if(approx_method == "TVBS"){
     out = TVBS_p(x_norm,Lambda_cor,1);
   }
   
   return out;
 }
 
 Eigen::VectorXd cal_choice_probs_1(Eigen::MatrixXd X, int alt, int y1, Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::MatrixXd Omega, Eigen::MatrixXd J_chol, Eigen::MatrixXd Hb, Eigen::MatrixXd HO, Eigen::MatrixXd HL,int lth_Omega_chol_vech,int lth_Sigma_chol_vech, std::string approx_method)
 {
   double out = 0;
   
   int alt2 = (alt-1);
   Eigen::MatrixXd M_Sig(alt2,alt2);
   M_Sig.setZero();
   
   Eigen::MatrixXd M1 = makeM((alt),(y1-1));
   Eigen::MatrixXd M1t = M1.transpose();
   Eigen::MatrixXd Sigma_norm_1 = M1*Sigma*M1t;
   
   int lthbb = X.cols(); 
   
   //............................................................//
   // Covariance Matrix
   //............................................................//
   
   M_Sig = Sigma_norm_1;
   
   // find the NA elements
   
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
   
   //
   // compute covariance matrix Omega of the random effects for this pair of observations using the observed effects
   int lRE = Omega.rows(); 
   
   Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
   Eigen::MatrixXd X_RE_t = X_RE.transpose();
   //
   
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
   
   int ltheta = Hb.cols() + HO.cols() + HL.cols()+1; 
   
   Eigen::VectorXd grad_output(ltheta);
   grad_output.setZero();
   
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
   
   int lth_output = grad_output.size();
   
   //
   //............................................................//
   // Gradient
   //............................................................//
   //
   
   Eigen::VectorXd grad_MVNCDF_approx;
   grad_MVNCDF_approx = grad_output.head(lth_output-1);
   out = grad_output[lth_output-1];
   //
   // initialize gradient contribution of the current pair
   Eigen::VectorXd grad_n_pair(Hb.cols()+HO.cols()+HL.cols());
   grad_n_pair.setZero();
   //
   //............................................................//
   // START: Manipulations to go from MVNCDF approximation gradient
   //        to gradient of the original problem
   //............................................................//
   //
   // // calculate the Jacobian matrix of w,vech(Lambda_cor) with respect to b and vech(Lambda)
   Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta(X_red.block(0,0,cur,X_red.cols()),b, Lambda);
   
   // // partial log prob / partial (b, vech(Lambda))
   Eigen::VectorXd g_bvechLambda = grad_MVNCDF_approx.transpose()*J_wLambda_cor_bvechLambda; //column to row
   
   // partial (b, vech(Lambda), vech(X Omega X')) / partial (b, vech(Omega), vech(Sigma))
   Eigen::MatrixXd L = elimmat(cur);
   Eigen::MatrixXd D_JOm  = duplmat(lRE);
   
   Eigen::VectorXd y12(1);
   y12(0)= y1;
   
   Eigen::MatrixXd J_dM_dbOmSig =  J_M_bOmSig_red(alt,ind, y12, lthbb, X_RE.block(0,0,cur,X_RE.cols()),L ,D_JOm);
   
   // // partial (b, vech(Omega), vech(Sigma)) / partial (b, vech(L_Omega),vech(L_Sigma))
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
   Eigen::VectorXd grad_n(Hb.cols()+HO.cols()+HL.cols());
   grad_n.setZero();
   
   Eigen::VectorXd grad(Hb.cols()+HO.cols()+HL.cols()+1);
   grad.setZero();
   
   if(HL.cols() == 0){
     if(HO.cols() == 0){
       grad_n << gradb;
     } else {
       grad_n << gradb,gradO;
     }
   } else {
     if(HO.cols() == 0){
       grad_n << gradb, gradL;
     } else {
       grad_n << gradb,gradO, gradL;
     }
   }
   
   grad.segment(0,Hb.cols()+HO.cols()+HL.cols()) = -grad_n; 
   grad[Hb.cols()+HO.cols()+HL.cols()] = out; 
   return grad;
 }
 
 double cal_choice_probs_2_nograd(Eigen::MatrixXd X1, Eigen::MatrixXd X2, int y1, int y2, Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::MatrixXd Omega, std::string approx_method)
 {
   
   double out = 0.0; 
   // Joint matrix for pairs
   int alt1 = X1.rows(); // Not needed yet but here we are flexible w.r.t to the number of choice alternatives in each time period
   int alt2 = X2.rows();
   int lthbb = X1.cols(); 
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
   int lRE = Omega.rows(); 
   
   Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
   Eigen::MatrixXd X_RE_t = X_RE.transpose();
   
   // Compute combined covariance matrix
   // Changed Sig to Lambda for combined Covariance
   Eigen::MatrixXd Lambda = M_Sig_red.block(0,0,cur,cur);
   if(lRE!=0){
     Lambda = Lambda + X_RE*Omega*X_RE_t;
   }
   
   Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
   Eigen::MatrixXd lambda_t = lambda.transpose();
   
   //............................................................//
   // Normalize upper limits and covariance matrix (to correlation matrix)
   //............................................................//
   
   Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
   Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
   Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
   
   if(approx_method == "SJ"){
     out = SJ(x_norm,Lambda_cor);
   }
   
   if(approx_method == "ME"){
     out = ME(x_norm,Lambda_cor);
   }
   
   if(approx_method == "TVBSv2"){
     out = TVBS_pv2(x_norm,Lambda_cor,1);
   }
   
   if(approx_method == "TVBS"){
     out = TVBS_p(x_norm,Lambda_cor,1);
   }
   
   return out;
   
 }
 
 
 
 Eigen::VectorXd cal_choice_probs_2(Eigen::MatrixXd X1, Eigen::MatrixXd X2, int y1, int y2, Eigen::VectorXd b, Eigen::MatrixXd Sigma, Eigen::MatrixXd Omega, Eigen::MatrixXd J_chol, Eigen::MatrixXd Hb, Eigen::MatrixXd HO, Eigen::MatrixXd HL,int lth_Omega_chol_vech,int lth_Sigma_chol_vech, std::string approx_method)
 {
   
   double out = 0.0; 
   // Joint matrix for pairs
   int alt1 = X1.rows(); // Not needed yet but here we are flexible w.r.t to the number of choice alternatives in each time period
   int alt2 = X2.rows();
   int lthbb = X1.cols(); 
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
   int lRE = Omega.rows(); 
   
   Eigen::MatrixXd X_RE = X_red.block(0,0,cur, lRE);
   Eigen::MatrixXd X_RE_t = X_RE.transpose();
   
   // Compute combined covariance matrix
   // Changed Sig to Lambda for combined Covariance
   Eigen::MatrixXd Lambda = M_Sig_red.block(0,0,cur,cur);
   if(lRE!=0){
     Lambda = Lambda + X_RE*Omega*X_RE_t;
   }
   
   Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
   Eigen::MatrixXd lambda_t = lambda.transpose();
   
   //............................................................//
   // Normalize upper limits and covariance matrix (to correlation matrix)
   //............................................................//
   
   Eigen::VectorXd x_norm = x_upper.array()/lambda.array();
   Eigen::MatrixXd lambda_li_lj = lambda*lambda_t;
   Eigen::MatrixXd Lambda_cor = (Lambda.array()/lambda_li_lj.array()).matrix();
   
   int ltheta = Hb.cols() + HO.cols() + HL.cols()+1; 
   
   Eigen::VectorXd grad_output(ltheta);
   grad_output.setZero();
   
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
   
   int lth_output = grad_output.size();
   
   //
   //............................................................//
   // Gradient
   //............................................................//
   //
   
   Eigen::VectorXd grad_MVNCDF_approx;
   grad_MVNCDF_approx = grad_output.head(lth_output-1);
   
   out = grad_output[lth_output-1];  
   //
   // initialize gradient contribution of the current pair
   Eigen::VectorXd grad_n_pair(Hb.cols()+HO.cols()+HL.cols());
   grad_n_pair.setZero();
   //
   //............................................................//
   // START: Manipulations to go from MVNCDF approximation gradient
   //        to gradient of the original problem
   //............................................................//
   //
   // // calculate the Jacobian matrix of w,vech(Lambda_cor) with respect to b and vech(Lambda)
   Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta(X_red.block(0,0,cur,X_red.cols()),b, Lambda);
   
   // // partial log prob / partial (b, vech(Lambda))
   Eigen::VectorXd g_bvechLambda = grad_MVNCDF_approx.transpose()*J_wLambda_cor_bvechLambda; //column to row
   
   // partial (b, vech(Lambda), vech(X Omega X')) / partial (b, vech(Omega), vech(Sigma))
   Eigen::MatrixXd L = elimmat(cur);
   Eigen::MatrixXd D_JOm  = duplmat(lRE);
   
   Eigen::VectorXd y12(2);
   y12(0)= y1;
   y12(1)= y2;
   
   Eigen::MatrixXd J_dM_dbOmSig =  J_M_bOmSig_red(alt1+1,ind, y12, lthbb, X_RE.block(0,0,cur,X_RE.cols()),L ,D_JOm);
   
   // // partial (b, vech(Omega), vech(Sigma)) / partial (b, vech(L_Omega),vech(L_Sigma))
   //Eigen::MatrixXd J_chol = J_bvOvSig_bLOLSig(Omega_chol_vech, Sigma_chol_vech, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);
   
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
   // fill in the gradients, depending on which is present.
   Eigen::VectorXd grad_n(Hb.cols()+HO.cols()+HL.cols());
   grad_n.setZero();
   
   Eigen::VectorXd grad(Hb.cols()+HO.cols()+HL.cols()+1);
   grad.setZero();
   
   if(HL.cols() == 0){
     if(HO.cols() == 0){
       grad_n << gradb;
     } else {
       grad_n << gradb,gradO;
     }
   } else {
     if(HO.cols() == 0){
       grad_n << gradb, gradL;
     } else {
       grad_n << gradb,gradO, gradL;
     }
   }
   
   grad.segment(0,Hb.cols()+HO.cols()+HL.cols()) = -grad_n; 
   grad[Hb.cols()+HO.cols()+HL.cols()] = out; 
   return grad;
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
 
//' Choice probabilities for probit models
//' @description
//' Computes the approximate choice porbabilities for the probit case.
//' @param theta
//' parameter vector
//' @param data
//' data list
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the probit estimation.
//' @return
//' A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
//' @keywords internal
//'
// [[Rcpp::export]]
Eigen::MatrixXd pred_probit_approx(Eigen::VectorXd theta, Rcpp::List data, Rcpp::List mod, Rcpp::List control)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   
   if(control.containsElementNamed("approx_method")){
     CharacterVector approx_method_temp = control["approx_method"];
     approx_method = Rcpp::as<std::string>(approx_method_temp);
   }
   
   int num_class = 1; 
   if (mod.containsElementNamed("num_class")){
     num_class = mod["num_class"];
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
   
   int lthb = mod["lthb"];
   int lthO = mod["lthO"];
   int lthL = mod["lthL"];
   int lRE = mod["lRE"];
   int alt = mod["alt"];
   int param_one = lthb + lthO +lthL; 
   Eigen::MatrixXd Hb = mod["Hb"];
   Eigen::VectorXd fb = mod["fb"];
   Eigen::VectorXd fL = mod["fL"];
   Eigen::MatrixXd HL = mod["HL"];
   Eigen::VectorXd fO = mod["fO"];
   Eigen::MatrixXd HO = mod["HO"];
   
   //............................................................//
   // extract data and calculate number of choice occasions      //
   //............................................................//
   
   // save the data of the current decision maker
   Rcpp::List data_n = data; //unpack the data
   Rcpp::List X_n = data_n["X"];
   Eigen::VectorXd y_n = data_n["y"];
   int Tp_n = y_n.rows();
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   Eigen::MatrixXd prob_predict(Tp_n,alt);
   prob_predict.setZero();
   
   //............................................................//
   // RECOVER PARAMETERS
   //............................................................//
   
   Eigen::MatrixXd b(Hb.rows(),num_class);
   b.setZero();
   
   Eigen::MatrixXd Sigma(alt*num_class,alt);
   Sigma.setZero();
   
   Eigen::MatrixXd Sigma_chol_vech;
   Eigen::MatrixXd Sigma_chol, Sigma_chol_t;
   Eigen::MatrixXd Omega(lRE*num_class,lRE);
   Omega.setZero();
   
   Eigen::MatrixXd ELt =  elimmat(alt).transpose();
   Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
   Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t;
   
   int lth_Omega_chol_vech;
   int lth_Sigma_chol_vech;
   
   for (int jb=0; jb<num_class; jb++){
     b.col(jb) =  ((Hb*theta.segment(param_one*jb,lthb)) + fb);
     
     // Recover the Sigma matrix
     // Without free parameters only fL is used
     if(lthL==0){
       Sigma_chol_vech = fL;
     }else{
       //Sigma_chol_vech = ((HL*l) + fL);
       Sigma_chol_vech = ((HL*theta.segment(param_one*jb+(lthO+lthb), lthL)) + fL);
     }
     
     lth_Sigma_chol_vech = Sigma_chol_vech.rows();
     Sigma_chol = ELt*Sigma_chol_vech;
     Sigma_chol.resize(alt,alt);
     Sigma_chol_t = Sigma_chol.transpose();
     
     Sigma.block(alt*jb,0,alt,alt) = Sigma_chol*Sigma_chol_t;
     
     // Recover the matrix Omega
     
     // check if we have mixed parameters
     if(lRE!=0){
       // Without free parameters only fO is used
       if(lthO==0){
         Omega_chol_vech = fO;
       }else{
         //Omega_chol_vech = ((HO*o) + fO);
         Omega_chol_vech = ((HO*theta.segment(param_one*jb+(lthb), lthO)) + fO);
       }
       lth_Omega_chol_vech = Omega_chol_vech.rows();
       Omega_chol = EOt*Omega_chol_vech;
       Omega_chol.resize(lRE,lRE);
       Omega_chol_t = Omega_chol.transpose();
       Omega.block(lRE*jb,0,lRE,lRE) = Omega_chol*Omega_chol_t;
     } 
   }
   
   Eigen::VectorXd pi_est(num_class),gamma_est(num_class);
   pi_est.setZero();
   gamma_est.setZero();
   
   gamma_est.segment(1,num_class-1) = theta.segment(param_one*num_class, num_class-1);
   
   double sum_pi = 0.0;
   
   for (int jg =0; jg < num_class; jg++){
     pi_est[jg]=  std::exp(gamma_est[jg]);
     sum_pi += pi_est[jg];
   }
   
   for (int jg =0; jg < num_class; jg++){
     pi_est[jg]=  pi_est[jg]/sum_pi;
   }
   
   //............................................................//
   // END of setup
   //............................................................//
   
   //............................................................//
   // CYCLE over choices to calculate the  probs.
   //............................................................//
   for(int i_choice_1 = 0; i_choice_1 < Tp_n; i_choice_1++){
     
     // extract the covariance matrix for the current pair of choice occasions
     Eigen::MatrixXd X1 = X_n[i_choice_1];
     Eigen::MatrixXd X_red = X1;
     
     
     // cycle over chosen alternative
     for (int ij = 0; ij< alt; ij++){
       
       // diff the matrix by the choosen alternative
       Eigen::MatrixXd X1work= X_red;
       Eigen::VectorXd X1_chosen = X1work.row((ij));
       
       removeRow(X1work, ij); //cpp indexing, choices start at 1
       
       X1work.rowwise() -= X1_chosen.transpose();
       
       // now data is prepared ->
       double ll_n_est =0.0;
       for (int jg = 0; jg< num_class; jg++){
         ll_n_est +=  std::exp(cal_choice_probs_1_nograd(X1work, alt, ij+1, b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE),approx_method))* pi_est[jg];
       }
       
       // take logarithm again
       prob_predict((i_choice_1),ij) = std::log(ll_n_est);
       
     } //end y_c loop over alternatives
     
   } //end i_choice_1 loop over choice decisions.
   
   // 
   // //............................................................//
   // // undo logarithm and normalize
   // //............................................................//
   prob_predict = prob_predict.array().exp();
   Eigen::VectorXd prob_sum = prob_predict.array().rowwise().sum();
   
   for(int y_c = 0; y_c < alt; y_c++){
     prob_predict.col(y_c) =  prob_predict.col(y_c).array() / prob_sum.array();
   }
   //............................................................//
   // Format the output
   //............................................................//
   
   return  prob_predict;
   
 }
 
 
 //////////////////////////////////////////////////////////////////
 /////// Latent Class MACML approach                     //////////
 //////////////////////////////////////////////////////////////////
 
 
 // input:
 //  - theta           VectorXd        parameter vector (to optimize)
 //  - data_obj        data_obj        observed data 
 //  - mod             mod_cl object   describes the model 
 //  - control         control_cl object   contains info on various parameter, e.g. 0=full pairwise, 1=adjacent pairs looped, 2=adjacent pairs chained (not yet correctly implemented -> needs weights)
 
//' Calculation of Approximated Composite Marginal Likelihood for latent class models. No gradient supplied!
//'
//' @description
//' This function computes the composite marginal likelihood for latent class model.
//'
//' @param theta
//' parameter vector
//' @param data_obj
//' data_cl object
//' @param mod
//' A \code{\link{mod_cl}} object.
//' @param control
//' Controls for the probit estimation.
//'
//' @return
//' A vector, containing the negative log-likelihood, no gradient and Hessian.
//' @export
//'
// [[Rcpp::export]]
NumericVector ll_macml_LC(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   int cml_pair_type = 0;
   int hess = 0;
   int el = 0;
   Rcpp::List pairs_list;
   
   //Rcout << "start " <<  std::endl; 
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
   
   // initialize output 
   NumericVector out;
   
   //Rcout << cml_pair_type << std::endl; 
   
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
   
   //Rcout << "pairs_list" << std::endl; 
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
   int num_class = mod["num_class"];
   Eigen::MatrixXd Hb = mod["Hb"];
   Eigen::VectorXd fb = mod["fb"];
   Eigen::VectorXd fL = mod["fL"];
   Eigen::MatrixXd HL = mod["HL"];
   Eigen::VectorXd fO = mod["fO"];
   Eigen::MatrixXd HO = mod["HO"];
   
   // number of parameters for one component
   
   int param_one = lthb+ lthO + lthL;
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   double ll = 0;
   
   //............................................................//
   // RECOVER PARAMETERS
   //............................................................//
   
   // Recover the vector b, matrices Omega and Sigma 
   Eigen::MatrixXd b(Hb.rows(),num_class);
   b.setZero();
   
   Eigen::MatrixXd Sigma(alt*num_class,alt);
   Sigma.setZero();
   
   Eigen::MatrixXd Sigma_chol_vech;
   Eigen::MatrixXd Sigma_chol, Sigma_chol_t;
   Eigen::MatrixXd Omega(lRE*num_class,lRE);
   Omega.setZero();
   
   Eigen::MatrixXd ELt =  elimmat(alt).transpose();
   Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
   Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t;
   
   int lth_Omega_chol_vech = 0;
   int lth_Sigma_chol_vech = 0;
   
   for (int jb=0; jb<num_class; jb++){
     b.col(jb) =  ((Hb*theta.segment(param_one*jb,lthb)) + fb);
     
     // Recover the Sigma matrix
     // Without free parameters only fL is used
     if(lthL==0){
       Sigma_chol_vech = fL;
     }else{
       //Sigma_chol_vech = ((HL*l) + fL);
       Sigma_chol_vech = ((HL*theta.segment(param_one*jb+(lthO+lthb), lthL)) + fL);
     }
     
     lth_Sigma_chol_vech = Sigma_chol_vech.rows();
     Sigma_chol = ELt*Sigma_chol_vech;
     Sigma_chol.resize(alt,alt);
     Sigma_chol_t = Sigma_chol.transpose();
     
     Sigma.block(alt*jb,0,alt,alt) = Sigma_chol*Sigma_chol_t;
     
     // Recover the matrix Omega
     
     // check if we have mixed parameters
     if(lRE!=0){
       // Without free parameters only fO is used
       if(lthO==0){
         Omega_chol_vech = fO;
       }else{
         //Omega_chol_vech = ((HO*o) + fO);
         Omega_chol_vech = ((HO*theta.segment(param_one*jb+(lthb), lthO)) + fO);
       }
       
       //Rcout << Omega_chol_vech << std::endl; 
       
       
       lth_Omega_chol_vech = Omega_chol_vech.rows();
       Omega_chol = EOt*Omega_chol_vech;
       Omega_chol.resize(lRE,lRE);
       Omega_chol_t = Omega_chol.transpose();
       Omega.block(lRE*jb,0,lRE,lRE) = Omega_chol*Omega_chol_t;
     } else{
       Omega_chol_vech = fO;
     }
   }
   
   
   
   int m_l =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech); //compute the dimension of the matrix from the length of the vector
   int m_o =  -0.5 + sqrt(0.25 + 2*lth_Omega_chol_vech); //compute the dimension of the matrix from the length of the vector
   
   Eigen::MatrixXd L_JcholL = elimmat(m_l);
   Eigen::MatrixXd D_JcholL  = duplmat(m_l);
   Eigen::MatrixXd K_JcholL  = commmat(m_l,m_l);
   
   Eigen::MatrixXd L_JcholO = elimmat(m_o);
   Eigen::MatrixXd D_JcholO  = duplmat(m_o);
   Eigen::MatrixXd K_JcholO  = commmat(m_o,m_o);
   
   int lthbb = b.rows();
   int lthJ = lthbb + lth_Omega_chol_vech + lth_Sigma_chol_vech;
   Eigen::MatrixXd J_chol(num_class*lthJ,lthJ);
   for (int jb=0;jb<num_class;jb++){
     if(lthL==0){
       Sigma_chol_vech = fL;
     }else{
       //Sigma_chol_vech = ((HL*l) + fL);
       Sigma_chol_vech = ((HL*theta.segment(param_one*jb+(lthO+lthb), lthL)) + fL);
     }
     if(lRE!=0){
       // Without free parameters only fO is used
       if(lthO==0){
         Omega_chol_vech = fO;
       }else{
         //Omega_chol_vech = ((HO*o) + fO);
         Omega_chol_vech = ((HO*theta.segment(param_one*jb+(lthb), lthO)) + fO);
       }
     }
     J_chol.block(lthJ*jb,0,lthJ,lthJ)= J_bvOvSig_bLOLSig(Omega_chol_vech, Sigma_chol_vech, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL);
   }
   
   //Rcout << J_chol << std::endl; 
   
   
   // mixing coefficients
   Eigen::VectorXd pi_est(num_class),gamma_est(num_class);
   pi_est.setZero();
   gamma_est.setZero();
   
   gamma_est.segment(1,num_class-1) = theta.segment(param_one*num_class, num_class-1);
   
   double sum_pi = 0;
   
   for (int jg =0; jg < num_class; jg++){
     pi_est[jg]=  std::exp(gamma_est[jg]);
     sum_pi += pi_est[jg];
   }
   
   for (int jg =0; jg < num_class; jg++){
     pi_est[jg]=  pi_est[jg]/sum_pi;
   }
   
   //............................................................//
   // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT AND HESSIAN
   //............................................................//
   
   int alt2 = alt+alt-2;
   int ltheta = theta.size();
   
   Eigen::MatrixXd Hess_approx(ltheta,ltheta);
   Hess_approx.setZero();
   Eigen::MatrixXd Hess_approx2(ltheta,ltheta);
   Hess_approx2.setZero();
   Eigen::VectorXd grad(ltheta);
   grad.setZero();
   
   //............................................................//
   // END of setup
   //............................................................//
   
   
   // freeth denotes the number of free parameters for b, LO and LSig.
   int freeth = Hb.rows() +HL.rows() + HO.rows(); // bsum, Osum, Lsum;
   
   int cur_lp =1;
   double ll_n = 0;
   //Rcout << "N:" << N << std::endl; 
   
   for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
     // initialize log-likelihood contribution, gradient contribution and hessian contribution for the current decision maker
     ll_n = 0.0;
     
     Eigen::MatrixXd Hess_approx_n(ltheta,ltheta);
     Hess_approx_n.setZero();
     
     Eigen::VectorXd grad_n(ltheta);
     grad_n.setZero();
     
     // save the data of the current decision maker
     Rcpp::List data = data_obj["data"]; // extract the data from the data object
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
     
     // if there are not random effects: calculate probit probs as product of marginals 
     if (lRE == 0) { 
       pairs_check = 1;
       lth_pairs_n = 0;
     }
     //Rcout << pairs_check << std::endl; 
     
     //pairs_check = 0; 
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
     
     // normalize to obtain similar likelihood as in full probit likelihood calculations. //
     // reweight by Tp/sum of weights*2. 
     
     // sum of weights 
     //double sum_of_weights = 0;
     //
     //for (int js=0;js<pairs_n.rows();js++){
     // sum_of_weights += pairs_n(js,2);
     //}
     //
     //sum_of_weights = sum_of_weights * 2; 
     //if (sum_of_weights == 0){ 
     // sum_of_weights = 1.0;
     //}
     //
     //for (int js=0;js<pairs_n.rows();js++){
     //   pairs_n(js,2) = pairs_n(js,2)*Tp_n / sum_of_weights;
     //}
     
     //............................................................//
     // deal with observations, where only one choice is observed  //
     //............................................................//
     
     if (lth_pairs_n == 0){
       
       // either only one observation or no random effects. 
       //  extract the chosen alternatives for the current pair of choice occasions
       int y1;
       int alt2 = (alt-1);
       Eigen::MatrixXd X(alt2,lthbb);
       X.setZero();
       
       for (int jp =0; jp< Tp_n; jp++){
         y1 = y_n[jp];
         X = X_n[jp];
         
         double ll_n_pair = 0.0;
         Eigen::VectorXd grad_n_pair(ltheta); 
         grad_n_pair.setZero();
         
         int lth_output;
         for (int jg = 0; jg< num_class; jg++){ 
           Eigen::VectorXd grad_output(ltheta+1); 
           grad_output  = cal_choice_probs_1(X,alt,y1,b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE), J_chol.block(jg*lthJ,0,lthJ,lthJ), Hb, HO, HL,lth_Omega_chol_vech,lth_Sigma_chol_vech, approx_method);
           
           Eigen::VectorXd dpi(num_class-1);
           dpi.setZero();
           
           for (int jk = 1;jk<num_class;jk++){
             if (jk == jg){ 
               dpi[jk-1] = pi_est[jk] * (1-pi_est[jk]); 
             } else {
               dpi[jk-1] = - pi_est[jg] * pi_est[jk];
             }
           }
           lth_output = grad_output.size();
           double pi_k = std::exp(grad_output[lth_output-1]);
           double pi_pij = pi_k* pi_est[jg];
           ll_n_pair +=  pi_pij;
           
           // gradient   // 
           grad_n_pair.segment(param_one*jg,param_one) = grad_output.segment(0,param_one)*pi_pij;
           grad_n_pair.segment(param_one*num_class,num_class-1) -= pi_k * dpi; 
         }
         //............................................................//
         // END of gradient manipulations
         //............................................................//
         
         
         // regularize for too small values 
         if (ll_n_pair < 0.0000000001){
           ll_n_pair = 0.0000000001;
         }
         // take logarithm again 
         ll_n -= std::log(ll_n_pair);
         
         // grad
         grad_n_pair = grad_n_pair / ll_n_pair; 
         grad_n += grad_n_pair;
         
         // Hess_approx
         //............................................................//
         //Hessian
         //............................................................//
         Hess_approx_n += grad_n_pair * grad_n_pair.transpose(); 
       }
       
     } // end marginal log likelihood
     
     // NOW LOOP OVER THE PAIRS MATRIX to compute the CML
     for(int i_pair = 0; i_pair < lth_pairs_n; i_pair ++){
       
       // only calculate the pair if the weight is !=0
       if(pairs_n(i_pair,2)!=0){
         double ll_n_pair = 0.0;
         Eigen::VectorXd grad_n_pair(ltheta); 
         grad_n_pair.setZero();
         
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
         
         for (int jg = 0; jg< num_class; jg++){ // cycle over class components
           int lth_output = (param_one+1)*num_class;
           
           //Eigen::VectorXd grad_output = cal_choice_probs_2(X1, X2, y1, y2, b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE),approx_method);
           Eigen::VectorXd grad_output  = cal_choice_probs_2(X1,X2,y1,y2,b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE), J_chol.block(jg*lthJ,0,lthJ,lthJ), Hb, HO, HL,lth_Omega_chol_vech,lth_Sigma_chol_vech, approx_method);
           
           Eigen::VectorXd dpi(num_class-1);
           dpi.setZero();
           
           //Rcout << pi_est << std::endl; 
           
           for (int jk = 1;jk<num_class;jk++){
             if (jk == jg){ 
               dpi[jk-1] = pi_est[jk] * (1-pi_est[jk]); 
             } else {
               dpi[jk-1] = - pi_est[jg] * pi_est[jk];
             }
           }
           lth_output = grad_output.size();
           double pi_k = std::exp(grad_output[lth_output-1]);
           double pi_pij = pi_k* pi_est[jg];
           ll_n_pair +=  pi_pij;
           
           // gradient   // 
           grad_n_pair.segment(param_one*jg,param_one) = grad_output.segment(0,param_one)*pi_pij;
           grad_n_pair.segment(param_one*num_class,num_class-1) -= pi_k * dpi; 
           //............................................................//
           // END of gradient manipulations
           //............................................................//
           
         }
         
         // regularize for too small values 
         if (ll_n_pair < 0.0000000001){
           ll_n_pair = 0.0000000001;
         }
         // take logarithm again 
         ll_n -= std::log(ll_n_pair)*pairs_n(i_pair,2);
         
         // grad
         grad_n_pair = grad_n_pair / ll_n_pair; 
         grad_n += grad_n_pair* pairs_n(i_pair,2);
         
         // Hess_approx
         //............................................................//
         //Hessian
         //............................................................//
         Hess_approx_n += grad_n_pair * grad_n_pair.transpose()* pairs_n(i_pair,2);         
       } // end if-weight!=0
     } //end pairs-loop
     
     
     ll += ll_n; //it is already negative in the pairs loop
     grad += grad_n; 
     Hess_approx2 += grad_n * grad_n.transpose();
     Hess_approx += Hess_approx_n;
   } //end N-loop
   
   
   //............................................................//
   // Format the output
   //............................................................//
   out = ll;
   out.attr("gradient") = (grad.array());
   out.attr("hessian2") =(Hess_approx2.array());
   out.attr("hessian1") =(Hess_approx.array());
   
   return out;
   
 }
 
 
 //////////////////////////////////////////////////////////////////
 /////// Non-parametric MACML approach                   //////////
 //////////////////////////////////////////////////////////////////
 
 
 // input:
 //  - theta           VectorXd        parameter vector (to optimize)
 //  - data_obj        data_obj        observed data 
 //  - mod             mod_cl object   describes the model 
 //  - control         control_cl object   contains info on various parameter, e.g. 0=full pairwise, 1=adjacent pairs looped, 2=adjacent pairs chained (not yet correctly implemented -> needs weights)
//' Calculation of choice probabilities for non-parametric models.
//'
//' @description
//' This function computes the choice probabilities for non-parametric models. Each component is not mixed and hence choice probabilities 
//' relate to each single choice, no parameters are needed, as the grid points are specified in the model object.  
//'
//' @param data_obj
//' data_cl object
//' @param mod
//' A mod_nonpara_cl object.
//' @param control
//' Controls for the probit estimation.
//' @param cml_pair_type
//' integer controls the probabilities calculated: -1: single choices, 0: all pairs, 1: adjacent pairs. 
//'
//' @return
//' A matrix, containing the choice probabilities for each choice in the rows; for each grid point in the columns. 
//' @export
//'
// [[Rcpp::export]]
Eigen::MatrixXd choice_probs_nonpara(Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control, int cml_pair_type)
{
   
   //............................................................//
   // UNPACK ESTIMATION CONTROLS
   //............................................................//
   
   std::string approx_method = "SJ";
   if (cml_pair_type>1){ 
     cml_pair_type = 1;
   }
   if (cml_pair_type<-1){ 
     cml_pair_type = -1;
   }
   
   //Rcout << "start " <<  std::endl; 
   if(control.containsElementNamed("approx_method")){
     CharacterVector approx_method_temp = control["approx_method"];
     approx_method = Rcpp::as<std::string>(approx_method_temp);
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
   int num_grid_points = mod["num_grid_points"];
   Eigen::MatrixXd Hb = mod["Hb"];
   Eigen::VectorXd fb = mod["fb"];
   Eigen::VectorXd fL = mod["fL"];
   Eigen::MatrixXd HL = mod["HL"];
   Eigen::VectorXd fO = mod["fO"];
   Eigen::MatrixXd HO = mod["HO"];
   
   // number of parameters for one component
   
   Eigen::MatrixXd params = mod["params"]; 
   
   //............................................................//
   // INITIALIZE OUTPUT
   //............................................................//
   
   // calculate the total number of choices 
   int number_choices = 0; 
   Rcpp::List data = data_obj["data"]; // extract the data from the data object
   
   for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
     Rcpp::List data_n = data[i_n]; //unpack the data
     Eigen::VectorXd y_n = data_n["y"];     
     int Tp_n = y_n.rows();
     if (cml_pair_type == -1){
       number_choices += Tp_n;
     }
     if (cml_pair_type == 0){
       number_choices += Tp_n*(Tp_n-1)*0.5;
     }
     if (cml_pair_type == 1){
       number_choices += Tp_n-1;
     }
   }
   
   int colout = params.cols()+1; 
   Eigen::MatrixXd out(number_choices,colout);
   out.setZero();
   
   //............................................................//
   // RECOVER PARAMETERS
   //............................................................//
   
   // Recover the vector b, matrices Omega and Sigma 
   Eigen::MatrixXd b(Hb.rows(),num_grid_points);
   b.setZero();
   
   int lthbb = Hb.rows(); 
   
   Eigen::MatrixXd Sigma(alt*num_grid_points,alt);
   Sigma.setZero();
   
   Eigen::MatrixXd Sigma_chol_vech;
   Eigen::MatrixXd Sigma_chol, Sigma_chol_t;
   Eigen::MatrixXd Omega(lRE*num_grid_points,lRE);
   Omega.setZero();
   
   Eigen::MatrixXd ELt =  elimmat(alt).transpose();
   Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
   Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t;
   
   int lth_Omega_chol_vech = 0;
   int lth_Sigma_chol_vech = 0;
   
   for (int jb=0; jb<num_grid_points; jb++){
     b.col(jb) =  ((Hb*params.block(0,jb,lthb,1)) + fb);
     
     // Recover the Sigma matrix
     // Without free parameters only fL is used
     if(lthL==0){
       Sigma_chol_vech = fL;
     }else{
       //Sigma_chol_vech = ((HL*l) + fL);
       Sigma_chol_vech = ((HL*params.block((lthO+lthb),jb,lthL,1)) + fL);
     }
     
     lth_Sigma_chol_vech = Sigma_chol_vech.rows();
     Sigma_chol = ELt*Sigma_chol_vech;
     Sigma_chol.resize(alt,alt);
     Sigma_chol_t = Sigma_chol.transpose();
     
     Sigma.block(alt*jb,0,alt,alt) = Sigma_chol*Sigma_chol_t;
     
     // Recover the matrix Omega
     
     // check if we have mixed parameters
     if(lRE!=0){
       // Without free parameters only fO is used
       if(lthO==0){
         Omega_chol_vech = fO;
       }else{
         //Omega_chol_vech = ((HO*o) + fO);
         Omega_chol_vech = ((HO*params.block((lthb),jb, lthO,1)) + fO);
       }
       
       //Rcout << Omega_chol_vech << std::endl; 
       
       
       lth_Omega_chol_vech = Omega_chol_vech.rows();
       Omega_chol = EOt*Omega_chol_vech;
       Omega_chol.resize(lRE,lRE);
       Omega_chol_t = Omega_chol.transpose();
       Omega.block(lRE*jb,0,lRE,lRE) = Omega_chol*Omega_chol_t;
     } else{
       Omega_chol_vech = fO;
     }
   }
   //int m_l =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech); //compute the dimension of the matrix from the length of the vector
   //int m_o =  -0.5 + sqrt(0.25 + 2*lth_Omega_chol_vech); //compute the dimension of the matrix from the length of the vector
   
   int alt2 = (alt-1);
   for (int jg=0; jg<num_grid_points; jg++){
     int cur_lp =0;
     for(int i_n = 0; i_n < N; i_n++){ //loop over the sample
       
       // save the data of the current decision maker
       Rcpp::List data = data_obj["data"]; // extract the data from the data object
       Rcpp::List data_n = data[i_n]; //unpack the data
       Rcpp::List X_n = data_n["X"];
       Eigen::VectorXd y_n = data_n["y"];
       
       // determine the panel length
       int Tp_n = y_n.rows();
       // split according to cml_pair_type 
       if ((cml_pair_type == -1)||(Tp_n == 1)){
         for (int jp =0; jp< Tp_n; jp++){
           
           int y1 = y_n[jp];
           Eigen::MatrixXd X(alt2,lthbb);
           X.setZero();
           X = X_n[jp];
           
           out(cur_lp,jg)  = std::exp(cal_choice_probs_1_nograd(X,alt,y1,b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE), approx_method));
           out(cur_lp,colout-1) = 1.0;
           cur_lp += 1;
         } // end probabilities calculation
       }
       if ((cml_pair_type == 0)&&(Tp_n > 1)){ // all pairs 
         for (int jpa =1; jpa< Tp_n; jpa++){
           int y1 = y_n[jpa];
           
           Eigen::MatrixXd X1(alt2,lthbb);
           X1.setZero();
           X1 = X_n[jpa];
           for (int jpb=0;jpb<jpa;jpb++){
             
             int y2 = y_n[jpb];
             
             Eigen::MatrixXd X2(alt2,lthbb);
             X2.setZero();
             X2 = X_n[jpb];
             
             out(cur_lp,jg)  = std::exp(cal_choice_probs_2_nograd(X1, X2, y1, y2, b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE), approx_method));
             out(cur_lp,colout-1) = 1.0/((double)Tp_n-1.0);
             
             //Rcout << cur_lp << colout << "out" << out(cur_lp,colout-1) << std::endl; 
             cur_lp += 1;
           } // end jpb loop 
         } // end probabilities calculation
       }
       if ((cml_pair_type == 1)&&(Tp_n > 1)){  // adjacent pairwise
         for (int jpa =1; jpa< Tp_n; jpa++){
           int y1 = y_n[jpa-1];
           
           Eigen::MatrixXd X1(alt2,lthbb);
           X1.setZero();
           X1 = X_n[jpa-1];
           
           int y2 = y_n[jpa];
           
           Eigen::MatrixXd X2(alt2,lthbb);
           X2.setZero();
           X2 = X_n[jpa];
           
           out(cur_lp,jg)  = std::exp(cal_choice_probs_2_nograd(X1, X2, y1, y2, b.col(jg), Sigma.block(jg*alt,0,alt,alt), Omega.block(jg*lRE,0,lRE,lRE), approx_method));
           out(cur_lp,colout-1) = (double)Tp_n/(2.0*((double)Tp_n-1.0));
           cur_lp += 1;
         } // end probabilities calculation
       }
       
     } //end N-loop
   } // end loop over grid points
   return out;
   
 }
 
 
