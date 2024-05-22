#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <omp.h>

#include "SJ_vdb.h"
#include "TVBS.h"
#include "ME.h"
#include "macml_grad_Hessian.h"
#include "lin_alg.h"
#include "toms462.h"                                         // allows bivariate normal cdf calculations
#include "cml_ordered_StateSpace_helper.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
//using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
//using Eigen::EigenSolver;    // one of the eigenvalue solvers


// // // // // // // // // // // // // // // // //
//                                              //
//    Ordered PROBIT CML WITH APPROXIMATED MVNCDF        //
//    includes state space model for error component.   //
// // // // // // // // // // // // // // // // //



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

//////////////////////////////////////////////////////
/// procedures to build the system             ///////
//////////////////////////////////////////////////////

///' build_system_from_model
///' @description
///' Calculates the model based on parameter vector theta and model structure mod
///' @param theta
///' parameter vector
///' @param mod
///' A \code{\link{mod_StSp_cl}} object.
///' @param time 
///' vector of observation times 
///' @return 
///' a structure containing beta, Omega, Sigma.
///' @keywords internal  
/// [[Rcpp::export]]

struct system build_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time)
{
  struct system syst; 
  
  int n = mod["dim_state"];
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
  
  bool stationary = mod["stationary"];
  
  //............................................................//
  // RECOVER PARAMETERS
  //............................................................//
  
  
  // Recover the vector b
  Eigen::VectorXd beta = ((Hb*theta.head(lthb)) + fb);
  
  syst.beta = beta; 
  
  // Recover the Sigma matrix
  // Without free parameters only fL is used
  Eigen::MatrixXd Sigma_chol_vech;
  if(lthL==0){
    Sigma_chol_vech = fL;
  }else{
    Sigma_chol_vech = ((HL*theta.segment((lthO+lthb), lthL)) + fL);
  }
  
  syst.Sigma_chol_vech = Sigma_chol_vech;
  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();
  
  int s =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech);
  Eigen::MatrixXd ELt =  elimmat(s).transpose();
  
  
  Eigen::MatrixXd Sigma_chol = ELt*Sigma_chol_vech;
  Sigma_chol.resize(s,s);
  Eigen::MatrixXd Sigma_chol_t = Sigma_chol.transpose();
  Eigen::MatrixXd Sigma = Sigma_chol*Sigma_chol_t;
  syst.Sigma = Sigma; 
  
  // Recover the matrix Omega
  Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
  Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t, Omega;
  
  
  // check if we have mixed parameters
  if(lRE!=0){
    // Without free parameters only fO is used
    if(lthO==0){
      Omega_chol_vech = fO;
    }else{
      //Omega_chol_vech = ((HO*o) + fO);
      Omega_chol_vech = ((HO*theta.segment((lthb), lthO)) + fO);
    }
    
    Omega_chol = EOt*Omega_chol_vech;
    Omega_chol.resize(lRE,lRE);
    Omega_chol_t = Omega_chol.transpose();
    Omega = Omega_chol*Omega_chol_t;
  }else{
    Omega_chol_vech = fO;
    Omega_chol = fO;
    Omega_chol_t = fO;
    Omega = fO;
  }
  
  syst.Omega = Omega; 
  syst.Omega_chol_vech = Omega_chol_vech;
  // get the cutoff points
  
  Eigen::VectorXd dtauk = theta.segment(lthb+lthO +lthL,alt-1);
  Eigen::VectorXd tauk(alt-1);
  tauk.setZero();
  
  tauk(0) = dtauk(0);
  for (int j=1;j<alt-1;j++){
    tauk(j)=tauk(j-1)+std::exp(dtauk(j));
  }
  
  syst.tauk = tauk; 
  
  // calculate the state space system 
  Eigen::VectorXd param(2*s*n);
  param = theta.tail(2*s*n);
  
  
  
  // calculate system 
  struct state_space_system StSp_syst;
  int grad_bool = 0;
  StSp_syst = param_to_system(param,s,n, grad_bool, stationary);
  
  syst.factor = StSp_syst.factor; 
  
  // calculate starting covariance matrix as stationary value. 
  Eigen::MatrixXd P0(n,n); 
  P0.setZero();
  
  if (stationary){
    Eigen::VectorXd vQ = vectorize(StSp_syst.K * Sigma *StSp_syst.K.transpose()); 
    P0 = solve_Lyapunov_equation(StSp_syst.A, vQ);
  } 
  // calculate variance matrix 
  int T = time.size();
  Eigen::MatrixXd GammaT(s*T,s*T);
  
  //Rcout << "T" << T  << std::endl;
  
  //Rcout << "A:" << StSp_syst.A << "K:" << StSp_syst.K << "C:" << StSp_syst.C << std::endl; 
  //Rcout << "P" << P0 << "Sig" << Sigma << "time" << time << std::endl; 
  
  GammaT = calculate_Cov_Seq(StSp_syst, P0, Sigma, time);
  syst.GammaT = GammaT; 
  
  // return system; 
  return syst;
}

//' build_system_from_model_R
//' @description
//' put evaluations back to R 
//' @param theta
//' parameter vector
//' @param mod
//' A \code{\link{mod_StSp_cl}} object.
//' @param time 
//' vector of observation times 
//' @return
//' a structure containing beta, Omega, Sigma.
//' @export
//'
// [[Rcpp::export]]
Rcpp::List build_system_from_model_R(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time)
{
   struct system syst = build_system_from_model(theta, mod, time);
   Rcpp::List out; 
   out["beta"] = syst.beta; 
   out["Omega"] = syst.Omega; 
   out["Sigma"] = syst.Sigma; 
   out["GammaT"] = syst.GammaT; 
   out["tauk"] = syst.tauk; 
   out["factor"] = syst.factor; 
   return out; 
}


//' build_derived_system_from_model
//' @description
//' Calculates the model based on parameter vector theta and model structure mod
//' @param theta
//' parameter vector
//' @param mod
//' A \code{\link{mod_StSp_cl}} object.
//' @param time 
//' vector of observation times 
//' @param coord 
//' integer, indicating with respect to which coordinate the derivative is taken
//' @return
//' matrix dGammaT of derivatives of covariances 
//' @export
//'
// [[Rcpp::export]]
Eigen::MatrixXd build_derived_system_from_model(Eigen::VectorXd theta, Rcpp::List mod, Eigen::VectorXd time, int coord)
{

   struct system syst;
   
   int n = mod["dim_state"];
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
   
   bool stationary = mod["stationary"];
   
   //............................................................//
   // RECOVER PARAMETERS
   //............................................................//
   
   
   // Recover the vector b
   Eigen::VectorXd beta = ((Hb*theta.head(lthb)) + fb);
   
   syst.beta = beta; 
   
   // Recover the Sigma matrix
   // Without free parameters only fL is used
   Eigen::MatrixXd Sigma_chol_vech;
   if(lthL==0){
     Sigma_chol_vech = fL;
   }else{
     Sigma_chol_vech = ((HL*theta.segment((lthO+lthb), lthL)) + fL);
   }
   
   int lth_Sigma_chol_vech = Sigma_chol_vech.rows();
   
   int s =  -0.5 + sqrt(0.25 + 2*lth_Sigma_chol_vech);
   Eigen::MatrixXd ELt =  elimmat(s).transpose();
   
   
   Eigen::MatrixXd Sigma_chol = ELt*Sigma_chol_vech;
   Sigma_chol.resize(s,s);
   Eigen::MatrixXd Sigma_chol_t = Sigma_chol.transpose();
   Eigen::MatrixXd Sigma = Sigma_chol*Sigma_chol_t;
   syst.Sigma = Sigma; 
   syst.Sigma_chol_vech = Sigma_chol_vech; 
   
   // Recover the matrix Omega
   Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
   Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t, Omega;
   
   
   // check if we have mixed parameters
   if(lRE!=0){
     // Without free parameters only fO is used
     if(lthO==0){
       Omega_chol_vech = fO;
     }else{
       //Omega_chol_vech = ((HO*o) + fO);
       Omega_chol_vech = ((HO*theta.segment((lthb), lthO)) + fO);
     }
     
     Omega_chol = EOt*Omega_chol_vech;
     Omega_chol.resize(lRE,lRE);
     Omega_chol_t = Omega_chol.transpose();
     Omega = Omega_chol*Omega_chol_t;
   }else{
     
     Omega_chol_vech = fO;
     Omega_chol = fO;
     Omega_chol_t = fO;
     Omega = fO;
   }
   
   syst.Omega = Omega; 
   syst.Omega_chol_vech = Omega_chol_vech; 
   // get the cutoff points
   
   Eigen::VectorXd dtauk = theta.segment(lthb+lthO +lthL,alt-1);
   Eigen::VectorXd tauk(alt-1);
   tauk.setZero();
   
   tauk(0) = dtauk(0);
   for (int j=1;j<alt-1;j++){
     tauk(j)=tauk(j-1)+std::exp(dtauk(j));
   }
   
   syst.tauk = tauk; 
   
   // calculate the state space system 
   Eigen::VectorXd param(2*s*n);
   param = theta.tail(2*s*n);
   
   // calculate system 
   struct state_space_system StSp_syst;
   int grad_bool = 0;
   StSp_syst = param_to_system(param,s,n, grad_bool,stationary);
   
   // calculate starting covariance matrix as stationary value. 
   Eigen::MatrixXd P0(n,n); 

   P0.setZero(); 


   if (stationary){
     Eigen::VectorXd vQ = vectorize(StSp_syst.K * Sigma *StSp_syst.K.transpose()); 
     Eigen::MatrixXd P0 = solve_Lyapunov_equation(StSp_syst.A, vQ);
   }
   
   // calculate variance matrix 
   int T = time.size();
   Eigen::MatrixXd GammaT(s*T,s*T);
   
   // //Rcout << "T" << T  << std::endl;
   

   GammaT = calculate_Cov_Seq(StSp_syst, P0, Sigma, time);
   syst.GammaT = GammaT; 
   
   // now turn to derivatives 
   Eigen::MatrixXd dGammaT = GammaT;
   dGammaT.setZero(); 
   
   struct state_space_system dStSp_syst;
   dStSp_syst = param_to_system_grad(2*s*n, s, n); // initialize with zeros. 
   
   // initialization using stationary state variance 
   // //P0(0,0)=-1; 
   Eigen::MatrixXd dP0 = 0*P0;
   dP0.setZero(); 
   
   // derivative of Sigma: initially set to zero.
   Eigen::MatrixXd dSigma = 0*Sigma;
   dSigma.setZero();
   
   // now determine derivative based on coordinate number:
   if (coord>lthb+lthO-1){ // not w.r.t. parameters for beta or Omega. 
     if (coord<lthb+lthO+lthL){
       int coord_sig = coord - lthb-lthO; 
       //Rcout << coord_sig << std::endl; 
       // parameter corresponds to Sigma = L_s L_s'.  
       Eigen::MatrixXd L_JcholL = elimmat(s);
       Eigen::MatrixXd D_JcholL  = duplmat(s);
       Eigen::MatrixXd K_JcholL  = commmat(s,s);
       
       Eigen::VectorXd HlL = Sigma_chol_vech.array();
       
       //Rcout << "s:" << s << "lth" << lth_Sigma_chol_vech << std::endl; 
       
       Eigen::MatrixXd JJ1(lth_Sigma_chol_vech,lth_Sigma_chol_vech);
       JJ1.setZero();
       JJ1  = J_chol1(HlL, L_JcholL, K_JcholL, D_JcholL);
       Eigen::VectorXd JJ2(lth_Sigma_chol_vech);
       JJ2.setZero();
       JJ2 = JJ1 * HL.col(coord_sig);
       Eigen::MatrixXd Sigma_chol = ELt*JJ2;
       Sigma_chol.resize(s,s);
       dSigma = Sigma_chol; 
       dStSp_syst = param_to_system_grad(2*s*n, s, n); // set to zero. 
       
       if (stationary){
         Eigen::MatrixXd dQ = 0*P0;
         Eigen::MatrixXd A = StSp_syst.A; 
         Eigen::MatrixXd K = dStSp_syst.K; 
         
         dQ = dQ + K * dSigma * K.transpose();
         Eigen::VectorXd dvQ = vectorize(dQ); 
         dP0 = solve_Lyapunov_equation(StSp_syst.A, dvQ);
       }
       //    //Rcout << P0 << "dP" << dP0 << "Sig" << Sigma << "dSig" << dSigma << std::endl; 
       //    
       //    //Rcout << "A:" << StSp_syst.A << "K:" << StSp_syst.K << "C:" << StSp_syst.C << std::endl; 
       //    //Rcout << "dA:" << dStSp_syst.A << "dK:" << dStSp_syst.K << "dC:" << dStSp_syst.C << std::endl; 
       //    
       dGammaT = calculate_dCov_Seq(StSp_syst, dStSp_syst,P0, dP0, Sigma, dSigma, time);
       //      
     } else {
       // parameter corresponds to either tauk or StSp system. 
       if (coord>lthb+lthO +lthL+alt-2){
         // not tauk -> state space system. 
         int coord_StSp = coord - lthb-lthO -lthL-alt+1; // coordinate number only counting the ones for the state space system. 
         dStSp_syst = param_to_system_grad(coord_StSp, s, n);
         dSigma.setZero(); 
         
         if (stationary){ 
           Eigen::MatrixXd dQ = 0*P0;
           Eigen::MatrixXd dA = dStSp_syst.A; 
           Eigen::MatrixXd A = StSp_syst.A; 
           Eigen::MatrixXd dK = dStSp_syst.K; 
           Eigen::MatrixXd K = dStSp_syst.K; 
           
           dQ = dQ + dA * P0 * A.transpose() + A * P0 * dA.transpose();
           dQ = dQ + dK * Sigma * K.transpose() + K * dSigma * K.transpose() + K * Sigma * dK.transpose();
           Eigen::VectorXd dvQ = vectorize(dQ); 
           dP0 = solve_Lyapunov_equation(StSp_syst.A, dvQ);
         } 
         dGammaT = calculate_dCov_Seq(StSp_syst, dStSp_syst,P0, dP0, Sigma, dSigma, time);
       }
     }
   }
   
   return dGammaT;

}







//////////////////////////////////////////////////////////////////////////////////////
////// calculate the model based on paramter vector theta and model structure mod  ///
//////////////////////////////////////////////////////////////////////////////////////



///' find_ind 
///' @description
///' finds an entry in a vector and returns its index.  
///' @param x
///' double 
///' @param vec
///' vector of doubles. 
///' @return
///' integer; index where vec(ind) = x.  
///' @keywords internal 

int find_ind(double x,Eigen::VectorXd vec)
{
  int ind =0;
  
  for (int jj=0;jj<vec.size();jj++){
    if (std::abs(vec(jj)-x)<0.00001){ 
      ind = jj; 
    }
  }
  return ind; 
}

//' get_observation_indices
//' @description
//' Compares the time and question numbers from an observation with the maximal ones and returns a vector containing the indices of observations contained for the 
//' individual. 
//' @param time 
//' vector of all observation times 
//' @param quest
//' integer of most questions asked 
//' @param time_n
//' vector of times when the n-th individual is observed. 
//' @param quest_n
//' vector of questions answered by the n-th individual. 
//' 
//' @return
//' a vector including all the observed questions. 
//' 
//' @export
//'
// [[Rcpp::export]] 
Eigen::VectorXd get_observation_indices(Eigen::VectorXd time,int quest,Eigen::VectorXd time_n, Eigen::VectorXd quest_n)
{
  int nn = time_n.size();
  int n = time.size() * quest; 
  
  Eigen::VectorXd ind_obs(nn); 
  ind_obs.setZero();
  
  Eigen::VectorXd time_q(n); 
  
  // generate vector of obs 
  for (int jj=0;jj<time.size();jj++){
    for (int jq=0;jq<quest;jq++){
      time_q(jj*quest+jq) = time(jj)*quest + (jq+1); 
    }
  }
  // generate vector of obs_n
  
  //Rcout << time_q << std::endl;
  for (int jj=0;jj<nn;jj++){
    // cycle over obs_n
    //Rcout << time_n(jj)*quest + quest_n(jj) << std::endl; 
    ind_obs(jj) = find_ind(time_n(jj)*quest + quest_n(jj),time_q); 
  }

  return ind_obs; 
}

//' take out rows and columns 
//' @description
//' Extracts the submatrix Gamma(vec,vec).  
//' @param Gam  
//' matrix.  
//' @param vec 
//' vector indices.
//'  
//' @return
//' a submatrix. 
//' @keywords internal 
// [[Rcpp::export]]
Eigen::MatrixXd subset(Eigen::MatrixXd Gam, Eigen::VectorXd vec)
{
  int n = vec.size(); 
  
  Eigen::MatrixXd Gam_ind(n,n);
  Gam_ind.setZero(); 
  
  for (int ir=0;ir<n;ir++){
    for (int ic=0;ic<n;ic++){
      Gam_ind(ir,ic) = Gam( (int)vec[ir], (int)vec[ic]);  
    }
  }

  return Gam_ind;
};


//' find indices for elimination of unwanted rows and columns
//' @description
//' Extracts the indices of the entries in the submatrix LSig (ind,ind).  
//' @param n  
//' integer; number of rows of large matrix   
//' @param ind 
//' vector indices.
//'  
//' @return
//' an indicator matrix selecting the smaller matrix elements.
//'  
//' @export
//'
// [[Rcpp::export]] 
Eigen::MatrixXd elim_ind(int n, Eigen::VectorXd ind)
{
    int m = ind.size(); 
    // write full matrix with all indices 
    Eigen::MatrixXd full_ind(n,n);
    full_ind.setZero();
    int cur =0;
    for (int ic=0;ic<n;ic++){
      for (int ir=ic;ir<n;ir++){
        full_ind(ir,ic) = cur;
        cur = cur+1; 
      }
    }
  
    //Rcout << full_ind << std::endl; 
    // reduce by subsetting only those rows and columns specified in ind. 
    Eigen:MatrixXd red_ind = subset(full_ind,ind);
  
    //Rcout << red_ind << std::endl;
  
    // define selection matrix. 
    int m_l = m*(m+1)/2;
    int n_l = n*(n+1)/2; 
  
  //Rcout << "m_l" << m_l << "n_l:" << n_l << std::endl; 
                        
    Eigen::MatrixXd Elim(m_l,n_l);
    Elim.setZero(); 
    cur = 0;
    for (int ic=0;ic<m;ic++){
      for (int ir=ic;ir<m;ir++){
        Elim(cur,(int)red_ind(ir,ic)) = 1;
        cur = cur+1; 
      }
    }
    
    // return result.
    return Elim;
 };


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
NumericVector ll_macml_o_StSp(Eigen::VectorXd theta, Rcpp::List data_obj, Rcpp::List mod, Rcpp::List control){
  //............................................................//
  // UNPACK ESTIMATION CONTROLS
  //............................................................//

  // fix choices  in this function 
  std::string approx_method = "SJ";
  int cml_pair_type = 1;
  int hess = 0;
  int el = 0;
  // 
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
   
  Eigen::VectorXd time = data_obj["time"];
  int quest = data_obj["quest"];
  struct system syst = build_system_from_model(theta, mod, time);
   
  //Rcout << "time:" << time << std::endl; 
  
  // Recover the vector b
  Eigen::VectorXd b = syst.beta;
  int lthbb = b.rows();
   
  //Rcout << "b" << b << std::endl; 
   
  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma = syst.Sigma;
  Eigen::MatrixXd Sigma_chol_vech = syst.Sigma_chol_vech;
   
  // number of questions per person 
  int s= Sigma.rows(); 
   
  //Rcout << "s:" << s << std::endl; 
   
  // state dimension
  int n = mod["dim_state"];
   
  // Recover GammaT matrix:
  Eigen::MatrixXd GammaT = syst.GammaT; 
   
  // dimension of this matrix 
  int m_l = GammaT.rows(); 
  int lth_Sigma_chol_vech = m_l* (m_l+1)/2; 
   
  //Rcout << "m_l" << m_l << std::endl;
   
  // Recover the matrix Omega
  Eigen::MatrixXd Omega = syst.Omega; 
  
  // get the cutoff points
  Eigen::VectorXd tauk = syst.tauk; 
   
  //Rcout << "end of syst preps" << std::endl; 
   
  // extract system components from mod object. 
   
  int N= data_obj["N"];
  int lthb = mod["lthb"];
  int lthO = mod["lthO"];
  int lthL = mod["lthL"];
  int lRE = mod["lRE"];
  int alt = mod["alt"];
  Eigen::VectorXd Tp = data_obj["Tp"];
  Eigen::MatrixXd Hb = mod["Hb"];
  Eigen::VectorXd fb = mod["fb"];
  Eigen::VectorXd fL = mod["fL"];
  Eigen::MatrixXd HL = mod["HL"];
  Eigen::VectorXd fO = mod["fO"];
  Eigen::MatrixXd HO = mod["HO"];
   
  bool stationary = mod["stationary"];
  
  // get the differences of the cutoff points
  Eigen::VectorXd dtauk(tauk.size());
  dtauk.setZero(); 
  dtauk(0)= tauk(0); 
  for(int i_t = 1; i_t < dtauk.size(); i_t++){
     dtauk(i_t) = std::log(tauk(i_t)-tauk(i_t-1));  
  }
   
  // freeth denotes the number of free parameters for b, LO and LSig, tauk and GammaT. 
  int freeth = theta.size(); // 
   
  //Rcout << "after freeth" << tauk << std::endl; 
   
  //............................................................//
  // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT
  //............................................................//
   
  // deriv of xb,vech(Omega),vech(Sigma) with respect to theta.
  Eigen::MatrixXd L_JcholL = elimmat(m_l);
  Eigen::MatrixXd D_JcholL  = duplmat(m_l);
  Eigen::MatrixXd K_JcholL  = commmat(m_l,m_l);
   
  Eigen::VectorXd HlL = L_JcholL * vectorize(GammaT); 
  Eigen::VectorXd HlL_n = HlL;
   
   // Recover the matrix Omega
  Eigen::MatrixXd EOt =  elimmat(lRE).transpose();
  Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t;
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
  }else{
     lth_Omega_chol_vech = 0;
     Omega_chol_vech = fO;
     Omega_chol = fO;
     Omega_chol_t = fO;
  }
   
   
  int m_o =  -0.5 + sqrt(0.25 + 2*lth_Omega_chol_vech); //compute the dimension of the matrix from the length of the vector
  Eigen::MatrixXd L_JcholO = elimmat(m_o);
  Eigen::MatrixXd D_JcholO  = duplmat(m_o);
  Eigen::MatrixXd K_JcholO  = commmat(m_o,m_o);
  //This is used later in the computation to compute the Hessian of the Cholesky
  Eigen::VectorXd HoO = Omega_chol_vech.array();
  Eigen::VectorXd HoO_n = Omega_chol_vech.array();
   
   
  // //Rcout << "before JJ1" << std::endl;
   
   
  //
  // set up the Jacobian matrices to relate (b,Lambda) to theta.
  //
  // This is new for state space systems.
  // // JJ1 converts matrices to Cholesky factors for Omega and Sigma.
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
   
  // change JJ1 for state space system!
  JJ1.block(lthbb+lth_Omega_chol_vech,lthbb+lth_Omega_chol_vech,lth_Sigma_chol_vech,lth_Sigma_chol_vech) = identity(lth_Sigma_chol_vech);
   
  // JJ2 converts Cholesky factors to parameters. 
  int dim_param = lthb+lthO+lthL+2*s*n;
   
  // //Rcout << dim_param << std::endl; 
   
  Eigen::MatrixXd JJ2(lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech,dim_param);
  JJ2.setZero();
  JJ2.block(0,0,lthbb,lthb) = Hb;
  if(lRE!=0){
    JJ2.block(lthbb,lthb,lth_Omega_chol_vech,lthO) = HO;
   }
   
    
  // // derivative of Lambda with respect to parameters for Sigma.
  int coord;
  Eigen::MatrixXd dGam(m_l,m_l);
  dGam.setZero();
  Eigen::VectorXd vdg(m_l*m_l);
  vdg.setZero();
  
  for (int js=0;js<lthL;js++){
    coord = lthb+lthO+js;
    dGam = build_derived_system_from_model(theta, mod, time, coord);
    vdg = vectorize(dGam);
    JJ2.block(lthbb+lth_Omega_chol_vech,coord,lth_Sigma_chol_vech,1) = L_JcholL * vdg;
  }
   
  // derivative of Lambda with respect to parameters in the transfer function
  for (int js=0;js<2*s*n;js++){
     coord = lthb+lthO+lthL+alt-1+js;
     dGam = build_derived_system_from_model(theta, mod, time, coord);
     vdg = vectorize(dGam);
     JJ2.block(lthbb+lth_Omega_chol_vech,coord-alt+1,lth_Sigma_chol_vech,1) = L_JcholL * vdg;
  }
   
   
  Eigen::MatrixXd IIl(lth_Sigma_chol_vech,lth_Sigma_chol_vech);
  IIl.setIdentity();
   
  //............................................................//
  // INITIALIZE GRADIENT OUTPUT
  //............................................................//
   
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
   
   
  // //Rcout << "JJ1" << JJ1 << "JJ2" << JJ2 << std::endl; 
   
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
     //m_l = X_n.rows();
     Eigen::VectorXd y_n = data_n["y"];
   
     Eigen::VectorXd time_n = data_n["time"];
     Eigen::VectorXd quest_n = data_n["quest"]; 
     
  //   //Rcout << "i_n" << i_n << "t" << time << "q" << quest << "tn" << time_n << "qn" << quest_n << std::endl; 
     Eigen::VectorXd ind_obs = get_observation_indices(time,quest,time_n, quest_n);
     
  //   //Rcout << "i_n" << i_n << "IO" << ind_obs << std::endl; 
  //   
  //   // compute X*b.
     Eigen::VectorXd xb;
   
     xb = X_n * b;
     int m_ln = xb.size();
     int lth_Sigma_chol_vechn = m_ln* (m_ln+1)/2;
     // adjust matrices 
     Eigen::MatrixXd L_JcholLn = elimmat(m_ln);
     
     // compute Lambda.
     Eigen::MatrixXd Lambda(m_ln,m_ln);
     Lambda.setZero();
     
     // extract GammaT parts 
     Lambda = subset(GammaT,ind_obs);
     //Lambda = GammaT; 
   
     Eigen::MatrixXd xRE(m_ln,lRE);
   
     if(lRE!=0){
       xRE = X_n.block(0,0,m_ln,lRE);
       Lambda+= xRE * Omega * xRE.transpose();
     }

  // //Rcout << Lambda << "xb" << xb << "y_n" << y_n << "dtauk" << dtauk << "cml" << cml_pair_type << std::endl; 
     
     ll_n = prob_ordered_CML(xb, y_n, Lambda, alt, dtauk, cml_pair_type);
     // Rcout << ll_n << "ll: " << ll << std::endl;
  
     Eigen::MatrixXd JJ0n(m_ln+lth_Sigma_chol_vechn,lthbb+lth_Omega_chol_vech+lth_Sigma_chol_vech);
     JJ0n.setZero();
   
     JJ0n.block(0,0,m_ln,lthbb)= X_n;
 
     if(lRE!=0){
       JJ0n.block(m_ln,lthbb,lth_Sigma_chol_vechn,lth_Omega_chol_vech) = L_JcholLn * kroneckerProduct(xRE,xRE) *  D_JcholO;
     }
   
     // calculate elimination matrix. 
  //   //Rcout << m_l << std::endl; 
     
     Eigen::MatrixXd IIln = elim_ind(m_l,ind_obs); 
   
  //   //Rcout << IIln.rows() << ",c:" << IIln.cols() << std::endl; 
  //   //Rcout << "lth_r:" << lth_Sigma_chol_vechn << ",c:" << lth_Sigma_chol_vech << std::endl; 
     
     JJ0n.block(m_ln,lthbb+lth_Omega_chol_vech,lth_Sigma_chol_vechn,lth_Sigma_chol_vech) = IIln;
   
  
     Eigen::MatrixXd Jn = JJ0n*JJ1*JJ2; // d(b,vech(Omega),vech(Sigma)/dtheta = d(b,vech(Omega),vech(Sigma))/d(thb,vech(LO),vech(LS))* d(thb...)/dtheta.
   
     Eigen::VectorXd gradCML = grad_ordered_CML(xb, y_n, Lambda, alt, dtauk, cml_pair_type);
     Eigen::MatrixXd hessCML = hess_ordered_CML_approx(xb, y_n, Lambda, alt, dtauk, cml_pair_type);
   
     // Rcout << "gradCML" << gradCML << std::endl;
     Eigen::MatrixXd J_full(Jn.cols()+alt-1,Jn.rows()+alt-1);
     J_full.setZero();
   
     J_full.block(0,0,Jn.cols(),Jn.rows()) = Jn.transpose();
     J_full.block(Jn.cols(),Jn.rows(),alt-1,alt-1) = I_alt;
     
     // reorder, since the tauk parameters occur in front of state space parameters. 
     Eigen::MatrixXd J_full2 = J_full;
     
     J_full2.block(lthb+lthO+lthL,0,alt-1,J_full.cols()) = J_full.block(lthb+lthO+lthL+ 2*s*n,0,alt-1,J_full.cols());
     J_full2.block(lthb+lthO+lthL+alt-1,0,2*s*n,J_full.cols()) = J_full.block(lthb+lthO+lthL,0,2*s*n,J_full.cols());
     
   
     grad_n = J_full2 * gradCML;
     Hess_approx_n2 = J_full2 * hessCML * J_full2.transpose();
   
     /////////////////////////////////////////////////////////
     //// Hessian     ///
     //////////////////////////////////////////////////////////
  //   /// TODO: adjust Hessian ///
  //   ////////////////////////////
  // 
     if (hess == 1){
       Hess_n = Hess_approx_n2;
     }
  // 
     // return value.
     ll += ll_n; //it is already negative in the pairs loop
   
     grad += grad_n;
     Hess_approx2 += Hess_approx_n2;
      
     if (hess == 1){
        Hess += Hess_n;//it is already negative in the pairs loop
     }
     Hess_approx += grad_n * grad_n.transpose();//is positive throughout because it is only used this way
   
  } //end N-loop
   
   

  // penalize unstable matrices A 
  if (syst.factor>0.99){
    ll = ll * std::exp((syst.factor-0.99)*100);
  }

  //............................................................//
  // Format the output
  //............................................................//
  out = -ll;
  
  // write gradient and Hessian into output. 
  out.attr("gradient") = -(grad.array());
   
  if(hess == 1){
     out.attr("hessian") =-(Hess.array());
  }
   
  out.attr("hessian1") =(Hess_approx.array());
  out.attr("hessian2") =(Hess_approx2.array());
  return out;

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
//' @param time 
//' vector of all observation times 
//' @param quest
//' integer of most questions asked 
//' @param timen
//' vector of times when the n-th individual is observed. 
//' @param questn
//' vector of questions answered by the n-th individual. 
//' 
//' @return
//' A matrix, containing the predicted choice probabilities for each choice marginally (ignoring correlations)
//' 
//' @export
//'
// [[Rcpp::export]]
Eigen::MatrixXd pred_probit_ordered_approx_StSp(Eigen::VectorXd theta, Eigen::MatrixXd Xn, Eigen::VectorXd yn, Rcpp::List mod, Eigen::VectorXd time, int quest, Eigen::VectorXd timen, Eigen::VectorXd questn)
{

  std::string approx_method = "SJ";

  //............................................................//
  // INITIALIZE OUTPUT
  //............................................................//
  int n = yn.rows();
  int alt = mod["alt"]; 
  
  Eigen::MatrixXd prob_predict(n,alt);
  prob_predict.setZero();
  
  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//

  struct system syst = build_system_from_model(theta, mod, time);
  // Recover the vector b
  Eigen::VectorXd b = syst.beta;
  int lthbb = b.rows();

  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma = syst.Sigma;

  // Recover GammaT matrix:
  Eigen::MatrixXd GammaT = syst.GammaT; 

  // check, which measurements are contained 
  Eigen::VectorXd ind_obs = get_observation_indices(time,quest,timen, questn);
  
  Eigen::MatrixXd Lambda = subset(GammaT,ind_obs); 
  
  // Recover the matrix Omega
  Eigen::MatrixXd Omega = syst.Omega; 

  // get the cutoff points
  Eigen::VectorXd tauk = syst.tauk; 


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
  //Eigen::MatrixXd Lambda(n,n);
  //Lambda.setZero();
  //Lambda = GammaT;
 
  int lRE = Omega.rows(); 
  Eigen::MatrixXd xRE(n,lRE);
 
  if(lRE!=0){
    xRE = Xn.block(0,0,n,lRE);
    Lambda+= xRE * Omega * xRE.transpose();
  }
   
  double xbn = 1.0;
   
  Eigen::MatrixXd lambda = Lambda.diagonal().array().sqrt();
   
  //Rcout << lambda << std::endl; 
   
  for(int i_choice_1 = 0; i_choice_1 < n; i_choice_1++){
    xbn = xb(i_choice_1) ;
    prob_predict(i_choice_1,0) = std_normal_cdf((tauk(0)-xbn)/ lambda(i_choice_1));
     for (int y=1; y < alt-1; y++){
      prob_predict(i_choice_1,y) = std_normal_cdf((tauk(y)-xbn)/ lambda(i_choice_1))-std_normal_cdf((tauk(y-1)-xbn)/ lambda(i_choice_1));
    }
  
    // last probability
    prob_predict(i_choice_1,alt-1) = 1.0- std_normal_cdf((tauk(alt-2)-xbn)/ lambda(i_choice_1));
   
  } //end i_choice_1 loop over choice decisions.
  

  //............................................................//
  // Format the output
  //............................................................//

  return  prob_predict;

}


