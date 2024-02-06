// 10-D states as: I,P,N,D, omega,delta,psi; Pp,omegap, irr
// means of initial states - mI0,mP0,mN0,mD0
// sd of initial states - sI0,sP0,sN0,sD0
// arrange as: mI0 - 0,mP0 - 1,mN0 - 2,mD0 - 3,momega - 4, mdelta - 5, mpsi - 6,
//      sI0 - 7,sP0 - 8,sN0 - 9,sD0 - 10,somega - 11,sdelta - 12,spsi - 13
// matrix known parms: the SD of the measurements

// first go at TB SSM
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// Unknown parameters theta:


// // a couple of utility functions
// double odds (double x){
//   return x/(1.0-x);
// }
// double iodds( double x){
//   return x/(1.0+x);
// }

// NOTE on odds
// iodds(R*odds(x)) = 
// [R *x/(1-x)]/[1+R *x/(1-x)] = [R*x]/[1-x + R*x]

// Function for the prior mean of alpha_1
// [[Rcpp::export]]
arma::vec a1_fn_ip(const arma::vec& theta, const arma::vec& known_params) {
  arma::vec a1(7);
  a1(0) = known_params(0) + theta(5);
  a1(1) = known_params(1);
  a1(2) = known_params(2) + theta(6);
  a1(3) = known_params(3) + theta(7);
  a1(4) = known_params(4);
  a1(5) = known_params(5);
  a1(6) = known_params(6);
  return a1;
}

// Function for the prior covariance matrix of alpha_1
// [[Rcpp::export]]
arma::mat P1_fn_ip(const arma::vec& theta, const arma::vec& known_params) {
  arma::mat P1(7, 7, arma::fill::zeros);
  P1(0,0) = std::pow(known_params(7),2);
  P1(1,1) = std::pow(known_params(8),2);
  P1(2,2) = std::pow(known_params(9),2);
  P1(3,3) = std::pow(known_params(10),2);
  P1(4,4) = std::pow(known_params(11),2);
  P1(5,5) = std::pow(known_params(12),2);
  P1(6,6) = std::pow(known_params(13),2);
  return P1;
}


// Function for the observational level standard deviation
// [[Rcpp::export]]
arma::mat H_fn_ip(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {
  arma::mat H(4,4);
  H(0,0) = known_tv_params(t,0);
  H(1,1) = known_tv_params(t,1);
  H(2,2) = known_tv_params(t,2);
  H(3,3) = known_tv_params(t,3);
  return H;
}

// NOTE chol
// Function for the Cholesky of state level covariance
// [[Rcpp::export]]
arma::mat R_fn_ip(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {

  double sI = exp(theta(0));
  // double sdelta = known_params(12);
  double sdelta = exp(theta(1)); // now unknown * inferred NOTE change
  // double somega = known_params(11);
  double somega = exp(theta(2)); // now also unknown
  double spsi = known_params(13);
  arma::mat R(7, 4, arma::fill::zeros);
  // I,P,N,D,  omega,delta,psi
  R(0, 0) = sI;                 // this is the unknown parm
  R(4, 1) = somega;
  R(5, 2) = sdelta;
  R(6, 3) = spsi;
  return R;
}


// Z function
// [[Rcpp::export]]
arma::vec Z_fn_ip(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {
  // I,P,N,D,  omega,delta,psi
  arma::vec Z(4);
  Z(0) = alpha(0);
  Z(1) = alpha(1);
  Z(2) = alpha(2);
  Z(3) = alpha(3);
  return Z;
}


// Jacobian of Z function
// [[Rcpp::export]]
arma::mat Z_gn_ip(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {

  arma::mat Z(4, 7, arma::fill::zeros);
  Z(0,0) = 1.0;
  Z(1,1) = 1.0;
  Z(2,2) = 1.0;
  Z(3,3) = 1.0;
  return Z;
}


// deltap / (Omegap + deltap) = delta / (Omega + delta) = CDR
// delta/Omega = iodds(CDR) = deltap/Omegap

// T function
// [[Rcpp::export]]
arma::vec T_fn_ipH(const unsigned int t, const arma::vec& alpha,
               const arma::vec& theta, const arma::vec& known_params,
               const arma::mat& known_tv_params) {
  // 10-D states as: I,P,N,D, omega,delta,psi; Pp,omegap, irr
  // known parms
  double momega = known_params(4);
  double mdelta = known_params(5);
  double mpsi = known_params(6);
  // I,P,N,D,  omega,delta,psi, omegap, irr
  double I = exp(alpha(0));     // now transformed
  double P = exp(alpha(1));     // now transformed
  double N = exp(alpha(2)); double D = exp(alpha(3)); // also transformed
  double omega = exp(alpha(4));
  double delta = exp(alpha(5) + known_tv_params(t,5)); // NOTE interventions
  double psi = 1.0 / (1.0 + exp(-alpha(6))); // expit
  double OM = omega + delta;
  double rho = ( 1-exp(-OM))/OM; // useful
  double Q = (1-rho) * I / OM + rho * P; // cumulative P for use in D & N
  double pr = 1.0 / (1.0 + exp(-theta(4)));
  // PLHIV
  double Pp = exp(alpha(7));
  double omegap = exp(alpha(8));
  double irr = exp(alpha(9));
  double h = known_tv_params(t,5);  // NOTE 5 is HIV prevalence
  double H = irr * h / (1-h+irr*h); // NOTE see NOTE on odds
  double In = I * H;
  double Ip = I * (1-H);
  double deltap = omegap * delta / omega; // see NOTE above; enforces same CDR in HIV+
  double OMp = omegap + deltap;
  double rhop = ( 1-exp(-OMp))/OMp; // useful
  double Qp = (1-rhop) * Ip / OMp + rhop * Pp; // cumulative P for use in D & N
  // dynamics
  arma::vec alpha_new(10);
  // I_{t+1} = exp(θ_3) * (I_{t}^{1-p} * P_t^p);   p=ilogit( θ_4 );
  // NOTE just use HIV-ve prevalence
  // alpha_new(0) = (1-pr) * alpha(0) + pr * logsumexp(alpha(1),alpha(7)) + theta(3) + known_tv_params(t,4); //
  alpha_new(0) = (1-pr) * alpha(0) + pr * alpha(1) + theta(3) + known_tv_params(t,4); //
  // log P: inflow only HIV-ve TB incidence
  alpha_new(1) = log(In * rho + P * exp(-OM) );
  // log N
  alpha_new(2) = logsumexp(log(delta) + log(Q), log(deltap) + log(Qp)); // 
  // log D
  alpha_new(3) = logsumexp(log(psi) + log(omega) + log(Q),
                           log(deltap) + log(Qp)); // NOTE assuming untreated CFR=1 -> log(psip) = 0
  // log-omega
  alpha_new(4) = alpha(4);
  // log-delta
  alpha_new(5) = alpha(5); // NOTE interventions not here
  // logit-psi
  alpha_new(6) = mpsi;
  // log-Pp
  alpha_new(7) = log(Ip * rhop + Pp * exp(-OMp) ); //
  // log-omegap
  alpha_new(8) = alpha(8);
  // log-irr
  alpha_new(9) = alpha(9);
  return alpha_new;
}

// Jacobian of T function
// [[Rcpp::export]]
arma::mat T_gn_ipH(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {
  // I,P,N,D,  omega,delta,psi, omegap, irr
  double I = exp(alpha(0));
  double P = exp(alpha(1));
  double N = exp(alpha(2)); double D = exp(alpha(3));
  double omega = exp(alpha(4));
  double delta = exp(alpha(5) + known_tv_params(t,5)); // NOTE interventions
  double psi = 1.0 / (1.0 + exp(-alpha(6))); // expit
  // P = I (1-exp(-O))/O + P exp(-O)
  // dP/domega = dP/ddelta = -P exp(-O) + I (exp(-O)-(1-exp(-O))/O^2)
  double OM = omega + delta;
  double rho = ( 1-exp(-OM))/OM; // useful
  double drdd = exp(-OM)/OM - ( 1-exp(-OM))/(OM*OM); // drho/ddelta
  double Q = (1-rho) * I / OM + rho * P; // cumulative P for use in D & N
  double dQdo = ((P-I/OM)*drdd-(1-rho)*I/(OM*OM));
  double dPt1do = -P * exp(-OM) + I * drdd;
  double Pt1 = I * ( 1-exp(-OM))/OM + P * exp(-OM); // P_{t+1} see above
  double Nt1 = delta * Q;
  double Dt1 = psi * omega * Q;
  double pr = 1.0 / (1.0 + exp(-theta(4)));
  // I,P,N,D,  omega,delta,psi
  arma::mat Tg(9, 9, arma::fill::zeros);
  // log I: I_{t+1} = exp(θ_3) ((1-p)I_{t} + p.P_t);   p=ilogit( θ_4 ) 
  // Tg(0,0) = (1-pr) * I / log( (1-pr) * I + pr * P );   // dlog I/d log I
  // Tg(0,1) = pr * P /log( (1-pr) * I + pr * P );   // dlog I/d log P
  Tg(0,0) = (1-pr);   // dlog I/d log I
  Tg(0,1) = pr;   // dlog I/d log P
  // log P
  Tg(1,0) =  rho*I/Pt1; // dlogP/dlogI
  Tg(1,1) = exp(-OM) *  P/ Pt1;                   // dlogP/dlogP
  Tg(1,4) = dPt1do * omega / Pt1;           // dP/domega * domega/dlomega
  Tg(1,5) = dPt1do * delta / Pt1;           // dP/ddelta * ddelta/dldelta
  // log N
  Tg(2,0) = ((1-rho) * I / OM)/Q; // dlogN/dlogI
  Tg(2,1) = (rho * P)/Q; // dlogN/dlogP
  Tg(2,4) = omega * dQdo/Q; // dlogN/dlogomega
  Tg(2,5) = 1.0 + delta * dQdo/Q; // dlogN/dlogdelta
  // log D
  Tg(3,0) = ((1-rho) * I / OM)/Q;// dlogD/logI
  Tg(3,1) = rho * P/Q;// dlogD/logP
  Tg(3,4) = 1.0 + omega * dQdo/Q;// dlogD/logomega
  Tg(3,5) = delta * dQdo/Q;// dlogD/logdelta
  Tg(3,6) = 1 - psi;// dlogD/logitpsi
  // NOTE log(psi) = -log(1+exp(-logitpsi))->d/d=exp(-lp)/(1+exp(-lp))=1-psi
  // log-omega
  Tg(4,4) = 1.0;
  // log-delta
  Tg(5,5) = 1.0;
  // log-omegap
  Tg(7,7) = 1.0;
  // log-irr
  Tg(8,8) = 1.0;

  return Tg;
}


double logsumexp(double logA, double logB){
  if(logA > logB){
    return logA + std::log1p(exp(logB-logA));
  } else{
    return logB + std::log1p(exp(logA-logB));
  }
}

// double logsumexp(double nums[], size_t ct) {
//   double max_exp = nums[0], sum = 0.0;
//   size_t i;

//   for (i = 1 ; i < ct ; i++)
//     if (nums[i] > max_exp)
//       max_exp = nums[i];

//   for (i = 0; i < ct ; i++)
//     sum += exp(nums[i] - max_exp);

//   return log(sum) + max_exp;
// }

// #include "logsumexp.h"

// #ifdef HAVE_LONG_DOUBLE
// #  define LDOUBLE long double
// #  define EXPL expl
// #else
// #  define LDOUBLE double
// #  define EXPL exp
// #endif

// // [[Rcpp::export]]
// double logSumExp(const arma::vec& x) {
//   unsigned int maxi = x.index_max();
//   LDOUBLE maxv = x(maxi);
//   if (!(maxv > -arma::datum::inf)) {
//     return -arma::datum::inf;
//   }
//   LDOUBLE cumsum = 0.0;
//   for (unsigned int i = 0; i < x.n_elem; i++) {
//     if ((i != maxi) && (x(i) > -arma::datum::inf)) {
//       cumsum += EXPL(x(i) - maxv);
//     }
//   }
  
//   return maxv + log1p(cumsum);
// }

// ================ hyperprior variants ==================

// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ip4(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), +4, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}


// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ip3(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), +3, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}


// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ip2(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), +2, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}



// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ip1(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), +1, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}


// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ip0(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), 0.0, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}



// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ipn1(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), -1, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}

// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ipn2(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), -2, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}

// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ipn3(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), -3, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}


// log-prior pdf for theta
// [[Rcpp::export]]
double log_prior_pdf_ipn4(const arma::vec& theta) {
  double ISW = 0.5;
  double noisehp = 1.0;
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), -1, 0.5, true)+   // phiIP
    R::dnorm(theta(4), -4, 0.1, true)+ 2 * log(1+exp(theta(4))) // logitN w/ below
    // start extra on IS
    +R::dnorm(theta(5), 0, ISW, true)+R::dnorm(theta(6), 0, ISW, true)+R::dnorm(theta(7), 0, ISW, true)
    // end extra on IS
    - arma::accu(theta); //jacobian term
  return log_pdf;
}

// [[Rcpp::export]]
Rcpp::List create_xptrs_ip_all() {

  // typedef for a pointer of nonlinear function returning vec (T, Z)
  typedef arma::vec (*nvec_fnPtr)(const unsigned int t, 
                                  const arma::vec& alpha, const arma::vec& theta, 
                                  const arma::vec& known_params, const arma::mat& known_tv_params);
  // for a pointer of nonlinear function returning mat (Tg, Zg, H, R)
  typedef arma::mat (*nmat_fnPtr)(const unsigned int t, 
                                  const arma::vec& alpha, const arma::vec& theta, 
                                  const arma::vec& known_params, const arma::mat& known_tv_params);

  // typedef for a pointer returning a1
  typedef arma::vec (*a1_fnPtr)(const arma::vec& theta, 
                                const arma::vec& known_params);
  // typedef for a pointer returning P1
  typedef arma::mat (*P1_fnPtr)(const arma::vec& theta, 
                                const arma::vec& known_params);
  // typedef for a pointer of log-prior function
  typedef double (*prior_fnPtr)(const arma::vec& theta);

  return Rcpp::List::create(
                            Rcpp::Named("a1_fn_ip") = Rcpp::XPtr<a1_fnPtr>(new a1_fnPtr(&a1_fn_ip)),
                            Rcpp::Named("P1_fn_ip") = Rcpp::XPtr<P1_fnPtr>(new P1_fnPtr(&P1_fn_ip)),
                            Rcpp::Named("Z_fn_ip") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&Z_fn_ip)),
                            Rcpp::Named("H_fn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&H_fn_ip)),
                            Rcpp::Named("T_fn_ip") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&T_fn_ip)),
                            Rcpp::Named("R_fn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R_fn_ip)),
                            Rcpp::Named("Z_gn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&Z_gn_ip)),
                            Rcpp::Named("T_gn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&T_gn_ip)),
                            Rcpp::Named("log_prior_pdf_ip4") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ip4)),
                            Rcpp::Named("log_prior_pdf_ip3") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ip3)),
                            Rcpp::Named("log_prior_pdf_ip2") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ip2)),
                            Rcpp::Named("log_prior_pdf_ip1") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ip1)),
                            Rcpp::Named("log_prior_pdf_ip0") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ip0)),
                            Rcpp::Named("log_prior_pdf_ipn1") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ipn1)),
                            Rcpp::Named("log_prior_pdf_ipn2") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ipn2)),
                            Rcpp::Named("log_prior_pdf_ipn3") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ipn3)),
                            Rcpp::Named("log_prior_pdf_ipn4") = 
                            Rcpp::XPtr<prior_fnPtr>(new prior_fnPtr(&log_prior_pdf_ipn4)));
}


