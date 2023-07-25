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


// Function for the prior mean of alpha_1
// [[Rcpp::export]]
arma::vec a1_fn_ip(const arma::vec& theta, const arma::vec& known_params) {
  arma::vec a1(7);
  a1(0) = known_params(0);
  a1(1) = known_params(1);
  a1(2) = known_params(2);
  a1(3) = known_params(3);
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


// T function
// [[Rcpp::export]]
arma::vec T_fn_ip(const unsigned int t, const arma::vec& alpha,
               const arma::vec& theta, const arma::vec& known_params,
               const arma::mat& known_tv_params) {

  // known parms
  double momega = known_params(4);
  double mdelta = known_params(5);
  double mpsi = known_params(6);

  // I,P,N,D,  omega,delta,psi
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
  // dynamics
  arma::vec alpha_new(7);
  // log I: I_{t+1} = exp(θ_3) ((1-p)I_{t} + p.P_t);   p=ilogit( θ_4 ) 
  // alpha_new(0) = log( (1-pr) * I + pr * P ) + theta(3); // see above
  alpha_new(0) = (1-pr) * alpha(0) + pr * alpha(1)  + theta(3); // log weighted average
  // log P
  alpha_new(1) = log(I * rho + P * exp(-OM) );
  // log N
  alpha_new(2) = log(delta) + log(Q); // change in this and below
  // log D
  alpha_new(3) = log(psi) + log(omega) + log(Q);
  // log-omega
  alpha_new(4) = alpha(4);
  // log-delta
  alpha_new(5) = alpha(5); // NOTE interventions not here
  // logit-psi
  alpha_new(6) = mpsi;

  return alpha_new;
}

// Jacobian of T function
// [[Rcpp::export]]
arma::mat T_gn_ip(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {

  // I,P,N,D,  omega,delta,psi
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
  arma::mat Tg(7, 7, arma::fill::zeros);
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
  return Tg;
}



// log-prior pdf for theta WITH hyperparams
// [[Rcpp::export]]
double log_prior_pdf_ip_HP(const arma::vec& theta,const arma::vec& hyperparms) {

  // weakly informative half-N(0, 4) priors. 
  // Note that the sampling is on log-scale, 
  // so we need to add jacobians of the corresponding transformations
  // we could also sample on natural scale with check such as
  // if(arma::any(theta < 0))
  //  return -std::numeric_limits<double>::infinity();
  // but this would be less efficient.

  // You can use R::dnorm and similar functions, see, e.g.
  // https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  double noisehp = hyperparms(0);
  double arhpmu = hyperparms(1);
  double arhpsg = hyperparms(2);
  double prhpmu = hyperparms(3);
  double prhpsg = hyperparms(4);
  double log_pdf =
    R::dnorm(theta(0), 0, noisehp, true) + // I noise
    R::dnorm(theta(1), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(2), 0, noisehp, true)+   // delta noise NOTE change
    R::dnorm(theta(3), arhpmu, arhpsg, true)+   // phiIP  AR
    R::dnorm(theta(4), prhpmu, prhpsg, true)+theta(4)*theta(4) // logitN w/ below
    - arma::accu(theta); //jacobian term
  // Rprintf("the value of log-likelihood : %f \n", log_pdf);
  return log_pdf;
}



// [[Rcpp::export]]
Rcpp::List create_xptrs_ip(const arma::vec& hyperparms) {

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

  // typedef
  typedef double PartialFn(const arma::vec& theta);
  typedef PartialFn* PartialFnPtr;
  // create lambda to prior function
  // auto log_prior_pdf_ipa = [](const arma::vec& theta) {
  //   return log_prior_pdf_ip_HP(std::forward<const arma::vec&>(theta),hyperparms);
  // };

  // new version
  PartialFn* log_prior_pdf_ipa = new auto([&](const arma::vec& theta) {
    return log_prior_pdf_ip_HP(std::forward<const arma::vec&>(theta), hyperparms);
  });


  // Pointer to prior fn with PartialFn Prototype
  PartialFn* log_prior_pdf_ipa_ptr = log_prior_pdf_ipa;
  // PartialFn* foo = log_prior_pdf_ip;

  return Rcpp::List::create(
                            Rcpp::Named("a1_fn_ip") = Rcpp::XPtr<a1_fnPtr>(new a1_fnPtr(&a1_fn_ip)),
                            Rcpp::Named("P1_fn_ip") = Rcpp::XPtr<P1_fnPtr>(new P1_fnPtr(&P1_fn_ip)),
                            Rcpp::Named("Z_fn_ip") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&Z_fn_ip)),
                            Rcpp::Named("H_fn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&H_fn_ip)),
                            Rcpp::Named("T_fn_ip") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&T_fn_ip)),
                            Rcpp::Named("R_fn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R_fn_ip)),
                            Rcpp::Named("Z_gn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&Z_gn_ip)),
                            Rcpp::Named("T_gn_ip") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&T_gn_ip)),
                            Rcpp::Named("log_prior_pdf_ip") =
                            // Rcpp::XPtr<PartialFnPtr>(new PartialFnPtr(&log_prior_pdf_ipa_ptr))
                            Rcpp::XPtr<PartialFnPtr>(new PartialFnPtr(log_prior_pdf_ipa))
                            // Rcpp::XPtr<PartialApp>(new PartialApp (&log_prior_pdf_ipa))
                            );
}

