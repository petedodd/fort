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
arma::vec a1_fn_ipH(const arma::vec& theta, const arma::vec& known_params) {
  arma::vec a1(10);
  a1(0) = known_params(0) + theta(5);
  a1(1) = known_params(1);
  a1(2) = known_params(2) + theta(6);
  a1(3) = known_params(3) + theta(7);
  a1(4) = known_params(4);
  a1(5) = known_params(5);
  a1(6) = known_params(6);
  a1(7) = known_params(7);
  a1(8) = known_params(8);
  a1(9) = known_params(9);
  return a1;
}

// Function for the prior covariance matrix of alpha_1
// [[Rcpp::export]]
arma::mat P1_fn_ipH(const arma::vec& theta, const arma::vec& known_params) {
  arma::mat P1(10, 10, arma::fill::zeros);
  P1(0,0) = std::pow(known_params(7),2);
  P1(1,1) = std::pow(known_params(8),2);
  P1(2,2) = std::pow(known_params(9),2);
  P1(3,3) = std::pow(known_params(10),2);
  P1(4,4) = std::pow(known_params(11),2);
  P1(5,5) = std::pow(known_params(12),2);
  P1(6,6) = std::pow(known_params(13),2);
  P1(7,7) = std::pow(known_params(14),2);
  P1(8,8) = std::pow(known_params(15),2);
  P1(9,9) = std::pow(known_params(16),2);
  return P1;
}

// NOTE chol
// Function for the Cholesky of state level covariance
// [[Rcpp::export]]
arma::mat R_fn_ipH(const unsigned int t, const arma::vec& alpha, 
               const arma::vec& theta, const arma::vec& known_params, 
               const arma::mat& known_tv_params) {

  double sI = exp(theta(0));
  // double sdelta = known_params(12);
  double sdelta = exp(theta(1)); // now unknown * inferred NOTE change
  // double somega = known_params(11);
  double somega = exp(theta(2)); // now also unknown
  double spsi = known_params(13);
  arma::mat R(10, 6, arma::fill::zeros);
  // I,P,N,D,  omega,delta,psi, Pp,omegap, irr
  R(0, 0) = sI;                 // this is the unknown parm
  R(4, 1) = somega;
  R(5, 2) = sdelta;
  R(6, 3) = spsi;
  R(8, 4) = spsi;
  R(9, 5) = spsi;
  // TODO consider extra parms
  return R;
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
  // I,P,N,D,  omega,delta,psi, Pp,omegap, irr
  double I = exp(alpha(0));     // now transformed
  double P = exp(alpha(1));     // now transformed
  double N = exp(alpha(2)); double D = exp(alpha(3)); // also transformed
  double omega = exp(alpha(4));
  double delta = exp(alpha(5) + known_tv_params(t,5)); // NOTE interventions
  double psi = 1.0 / (1.0 + exp(-alpha(6))); // expit
  double OM = omega + delta;
  double rho = ( 1-exp(-OM))/OM; // useful
  double Q = (1-rho) * In / OM + rho * P; // cumulative P for use in D & N
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
  // PLHIV
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
  double Q = (1-rho) * In / OM + rho * P; // cumulative P for use in D & N
  double dQdo = ((P-In/OM)*drdd-(1-rho)*In/(OM*OM));
  double dPt1do = -P * exp(-OM) + In * drdd;
  double Pt1 = In * ( 1-exp(-OM))/OM + P * exp(-OM); // P_{t+1} see above
  double Nt1 = delta * Q + deltap * Qp;
  double Dt1 = psi * omega * Q + omegap + Qp;
  double pr = 1.0 / (1.0 + exp(-theta(4)));
  // PLHIV
  double Pp = exp(alpha(7));
  double omegap = exp(alpha(8));
  double irr = exp(alpha(9));
  double h = known_tv_params(t,5);  // NOTE 5 is HIV prevalence
  double H = irr * h / (1-h+irr*h); // NOTE see NOTE on odds
  double dHdJ = (1-h) * h / ((1-h+irr*h)*(1-h+irr*h)); // dH/dIRR
  double In = I * H;
  double Ip = I * (1-H);
  double Ppt1 = Ip * ( 1-exp(-OMp))/OMp + Pp * exp(-OMp); // P^+_{t+1} see above
  double deltap = omegap * delta / omega; // see NOTE above; enforces same CDR in HIV+
  double OMp = omegap + deltap;
  double rhop = ( 1-exp(-OMp))/OMp; // useful
  double Qp = (1-rhop) * Ip / OMp + rhop * Pp; // cumulative P for use in D & N
  // double dQpdo = ((Pp-Ip/OMp)*drpdo*dopdo-(1-rhop)*Ip/(OMp*OMp));
  double dQpdOp = ((Pp-Ip/OMp)*drpdOp-(1-rhop)*Ip/(OMp*OMp)); // dQp/dOmegap
  double drpdOp = exp(-OMp)/OMp - ( 1-exp(-OMp))/(OMp*OMp); // drhop/dOmegap
  // double dopdo = - deltap / omega;    // dOMEGAp/domega [OMEGAp=omegap*(1+delt/om)]
  double dPpt1dOp = -Pp * exp(-OMp) + Ip * drpdOp; // dP^+_{t-1}/dOMEGAp
  // I,P,N,D,  omega,delta,psi, Pp,omegap, irr
  arma::mat Tg(10, 10, arma::fill::zeros);
  // log I: I_{t+1} = exp(θ_3) ((1-p)I_{t} + p.P_t);   p=ilogit( θ_4 ) 
  // Tg(0,0) = (1-pr) * I / log( (1-pr) * I + pr * P );   // dlog I/d log I
  // Tg(0,1) = pr * P /log( (1-pr) * I + pr * P );   // dlog I/d log P
  Tg(0,0) = (1-pr);   // dlog I/d log I
  Tg(0,1) = pr;   // dlog I/d log P
  // log P
  Tg(1,0) =  rho*In/Pt1; // dlogP/dlogI
  Tg(1,1) = exp(-OM) *  P/ Pt1;                   // dlogP/dlogP
  Tg(1,4) = dPt1do * omega / Pt1;           // dP/domega * domega/dlomega
  Tg(1,5) = dPt1do * delta / Pt1;           // dP/ddelta * ddelta/dldelta
  // log N
  //    dlogN/dlogI:
  Tg(2,0) = (delta * (1-rho) * In / OM + deltap * (1-rhop) * Ip / OMp) / Nt1;
  //    dlogN/dlogP
  Tg(2,1) = (delta * rho * P) / Nt1;
  //    dlogN/dlogomega
  Tg(2,4) = (delta * omega * dQdo  - deltap * Qp - deltap * deltap * dQpdOp)/ Nt1;
  //    dlogN/dlogdelta
  Tg(2,5) = (Nt + delta*delta*dQdo + deltap*dQpdOp ) / Nt;
  //    dlogN/dlogPp
  Tg(2,7) = (deltap * rhop * Pp) / Nt1;
  //    dlogN/dlogomegap
  Tg(2,8) = (delta *  Qp + OMp * dQpdOp)/ Nt1;
  //    dlogN/dlogirr
  Tg(2,9) = irr * deltap * (1-rhop) * I * dHdJ / (Nt * OMp);

  // log D: template=
  // (x/Dt1) * (Q * dompsidx + psi * omega * dQdo * dodx + Qp * ddpdx + deltap * dQpdOp * dOMpdx)
  //    dlogD/logI
  Tg(3,0) = (I/Dt1) * ((1-rho) * (1-H) / OM + (1-rhop) * H / OMp);
  //    dlogD/logP
  Tg(3,1) = (I/Dt1) * (rho * P);
  //    dlogD/logomega
  Tg(3,4) = (omega * Q * psi + psi*omega*omega * dQdo  - deltap*(Qp + deltap*dQpdOp) ) / Dt1;
  //    dlogD/logdelta
  Tg(3,5) = (deltap/Dt1) * ( Qp  + deltap * dQpdOp);
  //    dlogD/logitpsi
  Tg(3,6) = (Q * omega * psi * (1-psi) )/Dt1;
  // NOTE d psi/ dloig(psi) = 1/(logit psi)' = psi * (1-psi)
  //    dlogD/dlogPp
  Tg(3,7) = (deltap * rhop * Pp) / Dt1;
  //    dlogD/dlogomegap
  Tg(3,8) = (Qp  + OMp * dQpdOp ) * deltap / Dt1;
  //    dlogD/dlogirr
  Tg(3,9) = (irr/Dt1) * ( (1-rhop) / OMp - (1-rho) / OM ) * I * dHdJ;
  // log-omega
  Tg(4,4) = 1.0;
  // log-delta
  Tg(5,5) = 1.0;

  // PLHIV
  // log-Pp
  // alpha_new(7) = log(Ip * rhop + Pp * exp(-OMp) ); //TODO
  Tg(7,0) =  rhop*Ip/Ppt1; // dlogP/dlogI
  Tg(7,4) = -dPpt1dOp * deltap / Ppt1;          // dPp/dOmegap * dOmegap/dlomega
  Tg(7,5) = dPpt1dOp * deltap / Ppt1;           // dPp/dOmegap * dOmegap/dldelta
  Tg(7,7) = exp(-OMp) *  Pp/ Ppt1;              // dlogP/dlogPp
  Tg(7,8) = dPpt1dOp * OMp / Ppt1;              // dlogP/dlog omegap
  Tg(7,9) = I * rhop dHdJ/ Ppt1;              // dlogP/dlog IRR
  // log-omegap
  Tg(8,8) = 1.0;
  // log-irr
  Tg(9,9) = 1.0;


  return Tg;
}


double logsumexp(double logA, double logB){
  if(logA > logB){
    return logA + std::log1p(exp(logB-logA));
  } else{
    return logB + std::log1p(exp(logA-logB));
  }
}


// [[Rcpp::export]]
Rcpp::List create_xptrs_H_all() {

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
  return Rcpp::List::create(
                            Rcpp::Named("a1_fn_ipH") = Rcpp::XPtr<a1_fnPtr>(new a1_fnPtr(&a1_fn_ipH)),
                            Rcpp::Named("P1_fn_ipH") = Rcpp::XPtr<P1_fnPtr>(new P1_fnPtr(&P1_fn_ipH)),
                            Rcpp::Named("T_fn_ipH") = Rcpp::XPtr<nvec_fnPtr>(new nvec_fnPtr(&T_fn_ipH)),
                            Rcpp::Named("R_fn_ipH") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&R_fn_ipH)),
                            Rcpp::Named("T_gn_ipH") = Rcpp::XPtr<nmat_fnPtr>(new nmat_fnPtr(&T_gn_ipH)));
}


