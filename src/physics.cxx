/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "physics.h"
#include "TLorentzVector.h"

// Calcuating Q^2
//	Gotten from t channel
// -q^mu^2 = -(e^mu - e^mu')^2 = Q^2
double physics::Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  auto q_mu = (e_mu - e_mu_prime);
  return -q_mu.M2();
}
//	Calcualting W
//	Gotten from s channel [(gamma + P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double physics::W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  auto q_mu = (e_mu - e_mu_prime);
  auto p_mu = TLorentzVector(0, 0, 0, MASS_P);
  return (p_mu + q_mu).M();
}

double physics::Q2_calc(const Lorentz &gamma_mu) { return -gamma_mu->M2(); }

double physics::W_calc(const Lorentz &gamma_mu) {
  Lorentz p_mu = std::make_shared<TLorentzVector>(0, 0, 0, MASS_P);
  return (*p_mu + *gamma_mu).M();
}

double physics::xb_calc(const Lorentz &gamma_mu) {
  double Q2 = Q2_calc(gamma_mu);
  Lorentz p_mu = std::make_shared<TLorentzVector>(0, 0, 0, MASS_P);
  return (Q2 / (2 * (gamma_mu->Dot(*p_mu))));
}
