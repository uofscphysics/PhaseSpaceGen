/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include <memory>
#include "constants.h"

namespace physics {
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
double Q2_calc(const Lorentz &gamma_mu);
double W_calc(const Lorentz &gamma_mu);
double xb_calc(const Lorentz &gamma_mu);
}  // namespace physics
#endif
