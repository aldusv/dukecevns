#include "xscns.h"

#include <cmath>
#include <stdexcept>

namespace {

constexpr double kCrossSectionScaleCm2 = 8.43103e-45;

} // namespace

double diffxscnvec(
  double neutrinoEnergyMeV,
  double massMeV,
  double recoilMeV)
{
  return kCrossSectionScaleCm2 * massMeV *
    (2.0 - massMeV * recoilMeV / std::pow(neutrinoEnergyMeV, 2)
     - 2.0 * recoilMeV / neutrinoEnergyMeV
     + std::pow(recoilMeV / neutrinoEnergyMeV, 2));
}

double diffxscnaxial(
  double neutrinoEnergyMeV,
  double massMeV,
  double recoilMeV)
{
  return kCrossSectionScaleCm2 * massMeV *
    (2.0 + massMeV * recoilMeV / std::pow(neutrinoEnergyMeV, 2)
     - 2.0 * recoilMeV / neutrinoEnergyMeV
     + std::pow(recoilMeV / neutrinoEnergyMeV, 2));
}

double diffxscninterf(
  double neutrinoEnergyMeV,
  double massMeV,
  double recoilMeV)
{
  return kCrossSectionScaleCm2 * massMeV *
    (4.0 * recoilMeV / neutrinoEnergyMeV
     - 2.0 * std::pow(recoilMeV / neutrinoEnergyMeV, 2));
}

void sm_vector_couplings(int pdgYear, double* couplings)
{
  if (pdgYear == 0) {
    throw std::invalid_argument("custom vector couplings are not supported");
  }

  double proton = 0.0227;
  double neutron = -0.5117;
  if (pdgYear < 2004) {
    proton = 0.0152;
    neutron = -0.5122;
  } else if (pdgYear < 2011) {
    proton = 0.0304;
    neutron = -0.5122;
  } else if (pdgYear < 2012) {
    proton = 0.0306;
    neutron = -0.5120;
  } else if (pdgYear < 2014) {
    proton = 0.0307;
    neutron = -0.5120;
  } else if (pdgYear < 2020) {
    proton = 0.01836;
    neutron = -0.5117;
  }

  couplings[0] = proton;
  couplings[1] = neutron;
}

void sm_axial_couplings(int pdgYear, int flavorSign, double* couplings)
{
  if (pdgYear == 0) {
    throw std::invalid_argument("custom axial couplings are not supported");
  }
  if (flavorSign == 0) {
    throw std::invalid_argument("flavor sign must be nonzero");
  }

  const double sign = flavorSign > 0 ? 1.0 : -1.0;
  double proton = 0.4995 * sign;
  double neutron = -0.5121 * sign;
  if (pdgYear < 2004) {
    constexpr double strangeSpin = -0.15;
    proton = (1.27 - strangeSpin) / 2.0 * sign;
    neutron = -(1.27 - strangeSpin) / 2.0 * sign;
  } else if (pdgYear < 2011) {
    proton = 0.4955 * sign;
    neutron = -0.5125 * sign;
  } else if (pdgYear < 2012) {
    proton = 0.4942 * sign;
    neutron = -0.5123 * sign;
  } else if (pdgYear < 2014) {
    proton = 0.4953 * sign;
    neutron = -0.5124 * sign;
  }

  couplings[0] = proton;
  couplings[1] = neutron;
}
