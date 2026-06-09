#include "NuFlux.h"

#include <cmath>

void PiDAR::SetNorm(double norm)
{
  norm_ = norm;
}

double PiDAR::fluxval(
  double energyMeV,
  int flavor,
  double binWidthMeV) const
{
  constexpr double muonMassMeV = 105.66837;
  constexpr double promptNuMuEnergyMeV = 29.792;
  const double scaledEnergy = 2.0 * energyMeV / muonMassMeV;

  if (energyMeV < 0.0 || energyMeV > muonMassMeV / 2.0) {
    return 0.0;
  }

  if (flavor == 1) {
    return norm_ * 12.0 * std::pow(scaledEnergy, 2) *
           (1.0 - scaledEnergy) * (2.0 / muonMassMeV) * binWidthMeV;
  }
  if (flavor == 2) {
    return std::abs(energyMeV - promptNuMuEnergyMeV) < binWidthMeV / 2.0
      ? norm_
      : 0.0;
  }
  if (flavor == -2) {
    return norm_ * 2.0 * std::pow(scaledEnergy, 2) *
           (3.0 - 2.0 * scaledEnergy) * (2.0 / muonMassMeV) * binWidthMeV;
  }
  return 0.0;
}

double PiDAR::maxEnu() const
{
  return 105.66837 / 2.0;
}
