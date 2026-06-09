#include "FormFactor.h"

#include <cmath>

void FormFactor::SetA(int massNumber)
{
  massNumber_ = massNumber;
}

void FormFactor::SetRfac(double radiusScale)
{
  radiusScale_ = radiusScale;
}

void Helm::Setsval(double skinThicknessFm)
{
  skinThicknessFm_ = skinThicknessFm;
}

double Helm::FFval(double momentumTransferFmInverse) const
{
  const double q = momentumTransferFmInverse * radiusScale_;
  const double radiusFm = 1.2 * std::pow(massNumber_, 1.0 / 3.0);
  const double qr = q * radiusFm;
  const double formFactor =
    3.0 * (std::sin(qr) / (qr * qr) - std::cos(qr) / qr) / qr *
    std::exp(-q * q * skinThicknessFm_ * skinThicknessFm_ / 2.0);
  return std::isnan(formFactor) ? 1.0 : formFactor;
}
