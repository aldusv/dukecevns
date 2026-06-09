#ifndef DUKECEVNS_NUFLUX_H
#define DUKECEVNS_NUFLUX_H

class PiDAR {
 public:
  void SetNorm(double norm);
  double fluxval(double energyMeV, int flavor, double binWidthMeV) const;
  double maxEnu() const;

 private:
  double norm_ = 1.0;
};

#endif
