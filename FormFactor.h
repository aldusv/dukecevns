#ifndef DUKECEVNS_FORMFACTOR_H
#define DUKECEVNS_FORMFACTOR_H

class FormFactor {
 public:
  virtual ~FormFactor() = default;
  virtual double FFval(double momentumTransferFmInverse) const = 0;

  void SetA(int massNumber);
  void SetRfac(double radiusScale);

 protected:
  int massNumber_ = 0;
  double radiusScale_ = 1.0;
};

class Helm final : public FormFactor {
 public:
  void Setsval(double skinThicknessFm);
  double FFval(double momentumTransferFmInverse) const override;

 private:
  double skinThicknessFm_ = 0.9;
};

#endif
