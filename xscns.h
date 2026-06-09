#ifndef DUKECEVNS_XSCNS_H
#define DUKECEVNS_XSCNS_H

double diffxscnvec(double neutrinoEnergyMeV, double massMeV, double recoilMeV);
double diffxscnaxial(double neutrinoEnergyMeV, double massMeV, double recoilMeV);
double diffxscninterf(double neutrinoEnergyMeV, double massMeV, double recoilMeV);

void sm_vector_couplings(int pdgYear, double* couplings);
void sm_axial_couplings(int pdgYear, int flavorSign, double* couplings);

#endif
