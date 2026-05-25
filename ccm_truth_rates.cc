#include <nlohmann/json.hpp>

#include "FormFactor.h"
#include "NuFlux.h"
#include "xscns.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

constexpr double kAvogadro = 6.0221409e23;
constexpr double kHbarCMeVFm = 197.327;
constexpr double kAtomicMassUnitMeV = 931.4940955;
constexpr double kElectronMassMeV = 0.51099895;

struct Flavor {
  int id;
  const char* name;
  double fluxScale;
};

struct Isotope {
  std::string name;
  int z;
  int n;
  int a;
  int zSpinDiff;
  int nSpinDiff;
  double molarFraction;
  double massExcessMeV;
  double massMeV;
  double massFraction;
};

struct OutputPaths {
  std::string truthCsv;
  std::string samplerCsv;
  std::string metadataJson;
};

std::string ConfigPath(const std::string& baseName)
{
  if (baseName.find('/') != std::string::npos || baseName.size() >= 5) {
    const std::string suffix = ".json";
    if (baseName.size() >= suffix.size()
        && baseName.compare(baseName.size() - suffix.size(), suffix.size(), suffix) == 0) {
      return baseName;
    }
  }
  return "jsonfiles/" + baseName + ".json";
}

json ReadJson(const std::string& path)
{
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to open config '" + path + "'");
  }
  json config;
  input >> config;
  return config;
}

double RequiredDouble(const json& object, const char* key)
{
  if (!object.contains(key)) {
    throw std::runtime_error(std::string("missing required key '") + key + "'");
  }
  return object.at(key).get<double>();
}

std::unique_ptr<FormFactor> MakeFormFactor(const json& config, int z, int a)
{
  const std::string type = config.value("type", "helm");
  std::unique_ptr<FormFactor> formFactor;
  if (type == "unity") {
    formFactor = std::make_unique<UnityFF>();
  } else if (type == "helm") {
    auto helm = std::make_unique<Helm>();
    helm->Setsval(config.value("sval", 0.9));
    formFactor = std::move(helm);
  } else {
    throw std::runtime_error("unsupported formfactor.type '" + type + "'");
  }

  formFactor->SetZ(z);
  formFactor->SetA(a);
  formFactor->SetRfac(config.value("rfac", 1.0));
  return formFactor;
}

std::vector<Isotope> BuildIsotopes(const std::string& material)
{
  if (material != "CsI") {
    throw std::runtime_error(
      "ccm_truth_rates currently supports material='CsI' only");
  }

  std::vector<Isotope> result = {
    {"Cs133", 55, 78, 133, 1, 0, 0.5, -88.070, 0.0, 0.0},
    {"I127", 53, 74, 127, 1, 0, 0.5, -88.984, 0.0, 0.0},
  };

  double averageMolarMassMeV = 0.0;
  for (auto& isotope : result) {
    isotope.massMeV =
      isotope.a * kAtomicMassUnitMeV - isotope.z * kElectronMassMeV
      + isotope.massExcessMeV;
    averageMolarMassMeV += isotope.massMeV * isotope.molarFraction;
  }
  for (auto& isotope : result) {
    isotope.massFraction =
      isotope.massMeV / averageMolarMassMeV * isotope.molarFraction;
  }
  return result;
}

double TargetsPerTon(const Isotope& isotope)
{
  return 1.0e6 / (isotope.massMeV / kAtomicMassUnitMeV) * kAvogadro;
}

double MaximumRecoilMeV(double maxNeutrinoEnergyMeV, double massMeV)
{
  return 2.0 * maxNeutrinoEnergyMeV * maxNeutrinoEnergyMeV
       / (massMeV + 2.0 * maxNeutrinoEnergyMeV);
}

double MinimumNeutrinoEnergyMeV(double recoilEnergyMeV, double massMeV)
{
  return 0.5 * (recoilEnergyMeV
                + std::sqrt(recoilEnergyMeV * recoilEnergyMeV
                            + 2.0 * massMeV * recoilEnergyMeV));
}

double RecoilRatePerTonSecondMeV(
  const Isotope& isotope,
  FormFactor& formFactor,
  NuFlux& flux,
  const Flavor& flavor,
  double recoilEnergyMeV,
  double neutrinoStepMeV,
  int pdgYear)
{
  double gv[2];
  double ga[2];
  sm_vector_couplings(pdgYear, gv);
  sm_axial_couplings(pdgYear, flavor.id > 0 ? 1 : -1, ga);

  const double qMeV = std::sqrt(2.0 * isotope.massMeV * recoilEnergyMeV);
  const double qFmInv = qMeV / kHbarCMeVFm;
  const double ff = formFactor.FFval(qFmInv);

  const int zdiff = isotope.zSpinDiff;
  const int ndiff = isotope.nSpinDiff;
  const double gvWff = (isotope.z * gv[0] + isotope.n * gv[1]) * ff;
  const double gaWff = (zdiff * ga[0] + ndiff * ga[1]) * ff;

  const double minNu = MinimumNeutrinoEnergyMeV(recoilEnergyMeV, isotope.massMeV);
  const double maxNu = flux.maxEnu();
  double rate = 0.0;
  for (double enu = minNu; enu <= maxNu; enu += neutrinoStepMeV) {
    const double fluxBin = flux.fluxval(enu, flavor.id, neutrinoStepMeV) * flavor.fluxScale;
    rate += std::pow(gvWff, 2) * diffxscnvec(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
    rate += std::pow(gaWff, 2) * diffxscnaxial(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
    rate += gvWff * gaWff * diffxscninterf(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
  }

  return rate * TargetsPerTon(isotope) * isotope.massFraction;
}

std::vector<Flavor> ReadFlavorWeights(const json& config)
{
  const json weights = config.value("flavor_weights", json::object());
  return {
    {1, "nue", weights.value("nue", 1.0)},
    {2, "numu", weights.value("numu", 1.0)},
    {-2, "numubar", weights.value("numubar", 1.0)},
  };
}

OutputPaths ReadOutputPaths(const json& config)
{
  return {
    config.value("output_csv", "out/ccm_truth_recoil_rates.csv"),
    config.value("sampler_csv", "out/ccm_recoil_sampler.csv"),
    config.value("metadata_json", "out/ccm_truth_recoil_rates.meta.json"),
  };
}

} // namespace

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: ./ccm_truth_rates [jsonfile-or-name]\n";
    return 2;
  }

  try {
    const std::string configPath = ConfigPath(argv[1]);
    const json config = ReadJson(configPath);

    const std::string material = config.value("material", "CsI");
    const double detectorMassTons = RequiredDouble(config, "mass_tons");
    const double exposureSeconds = RequiredDouble(config, "exposure_seconds");
    const double fluxPerFlavorCm2Second =
      RequiredDouble(config.at("flux"), "per_flavor_per_cm2_s");
    const double recoilStepMeV = config.value("recoil_step_mev", 0.0001);
    const double recoilMinMeV = config.value("recoil_min_mev", 0.0);
    const double recoilMaxOverrideMeV = config.value("recoil_max_mev", 0.0);
    const double neutrinoStepMeV = config.value("neutrino_step_mev", 0.0005);
    const int pdgYear = config.value("pdg_year", 2024);

    PiDAR flux;
    flux.SetNorm(fluxPerFlavorCm2Second);
    const double maxNu = flux.maxEnu();

    const std::vector<Flavor> flavors = ReadFlavorWeights(config.at("flux"));
    const std::vector<Isotope> isotopeList = BuildIsotopes(material);

    double minMass = isotopeList.front().massMeV;
    for (const auto& isotope : isotopeList) {
      minMass = std::min(minMass, isotope.massMeV);
    }
    const double recoilMaxMeV =
      recoilMaxOverrideMeV > recoilMinMeV
        ? recoilMaxOverrideMeV
        : MaximumRecoilMeV(maxNu, minMass);

    const OutputPaths paths = ReadOutputPaths(config);
    std::ofstream output(paths.truthCsv);
    if (!output) {
      throw std::runtime_error("failed to open output '" + paths.truthCsv + "'");
    }

    std::ofstream sampler(paths.samplerCsv);
    if (!sampler) {
      throw std::runtime_error("failed to open output '" + paths.samplerCsv + "'");
    }

    std::ofstream metadata(paths.metadataJson);
    if (!metadata) {
      throw std::runtime_error("failed to open metadata '" + paths.metadataJson + "'");
    }

    output << "material,isotope,z,a,recoil_energy_mev,recoil_bin_width_mev,"
           << "dnde_events_per_mev,expected_events_bin";
    for (const auto& flavor : flavors) {
      output << ",dnde_" << flavor.name << "_events_per_mev";
    }
    output << '\n';

    sampler << "target_z,target_a,recoil_energy_mev,probability,"
            << "rate_density_events_per_mev,expected_events_bin,"
            << "recoil_bin_width_mev,nu_energy_mev,source\n";

    double totalEvents = 0.0;
    for (const auto& isotope : isotopeList) {
      const auto formFactor = MakeFormFactor(config.value("formfactor", json::object()),
                                             isotope.z,
                                             isotope.a);
      for (double recoil = recoilMinMeV + recoilStepMeV; recoil <= recoilMaxMeV;
           recoil += recoilStepMeV) {
        std::vector<double> flavorRates;
        double ratePerTonSecondMeV = 0.0;
        for (const auto& flavor : flavors) {
          const double rate = RecoilRatePerTonSecondMeV(
            isotope, *formFactor, flux, flavor, recoil, neutrinoStepMeV, pdgYear);
          flavorRates.push_back(rate * detectorMassTons * exposureSeconds);
          ratePerTonSecondMeV += rate;
        }

        const double dnde = ratePerTonSecondMeV * detectorMassTons * exposureSeconds;
        const double binEvents = dnde * recoilStepMeV;
        totalEvents += binEvents;

        output << material << ',' << isotope.name << ',' << isotope.z << ','
               << isotope.a << ',' << std::setprecision(12) << recoil << ','
               << recoilStepMeV << ',' << std::scientific << dnde << ','
               << binEvents;
        for (double flavorRate : flavorRates) {
          output << ',' << flavorRate;
        }
        output << std::defaultfloat << '\n';

        sampler << isotope.z << ',' << isotope.a << ','
                << std::setprecision(12) << recoil << ','
                << std::scientific << binEvents << ','
                << dnde << ','
                << binEvents << ','
                << std::defaultfloat << recoilStepMeV << ','
                << 0.0 << ','
                << isotope.name << '\n';
      }
    }

    json meta = {
      {"config", config},
      {"config_path", configPath},
      {"flux_model", "pi_decay_at_rest"},
      {"outputs", {
        {"truth_csv", paths.truthCsv},
        {"sampler_csv", paths.samplerCsv},
        {"metadata_json", paths.metadataJson},
      }},
      {"output_units", {
        {"recoil_energy", "MeVnr"},
        {"dnde_events_per_mev", "expected events / MeVnr over configured exposure"},
        {"expected_events_bin", "expected events in recoil bin over configured exposure"},
        {"sampler_probability", "expected events in recoil bin; current CCI sampler normalizes this"},
        {"rate_density_events_per_mev", "expected events / MeVnr over configured exposure"},
      }},
      {"total_expected_events", totalEvents},
    };
    metadata << std::setw(2) << meta << '\n';

    std::cout << "Wrote " << paths.truthCsv << '\n';
    std::cout << "Wrote " << paths.samplerCsv << '\n';
    std::cout << "Wrote " << paths.metadataJson << '\n';
    std::cout << "Total expected events: " << std::setprecision(8)
              << totalEvents << '\n';
  } catch (const std::exception& error) {
    std::cerr << "ccm_truth_rates: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
