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
constexpr double kPi = 3.14159265358979323846;
constexpr double kPromptNuMuEnergyMeV = 29.792;
constexpr double kCentimetersPerInch = 2.54;

struct SourceComponent {
  std::string component;
  int id;
  std::string flavor;
  double yieldPerStoppedPiPlus;
  std::string timeProfile;
  double decayLifetimeNs;
  std::string beamPulseShape;
  double beamPulseWidthNs;
  double beamPulsePeakNs;
  double beamPeriodNs;
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
  std::string recoilByComponentCsv;
  std::string timeByComponentCsv;
  std::string recoilTimeCsv;
  std::string nuRecoilCsv;
};

struct TimeWindow {
  double startNs;
  double endNs;
  std::string mode;
};

struct Normalization {
  std::string angularModel;
  double baselineM;
  double sphereAreaCm2;
  double potPerSecond;
  double stoppedPiPlusPerPot;
  double baseFluxCm2Second;
  double referenceTotalFluxCm2Second;
};

struct RateLedger {
  std::map<std::string, double> expectedEventsByComponent;
  std::map<std::string, double> expectedEventsByIsotope;
  std::map<std::string, double> timeAcceptanceByComponent;
  double totalUnwindowedExpectedEvents = 0.0;
  double totalExpectedEvents = 0.0;
};

using RecoilComponentEvents =
  std::map<double, std::map<std::string, double>>;
using NuRecoilEvents =
  std::map<double, std::map<std::string, std::map<double, double>>>;

struct OutputFiles {
  std::ofstream truthRateAudit;
  std::ofstream cciSampler;
  std::ofstream metadata;
  std::ofstream recoilByComponentDiagnostic;
  std::ofstream timeByComponentDiagnostic;
  std::ofstream recoilTimeDiagnostic;
  std::ofstream nuRecoilDiagnostic;
};

double RequiredDouble(const json& object, const char* key);

double ReadDetectorMassTons(const json& config)
{
  const json& detector = config.at("detector");
  const json& dimensions = detector.at("crystal_dimensions_in");
  if (!dimensions.is_array() || dimensions.size() != 3) {
    throw std::runtime_error(
      "detector.crystal_dimensions_in must contain exactly three values");
  }

  double crystalVolumeCm3 = 1.0;
  for (const json& dimension : dimensions) {
    crystalVolumeCm3 *= dimension.get<double>() * kCentimetersPerInch;
  }
  const double crystalMassKg =
    crystalVolumeCm3 * RequiredDouble(detector, "density_g_cm3") / 1000.0;
  return RequiredDouble(detector, "crystal_count") * crystalMassKg / 1000.0;
}

struct NuEnergyContribution {
  double nuEnergyMeV = 0.0;
  double ratePerTonSecondMeV = 0.0;
};

struct RateBreakdown {
  double ratePerTonSecondMeV = 0.0;
  std::vector<NuEnergyContribution> nuEnergyContributions;
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
  if (!object.at(key).is_number()) {
    throw std::runtime_error(std::string("key '") + key + "' must be numeric");
  }
  return object.at(key).get<double>();
}

json BookkeepingEntry(const json& config, const std::string& name)
{
  const json bookkeeping = config.value("bookkeeping", json::array());
  if (!bookkeeping.is_array()) {
    throw std::runtime_error("bookkeeping must be an array of named entries");
  }
  for (const json& entry : bookkeeping) {
    if (entry.value("name", "") == name) {
      return entry;
    }
  }
  return json::object();
}

TimeWindow ReadTimeWindow(const json& config)
{
  const json window = config.value("time_window", json::object());
  return {
    window.value("start_ns", 0.0),
    window.value("end_ns", 0.0),
    window.value("mode", "annotate_only"),
  };
}

Normalization ReadNormalization(const json& config)
{
  const json& settings = config.at("normalization");
  const json& beam = config.at("beam");
  const json referenceFlux = BookkeepingEntry(config, "reference_flux");
  if (settings.value("mode", "") != "source_yield") {
    throw std::runtime_error(
      "normalization.mode must be 'source_yield'");
  }

  const double baselineM = RequiredDouble(beam, "baseline_m");
  const double sphereAreaCm2 =
    baselineM > 0.0 ? 4.0 * kPi * std::pow(baselineM * 100.0, 2) : 0.0;
  const double potPerSecond = RequiredDouble(settings, "pot_per_second");
  const double referenceTotalFluxCm2Second =
    referenceFlux.value("total_per_cm2_s", 0.0);
  const json& sourceYield = settings.at("source_yield");
  const std::string angularModel =
    sourceYield.value("angular_model", "isotropic");
  if (angularModel != "isotropic") {
    throw std::runtime_error(
      "normalization.source_yield.angular_model must be 'isotropic'");
  }
  if (sphereAreaCm2 <= 0.0 || potPerSecond <= 0.0) {
    throw std::runtime_error(
      "baseline_m and pot_per_second must be positive");
  }
  const double stoppedPiPlusPerPot =
    RequiredDouble(sourceYield, "stopped_piplus_per_pot");
  const double baseFluxCm2Second =
    potPerSecond * stoppedPiPlusPerPot / sphereAreaCm2;
  return {
    angularModel,
    baselineM,
    sphereAreaCm2,
    potPerSecond,
    stoppedPiPlusPerPot,
    baseFluxCm2Second,
    referenceTotalFluxCm2Second,
  };
}

std::unique_ptr<FormFactor> MakeFormFactor(const json& config, int a)
{
  if (config.value("type", "") != "helm") {
    throw std::runtime_error("formfactor.type must be 'helm'");
  }

  auto formFactor = std::make_unique<Helm>();
  formFactor->Setsval(config.value("sval", 0.9));
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

double BeamPulseCdf(double timeNs, const SourceComponent& component)
{
  if (timeNs <= 0.0) {
    return 0.0;
  }
  if (component.beamPulseWidthNs <= 0.0) {
    return 1.0;
  }
  if (timeNs >= component.beamPulseWidthNs) {
    return 1.0;
  }

  if (component.beamPulseShape == "triangular") {
    const double width = component.beamPulseWidthNs;
    const double peak = std::clamp(component.beamPulsePeakNs, 0.0, width);
    if (peak <= 0.0) {
      return 1.0 - std::pow((width - timeNs) / width, 2);
    }
    if (peak >= width) {
      return std::pow(timeNs / width, 2);
    }
    if (timeNs <= peak) {
      return timeNs * timeNs / (width * peak);
    }
    return 1.0 - std::pow(width - timeNs, 2) / (width * (width - peak));
  }

  return std::clamp(timeNs / component.beamPulseWidthNs, 0.0, 1.0);
}

double BeamPulsePdf(double timeNs, const SourceComponent& component)
{
  if (component.beamPulseWidthNs <= 0.0 || timeNs < 0.0 ||
      timeNs > component.beamPulseWidthNs) {
    return 0.0;
  }

  if (component.beamPulseShape == "triangular") {
    const double width = component.beamPulseWidthNs;
    const double peak = std::clamp(component.beamPulsePeakNs, 0.0, width);
    if (peak <= 0.0) {
      return 2.0 * (width - timeNs) / (width * width);
    }
    if (peak >= width) {
      return 2.0 * timeNs / (width * width);
    }
    if (timeNs <= peak) {
      return 2.0 * timeNs / (width * peak);
    }
    return 2.0 * (width - timeNs) / (width * (width - peak));
  }

  return 1.0 / component.beamPulseWidthNs;
}

double DecayCdf(double timeNs, double decayLifetimeNs)
{
  if (timeNs <= 0.0) {
    return 0.0;
  }
  if (decayLifetimeNs <= 0.0) {
    return 1.0;
  }
  return 1.0 - std::exp(-timeNs / decayLifetimeNs);
}

double TimeCdf(double timeNs, const SourceComponent& component)
{
  if (timeNs <= 0.0) {
    return 0.0;
  }
  if (component.decayLifetimeNs <= 0.0) {
    return BeamPulseCdf(timeNs, component);
  }
  if (component.beamPulseWidthNs <= 0.0) {
    return DecayCdf(timeNs, component.decayLifetimeNs);
  }

  const int steps = 512;
  const double width = component.beamPulseWidthNs;
  const double step = width / static_cast<double>(steps);
  double cdf = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double pulseTime = (static_cast<double>(i) + 0.5) * step;
    cdf += BeamPulsePdf(pulseTime, component) *
           DecayCdf(timeNs - pulseTime, component.decayLifetimeNs) *
           step;
  }
  return std::clamp(cdf, 0.0, 1.0);
}

double TimeAcceptance(const SourceComponent& component, const TimeWindow& window)
{
  if (window.endNs <= window.startNs) {
    return 1.0;
  }
  const double startCdf =
    TimeCdf(window.startNs, component);
  const double endCdf =
    TimeCdf(window.endNs, component);
  return std::clamp(endCdf - startCdf, 0.0, 1.0);
}

RateBreakdown RecoilRateBreakdownPerTonSecondMeV(
  const Isotope& isotope,
  FormFactor& formFactor,
  PiDAR& flux,
  const SourceComponent& component,
  double recoilEnergyMeV,
  double neutrinoStepMeV,
  int pdgYear)
{
  double gv[2];
  double ga[2];
  sm_vector_couplings(pdgYear, gv);
  sm_axial_couplings(pdgYear, component.id > 0 ? 1 : -1, ga);

  const double qMeV = std::sqrt(2.0 * isotope.massMeV * recoilEnergyMeV);
  const double qFmInv = qMeV / kHbarCMeVFm;
  const double ff = formFactor.FFval(qFmInv);

  const int zdiff = isotope.zSpinDiff;
  const int ndiff = isotope.nSpinDiff;
  const double gvWff = (isotope.z * gv[0] + isotope.n * gv[1]) * ff;
  const double gaWff = (zdiff * ga[0] + ndiff * ga[1]) * ff;

  const double minNu = MinimumNeutrinoEnergyMeV(recoilEnergyMeV, isotope.massMeV);
  const double maxNu = flux.maxEnu();
  RateBreakdown result;
  const double targetScale = TargetsPerTon(isotope) * isotope.massFraction;
  for (double enu = minNu; enu <= maxNu; enu += neutrinoStepMeV) {
    const double fluxBin =
      flux.fluxval(enu, component.id, neutrinoStepMeV) *
      component.yieldPerStoppedPiPlus;
    double rateContribution = 0.0;
    rateContribution +=
      std::pow(gvWff, 2) *
      diffxscnvec(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
    rateContribution +=
      std::pow(gaWff, 2) *
      diffxscnaxial(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
    rateContribution +=
      gvWff * gaWff *
      diffxscninterf(enu, isotope.massMeV, recoilEnergyMeV) * fluxBin;
    rateContribution *= targetScale;
    result.ratePerTonSecondMeV += rateContribution;
    if (rateContribution > 0.0) {
      result.nuEnergyContributions.push_back({enu, rateContribution});
    }
  }

  return result;
}

int FlavorId(const std::string& flavor)
{
  if (flavor == "nue") {
    return 1;
  }
  if (flavor == "numu") {
    return 2;
  }
  if (flavor == "numubar") {
    return -2;
  }
  throw std::runtime_error("unsupported neutrino flavor '" + flavor + "'");
}

std::vector<SourceComponent> ReadSourceComponents(
  const json& config,
  double beamPulseWidthNs,
  const std::string& beamPulseShape,
  double beamPulsePeakNs,
  double beamPeriodNs)
{
  if (!config.contains("source_components")) {
    throw std::runtime_error("flux.source_components is required");
  }

  std::vector<SourceComponent> components;
  for (const auto& raw : config.at("source_components")) {
    const std::string flavor = raw.at("flavor").get<std::string>();
    components.push_back({
      raw.at("name").get<std::string>(),
      FlavorId(flavor),
      flavor,
      raw.value("yield_per_stopped_piplus", 1.0),
      raw.at("time_profile").get<std::string>(),
      raw.value("decay_lifetime_ns", 0.0),
      raw.value("beam_pulse_shape", beamPulseShape),
      beamPulseWidthNs,
      raw.value("beam_pulse_peak_ns", beamPulsePeakNs),
      beamPeriodNs,
    });
  }
  return components;
}

OutputPaths ReadOutputPaths(const json& config)
{
  const json diagnostics = config.value("diagnostics", json::object());
  return {
    config.value("output_csv", "out/ccm_truth_recoil_rates.csv"),
    config.value("sampler_csv", "out/ccm_recoil_sampler.csv"),
    config.value("metadata_json", "out/ccm_truth_recoil_rates.meta.json"),
    diagnostics.value("recoil_by_component_csv", "out/ccm_recoil_by_component.csv"),
    diagnostics.value("time_by_component_csv", "out/ccm_time_by_component.csv"),
    diagnostics.value("recoil_time_csv", "out/ccm_recoil_time_2d.csv"),
    diagnostics.value("nu_recoil_csv", "out/ccm_nu_recoil_2d.csv"),
  };
}

std::ofstream OpenOutputFile(const std::string& path, const std::string& label)
{
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("failed to open " + label + " '" + path + "'");
  }
  return output;
}

OutputFiles OpenOutputFiles(const OutputPaths& paths)
{
  return {
    OpenOutputFile(paths.truthCsv, "output"),
    OpenOutputFile(paths.samplerCsv, "output"),
    OpenOutputFile(paths.metadataJson, "metadata"),
    OpenOutputFile(paths.recoilByComponentCsv, "output"),
    OpenOutputFile(paths.timeByComponentCsv, "output"),
    OpenOutputFile(paths.recoilTimeCsv, "output"),
    OpenOutputFile(paths.nuRecoilCsv, "output"),
  };
}

void WriteTruthRateHeader(std::ostream& output)
{
  output << "material,isotope,z,a,recoil_energy_mev,recoil_bin_width_mev,"
         << "source_component,neutrino_flavor,time_profile,decay_lifetime_ns,"
         << "beam_pulse_shape,beam_pulse_width_ns,beam_pulse_peak_ns,"
         << "beam_period_ns,time_window_start_ns,time_window_end_ns,"
         << "time_window_mode,time_acceptance,exposure_seconds,"
         << "dnde_events_per_mev,"
         << "expected_events_bin\n";
}

void WriteCciSamplerHeader(std::ostream& sampler)
{
  sampler << "target_z,target_a,recoil_energy_mev,probability,"
          << "rate_density_events_per_mev,expected_events_bin,"
          << "recoil_bin_width_mev,nu_energy_mev,source,source_component,"
          << "neutrino_flavor,time_profile,decay_lifetime_ns,"
          << "beam_pulse_shape,beam_pulse_width_ns,beam_pulse_peak_ns,"
          << "beam_period_ns,time_window_start_ns,time_window_end_ns,"
          << "time_window_mode,time_acceptance,exposure_seconds\n";
}

void WriteNuRecoilDiagnostic(
  std::ostream& nuRecoil,
  const NuRecoilEvents& nuRecoilEventsByComponent,
  const std::vector<SourceComponent>& components,
  double recoilStepMeV,
  double nuDiagnosticStepMeV)
{
  nuRecoil << "nu_energy_mev,nu_energy_bin_width_mev,recoil_energy_mev,"
           << "recoil_bin_width_mev,source_component,neutrino_flavor,"
           << "expected_events_bin\n";
  for (const auto& recoilBin : nuRecoilEventsByComponent) {
    for (const auto& component : components) {
      const auto componentIt = recoilBin.second.find(component.component);
      if (componentIt == recoilBin.second.end()) {
        continue;
      }
      for (const auto& nuBin : componentIt->second) {
        nuRecoil << std::setprecision(12)
                 << nuBin.first << ',' << nuDiagnosticStepMeV << ','
                 << recoilBin.first << ',' << recoilStepMeV << ','
                 << component.component << ',' << component.flavor << ','
                 << std::scientific << nuBin.second << std::defaultfloat << '\n';
      }
    }
  }
}

void WriteRecoilByComponentDiagnostic(
  std::ostream& recoilByComponent,
  const RecoilComponentEvents& recoilEventsByComponent,
  const std::vector<SourceComponent>& components,
  double recoilStepMeV)
{
  recoilByComponent << "recoil_energy_mev,recoil_bin_width_mev";
  for (const auto& component : components) {
    recoilByComponent << ',' << component.component << "_events";
  }
  recoilByComponent << ",total_events\n";

  for (const auto& recoilBin : recoilEventsByComponent) {
    double total = 0.0;
    recoilByComponent << std::setprecision(12) << recoilBin.first << ','
                      << recoilStepMeV;
    for (const auto& component : components) {
      const double events =
        recoilBin.second.count(component.component)
          ? recoilBin.second.at(component.component)
          : 0.0;
      total += events;
      recoilByComponent << ',' << std::scientific << events << std::defaultfloat;
    }
    recoilByComponent << ',' << std::scientific << total << std::defaultfloat << '\n';
  }
}

void WriteTimeDiagnostics(
  std::ostream& timeByComponent,
  std::ostream& recoilTime,
  const RateLedger& ledger,
  const RecoilComponentEvents& recoilEventsByComponent,
  const std::vector<SourceComponent>& components,
  double recoilStepMeV,
  double timeMinNs,
  double timeMaxNs,
  double timeStepNs)
{
  std::map<std::string, std::vector<double>> timeProbByComponent;
  timeByComponent << "time_start_ns,time_end_ns,time_center_ns,time_bin_width_ns,"
                  << "source_component,neutrino_flavor,probability,"
                  << "expected_events_bin\n";
  recoilTime << "recoil_energy_mev,recoil_bin_width_mev,time_start_ns,"
             << "time_end_ns,time_center_ns,time_bin_width_ns,"
             << "source_component,neutrino_flavor,expected_events_bin\n";

  for (const auto& component : components) {
    std::vector<double> probabilities;
    for (double timeStart = timeMinNs; timeStart < timeMaxNs; timeStart += timeStepNs) {
      const double timeEnd = std::min(timeStart + timeStepNs, timeMaxNs);
      const double probability =
        std::clamp(TimeCdf(timeEnd, component) - TimeCdf(timeStart, component),
                   0.0, 1.0);
      probabilities.push_back(probability);
      const double expectedEvents =
        ledger.expectedEventsByComponent.at(component.component) * probability;
      timeByComponent << std::setprecision(12)
                      << timeStart << ',' << timeEnd << ','
                      << 0.5 * (timeStart + timeEnd) << ','
                      << timeEnd - timeStart << ','
                      << component.component << ',' << component.flavor << ','
                      << std::scientific << probability << ','
                      << expectedEvents << std::defaultfloat << '\n';
    }
    timeProbByComponent[component.component] = probabilities;
  }

  for (const auto& recoilBin : recoilEventsByComponent) {
    for (const auto& component : components) {
      const double recoilComponentEvents =
        recoilBin.second.count(component.component)
          ? recoilBin.second.at(component.component)
          : 0.0;
      const std::vector<double>& probabilities =
        timeProbByComponent.at(component.component);
      std::size_t timeIndex = 0;
      for (double timeStart = timeMinNs; timeStart < timeMaxNs; timeStart += timeStepNs) {
        const double timeEnd = std::min(timeStart + timeStepNs, timeMaxNs);
        const double expectedEvents =
          recoilComponentEvents * probabilities.at(timeIndex++);
        recoilTime << std::setprecision(12) << recoilBin.first << ','
                   << recoilStepMeV << ',' << timeStart << ',' << timeEnd
                   << ',' << 0.5 * (timeStart + timeEnd) << ','
                   << timeEnd - timeStart << ','
                   << component.component << ',' << component.flavor << ','
                   << std::scientific << expectedEvents << std::defaultfloat << '\n';
      }
    }
  }
}

json BuildMetadata(
  const json& config,
  const std::string& configPath,
  const std::string& material,
  const json& detectorConfig,
  const OutputPaths& paths,
  const Normalization& normalization,
  const RateLedger& ledger,
  const std::vector<Isotope>& isotopeList,
  const std::vector<SourceComponent>& components,
  double detectorMassTons,
  double exposureSeconds,
  const std::string& beamPulseShape,
  double beamPulseWidthNs,
  double beamPulsePeakNs,
  double pulseRateHz,
  double beamPeriodNs,
  double recoilMinMeV,
  double recoilMaxMeV,
  double recoilStepMeV,
  double neutrinoStepMeV,
  double maxNu,
  int pdgYear,
  const TimeWindow& timeWindow,
  double timeMinNs,
  double timeMaxNs,
  double timeStepNs,
  double nuDiagnosticStepMeV)
{
  json isotopeLedger = json::array();
  for (const auto& isotope : isotopeList) {
    isotopeLedger.push_back({
      {"name", isotope.name},
      {"z", isotope.z},
      {"a", isotope.a},
      {"molar_fraction", isotope.molarFraction},
      {"mass_fraction", isotope.massFraction},
      {"targets_per_ton", TargetsPerTon(isotope) * isotope.massFraction},
      {"targets_in_detector", TargetsPerTon(isotope) * isotope.massFraction * detectorMassTons},
      {"expected_events", ledger.expectedEventsByIsotope.at(isotope.name)},
    });
  }

  json componentLedger = json::array();
  double configuredFluxCm2Second = 0.0;
  for (const auto& component : components) {
    const double componentFlux =
      normalization.baseFluxCm2Second * component.yieldPerStoppedPiPlus;
    configuredFluxCm2Second += componentFlux;
    componentLedger.push_back({
      {"name", component.component},
      {"flavor", component.flavor},
      {"flavor_id", component.id},
      {"yield_per_stopped_piplus", component.yieldPerStoppedPiPlus},
      {"flux_per_cm2_s", componentFlux},
      {"time_profile", component.timeProfile},
      {"decay_lifetime_ns", component.decayLifetimeNs},
      {"beam_pulse_shape", component.beamPulseShape},
      {"beam_pulse_width_ns", component.beamPulseWidthNs},
      {"beam_pulse_peak_ns", component.beamPulsePeakNs},
      {"time_acceptance", ledger.timeAcceptanceByComponent.at(component.component)},
      {"expected_events", ledger.expectedEventsByComponent.at(component.component)},
    });
  }

  const double inferredNeutrinosPerSecond =
    normalization.sphereAreaCm2 > 0.0
      ? configuredFluxCm2Second * normalization.sphereAreaCm2
      : 0.0;
  const double inferredNeutrinosPerPot =
    normalization.potPerSecond > 0.0
      ? inferredNeutrinosPerSecond / normalization.potPerSecond
      : 0.0;
  const double referenceToConfiguredFluxRatio =
    normalization.referenceTotalFluxCm2Second > 0.0
      ? configuredFluxCm2Second / normalization.referenceTotalFluxCm2Second
      : 0.0;

  return {
    {"config", config},
    {"config_path", configPath},
    {"flux_model", "pi_decay_at_rest"},
    {"normalization", {
      {"mode", "source_yield"},
      {"angular_model", normalization.angularModel},
      {"base_flux_cm2_s", normalization.baseFluxCm2Second},
      {"total_configured_flux_cm2_s", configuredFluxCm2Second},
      {"reference_total_flux_cm2_s", normalization.referenceTotalFluxCm2Second},
      {"configured_to_reference_flux_ratio", referenceToConfiguredFluxRatio},
      {"baseline_m", normalization.baselineM},
      {"sphere_area_cm2_at_baseline", normalization.sphereAreaCm2},
      {"exposure_seconds", exposureSeconds},
      {"mass_tons", detectorMassTons},
      {"pot_per_second", normalization.potPerSecond},
      {"stopped_piplus_per_pot", normalization.stoppedPiPlusPerPot},
      {"inferred_total_neutrinos_per_second_4pi", inferredNeutrinosPerSecond},
      {"inferred_total_neutrinos_per_pot_4pi", inferredNeutrinosPerPot},
      {"note", "Rates are normalized from POT, stopped pi+ yield, component yields, and isotropic geometry; reference flux is an audit quantity."},
    }},
    {"beam", {
      {"pulse_shape", beamPulseShape},
      {"pulse_width_ns", beamPulseWidthNs},
      {"pulse_peak_ns", beamPulsePeakNs},
      {"pulse_rate_hz", pulseRateHz},
      {"period_ns", beamPeriodNs},
    }},
    {"target", {
      {"material", material},
      {"mass_tons", detectorMassTons},
      {"detector", detectorConfig},
      {"isotopes", isotopeLedger},
    }},
    {"integration", {
      {"recoil_min_mev", recoilMinMeV},
      {"recoil_max_mev", recoilMaxMeV},
      {"recoil_step_mev", recoilStepMeV},
      {"neutrino_step_mev", neutrinoStepMeV},
      {"max_neutrino_energy_mev", maxNu},
      {"pdg_year", pdgYear},
    }},
    {"source_components", componentLedger},
    {"time_window", {
      {"start_ns", timeWindow.startNs},
      {"end_ns", timeWindow.endNs},
      {"mode", timeWindow.mode},
    }},
    {"diagnostics", {
      {"time_min_ns", timeMinNs},
      {"time_max_ns", timeMaxNs},
      {"time_step_ns", timeStepNs},
      {"nu_energy_step_mev", nuDiagnosticStepMeV},
    }},
    {"outputs", {
      {"truth_csv", paths.truthCsv},
      {"sampler_csv", paths.samplerCsv},
      {"metadata_json", paths.metadataJson},
      {"recoil_by_component_csv", paths.recoilByComponentCsv},
      {"time_by_component_csv", paths.timeByComponentCsv},
      {"recoil_time_csv", paths.recoilTimeCsv},
      {"nu_recoil_csv", paths.nuRecoilCsv},
    }},
    {"output_units", {
      {"recoil_energy", "MeVnr"},
      {"dnde_events_per_mev", "expected events / MeVnr over configured exposure"},
      {"expected_events_bin", "expected events in recoil bin over configured exposure"},
      {"sampler_probability", "expected events in recoil bin; current CCI sampler normalizes this"},
      {"rate_density_events_per_mev", "expected events / MeVnr over configured exposure"},
    }},
    {"rate_ledger", {
      {"total_unwindowed_expected_events", ledger.totalUnwindowedExpectedEvents},
      {"total_expected_events", ledger.totalExpectedEvents},
      {"expected_events_by_component", ledger.expectedEventsByComponent},
      {"expected_events_by_isotope", ledger.expectedEventsByIsotope},
      {"time_acceptance_by_component", ledger.timeAcceptanceByComponent},
    }},
    {"total_expected_events", ledger.totalExpectedEvents},
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
    const json detectorConfig = config.value("detector", json::object());
    const double detectorMassTons = ReadDetectorMassTons(config);
    const double exposureSeconds = RequiredDouble(config, "exposure_seconds");
    const double beamPulseWidthNs = RequiredDouble(config.at("beam"), "pulse_width_ns");
    const std::string beamPulseShape =
      config.at("beam").value("pulse_shape", "uniform");
    const double beamPulsePeakNs =
      config.at("beam").value("pulse_peak_ns", 0.5 * beamPulseWidthNs);
    const double pulseRateHz = RequiredDouble(config.at("beam"), "pulse_rate_hz");
    const double beamPeriodNs = 1.0e9 / pulseRateHz;
    const double recoilStepMeV = config.value("recoil_step_mev", 0.0001);
    const double recoilMinMeV = config.value("recoil_min_mev", 0.0);
    const double recoilMaxOverrideMeV = config.value("recoil_max_mev", 0.0);
    const double neutrinoStepMeV = config.value("neutrino_step_mev", 0.0005);
    const int pdgYear = config.value("pdg_year", 2024);
    const TimeWindow timeWindow = ReadTimeWindow(config);
    const json diagnostics = config.value("diagnostics", json::object());
    const double timeMinNs = diagnostics.value("time_min_ns", 0.0);
    const double timeMaxNs = diagnostics.value("time_max_ns", 10000.0);
    const double timeStepNs = diagnostics.value("time_step_ns", 100.0);
    const double nuDiagnosticStepMeV =
      diagnostics.value("nu_energy_step_mev", 0.25);
    if (timeMaxNs <= timeMinNs || timeStepNs <= 0.0) {
      throw std::runtime_error(
        "diagnostics requires time_max_ns > time_min_ns and time_step_ns > 0");
    }
    if (nuDiagnosticStepMeV <= 0.0) {
      throw std::runtime_error("diagnostics.nu_energy_step_mev must be > 0");
    }

    const std::vector<SourceComponent> components =
      ReadSourceComponents(config.at("flux"), beamPulseWidthNs, beamPulseShape,
        beamPulsePeakNs, beamPeriodNs);
    const Normalization normalization = ReadNormalization(config);

    PiDAR flux;
    flux.SetNorm(normalization.baseFluxCm2Second);
    const double maxNu = flux.maxEnu();
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
    OutputFiles files = OpenOutputFiles(paths);
    WriteTruthRateHeader(files.truthRateAudit);
    WriteCciSamplerHeader(files.cciSampler);

    RateLedger ledger;
    RecoilComponentEvents recoilEventsByComponent;
    NuRecoilEvents nuRecoilEventsByComponent;
    for (const auto& isotope : isotopeList) {
      const auto formFactor =
        MakeFormFactor(config.value("formfactor", json::object()), isotope.a);
      for (double recoil = recoilMinMeV + recoilStepMeV; recoil <= recoilMaxMeV;
           recoil += recoilStepMeV) {
        for (const auto& component : components) {
          const RateBreakdown rateBreakdown = RecoilRateBreakdownPerTonSecondMeV(
            isotope, *formFactor, flux, component, recoil, neutrinoStepMeV, pdgYear);
          const double ratePerTonSecondMeV = rateBreakdown.ratePerTonSecondMeV;
          const double timeAcceptance = TimeAcceptance(component, timeWindow);
          const double windowScale =
            timeWindow.mode == "rate_windowed" ? timeAcceptance : 1.0;
          const double unwindowedDnde =
            ratePerTonSecondMeV * detectorMassTons * exposureSeconds;
          const double dnde =
            ratePerTonSecondMeV * detectorMassTons * exposureSeconds *
            windowScale;
          ledger.totalUnwindowedExpectedEvents += unwindowedDnde * recoilStepMeV;
          const double binEvents = dnde * recoilStepMeV;
          ledger.totalExpectedEvents += binEvents;
          ledger.expectedEventsByComponent[component.component] += binEvents;
          ledger.expectedEventsByIsotope[isotope.name] += binEvents;
          ledger.timeAcceptanceByComponent[component.component] = timeAcceptance;
          recoilEventsByComponent[recoil][component.component] += binEvents;

          for (const auto& contribution : rateBreakdown.nuEnergyContributions) {
            const double contributionEvents =
              contribution.ratePerTonSecondMeV * detectorMassTons *
              exposureSeconds * windowScale * recoilStepMeV;
            const double nuBinCenter =
              (std::floor(contribution.nuEnergyMeV / nuDiagnosticStepMeV) + 0.5) *
              nuDiagnosticStepMeV;
            nuRecoilEventsByComponent[recoil][component.component][nuBinCenter] +=
              contributionEvents;
          }

          files.truthRateAudit
            << material << ',' << isotope.name << ',' << isotope.z << ','
            << isotope.a << ',' << std::setprecision(12) << recoil << ','
            << recoilStepMeV << ',' << component.component << ','
            << component.flavor << ',' << component.timeProfile << ','
            << component.decayLifetimeNs << ',' << component.beamPulseShape
            << ',' << component.beamPulseWidthNs << ','
            << component.beamPulsePeakNs << ',' << component.beamPeriodNs << ','
            << timeWindow.startNs << ',' << timeWindow.endNs << ','
            << timeWindow.mode << ',' << timeAcceptance << ','
            << exposureSeconds << ','
            << std::scientific
            << dnde << ',' << binEvents << std::defaultfloat << '\n';

          std::map<double, std::pair<double, double>> samplerNuBins;
          for (const auto& contribution : rateBreakdown.nuEnergyContributions) {
            const double nuBin =
              std::floor(contribution.nuEnergyMeV / nuDiagnosticStepMeV);
            auto& samplerBin = samplerNuBins[nuBin];
            samplerBin.first += contribution.ratePerTonSecondMeV;
            samplerBin.second +=
              contribution.nuEnergyMeV * contribution.ratePerTonSecondMeV;
          }
          for (const auto& nuBin : samplerNuBins) {
            const double contributionRate = nuBin.second.first;
            if (contributionRate <= 0.0) {
              continue;
            }
            const double contributionEvents =
              contributionRate * detectorMassTons * exposureSeconds *
              windowScale * recoilStepMeV;
            const double meanNuEnergyMeV =
              component.flavor == "numu"
                ? kPromptNuMuEnergyMeV
                : nuBin.second.second / contributionRate;
            files.cciSampler << isotope.z << ',' << isotope.a << ','
                              << std::setprecision(12) << recoil << ','
                              << std::scientific << contributionEvents << ','
                              << contributionEvents / recoilStepMeV << ','
                              << contributionEvents << ','
                              << std::defaultfloat << recoilStepMeV << ','
                              << meanNuEnergyMeV << ','
                              << isotope.name << ','
                              << component.component << ','
                              << component.flavor << ','
                              << component.timeProfile << ','
                              << component.decayLifetimeNs << ','
                              << component.beamPulseShape << ','
                              << component.beamPulseWidthNs << ','
                              << component.beamPulsePeakNs << ','
                              << component.beamPeriodNs << ','
                              << timeWindow.startNs << ','
                              << timeWindow.endNs << ','
                              << timeWindow.mode << ','
                              << timeAcceptance << ','
                              << exposureSeconds << '\n';
          }
        }
      }
    }

    WriteNuRecoilDiagnostic(files.nuRecoilDiagnostic, nuRecoilEventsByComponent,
      components, recoilStepMeV, nuDiagnosticStepMeV);
    WriteRecoilByComponentDiagnostic(files.recoilByComponentDiagnostic,
      recoilEventsByComponent, components, recoilStepMeV);
    WriteTimeDiagnostics(files.timeByComponentDiagnostic,
      files.recoilTimeDiagnostic, ledger, recoilEventsByComponent, components,
      recoilStepMeV, timeMinNs, timeMaxNs, timeStepNs);

    const json meta = BuildMetadata(config, configPath, material, detectorConfig,
      paths, normalization, ledger, isotopeList, components, detectorMassTons,
      exposureSeconds, beamPulseShape, beamPulseWidthNs, beamPulsePeakNs,
      pulseRateHz, beamPeriodNs, recoilMinMeV, recoilMaxMeV, recoilStepMeV,
      neutrinoStepMeV, maxNu, pdgYear, timeWindow, timeMinNs, timeMaxNs,
      timeStepNs, nuDiagnosticStepMeV);
    files.metadata << std::setw(2) << meta << '\n';

    std::cout << "Wrote " << paths.truthCsv << '\n';
    std::cout << "Wrote " << paths.samplerCsv << '\n';
    std::cout << "Wrote " << paths.metadataJson << '\n';
    std::cout << "Wrote " << paths.recoilByComponentCsv << '\n';
    std::cout << "Wrote " << paths.timeByComponentCsv << '\n';
    std::cout << "Wrote " << paths.recoilTimeCsv << '\n';
    std::cout << "Wrote " << paths.nuRecoilCsv << '\n';
    std::cout << "Total expected events: " << std::setprecision(8)
              << ledger.totalExpectedEvents << '\n';
  } catch (const std::exception& error) {
    std::cerr << "ccm_truth_rates: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
