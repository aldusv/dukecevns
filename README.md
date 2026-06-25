CCM CsI CEvNS truth-rate and recoil-sampling tools for CCI.

## Run

```sh
make ccm_truth_rates
./ccm_truth_rates ccm_csi_truth
uv run python plot_ccm_diagnostics.py
```

The calculation is configured by `jsonfiles/ccm_csi_truth.json`.

## Calculation

Rates are normalized from the source yield:

```text
stopped pi+/s = POT/s * stopped pi+/POT
component flux = stopped pi+/s * neutrinos/stopped pi+ / (4 pi R^2)
dN/dT = exposure * target count * integral[dPhi/dEnu * dSigma/dT dEnu]
```

The neutrino-energy and recoil-energy integrals use fixed-width rectangular
Riemann sums. The neutrino grid begins at the recoil-dependent kinematic
minimum; recoil bins are evaluated at their upper-edge grid value. CEvNS
includes Standard Model vector, axial, and interference terms with a Helm
nuclear form factor. Cs-133 and I-127 target counts are calculated from their
mass fractions, molar masses, detector mass, and Avogadro's constant.

## Assumptions

- The source emits stopped-pion decay-at-rest neutrinos isotropically.
- Each stopped pi+ produces one prompt nu_mu, one delayed nu_e, and one delayed
  anti-nu_mu. Pi- and decay-in-flight contributions are omitted.
- Prompt nu_mu energy is 29.792 MeV. Delayed flavors use the standard muon-DAR
  Michel spectra up to 52.834185 MeV.
- The proton pulse is triangular, 290 ns wide, centered at 145 ns, and repeats
  at 20 Hz.
- Production time is the beam-pulse PDF convolved with exponential pion
  (26.03 ns) or muon (2196.9811 ns) decay.
- The configured timing window is finite. `annotate_only` records its
  acceptance without changing rates; `rate_windowed` applies the acceptance.
- The baseline is 23 m. The `1/(4 pi R^2)` term gives flux density at that
  distance; no detector projected-area or orientation factor is applied.
- The detector contains 25 full-crystal equivalents at `2 x 2 x 12 in`.
  DukeCEvNS calculates their mass from `4.51 g/cm3` CsI: `3.547471615 kg`
  each and `0.088686790368` tonnes total.
- Normalization uses `5.0e14 POT/s`, `0.099774 stopped pi+/POT`, and
  94,608,000 s exposure. The historical LAr flux is comparison metadata only.
- Recoils are unquenched nuclear-recoil energies. QF, CsI[Na/Tl] light yield,
  thresholds, efficiency, resolution, transport, electronics, and waveform
  response are deliberately excluded.
- Calculated event counts are expectations. The sampler output is used by the
  downstream Monte Carlo to generate discrete events.

## Outputs

```text
out/ccm_csi_truth_recoil_rates.csv
out/ccm_csi_recoil_sampler.csv
out/ccm_csi_truth_recoil_rates.meta.json
out/ccm_csi_recoil_by_component.csv
out/ccm_csi_time_by_component.csv
out/ccm_csi_recoil_time_2d.csv
out/ccm_csi_nu_recoil_2d.csv
out/plots/
```

The metadata JSON records the complete input configuration, normalization,
integration settings, timing acceptance, isotope totals, and flavor totals.
