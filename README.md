Simple CEvNS rate tools, with a CCM CsI truth-recoil path for CCI.

## CCM Truth Recoils

Build and run:

```sh
make ccm_truth_rates
./ccm_truth_rates ccm_csi_truth
python3 plot_ccm_diagnostics.py
```

The CCM path writes unquenched nuclear recoil rates only. It does not apply QF,
thresholds, smearing, trigger efficiency, or detector response.

Detector mass is derived from crystal count and per-crystal mass. The current
config uses 16 crystals at 4 kg each, with 2 x 2 x 12 inch dimensions.

Normalization is:

```text
POT/s * stopped pi+ / POT * neutrinos / stopped pi+ * 1 / (4 pi R^2)
```

Timing uses a triangular 290 ns beam pulse convolved with prompt pion or delayed
muon decay. Timing windows are annotated by default; `rate_windowed` scales rates
by the finite window acceptance.

Main outputs:

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
