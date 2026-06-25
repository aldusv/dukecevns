#!/usr/bin/env python3
import csv
import json
import math
import subprocess
import tempfile
from pathlib import Path


DUKE_DIR = Path(__file__).resolve().parent
EXE = DUKE_DIR / "ccm_truth_rates"
CONFIG = DUKE_DIR / "jsonfiles" / "ccm_csi_truth.json"


EXPECTED_TOTAL = 1667.6856635745974
EXPECTED_COMPONENTS = {
    "prompt_numu": 484.1955703106046,
    "delayed_nue": 532.0977293661168,
    "delayed_numubar": 651.3923638978722,
}
EXPECTED_ISOTOPES = {
    "Cs133": 874.1055182750097,
    "I127": 793.5801452995828,
}


def close(actual: float, expected: float, rel: float = 1e-12) -> None:
    if not math.isclose(actual, expected, rel_tol=rel, abs_tol=rel):
        raise AssertionError(f"{actual} != {expected}")


def make_temp_config(out_dir: Path) -> Path:
    with CONFIG.open() as stream:
        config = json.load(stream)

    config["output_csv"] = str(out_dir / "truth.csv")
    config["sampler_csv"] = str(out_dir / "sampler.csv")
    config["metadata_json"] = str(out_dir / "metadata.json")
    config["diagnostics"]["recoil_by_component_csv"] = str(out_dir / "recoil_by_component.csv")
    config["diagnostics"]["time_by_component_csv"] = str(out_dir / "time_by_component.csv")
    config["diagnostics"]["recoil_time_csv"] = str(out_dir / "recoil_time.csv")
    config["diagnostics"]["nu_recoil_csv"] = str(out_dir / "nu_recoil.csv")

    path = out_dir / "ccm_csi_truth_temp.json"
    with path.open("w") as stream:
        json.dump(config, stream, indent=2)
        stream.write("\n")
    return path


def require_outputs(out_dir: Path) -> None:
    for name in (
        "truth.csv",
        "sampler.csv",
        "metadata.json",
        "recoil_by_component.csv",
        "time_by_component.csv",
        "recoil_time.csv",
        "nu_recoil.csv",
    ):
        path = out_dir / name
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing nonempty output {path}")


def check_metadata(path: Path) -> None:
    with path.open() as stream:
        meta = json.load(stream)

    close(meta["total_expected_events"], EXPECTED_TOTAL)
    close(meta["rate_ledger"]["total_expected_events"], EXPECTED_TOTAL)
    close(meta["rate_ledger"]["total_unwindowed_expected_events"], EXPECTED_TOTAL)
    for name, expected in EXPECTED_COMPONENTS.items():
      close(meta["rate_ledger"]["expected_events_by_component"][name], expected)
    for name, expected in EXPECTED_ISOTOPES.items():
      close(meta["rate_ledger"]["expected_events_by_isotope"][name], expected)

    if meta["cci_consumed_outputs"] != ["sampler_csv", "metadata_json"]:
        raise AssertionError("unexpected CCI consumed output list")
    if meta["output_roles"]["sampler_csv"] != "cci_recoil_sampler":
        raise AssertionError("sampler output role missing")
    if meta["output_roles"]["recoil_time_csv"] != "diagnostic":
        raise AssertionError("diagnostic output role missing")
    if meta["unsupported_physics"]["decay_in_flight_contributions"] != "omitted":
        raise AssertionError("unsupported physics metadata missing")


def check_sampler(path: Path) -> None:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)

    if not rows:
        raise AssertionError("sampler has no rows")
    first = rows[0]
    if first["target_z"] != "55" or first["target_a"] != "133":
        raise AssertionError("unexpected first sampler target")
    close(float(first["recoil_energy_mev"]), 0.0001)
    close(float(first["expected_events_bin"]), 4.499764289665)
    close(float(first["nu_energy_mev"]), 29.792)

    total = sum(float(row["probability"]) for row in rows)
    close(total, EXPECTED_TOTAL, rel=1e-11)


def check_truth(path: Path) -> None:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        first = next(reader)

    if first["material"] != "CsI" or first["isotope"] != "Cs133":
        raise AssertionError("unexpected first truth row")
    close(float(first["recoil_energy_mev"]), 0.0001)
    close(float(first["expected_events_bin"]), 4.499764289665)


def check_summary_mode(config: Path, out_dir: Path) -> None:
    summary_dir = out_dir / "summary_only"
    summary_dir.mkdir()
    summary_config = make_temp_config(summary_dir)
    result = subprocess.run(
        [str(EXE), "--summary", str(summary_config)],
        cwd=DUKE_DIR,
        check=True,
        text=True,
        capture_output=True,
    )
    if "total_expected_events: 1667.68566357" not in result.stdout:
        raise AssertionError("summary output does not include expected total")
    if any(path.name != "ccm_csi_truth_temp.json" for path in summary_dir.iterdir()):
        raise AssertionError("summary mode wrote output files")


def main() -> int:
    if not EXE.exists():
        raise SystemExit("ccm_truth_rates executable is missing; run `make` first")

    with tempfile.TemporaryDirectory(prefix="dukecevns_test_") as temp:
        out_dir = Path(temp)
        config = make_temp_config(out_dir)
        subprocess.run([str(EXE), str(config)], cwd=DUKE_DIR, check=True)
        require_outputs(out_dir)
        check_metadata(out_dir / "metadata.json")
        check_sampler(out_dir / "sampler.csv")
        check_truth(out_dir / "truth.csv")
        check_summary_mode(config, out_dir)

    print("DukeCEvNS regression checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
