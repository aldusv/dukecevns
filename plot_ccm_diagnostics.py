#!/usr/bin/env python3
"""Plot CCM truth-rate diagnostics produced by ccm_truth_rates."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="ccm-mpl-"))

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def padded_extent(values: list[float]) -> tuple[float, float]:
    low = min(values)
    high = max(values)
    if low == high:
        pad = max(abs(low) * 0.01, 0.5)
        return low - pad, high + pad
    return low, high


def binned_extent(values: list[float], bin_width: float | None = None) -> tuple[float, float]:
    low = min(values)
    high = max(values)
    if bin_width is not None and bin_width > 0.0:
        return low - 0.5 * bin_width, high + 0.5 * bin_width
    if low == high:
        return padded_extent(values)
    sorted_values = sorted(values)
    low_step = sorted_values[1] - sorted_values[0]
    high_step = sorted_values[-1] - sorted_values[-2]
    return sorted_values[0] - 0.5 * low_step, sorted_values[-1] + 0.5 * high_step


def plot_recoil_by_component(rows: list[dict[str, str]], output_dir: Path) -> None:
    energy_kev = np.array([float(row["recoil_energy_mev"]) * 1000.0 for row in rows])
    component_columns = [
        key for key in rows[0].keys()
        if key.endswith("_events") and key != "total_events"
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    for column in component_columns:
        label = column.removesuffix("_events")
        values = np.array([float(row[column]) for row in rows])
        ax.step(energy_kev, values, where="mid", label=label)
    total = np.array([float(row["total_events"]) for row in rows])
    ax.step(energy_kev, total, where="mid", color="black", linewidth=1.2, label="total")
    ax.set_xlabel("Recoil energy [keVnr]")
    ax.set_ylabel("Expected events / recoil bin")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "ccm_recoil_by_component.png", dpi=180)
    plt.close(fig)


def plot_time_by_component(rows: list[dict[str, str]], output_dir: Path) -> None:
    by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_component[row["source_component"]].append(row)

    fig, ax = plt.subplots(figsize=(8, 5))
    for component, component_rows in by_component.items():
        time_us = np.array([float(row["time_center_ns"]) / 1000.0 for row in component_rows])
        events = np.array([float(row["expected_events_bin"]) for row in component_rows])
        ax.step(time_us, events, where="mid", label=component)
    ax.set_xlabel("Time since pulse start [us]")
    ax.set_ylabel("Expected events / time bin")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "ccm_time_by_component.png", dpi=180)
    plt.close(fig)


def plot_recoil_time(rows: list[dict[str, str]], output_dir: Path) -> None:
    by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_component[row["source_component"]].append(row)

    for component, component_rows in by_component.items():
        recoil_values = sorted({float(row["recoil_energy_mev"]) * 1000.0 for row in component_rows})
        time_values = sorted({float(row["time_center_ns"]) / 1000.0 for row in component_rows})
        recoil_bin_width = float(component_rows[0]["recoil_bin_width_mev"]) * 1000.0
        time_bin_width = float(component_rows[0]["time_bin_width_ns"]) / 1000.0
        recoil_index = {value: index for index, value in enumerate(recoil_values)}
        time_index = {value: index for index, value in enumerate(time_values)}
        grid = np.zeros((len(time_values), len(recoil_values)))
        for row in component_rows:
            recoil = float(row["recoil_energy_mev"]) * 1000.0
            time = float(row["time_center_ns"]) / 1000.0
            grid[time_index[time], recoil_index[recoil]] += float(row["expected_events_bin"])

        fig, ax = plt.subplots(figsize=(8, 5))
        mesh = ax.imshow(
            grid,
            origin="lower",
            aspect="auto",
            extent=[
                *binned_extent(recoil_values, recoil_bin_width),
                *binned_extent(time_values, time_bin_width),
            ],
        )
        ax.set_xlabel("Recoil energy [keVnr]")
        ax.set_ylabel("Time since pulse start [us]")
        ax.set_title(component)
        fig.colorbar(mesh, ax=ax, label="Expected events / 2D bin")
        fig.tight_layout()
        fig.savefig(output_dir / f"ccm_recoil_time_{component}.png", dpi=180)
        plt.close(fig)


def plot_nu_recoil(rows: list[dict[str, str]], output_dir: Path) -> None:
    by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_component[row["source_component"]].append(row)

    for component, component_rows in by_component.items():
        recoil_values = sorted({float(row["recoil_energy_mev"]) * 1000.0 for row in component_rows})
        nu_values = sorted({float(row["nu_energy_mev"]) for row in component_rows})
        recoil_bin_width = float(component_rows[0]["recoil_bin_width_mev"]) * 1000.0
        nu_bin_width = float(component_rows[0]["nu_energy_bin_width_mev"])
        recoil_index = {value: index for index, value in enumerate(recoil_values)}
        nu_index = {value: index for index, value in enumerate(nu_values)}
        grid = np.zeros((len(nu_values), len(recoil_values)))
        for row in component_rows:
            recoil = float(row["recoil_energy_mev"]) * 1000.0
            nu_energy = float(row["nu_energy_mev"])
            grid[nu_index[nu_energy], recoil_index[recoil]] += float(row["expected_events_bin"])

        fig, ax = plt.subplots(figsize=(8, 5))
        mesh = ax.imshow(
            grid,
            origin="lower",
            aspect="auto",
            extent=[
                *binned_extent(recoil_values, recoil_bin_width),
                *binned_extent(nu_values, nu_bin_width),
            ],
        )
        ax.set_xlabel("Recoil energy [keVnr]")
        ax.set_ylabel("Neutrino energy [MeV]")
        ax.set_title(component)
        fig.colorbar(mesh, ax=ax, label="Expected events / 2D bin")
        fig.tight_layout()
        fig.savefig(output_dir / f"ccm_nu_recoil_{component}.png", dpi=180)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("out"))
    parser.add_argument("--output-dir", type=Path, default=Path("out/plots"))
    parser.add_argument(
        "--prefix",
        default="ccm_csi",
        help="Diagnostic filename prefix, e.g. ccm_csi for ccm_csi_recoil_by_component.csv.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    recoil_rows = read_csv(args.input_dir / f"{args.prefix}_recoil_by_component.csv")
    time_rows = read_csv(args.input_dir / f"{args.prefix}_time_by_component.csv")
    recoil_time_rows = read_csv(args.input_dir / f"{args.prefix}_recoil_time_2d.csv")
    nu_recoil_rows = read_csv(args.input_dir / f"{args.prefix}_nu_recoil_2d.csv")

    plot_recoil_by_component(recoil_rows, args.output_dir)
    plot_time_by_component(time_rows, args.output_dir)
    plot_recoil_time(recoil_time_rows, args.output_dir)
    plot_nu_recoil(nu_recoil_rows, args.output_dir)
    print(f"Wrote CCM diagnostic plots to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
