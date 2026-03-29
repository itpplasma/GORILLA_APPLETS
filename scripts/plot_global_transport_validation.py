#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_table(path: pathlib.Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="ascii") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    header = lines[0].split()
    columns: dict[str, list[float] | list[bool]] = {name: [] for name in header}
    for line in lines[1:]:
        fields = line.split()
        for name, value in zip(header, fields):
            if value == "T":
                columns[name].append(True)  # type: ignore[arg-type]
            elif value == "F":
                columns[name].append(False)  # type: ignore[arg-type]
            else:
                columns[name].append(float(value))  # type: ignore[arg-type]
    return {name: np.asarray(values) for name, values in columns.items()}


def parse_experiment(path: pathlib.Path) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    with path.open("r", encoding="ascii") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    shell_path = path.with_suffix(".shell.tmp")
    shell_path.write_text("\n".join(lines[:-2]) + "\n", encoding="ascii")
    try:
        shell_data = parse_table(shell_path)
    finally:
        shell_path.unlink(missing_ok=True)
    summary_header = lines[-2].split()
    summary_values = lines[-1].split()
    return shell_data, dict(zip(summary_header, summary_values))


def save_profile_plot(
    outpath: pathlib.Path,
    title: str,
    ylabel: str,
    fit_s: np.ndarray,
    fit_mean: np.ndarray,
    fit_2sigma: np.ndarray,
    local_s: np.ndarray,
    local_mean: np.ndarray,
    local_2sigma: np.ndarray,
    local_valid: np.ndarray,
    s_max: float = 1.0,
) -> None:
    mask = fit_s <= s_max
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.fill_between(fit_s[mask], (fit_mean - fit_2sigma)[mask], (fit_mean + fit_2sigma)[mask], color="tab:blue", alpha=0.2)
    ax.plot(fit_s[mask], fit_mean[mask], color="tab:blue", linewidth=2.0, label="global fit")
    local_mask = local_s <= s_max
    valid = local_valid & local_mask
    invalid = ~local_valid & local_mask
    ax.errorbar(
        local_s[valid],
        local_mean[valid],
        yerr=local_2sigma[valid],
        fmt="o",
        color="tab:orange",
        capsize=3.0,
        label="local reference",
    )
    if np.any(invalid):
        ax.plot(local_s[invalid], np.zeros(np.count_nonzero(invalid)), "x", color="tab:red", label="invalid local")
    ax.set_xlabel("s")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def save_signal_plot(
    outpath: pathlib.Path,
    experiments: list[tuple[dict[str, np.ndarray], dict[str, str]]],
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
    colors = ["tab:green", "tab:purple"]
    labels = ["experiment 1", "experiment 2"]

    for (rows, summary), color, label in zip(experiments, colors, labels):
        s = rows["shell_center_s"].astype(float)
        source = rows["source_density"].astype(float)
        density = rows["density"].astype(float)
        density_2sigma = rows["density_2sigma"].astype(float)
        flux = rows["right_boundary_integrated_flux"].astype(float)
        flux_2sigma = rows["right_boundary_integrated_flux_2sigma"].astype(float)
        supported = rows["right_boundary_supported"].astype(bool)
        supported_count = int(summary["supported_flux_boundaries"])

        axes[0].plot(s, source, "--", color=color, linewidth=1.5, label=f"{label} source")
        axes[0].fill_between(s, density - density_2sigma, density + density_2sigma, color=color, alpha=0.15)
        axes[0].plot(s, density, color=color, linewidth=2.0, label=f"{label} density ({supported_count} flux boundaries)")

        axes[1].fill_between(s, flux - flux_2sigma, flux + flux_2sigma, color=color, alpha=0.15)
        axes[1].plot(s, flux, color=color, linewidth=2.0, label=f"{label} flux")
        axes[1].scatter(s[supported], flux[supported], color=color, marker="o", s=20)

    axes[0].set_ylabel("density / source")
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].set_xlabel("s")
    axes[1].set_ylabel("integrated boundary flux")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: plot_global_transport_validation.py <case_dir> <output_dir>", file=sys.stderr)
        return 1

    case_dir = pathlib.Path(sys.argv[1]).resolve()
    output_dir = pathlib.Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    experiments: list[tuple[dict[str, np.ndarray], dict[str, str]]] = []
    for name in ("global_transport_experiment_1.dat", "global_transport_experiment_2.dat"):
        path = case_dir / name
        if path.exists():
            experiments.append(parse_experiment(path))
    if experiments:
        save_signal_plot(output_dir / "global_transport_signal_quality.png", experiments)

    history_path = case_dir / "global_transport_convergence_history.dat"
    if history_path.exists():
        hist = parse_table(history_path)
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.0), sharex=True)
        iters = hist["iteration"].astype(int)
        accepted = hist["accepted"].astype(bool)

        axes[0].semilogy(iters, hist["objective"], "o-", markersize=4)
        axes[0].set_ylabel("objective")
        axes[0].grid(True, alpha=0.3)

        axes[1].semilogy(iters, hist["gradient_norm"], "o-", markersize=4, color="tab:orange")
        axes[1].set_ylabel("gradient norm")
        axes[1].grid(True, alpha=0.3)

        axes[2].semilogy(iters, hist["damping"], "o-", markersize=4, color="tab:green")
        if any(~accepted):
            axes[2].scatter(iters[~accepted], hist["damping"][~accepted], color="tab:red", marker="x", s=60, zorder=5, label="rejected")
            axes[2].legend()
        axes[2].set_ylabel("LM damping")
        axes[2].set_xlabel("iteration")
        axes[2].grid(True, alpha=0.3)

        fig.suptitle("LM convergence history")
        fig.tight_layout()
        fig.savefig(output_dir / "global_transport_convergence_history.png", dpi=200)
        plt.close(fig)

    fit_path = case_dir / "global_transport_fit_profiles.dat"
    comp_path = case_dir / "global_transport_comparison_profiles.dat"
    if not fit_path.exists() or not comp_path.exists():
        return 0

    fit_data = parse_table(fit_path)
    comp_data = parse_table(comp_path)

    s_max = 1.0
    for exp_data, _ in experiments:
        supported = exp_data.get("right_boundary_supported", np.array([])).astype(bool)
        exp_s = exp_data.get("shell_center_s", np.array([])).astype(float)
        if len(supported) > 0 and np.any(supported):
            s_max = min(s_max, float(exp_s[supported].max()) + 0.05)
    s_max = min(s_max, 0.95)

    save_profile_plot(
        output_dir / "global_vs_local_A_2sigma.png",
        "Global vs local convection coefficient",
        "A(s)",
        fit_data["boundary_s"].astype(float),
        fit_data["A_s"].astype(float),
        fit_data["A_s_2sigma"].astype(float),
        comp_data["boundary_s"].astype(float),
        comp_data["local_A"].astype(float),
        comp_data["local_A_2sigma"].astype(float),
        comp_data["local_valid"].astype(bool),
        s_max=s_max,
    )
    save_profile_plot(
        output_dir / "global_vs_local_B_2sigma.png",
        "Global vs local diffusion coefficient",
        "B(s)",
        fit_data["boundary_s"].astype(float),
        fit_data["B_s"].astype(float),
        fit_data["B_s_2sigma"].astype(float),
        comp_data["boundary_s"].astype(float),
        comp_data["local_B"].astype(float),
        comp_data["local_B_2sigma"].astype(float),
        comp_data["local_valid"].astype(bool),
        s_max=s_max,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
