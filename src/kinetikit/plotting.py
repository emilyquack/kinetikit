from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .models import michaelis_menten, substrate_inhibition
from .results import FitResult


def save_enzyme_plot(result: FitResult, output: str | Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    substrate = np.asarray(result.data["substrate"], dtype=float)
    rate = np.asarray(result.data["rate"], dtype=float)
    x_smooth = np.linspace(0, float(substrate.max()) * 1.08, 300)
    if result.model == "michaelis_menten":
        y_smooth = michaelis_menten(x_smooth, result.parameters["vmax"], result.parameters["km"])
    elif result.model == "substrate_inhibition":
        y_smooth = substrate_inhibition(x_smooth, result.parameters["vmax"], result.parameters["km"], result.parameters["ki"])
    else:
        raise ValueError(f"Unsupported enzyme plot model: {result.model}")

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(substrate, rate, s=70, color="#2f2a27", edgecolor="white", linewidth=1.5, label="data")
    ax.plot(x_smooth, y_smooth, color="#4472c4", linewidth=2.5, label="fit")
    ax.set_xlabel("Substrate concentration [S]")
    ax.set_ylabel("Initial rate v")
    ax.set_title("KinetiKit enzyme kinetics fit")
    ax.grid(True, alpha=0.25)
    ax.legend()
    text = [f"{k} = {v:.4g}" for k, v in result.parameters.items() if not k.endswith("_stderr")]
    text.append(f"R² = {result.goodness['r_squared']:.4f}")
    ax.text(0.04, 0.96, "\n".join(text), transform=ax.transAxes, va="top", bbox={"boxstyle": "round", "facecolor": "#fffaf4", "edgecolor": "#d9b99b"})
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def save_rate_plot(result: FitResult, output: str | Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    time = np.asarray(result.data["time"], dtype=float)
    concentration = np.asarray(result.data["concentration"], dtype=float)
    pred = np.asarray(result.predictions, dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(time, concentration, s=65, color="#2f2a27", edgecolor="white", linewidth=1.5, label="data")
    ax.plot(time, pred, color="#c26d67", linewidth=2.5, label=result.model.replace("_", " "))
    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration [A]")
    ax.set_title("KinetiKit rate-law fit")
    ax.grid(True, alpha=0.25)
    ax.legend()
    text = [f"{k} = {v:.4g}" for k, v in result.parameters.items()]
    text.append(f"R² = {result.goodness['r_squared']:.4f}")
    ax.text(0.04, 0.96, "\n".join(text), transform=ax.transAxes, va="top", bbox={"boxstyle": "round", "facecolor": "#fffaf4", "edgecolor": "#d9b99b"})
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def save_json(result: FitResult, output: str | Path) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return output
