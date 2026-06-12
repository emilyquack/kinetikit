from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(slots=True)
class FitResult:
    """Serializable result from a KinetiKit analysis."""

    model: str
    equation: str
    parameters: dict[str, float]
    goodness: dict[str, float]
    n_points: int
    derived: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    data: dict[str, list[float]] = field(default_factory=dict)
    predictions: list[float] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "model": self.model,
            "equation": self.equation,
            "parameters": self.parameters,
            "derived": self.derived,
            "goodness": self.goodness,
            "n_points": self.n_points,
            "warnings": self.warnings,
            "data": self.data,
            "predictions": self.predictions,
        }

    def summary_lines(self) -> list[str]:
        lines = [f"Model: {self.model}", f"Equation: {self.equation}", f"n = {self.n_points}"]
        if self.parameters:
            lines.append("Parameters:")
            lines.extend(f"  {name}: {value:.6g}" for name, value in self.parameters.items())
        if self.derived:
            lines.append("Derived:")
            lines.extend(f"  {name}: {value:.6g}" for name, value in self.derived.items())
        if self.goodness:
            lines.append("Goodness of fit:")
            lines.extend(f"  {name}: {value:.6g}" for name, value in self.goodness.items())
        if self.warnings:
            lines.append("Warnings:")
            lines.extend(f"  - {warning}" for warning in self.warnings)
        return lines

    def summary(self) -> str:
        return "\n".join(self.summary_lines())
