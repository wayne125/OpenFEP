from __future__ import annotations

from openfep.constants import RBFE_DEFAULT_LAMBDAS


def rbfe_lambdas(n: int = RBFE_DEFAULT_LAMBDAS) -> list[float]:
    """Return n evenly-spaced lambda values from 0.0 to 1.0 inclusive."""
    if n < 2:
        raise ValueError(f"n must be >= 2, got {n}")
    return [round(i / (n - 1), 6) for i in range(n)]
