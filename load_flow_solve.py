"""
Compute bus voltages for the 4-bus load-flow exercise using nodal analysis.
All reactances are purely imaginary (inductive), so admittances are purely
imaginary and negative.

There seem to be some errors in the solutions: 
1. The current values aren't correct, they do not equate to the results from the computations on paper
2. The bus admittance matrix is not symmetrical

This needs to be checked with Fabrizio
"""

from __future__ import annotations

import cmath
import numpy as np


def polar(mag: float, deg: float) -> complex:
    """Return a complex number from magnitude and degrees."""
    return cmath.rect(mag, np.deg2rad(deg))


def build_ybus(z_branches: dict[tuple[int, int], complex], y_shunts: dict[int, complex]) -> np.ndarray:
    """Assemble the bus admittance matrix."""
    bus_ids = sorted({bus for pair in z_branches for bus in pair} | set(y_shunts))
    n = len(bus_ids)
    ybus = np.zeros((n, n), dtype=complex)
    index = {bus: i for i, bus in enumerate(bus_ids)}

    for (a, b), z in z_branches.items():
        y = 1 / z
        i, j = index[a], index[b]
        ybus[i, i] += y
        ybus[j, j] += y
        ybus[i, j] -= y
        ybus[j, i] -= y

    for bus, y in y_shunts.items():
        ybus[index[bus], index[bus]] += y

    return ybus


def solve_and_report(name: str, ybus: np.ndarray, i_vec: np.ndarray, bus_ids: list[int]) -> None:
    """Solve Y_bus V = I and print voltages, magnitudes, and angles."""
    v = np.linalg.solve(ybus, i_vec)

    print(f"\n=== {name} ===")
    with np.printoptions(precision=4, suppress=True):
        print("Y_bus:\n", ybus)
        print("\nInjected currents I:\n", i_vec)
        print("\nBus voltages:\n", v)

    mags = np.abs(v)
    angles_deg = np.rad2deg(np.angle(v))
    print("\nVoltage magnitudes (pu):")
    for bus, mag in zip(bus_ids, mags):
        print(f"Bus {bus}: {mag:.4f}")
    print("\nVoltage angles (deg):")
    for bus, ang in zip(bus_ids, angles_deg):
        print(f"Bus {bus}: {ang:.2f}")


def main() -> None:
    # Source definitions
    v3 = polar(1.25, 0)       # V3 magnitude and angle (deg)
    v4 = polar(0.85, -45)     # V4 magnitude and angle (deg)

    # Line/source reactances (jX)
    z = {
        "a": 1j * 1.25,
        "b": 1j * 0.25,
        "c": 1j * 0.25,
        "d": 1j * 0.125,
        "e": 1j * 0.2,
        "f": 1j * 0.4,
        "g": 1j * 1.25,
    }

    # Branches between buses (bus numbering matches the diagram: 1–4 internal nodes)
    z_branches = {
        (3, 2): z["b"],
        (3, 1): z["c"],
        (2, 1): z["d"],
        (2, 4): z["e"],
        (1, 4): z["f"],
    }

    # Thevenin source admittances (shunts at buses 3 and 4)
    y_shunts = {
        3: 1 / z["a"],
        4: 1 / z["g"],
    }

    # Current injections from sources
    injections = {
        3: v3 / z["a"],
        4: v4 / z["g"],
    }

    # Build Ybus and I vector
    ybus = build_ybus(z_branches, y_shunts)
    bus_ids = sorted({bus for pair in z_branches for bus in pair} | set(y_shunts))
    i_vec = np.zeros(len(bus_ids), dtype=complex)
    for bus, current in injections.items():
        idx = bus_ids.index(bus)
        i_vec[idx] = current

    # Solve using the full complex admittance matrix (physically faithful).
    solve_and_report("Full complex Y_bus (direct nodal solution)", ybus, i_vec, bus_ids)

    # Solve using the susceptance-only real matrix (factor out -j) as in the slide deck.
    # This rotates both Y_bus and I by -j so the linear system is purely real.
    ybus_real = (-1j) * ybus  # removes the global j factor; should be real-valued
    i_vec_rot = (-1j) * i_vec
    solve_and_report("Real-valued susceptance matrix (-j factored out)", ybus_real, i_vec_rot, bus_ids)

    # If you need to reproduce the exact slide numerics (including their asymmetric Y entry),
    # uncomment below. The asymmetry (8 vs 4) likely comes from a transcription error.
    ybus_slide = np.array([
        [-14.5, 8.0, 8.0, 2.5],
        [8.0, -17.0, 4.0, 5.0],
        [4.0, 4.0, -8.8, 0.0],
        [2.5, 5.0, 0.0, -8.3],
    ], dtype=float)
    i_slide = np.array([0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 1.0j, -0.4808 + 0.4808j])
    solve_and_report("Slide's printed Y_bus/I (asymmetric reproduction)", ybus_slide, i_slide, bus_ids)


if __name__ == "__main__":
    main()
