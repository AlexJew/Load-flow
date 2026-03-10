# Load Flow Solver (4-Bus Example)

This repository contains a small, self-contained Python script (`load_flow_solve.py`) that computes bus voltages for the 4-bus load-flow exercise shown in the provided slides.

## What the script does
- Defines source phasors `V3`
- Builds branch admittances from the given purely reactive impedances `Za...Zg`.
- Converts the voltage sources to Thevenin equivalents: shunt admittances at buses 3 and 4 with injected currents `I = V/Z`.
- Assembles the bus admittance matrix `Y_bus` for the 4 internal buses.
- Solves the nodal equation `Y_bus · V = I` for complex bus voltages, then reports magnitudes and angles.
- Also shows the “susceptance-only” variant (factor out the common `-j`) so you can match the real-valued matrix used in the slides.

## How to run
```bash
python3 load_flow_solve.py
```

The script prints two sections:
1) **Full complex Y_bus**: direct solution with the complex admittance matrix.
2) **Real-valued susceptance matrix**: same physics with `-j` factored out from `Y_bus` and `I`, yielding identical voltages.

If you need to mirror the exact numeric table from the slides (including their asymmetric entry), uncomment the “Slide's printed Y_bus/I” block in `load_flow_solve.py`.

## Key outputs (default setup)
- Bus order: `[1, 2, 3, 4]`
- Voltage magnitudes (pu): approximately `[0.9750, 0.9728, 0.9941, 0.9534]`
- Voltage angles (deg): approximately `[-17.78, -18.02, -15.89, -20.18]`
