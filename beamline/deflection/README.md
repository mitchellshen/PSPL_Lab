# Beam deflection in the beamline's magnetic fields

`gyroradiusCalculator.m` estimates how far the ion beam is pushed sideways by magnetic fields along the beamline. It is split into cells (`%%`). Run them in order: the detailed-deflection cell reuses `E` from the first cell.

| Cell | What it computes | Inputs to edit |
|---|---|---|
| Simple Gyroradius calculator and simple deflection | Gyroradius `rg` (m, km, AU) and deflection `dX` (m, cm) for one case, using the simple geometry `dX = rg (1 − cos(asin(L/rg)))` | `E` (eV), `mp` (particle mass in kg; default O⁺), `L` (field length, m), `BE` (field, T) |
| A loop of Gyroradius calculator | The same over a range of field strengths | `Brange`, `L` |
| Detailed deflection | Deflection `dy` (m, cm) along the whole beamline: ambient field outside, Helmholtz-coil field (4/5)^(3/2) μ₀ n I / R plus ambient in the middle section | `I` (coil current, A), `amu`, `L` (coil section length, m), `Lfull`, `n`, `R` |
| Plot | `rg` (log scale) and `dX` against field, with a dotted line at 8 cm | — |

The helper `computeRg_dX(L, B, m, E)` used by the first two cells is defined at the bottom of the script.

## Notes in this folder

- `particle_deflection_Part_a.pptx` — sketch of the simple geometry used in the first two cells.
- `particle_deflection_Bfield.pdf` — derivation of the three-region `Δy` used in the detailed cell. It follows Eq. A14 of Khoo et al. (2023), expanded without assuming Δy = 0, with constant v_x and v_y,0 ≈ 0.
