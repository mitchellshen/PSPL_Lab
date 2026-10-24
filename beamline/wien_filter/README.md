# Wien filter

`Ion_source_Wien.m` has two parts.

## 1. Plate-voltage calculator

For plate gap `d` (m), Wien filter magnetic field `B` (T) and beam energy `E`, the script prints the plate voltage that passes each species in its list (H, H₂, mass 3, He, N, Ne, N₂, Ar, Xe):

```
U = d · sqrt(2 B² E q / m)
```

with `E` in volts (energy per charge). `d` is an effective gap fitted to measured peaks; the commented lines under the calculation show how it was back-computed from measured H₂⁺ and He⁺ peak voltages.

## 2. Scan plot

Reads a Wien filter plate-potential scan and plots beam current against plate voltage and against mass.

- **Input file:** `Wien-Filter-Plate-Potential-Scan__YYYY.MM.DD_N.txt`, two numeric columns: plate voltage (V) and picoammeter current (A). Set `path`, `yyyy`, `mm`, `dd` and `num` to pick the file.
- **Processing:** the current is averaged at each voltage in `V_wien_seq` (0–256 V in 1 V steps, then 258–300 V in 2 V steps) and shifted so its minimum is zero. Change `V_wien_seq` if the scan used different steps.
- **Output:** a two-panel figure (current in pA vs plate voltage, and vs mass in amu using the same `d`, `B` and `E`), saved at 300 dpi as `<scan file name>_v2.png` next to the data file.
