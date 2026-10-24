# C-foil absolute beam monitor (ABM)

Tools for the carbon-foil ABM's counts. Two CEMs (Aft and Fore) and their coincidence rate give the absolute detection efficiency (ADE) and the absolute beam flux.

## Data files

Raw data are not kept in this repo; they live on OneDrive (for example under `Beam_CAL_database/` or `IMAP-Lo_CAL_database/.../ABM_reference/`). Each acquisition day is a folder of numbered CSV files:

```
HVPS-ABM-Data-YYYY.MM.DD/
    ABM-Counts_YYYY.MM.DD_1.csv
    ABM-Counts_YYYY.MM.DD_2.csv
    ...
```

| Column | Name | Notes |
|---|---|---|
| 1 | `t` | timestamp, `MM/dd/yyyy/HH:mm:ss.S` |
| 2 | `IN_CEM` | |
| 3 | `OUT_CEM` | |
| 4 | `IN_CEM_rate` | Aft CEM |
| 5 | `OUT_CEM_rate` | Fore CEM |
| 6 | `COIN_rate` | coincidences |

In the current setup `IN` is the Aft CEM and `OUT` is the Fore CEM. The analysis uses columns 4–6: each row is treated as the counts in one sample, and dividing by the elapsed time gives counts/s.

## Functions

**`import_ABM_file`** reads one CSV into column vectors. The optional second input selects data rows (default `[2, Inf]`, which skips the header).

```matlab
[t, IN_CEM, OUT_CEM, IN_CEM_rate, OUT_CEM_rate, COIN_rate] = import_ABM_file(filename);
```

**`determine_avg_ABM_rate`** averages over a time window `[t1, t2]` (datetimes). Rates are in counts/s with Poisson uncertainties; `ADE = COIN² / (IN · OUT)`, with its uncertainty propagated from the three rates.

```matlab
[avg_IN, avg_OUT, avg_COIN, avg_ADE, avg_IN_err, avg_OUT_err, avg_COIN_err, avg_ADE_err, t_desired, t_elapsed] = ...
    determine_avg_ABM_rate(t1, t2, t, IN_CEM_rate, OUT_CEM_rate, COIN_rate);
```

Example:

```matlab
f = 'HVPS-ABM-Data-2024.12.06/ABM-Counts_2024.12.06_6.csv';
[t, ~, ~, aft, fore, coin] = import_ABM_file(f);

t1 = datetime(2024, 12, 6, 14, 30, 0);
t2 = datetime(2024, 12, 6, 14, 40, 0);
[A, F, C, ADE] = determine_avg_ABM_rate(t1, t2, t, aft, fore, coin);

beam_flux = A * F / C;   % absolute beam flux, counts/s
```

## Batch script: `ABM_process.m`

Processes every file from one acquisition day. Edit the inputs at the top: `yr` (two-digit year, e.g. `24`), `mon`, `day`, `filenum_max` (number of files) and `dir` (the folder that contains `HVPS-ABM-Data-...`).

It writes `tbl_YYYY.MM.DD.txt` into `dir`: tab-delimited, no header, one row per file.

| Column | Quantity |
|---|---|
| 1 | file number |
| 2 | time span (s) |
| 3, 4, 5 | mean Fore, Aft and coincidence rates (counts/s) |
| 6, 7, 8 | uncertainties of columns 3–5 |
| 9, 10 | ADE as the mean of per-row values, and its uncertainty |
| 11, 12 | ADE from summed counts, and uncertainty (same value as column 10) |
| 13, 14 | beam flux Fore · Aft / Coin as the mean of per-row values (counts/s), and its uncertainty |
| 15, 16 | beam flux from summed counts (counts/s), and uncertainty (same value as column 14) |

Rows that would divide by zero are skipped in the per-row means.
