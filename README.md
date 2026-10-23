# SSA gap filling for GRACE / GRACE-FO

This package adapts your `ssa_missing_iterative.m`, `fun_SSA_filling_a.m`, and
`fun_SSA_filling_b.m` to Python. It processes one monthly time series at each
spatial location, with bounded spatial chunks. It works with grids or mascon
dimensions, provided the chosen variable contains a decoded `time` dimension.
Python 3.10+ is recommended.

## Install on Windows

Download and extract this package. Open its folder in File Explorer, click the
address bar, type `cmd`, and press Enter. In that Command Prompt run:

```bat
py -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
python demo.py
```

If `py` is not found, install Python from https://www.python.org/downloads/ first.
The demo creates `demo_input.nc` and `demo_filled.nc` and prints a reconstruction
error for a known synthetic seasonal signal. These are not real GRACE data.
Rename/delete the demo files if you want to repeat the demo.

## Install on Derecho / Linux

Use your site's supported Python module or an existing Python environment, then:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python demo.py
```

For production runs, use a compute allocation following your site's scheduler
policy. No site-specific module names or account/queue settings are assumed.

## Inspect the real file first

The command you supplied printed help because `-g` requires a group argument.
Use `-h` for the header:

```bash
ncdump -h CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections_FillNaN.nc
```

Without `ncdump`, Python can show the same essential information:

```bash
python -c "import xarray as xr; ds=xr.open_dataset('CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections_FillNaN.nc'); print(ds); print(ds.time.attrs); print(ds.time.encoding); ds.close()"
```

Your supplied header confirms `lwe_thickness(time,lat,lon)`, `Units="cm"`,
NaN `_FillValue`, and time in days since 2002-04-01 with a proleptic Gregorian
calendar. Dimensions are time=293, lat=720, lon=1440. Time values themselves
were not supplied, so the program still checks monthly continuity when run.
Only root-group NetCDF variables are supported by this command-line interface.

## Run on a NetCDF file

Place the Python script alongside your file, or use full paths. This single-line
command works in Windows Command Prompt and Linux after environment activation:

```bash
python ssa_gapfill.py CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections_FillNaN.nc CSR_SSA_filled.nc --variable lwe_thickness --window 60 --rank 10
```

For your file, `lwe_thickness` is the correct variable name. The output
must be a NEW filename. Input and pre-existing output files are never overwritten.
Storage is needed for the original dataset plus the filled variable and flags.

By default, all monthly records must exist already, including records whose data
are NaN. If missing months are entirely absent, add `--regularize-months`:

```bash
python ssa_gapfill.py input.nc output.nc --variable lwe_thickness --regularize-months
```

This option assigns records to their calendar months, labels all time coordinates
on the first day of the month, and inserts missing months throughout the ENTIRE
dataset. Existing time-dependent ancillary data, including time bounds, will be
missing in inserted months; no ancillary data are interpolated. Original dates
are retained only in the unchanged source file. Inspect this choice if your time
values denote interval endpoints rather than representative monthly dates.
Multiple records in a month, unsorted dates, and undecodable time axes are rejected.
CF standard and cftime calendars are supported; numeric decimal-year axes without
CF date units require conversion before using this NetCDF interface.

CF `_FillValue` and `missing_value` are decoded automatically. For an undocumented
sentinel such as -9999999, optionally use `--missing-below -1000000`. This treats
values AT OR BELOW the threshold as missing in the selected variable. Do not use
this flag unless the threshold is valid for your dataset. Infinity is also treated
as missing by the NetCDF workflow. Zero is a valid observation, never a gap.

Your full grid contains 1,036,800 locations. Start with a 2x2 pilot subset:

```bash
python ssa_gapfill.py CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections_FillNaN.nc CSR_SSA_pilot.nc --variable lwe_thickness --isel lat=300:302 --isel lon=300:302
```

Indices are zero-based and the stop is excluded. The example indices are only
for testing; choose indices over your study area with valid observations.
Subset outputs contain only selected coordinates; inherited global coverage
attributes still describe the parent file and a history note records the subset.
Run the full command without `--isel` after checking the pilot.

You can also create a separate spatial subset beforehand:

```python
import xarray as xr
with xr.open_dataset('input.nc') as ds:
    ds.isel(lat=slice(0, 2), lon=slice(0, 2)).to_netcdf('small_input.nc')
```

Use dimension names appropriate to your file. `--spatial-chunk 8` (default)
limits each spatial dimension's chunk size, while retaining the entire time
series per task. Default `--workers 1` avoids oversubscribing CPU libraries.
More workers can help, but the full MATLAB parameter search at every grid cell
can be very expensive. Test timing and reconstruction quality on representative
land cells before scaling up. Permanently masked all-NaN cells remain NaN.
The filled float64 array alone would occupy about 2.43 GB uncompressed for your
293x720x1440 grid, in addition to the retained input variable and diagnostics.
The actual NetCDF size depends on compression. Many displayed grid pixels may
represent the same native mascon; this code processes grid cells independently.

## Method and differences from MATLAB

1. Embed the series as overlapping windows, compute an SVD, reconstruct its first
   K components, and average anti-diagonals back to a time series.
2. Initialize missing samples to zero. Starting at K=1, update missing values
   until the relative change is below the tolerance, then increase K toward the
   requested rank. Observed values stay fixed throughout.
3. Stage A uses M=24 and K=10 and fills pre-2017.5 gaps after excluding post-split
   missing samples, following your MATLAB method. This stage intentionally uses
   the compressed time sequence from the supplied routine.
4. Stage B uses M=60 and K=10 by default for remaining gaps. Use `--single-stage`
   to skip stage A. M=60/K=10 is a starting configuration, not a validated optimum.
5. `--cv` enables the MATLAB-style search over M=24:12:96 and
   K=[1,2,4,6,8,10,12], with held-out years 2004, 2007, 2010 and 2013 within
   2003–2016. Parameters are selected separately for each series needing stage B.
   Unlike the supplied MATLAB code, only genuinely observed values are scored:
   stage-A estimates are not used as validation truth. Candidates with failed
   convergence are excluded. No viable candidate leaves remaining gaps unfilled.

The default tolerance is 0.01 and the total iteration budget is 10000, matching
the MATLAB settings. Values are scaled internally to avoid problems with very
small magnitudes. This implementation handles a zero convergence denominator,
uses an economy SVD, and explicitly reports convergence. A stage A failure falls
back to stage B. Nonconverged stage B estimates are not written as filled data.
Unknown leading/trailing values are deliberately left missing; no extrapolation
is attempted. At least 24 observed samples are required by default.

This is a methodological adaptation, not a numerical-equivalence certification.
The MATLAB EOF periodogram flag and residual-based `verror` are not reproduced;
they do not select the components used for gap filling in your supplied routine.
CV RMSE is provided when CV runs; it is not a confidence interval. No uncertainty
is invented for a fixed-parameter fill. SSA can reconstruct plausible cycles
without recovering the true missing observations. Long mission-transition gaps
need independent scientific validation, including trend and amplitude checks.

## Output variables

The original variable and other dataset variables are retained (with monthly
reindexing if requested). NetCDF container layout and packing are not retained
byte for byte; output is decoded/re-encoded and compressed. The new filled data
are float64, retaining the original physical units.

For an input variable `lwe_thickness`, the additions are:

| Variable | Meaning |
| --- | --- |
| `lwe_thickness_ssa` | Filled data; finite observed values unchanged |
| `lwe_thickness_ssa_flag` | 0 observed, 1 stage A fill, 2 stage B fill, 3 still missing |
| `lwe_thickness_ssa_status` | Per-location status: 0 no gaps, 1 success, 2 insufficient data or unknown endpoints, 3 nonconvergence, 4 CV failure |
| `lwe_thickness_ssa_window` | Stage B window, 0 if not selected |
| `lwe_thickness_ssa_rank` | Stage B rank, 0 if not selected |
| `lwe_thickness_ssa_cv_rmse` | CV RMSE in original units; NaN when CV was not run |

Window/rank describe the configured stage B values even when stage A filled all
interior gaps. Read flags to identify which stage actually supplied a value.
An exit code of zero means the file was written, not that every series could be
filled. The command prints status counts; inspect statuses 2, 3, and 4 and the
remaining-NaN mask before analysis. File history and attributes record settings.

## Use from Python

```python
import numpy as np
from ssa_gapfill import ssa_fill

# x: one regularly spaced monthly series, with gaps encoded as NaN
x = np.array([ ... ], dtype=float)  # supply your data
result = ssa_fill(x, window=60, rank=10)
if result.converged:
    filled = result.values
```

`ssa_fill` returns estimates even if convergence fails so API users can diagnose
them; always check `result.converged`. The NetCDF workflow leaves failed gaps
missing. `fill_series()` supplies the two-stage workflow and flags.

## Validation

```bash
python -m pip install pytest
python -m pytest -q
```

Tests exercise known seasonal signals, short and 11-month gaps, two-stage flags,
observed-value preservation, all-NaN masks, zeros, endpoints, failed convergence,
blocked-year CV, calendar regularization, nonstandard calendars, and an actual
NetCDF CLI round trip. They do not establish accuracy on your GRACE dataset.

Your `example_C305.mat` has only a 300x1 `ser` time array, rather than a two-column
time/value series. Your two uploaded text files contain year/month pairs only.
The MATLAB main script refers to `example_C30.mat`, which was not supplied.
Consequently, MATLAB-output comparisons and real-data accuracy validation remain
pending until a measured series or your NetCDF file is available.
