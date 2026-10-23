# SSA gap filling for GRACE / GRACE-FO

This package adapts your `ssa_missing_iterative.m`, `fun_SSA_filling_a.m`, and
`fun_SSA_filling_b.m` to Python. It processes one monthly time series at each
spatial location, with bounded spatial chunks. It works with grids or mascon
dimensions, provided the chosen variable contains a decoded `time` dimension.
Python 3.12 is recommended, particularly on Windows where newer Python
versions may lack precompiled dependency wheels.

## Install on Windows

Download and extract this package. Open its folder in File Explorer, click the
address bar, type `cmd`, and press Enter. In that Command Prompt run:

```bat
py -3.12 -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
python demo.py
```

If `py` is not found, install Python from https://www.python.org/downloads/ first.
If the launcher reports no Python 3.12 runtime, install Python 3.12 first; the
Windows Python install manager supports `py install 3.12`. Stop after any failed
command rather than continuing with a different Python version.
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
Unknown leading/trailing values are left missing by default. The optional
endpoint procedure below estimates them explicitly. At least 24 observed
samples are required by default.

## Estimate exterior months: optional SSA extrapolation

First prepare an input file with one record for each month from January 2002 to
December 2026, with NaNs at absent months. The input must preserve the separate
CSR special-month solutions (November 2011 and May 2015), rather than averaging
their representative timestamps into October/April. This script cannot recover
observations previously averaged together.

Update `ssa_gapfill.py` from this repository, then run a pilot:

```bash
python ssa_gapfill.py "../CSR_GRACE_GRACE-FO_RL0603_Mascons_Jan2002_Dec2026_FillNaN.nc" CSR_SSA_endpoint_pilot.nc --variable lwe_thickness --extrapolate --max-extrapolation 4 --validate-endpoints --isel lat=300:302 --isel lon=300:302
```

Adjust the input path to the actual file location. After examining the pilot,
remove the two `--isel` arguments and use a new output filename for the full grid:

```bash
python ssa_gapfill.py "../CSR_GRACE_GRACE-FO_RL0603_Mascons_Jan2002_Dec2026_FillNaN.nc" CSR_SSA_Jan2002_Dec2026_full.nc --variable lwe_thickness --extrapolate --max-extrapolation 4 --validate-endpoints
```

`--extrapolate` runs an additional iterative low-rank SSA completion after
interior filling. Only unknown endpoints are assigned; interior results and
observations stay fixed. This is iterative SSA extrapolation, not a recurrent
SSA forecasting algorithm. Backward estimates get flag 4 and forward estimates
get flag 5. For a series observed from April 2002 to August 2026, it estimates
the three preceding and four following months. All-NaN locations remain missing.

`--max-extrapolation 4` permits at most four exterior months on each side. Each
side is checked independently. Longer exterior masks remain NaN instead of
being extended arbitrarily. The default limit is 12 months. This limit applies
to endpoint gaps only, not the internal GRACE/GRACE-FO mission gap.

`--validate-endpoints` runs additional independent, single-stage SSA fits to
the raw series with the same number of months withheld at either end of its
observed span. Predictions are compared only with genuinely measured holdout
values. This checks the extrapolation method without using previously filled
values as validation truth. Parameters are the selected/configured M/K; this is
a diagnostic, not an independent parameter search or a confidence interval.
The two saved RMSE variables are in the data's physical units. NaN means the
diagnostic was not requested, insufficient observations were available, or its
fit failed to converge. Validation failure does not automatically reject an
otherwise converged extrapolation; inspect the diagnostics before analysis.

The extra endpoint fits increase runtime. A written file or successful
convergence does not establish accurate extrapolation. Check endpoint holdout
RMSE, trend/seasonal plausibility, and the flags on representative land cells.
Original values remain in the original variable. Missing values elsewhere in
the dataset, including original-observation timestamp variables, are retained.

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
| `lwe_thickness_ssa_flag` | 0 observed, 1 stage A fill, 2 stage B fill, 3 still missing, 4 backward extrapolation, 5 forward extrapolation |
| `lwe_thickness_ssa_status` | Per-location status: 0 no gaps, 1 success, 2 insufficient data or unknown endpoints, 3 nonconvergence, 4 CV failure, 5 endpoint nonconvergence |
| `lwe_thickness_ssa_window` | Stage B window, 0 if not selected |
| `lwe_thickness_ssa_rank` | Stage B rank, 0 if not selected |
| `lwe_thickness_ssa_cv_rmse` | CV RMSE in original units; NaN when CV was not run |
| `lwe_thickness_ssa_backward_holdout_rmse` | Withheld-start diagnostic RMSE; NaN if unavailable/not requested |
| `lwe_thickness_ssa_forward_holdout_rmse` | Withheld-end diagnostic RMSE; NaN if unavailable/not requested |

Window/rank describe the configured stage B values even when stage A filled all
interior gaps. Read flags to identify which stage actually supplied a value.
An exit code of zero means the file was written, not that every series could be
filled. The command prints status counts; inspect statuses 2, 3, 4, and 5 and the
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
`process_series()` adds optional endpoint estimation and returns the six
`fill_series()` outputs followed by backward and forward holdout RMSE.

## Validation

```bash
python -m pip install pytest
python -m pytest -q
```

Tests exercise known seasonal signals, short and 11-month gaps, two-stage flags,
observed-value preservation, all-NaN masks, zeros, endpoints, failed convergence,
blocked-year CV, calendar regularization, nonstandard calendars, and actual
NetCDF CLI round trips. Endpoint tests include a known seasonal signal plus
trend, independent holdout refits, horizon limits, nonconvergence, and
separate backward/forward flags. They do not establish accuracy on your GRACE dataset.

Your `example_C305.mat` has only a 300x1 `ser` time array, rather than a two-column
time/value series. Your two uploaded text files contain year/month pairs only.
The MATLAB main script refers to `example_C30.mat`, which was not supplied.
Consequently, MATLAB-output comparisons and real-data accuracy validation remain
pending until a measured series or your NetCDF file is available.
