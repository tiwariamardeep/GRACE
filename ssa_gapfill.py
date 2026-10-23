"""Iterative SSA gap filling adapted from the MATLAB routines supplied by the user.

This is a methodological adaptation, not a claim of bitwise MATLAB equivalence.
Only missing samples are replaced. No spatial interpolation is performed.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import svd


@dataclass
class SSAResult:
    values: np.ndarray
    converged: bool
    iterations: int


def ssa_fill(x, window=60, rank=10, tolerance=0.01, max_iterations=10000):
    """Zero-initialized SSA, increasing rank from 1 to `rank` as in MATLAB.

    Convergence uses relative RMS change at missing positions. The iteration
    budget is shared across ranks. Observations are fixed at every iteration.
    An all-missing series is returned unchanged, with converged=False.
    Nonconverged estimates are returned to API callers for diagnosis; the
    NetCDF/two-stage workflow below leaves those gaps unfilled instead.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1:
        raise ValueError('Input must be one-dimensional')
    n = x.size
    if not 2 <= window <= n:
        raise ValueError(f'window must be between 2 and series length ({n})')
    if not 1 <= rank <= min(window, n - window + 1):
        raise ValueError('rank exceeds the trajectory matrix rank limit')
    if tolerance <= 0 or max_iterations < 1:
        raise ValueError('tolerance and max_iterations must be positive')
    if np.isinf(x).any():
        raise ValueError('Convert infinite values to NaN before filling')
    missing = np.isnan(x)
    if not missing.any():
        return SSAResult(x.copy(), True, 0)
    if missing.all():
        return SSAResult(x.copy(), False, 0)
    # Scaling improves numerical conditioning for tiny spherical harmonics.
    scale = max(float(np.median(np.abs(x[~missing]))), 1e-14)
    original = x / scale
    working = np.where(missing, 0.0, original)
    k = 1
    rows = n - window + 1
    counts = np.convolve(np.ones(rows), np.ones(window))
    for iteration in range(1, max_iterations + 1):
        trajectory = np.lib.stride_tricks.sliding_window_view(working, window)
        u, s, vh = svd(trajectory, full_matrices=False, check_finite=False)
        low_rank = (u[:, :k] * s[:k]) @ vh[:k]
        reconstructed = np.zeros(n)
        # Average anti-diagonals (Hankelization), matching the MATLAB RC sum.
        for j in range(window):
            reconstructed[j:j + rows] += low_rank[:, j]
        reconstructed /= counts
        old = working[missing].copy()
        working[missing] = reconstructed[missing]
        denominator = max(np.linalg.norm(old), np.finfo(float).eps)
        change = np.linalg.norm(working[missing] - old) / denominator
        if change < tolerance:
            if k == rank:
                filled = x.copy()
                filled[missing] = working[missing] * scale
                return SSAResult(filled, True, iteration)
            k += 1
    filled = x.copy()
    filled[missing] = working[missing] * scale
    return SSAResult(filled, False, max_iterations)


def select_parameters(time, x, windows=(24, 36, 48, 60, 72, 84, 96),
                      ranks=(1, 2, 4, 6, 8, 10, 12), **solver_options):
    """Blocked-year CV on 2003 <= t < 2017; hold out 2004, 2007, 2010, 2013.

    Unlike the MATLAB routine, score *only originally observed* withheld data.
    The mean of per-year MSEs weights each year equally. Return (M,K,RMSE).
    RMSE is a validation statistic, not a formal uncertainty interval.
    """
    time, x = np.asarray(time), np.asarray(x)
    use = (time >= 2003) & (time < 2017)
    t, y = time[use], x[use]
    folds = []
    for year in (2004, 2007, 2010, 2013):
        holdout = (np.floor(t) == year) & np.isfinite(y)
        if holdout.any():
            folds.append(holdout)
    if not folds:
        raise ValueError('No observed holdout samples in the CV years')
    best = None
    # K outermost retains MATLAB's first-minimum tie ordering.
    for k in ranks:
        for m in windows:
            if not 2 <= m <= len(y) or not 1 <= k <= min(m, len(y)-m+1):
                continue
            errors = []
            for holdout in folds:
                trial = y.copy()
                trial[holdout] = np.nan
                fit = ssa_fill(trial, m, k, **solver_options)
                if not fit.converged:
                    break
                errors.append(np.mean((fit.values[holdout] - y[holdout])**2))
            if len(errors) == len(folds):
                score = float(np.sqrt(np.mean(errors)))
                if best is None or score < best[2]:
                    best = (m, k, score)
    if best is None:
        raise ValueError('No valid, converged cross-validation candidate')
    return best


def fill_series(time, x, window=60, rank=10, two_stage=True,
                stage_a_window=24, stage_a_rank=10, split=2017.5,
                min_observations=24, cv=False, cv_windows=(24,36,48,60,72,84,96),
                cv_ranks=(1,2,4,6,8,10,12), tolerance=0.01,
                max_iterations=10000):
    """Return values, fill flags, convergence status, chosen M/K, and CV RMSE.

    Flags: 0=observed, 1=stage A, 2=stage B, 3=unfilled.
    Status: 0=no gaps, 1=success, 2=insufficient data,
    3=nonconvergence, 4=no valid CV candidate. No endpoint extrapolation.
    """
    t, x = np.asarray(time), np.asarray(x, dtype=float)
    if t.shape != x.shape or x.ndim != 1:
        raise ValueError('Time and values must be equal-length 1D arrays')
    good = np.isfinite(x)
    if np.isinf(x).any():
        raise ValueError('Infinite values are not supported')
    output = x.copy()
    flags = np.where(good, 0, 3).astype(np.int8)
    if good.all():
        return output, flags, np.int8(0), np.int16(0), np.int16(0), np.nan
    if good.sum() < min_observations:
        return output, flags, np.int8(2), np.int16(0), np.int16(0), np.nan
    indices = np.flatnonzero(good)
    interior = np.zeros(len(x), dtype=bool)
    interior[indices[0]:indices[-1]+1] = True
    target = ~good & interior
    options = dict(tolerance=tolerance, max_iterations=max_iterations)
    if two_stage and np.any(target & (t < split)):
        # MATLAB excludes post-split gaps before stage A embedding.
        keep = good | ((t < split) & interior)
        fit = None
        if stage_a_window <= keep.sum() and stage_a_rank <= min(stage_a_window, keep.sum()-stage_a_window+1):
            fit = ssa_fill(output[keep], stage_a_window, stage_a_rank, **options)
        if fit is not None and fit.converged:
            locations = np.flatnonzero(keep)
            assign = target[locations] & (t[locations] < split)
            output[locations[assign]] = fit.values[assign]
            flags[locations[assign]] = 1
    score = np.nan
    if cv and np.any(np.isnan(output) & target):
        # Do not treat stage-A estimates as validation truth.
        try:
            window, rank, score = select_parameters(t, x, cv_windows, cv_ranks, **options)
        except ValueError:
            return output, flags, np.int8(4), np.int16(0), np.int16(0), score
    remaining = np.isnan(output) & target
    status = 1
    if remaining.any():
        # Trim unknown endpoints; SSA here is for interpolation only.
        start, stop = indices[0], indices[-1]+1
        segment = output[start:stop]
        if window > len(segment) or rank > min(window, len(segment)-window+1):
            return output, flags, np.int8(2), np.int16(window), np.int16(rank), score
        fit = ssa_fill(segment, window, rank, **options)
        if fit.converged:
            local = remaining[start:stop]
            locations = np.arange(start, stop)[local]
            output[locations] = fit.values[local]
            flags[locations] = 2
        else:
            status = 3
    if status == 1 and np.any(flags == 3):
        status = 2  # endpoints deliberately remain missing
    return output, flags, np.int8(status), np.int16(window), np.int16(rank), score


def process_series(time, x, extrapolate=False, max_extrapolation=12,
                   validate_endpoints=False, **options):
    """Optional iterative low-rank SSA completion of leading/trailing gaps.

    Returns the six fill_series outputs plus backward/forward holdout RMSE.
    Extrapolated values are estimates, flagged 4 (backward) or 5 (forward).
    Holdout diagnostics refit raw data without the withheld endpoint samples.
    """
    values, flags, status, m, k, score = fill_series(time, x, **options)
    errors = [np.nan, np.nan]
    raw = np.asarray(x, dtype=float)
    known = np.flatnonzero(np.isfinite(raw))
    if not extrapolate or known.size < options.get('min_observations', 24):
        return values, flags, status, m, k, score, *errors
    start, stop = int(known[0]), int(known[-1])+1
    left, right = start, len(raw)-stop
    if not (left or right) or not np.isfinite(values[start:stop]).all():
        return values, flags, status, m, k, score, *errors
    window = int(m) if m else options.get('window', 60)
    rank = int(k) if k else options.get('rank', 10)
    solver = dict(tolerance=options.get('tolerance', 0.01),
                  max_iterations=options.get('max_iterations', 10000))
    # Each side is bounded independently; long masked margins remain missing.
    extend_left = left if left <= max_extrapolation else 0
    extend_right = right if right <= max_extrapolation else 0
    if not (extend_left or extend_right):
        return values, flags, status, m, k, score, *errors
    a, b = start-extend_left, stop+extend_right
    if window > b-a or rank > min(window, b-a-window+1):
        return values, flags, np.int8(2), m, k, score, *errors
    fit = ssa_fill(values[a:b], window, rank, **solver)
    if not fit.converged:
        return values, flags, np.int8(5), m, k, score, *errors
    if extend_left:
        values[:start] = fit.values[:extend_left]
        flags[:start] = 4
    if extend_right:
        values[stop:] = fit.values[-extend_right:]
        flags[stop:] = 5
    status = np.int8(2 if np.any(flags == 3) else 1)
    if validate_endpoints:
        # Independent tests: withhold h calendar months from either end of
        # the original observed span. Score only measured holdout values.
        segment = raw[start:stop]
        for side, horizon in enumerate((extend_left, extend_right)):
            if not horizon:
                continue
            mask = np.zeros(len(segment), dtype=bool)
            mask[:horizon] = side == 0
            if side == 1:
                mask[-horizon:] = True
            scored = mask & np.isfinite(segment)
            trial = segment.copy()
            trial[mask] = np.nan
            if np.isfinite(trial).sum() < options.get('min_observations', 24):
                continue
            if window > len(trial) or rank > min(window, len(trial)-window+1):
                continue
            check = ssa_fill(trial, window, rank, **solver)
            if check.converged and scored.any():
                errors[side] = float(np.sqrt(np.mean((check.values[scored]-segment[scored])**2)))
    return values, flags, status, m, k, score, *errors


def monthly_axis(ds, time_name, regularize=False):
    """Validate monthly sampling; optionally insert absent calendar months.

    Regularization labels each month on day 1 and reindexes the entire Dataset.
    Existing samples are assigned by their calendar month, never interpolated.
    """
    import xarray as xr
    time = ds[time_name]
    if time.dims != (time_name,):
        raise ValueError('Time must be a one-dimensional coordinate')
    try:
        years, months = time.dt.year.values, time.dt.month.values
    except (AttributeError, TypeError) as exc:
        raise ValueError('Time must decode as CF dates; inspect time units/calendar') from exc
    ids = years * 12 + months - 1
    if len(ids) < 2 or not np.all(np.diff(ids) > 0):
        raise ValueError('Time must be sorted with exactly one sample per calendar month')
    if regularize:
        first, last = time.values[0], time.values[-1]
        encoding, attrs = dict(time.encoding), dict(time.attrs)
        if np.issubdtype(time.dtype, np.datetime64):
            labels = np.array([f'{y:04d}-{m:02d}-01' for y,m in zip(years,months)], dtype='datetime64[ns]')
            full = xr.date_range(str(labels[0])[:10], str(labels[-1])[:10], freq='MS')
        else:
            cls = type(first)
            labels = [cls(int(y), int(m), 1) for y,m in zip(years,months)]
            full = xr.date_range(labels[0], cls(last.year,last.month,1), freq='MS',
                                 calendar=time.dt.calendar, use_cftime=True)
        ds = ds.assign_coords({time_name: labels}).reindex({time_name: full})
        ds[time_name].attrs = attrs
        ds[time_name].encoding = {k:v for k,v in encoding.items() if k in ('units','calendar','dtype')}
        years, months = ds[time_name].dt.year.values, ds[time_name].dt.month.values
    elif np.any(np.diff(ids) != 1):
        raise ValueError('Absent monthly records detected: rerun with --regularize-months')
    return ds, years + (months - 0.5)/12.0


def run_netcdf(args):
    import dask
    import xarray as xr
    from dask.diagnostics import ProgressBar
    from threadpoolctl import threadpool_limits
    source, destination = Path(args.input), Path(args.output)
    if source.resolve() == destination.resolve() or destination.exists():
        raise ValueError('Choose a new output filename; existing files are never overwritten')
    with xr.open_dataset(source, engine='netcdf4', chunks={}, decode_times=True) as original:
        if args.variable not in original.data_vars:
            raise ValueError(f'Unknown variable {args.variable!r}; available: {list(original.data_vars)}')
        subset = {}
        for selection in args.isel:
            try:
                dim, interval = selection.split('=')
                start, stop = map(int, interval.split(':'))
            except ValueError as exc:
                raise ValueError('--isel must use DIM=START:STOP, e.g. lat=300:302') from exc
            if dim == args.time or dim not in original.dims or not 0 <= start < stop <= original.sizes[dim]:
                raise ValueError('Subset must be a valid non-time dimension index range')
            subset[dim] = slice(start, stop)
        ds, t = monthly_axis(original.isel(subset), args.time, args.regularize_months)
        data = ds[args.variable]
        if args.time not in data.dims or not np.issubdtype(data.dtype, np.number):
            raise ValueError('Selected variable must be numeric and contain the time dimension')
        # CF _FillValue/missing_value are decoded automatically by xarray.
        data = data.astype(np.float64)
        if args.missing_below is not None:
            data = data.where(data > args.missing_below)
        data = data.where(np.isfinite(data))
        data = data.chunk({dim: (-1 if dim == args.time else args.spatial_chunk) for dim in data.dims})
        time = xr.DataArray(t, dims=[args.time], coords={args.time: ds[args.time]})
        options = dict(window=args.window, rank=args.rank, two_stage=not args.single_stage,
                       stage_a_window=args.stage_a_window, stage_a_rank=args.stage_a_rank,
                       min_observations=args.min_observations, cv=args.cv,
                       tolerance=args.tolerance, max_iterations=args.max_iterations,
                       extrapolate=args.extrapolate, max_extrapolation=args.max_extrapolation,
                       validate_endpoints=args.validate_endpoints)
        values, flags, status, m, k, rmse, backward_rmse, forward_rmse = xr.apply_ufunc(
            process_series, time, data,
            input_core_dims=[[args.time], [args.time]],
            output_core_dims=[[args.time], [args.time], [], [], [], [], [], []],
            vectorize=True, dask='parallelized', kwargs=options,
            output_dtypes=[np.float64, np.int8, np.int8, np.int16, np.int16, np.float64, np.float64, np.float64])
        prefix = args.variable + '_ssa'
        names = [prefix, prefix+'_flag', prefix+'_status', prefix+'_window', prefix+'_rank', prefix+'_cv_rmse',
                 prefix+'_backward_holdout_rmse', prefix+'_forward_holdout_rmse']
        if any(name in ds.variables for name in names):
            raise ValueError('Output SSA variable names already exist in the input')
        out = ds.copy()
        out[names[0]] = values.transpose(*data.dims)
        out[names[0]].attrs = {key:value for key,value in ds[args.variable].attrs.items()
                              if key not in ('_FillValue','missing_value','valid_min','valid_max','valid_range','actual_range','ancillary_variables')}
        out[names[0]].attrs.update(long_name=f'SSA gap-filled {args.variable}',
            ancillary_variables=' '.join(names[1:]), ssa_window=args.window, ssa_rank=args.rank,
            ssa_two_stage=int(not args.single_stage), ssa_cv=int(args.cv),
            ssa_tolerance=args.tolerance, ssa_max_iterations=args.max_iterations,
            ssa_stage_a_window=args.stage_a_window, ssa_stage_a_rank=args.stage_a_rank,
            ssa_min_observations=args.min_observations,
            ssa_extrapolate=int(args.extrapolate), ssa_max_extrapolation=args.max_extrapolation,
            ssa_validate_endpoints=int(args.validate_endpoints),
            comment='Observed samples preserved. Flags 4/5 are endpoint estimates, not observations; inspect holdout diagnostics.')
        out[names[1]] = flags.transpose(*data.dims)
        out[names[1]].attrs.update(flag_values=np.arange(6,dtype=np.int8),
            flag_meanings='observed filled_stage_a filled_stage_b unfilled backward_extrapolation forward_extrapolation')
        out[names[2]] = status
        out[names[2]].attrs.update(flag_values=np.arange(6,dtype=np.int8),
            flag_meanings='no_gaps success insufficient_data nonconvergence no_cv_candidate endpoint_nonconvergence')
        out[names[3]], out[names[4]], out[names[5]] = m, k, rmse
        out[names[3]].attrs['long_name'] = 'Stage B SSA window (0 if not selected)'
        out[names[4]].attrs['long_name'] = 'Stage B SSA rank (0 if not selected)'
        out[names[5]].attrs.update(long_name='Blocked-year cross-validation RMSE; NaN when not run',
                                   units=ds[args.variable].attrs.get('units',ds[args.variable].attrs.get('Units','')))
        for name, diagnostic, direction in zip(names[6:], [backward_rmse, forward_rmse], ['Backward','Forward']):
            out[name] = diagnostic
            out[name].attrs.update(long_name=f'{direction} endpoint holdout RMSE',
                units=ds[args.variable].attrs.get('units',ds[args.variable].attrs.get('Units','')),
                comment='Single-stage SSA refit withholding the same number of endpoint calendar months as extrapolated. NaN means not run or unavailable. Not a confidence interval.')
        history = ds.attrs.get('history', '')
        out.attrs['history'] = history + '\n' + datetime.now(timezone.utc).isoformat() + ' SSA gap filling; original variable retained'
        if subset:
            out.attrs['history'] += '; index subset: ' + ', '.join(args.isel)
            out.attrs['ssa_subset_note'] = 'Global source coverage attributes describe the parent dataset; consult output coordinates for this subset.'
        if args.regularize_months:
            out.attrs['history'] += '; time labels normalized to month start, absent months inserted in all time-dependent variables'
            # Reindexing invalidates time-bounds references at inserted months.
            # Retain those bounds as missing rather than inventing intervals.
        # Avoid restoring integer packing after reindexing introduces NaNs.
        for name in out.variables:
            out[name].encoding = {}
        encoding = {name: dict(zlib=True, complevel=4) for name in out.data_vars}
        encoding[names[0]].update(dtype='float64', _FillValue=np.nan)
        temporary = destination.with_name(destination.name + '.partial')
        if temporary.exists():
            raise ValueError(f'Remove or rename stale temporary output first: {temporary}')
        print(f'Processing {dict(data.sizes)}; output: {destination}', flush=True)
        try:
            with threadpool_limits(limits=1), dask.config.set(scheduler='threads', num_workers=args.workers), ProgressBar():
                out.to_netcdf(temporary, engine='netcdf4', encoding=encoding)
            temporary.rename(destination)
        except BaseException:
            if temporary.exists():
                temporary.unlink()
            raise
    with xr.open_dataset(destination) as check:
        counts = check[names[2]].values
        unique, number = np.unique(counts, return_counts=True)
        print('Series status counts:', dict(zip(unique.tolist(), number.tolist())))
        print('0=no gaps, 1=success, 2=insufficient data/endpoints, 3=nonconvergence, 4=CV failure, 5=endpoint nonconvergence')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', help='Input NetCDF file')
    p.add_argument('output', help='New output NetCDF file')
    p.add_argument('--variable', required=True, help='Data variable, e.g. lwe_thickness; confirm with ncdump -h')
    p.add_argument('--time', default='time', help='Time dimension and coordinate name')
    p.add_argument('--window', type=int, default=60)
    p.add_argument('--rank', type=int, default=10)
    p.add_argument('--stage-a-window', type=int, default=24)
    p.add_argument('--stage-a-rank', type=int, default=10)
    p.add_argument('--single-stage', action='store_true')
    p.add_argument('--cv', action='store_true', help='Expensive per-series blocked-year parameter selection')
    p.add_argument('--extrapolate', action='store_true', help='Estimate leading/trailing NaNs using iterative SSA; flags 4/5')
    p.add_argument('--max-extrapolation', type=int, default=12, help='Maximum endpoint gap per side in months (default 12)')
    p.add_argument('--validate-endpoints', action='store_true', help='With --extrapolate, save independent endpoint holdout RMSE diagnostics')
    p.add_argument('--regularize-months', action='store_true', help='Normalize dates to month start and insert absent months across the dataset')
    p.add_argument('--missing-below', type=float, help='Also treat selected-variable values <= this threshold as missing')
    p.add_argument('--min-observations', type=int, default=24)
    p.add_argument('--tolerance', type=float, default=0.01)
    p.add_argument('--max-iterations', type=int, default=10000)
    p.add_argument('--spatial-chunk', type=int, default=8)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--isel', action='append', default=[], metavar='DIM=START:STOP',
                   help='Process only an index subset; repeat for multiple spatial dimensions (stop excluded)')
    args = p.parse_args()
    if args.max_extrapolation < 1:
        p.error('max-extrapolation must be positive')
    if args.validate_endpoints and not args.extrapolate:
        p.error('--validate-endpoints requires --extrapolate')
    if args.workers < 1 or args.spatial_chunk < 1 or args.min_observations < 2:
        p.error('workers/chunks must be positive and min-observations >= 2')
    if args.window < 2 or args.rank < 1 or args.rank > args.window:
        p.error('Require window >= 2 and 1 <= rank <= window')
    if args.stage_a_window < 2 or not 1 <= args.stage_a_rank <= args.stage_a_window:
        p.error('Invalid stage A window/rank')
    if args.tolerance <= 0 or args.max_iterations < 1:
        p.error('tolerance/max-iterations must be positive')
    run_netcdf(args)


if __name__ == '__main__':
    main()
