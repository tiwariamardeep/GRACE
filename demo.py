"""Run with: python demo.py. Creates a small synthetic NetCDF example."""
from pathlib import Path
import subprocess
import sys

import numpy as np
import xarray as xr
from threadpoolctl import threadpool_limits
from ssa_gapfill import ssa_fill


def main():
    root = Path(__file__).resolve().parent
    input_file, output_file = root/'demo_input.nc', root/'demo_filled.nc'
    if input_file.exists() or output_file.exists():
        raise SystemExit('Demo files already exist; rename/remove them before rerunning.')
    i = np.arange(300)
    truth = 2 + 3*np.sin(2*np.pi*i/12) + np.cos(2*np.pi*i/6)
    missing = truth.copy()
    missing[50:53] = np.nan
    missing[186:197] = np.nan
    with threadpool_limits(limits=1):
        fit = ssa_fill(missing,48,5,tolerance=1e-6)
    mask = np.isnan(missing)
    print(f'Synthetic gap RMSE: {np.sqrt(np.mean((fit.values[mask]-truth[mask])**2)):.6g}')
    print(f'Converged: {fit.converged}; iterations: {fit.iterations}')
    ds = xr.Dataset({'lwe_thickness': (('time','lat','lon'), missing[:,None,None], {'units':'cm'}),
                     'synthetic_truth': ('time',truth)},
                    coords={'time':xr.date_range('2002-01-01',periods=300,freq='MS'),
                            'lat':[25.0], 'lon':[80.0]})
    ds.to_netcdf(input_file)
    subprocess.run([sys.executable,str(root/'ssa_gapfill.py'),str(input_file),str(output_file),
                    '--variable','lwe_thickness','--window','48','--rank','5',
                    '--stage-a-rank','5','--tolerance','0.000001'],check=True)


if __name__ == '__main__':
    main()
