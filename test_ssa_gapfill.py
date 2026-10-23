import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from threadpoolctl import threadpool_limits

from ssa_gapfill import ssa_fill, fill_series, monthly_axis, select_parameters


@pytest.fixture(autouse=True)
def single_thread():
    with threadpool_limits(limits=1):
        yield


def signal(n=240):
    i = np.arange(n)
    return 2 + 3*np.sin(2*np.pi*i/12) + np.cos(2*np.pi*i/6)


def test_long_and_short_gaps():
    truth = signal()
    x = truth.copy()
    x[45:48] = np.nan
    x[150:161] = np.nan
    result = ssa_fill(x, 48, 5, tolerance=1e-6)
    assert result.converged
    np.testing.assert_array_equal(result.values[np.isfinite(x)], x[np.isfinite(x)])
    assert np.sqrt(np.mean((result.values[np.isnan(x)]-truth[np.isnan(x)])**2)) < 0.001


def test_zero_series_and_all_missing():
    x = np.zeros(100)
    x[40:50] = np.nan
    fit = ssa_fill(x,24,4)
    assert fit.converged
    np.testing.assert_array_equal(fit.values, np.zeros(100))
    fit = ssa_fill(np.full(100,np.nan),24,4)
    assert not fit.converged and np.isnan(fit.values).all()


def test_endpoints_and_failure_not_filled():
    x = signal()
    x[:3] = np.nan
    x[80:90] = np.nan
    t = 2002 + (np.arange(240)+0.5)/12
    out, flags, status, *_ = fill_series(t,x,window=48,rank=5,two_stage=False,tolerance=1e-6)
    assert np.isnan(out[:3]).all() and np.all(flags[:3] == 3)
    assert np.isfinite(out[80:90]).all() and status == 2
    out, flags, status, *_ = fill_series(t,x,window=48,rank=5,two_stage=False,max_iterations=1)
    assert np.isnan(out[80:90]).all() and status == 3


def test_two_stage_flags():
    t = 2002 + (np.arange(300)+0.5)/12
    x = signal(300)
    x[30:32] = np.nan
    x[186:197] = np.nan
    out, flags, status, *_ = fill_series(t,x,window=48,rank=5,stage_a_rank=5,tolerance=1e-5)
    assert status == 1
    assert (flags[30:32] == 1).all() and (flags[186:197] == 2).all()
    assert np.isfinite(out).all()


def test_cv():
    t = 2002 + (np.arange(240)+0.5)/12
    x = signal()
    x[70] = np.nan
    m,k,score = select_parameters(t,x,windows=[24,36],ranks=[5],tolerance=1e-5)
    assert m in (24,36) and k == 5 and score < 0.01


def test_month_validation_and_regularize():
    time = xr.date_range('2002-01-15', periods=6, freq='MS').delete(2)
    ds = xr.Dataset({'v': ('time',np.arange(5.))},coords={'time':time})
    with pytest.raises(ValueError,match='Absent monthly'):
        monthly_axis(ds,'time')
    out,t = monthly_axis(ds,'time',True)
    assert out.sizes['time'] == 6 and np.isnan(out.v.values[2])
    np.testing.assert_allclose(np.diff(t), 1/12)


def test_netcdf_end_to_end(tmp_path):
    truth = signal(240)
    data = np.tile(truth[:,None,None],(1,2,2))
    data[60:63,0,0] = np.nan
    data[185:196,0,0] = np.nan
    data[:,:,1] = np.nan  # permanently masked domain
    ds = xr.Dataset({'lwe_thickness': (('time','lat','lon'),data,{'units':'cm'}),
                     'land_mask': (('lat','lon'),[[1,0],[1,0]])},
                    coords={'time':xr.date_range('2002-01-01',periods=240,freq='MS'),
                            'lat':[10.,11.], 'lon':[70.,71.]})
    source, dest = tmp_path/'in.nc', tmp_path/'out.nc'
    ds.to_netcdf(source,encoding={'lwe_thickness':{'_FillValue':-9999.}})
    cli = Path(__file__).with_name('ssa_gapfill.py')
    result = subprocess.run([sys.executable,str(cli),str(source),str(dest),
                             '--variable','lwe_thickness','--window','48','--rank','5',
                             '--stage-a-rank','5','--tolerance','0.00001'],
                            capture_output=True,text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    with xr.open_dataset(dest) as out:
        np.testing.assert_array_equal(out.lwe_thickness.values,data)
        filled = out.lwe_thickness_ssa.values
        observed = np.isfinite(data)
        np.testing.assert_array_equal(filled[observed],data[observed])
        assert np.isfinite(filled[:,0,0]).all()
        assert np.isnan(filled[:,:,1]).all()
        assert out.lwe_thickness_ssa_status.values[0,0] == 1
        np.testing.assert_array_equal(out.land_mask.values,ds.land_mask.values)


def test_nonstandard_calendar():
    dates = xr.date_range('2002-01-15', periods=6, freq='MS', calendar='360_day', use_cftime=True).delete(2)
    ds = xr.Dataset({'v':('time',np.arange(5.))},coords={'time':dates})
    out,t = monthly_axis(ds,'time',True)
    assert out.sizes['time'] == 6 and out.time.dt.calendar == '360_day'


def test_invalid_parameters():
    with pytest.raises(ValueError):
        ssa_fill(signal(24),24,10)


def test_csr_header_subset_and_absent_month(tmp_path):
    time = xr.date_range('2002-04-01',periods=293,freq='MS')
    truth = signal(293)
    data = np.tile(truth[:,None,None],(1,3,3)).astype('float32')
    data[30:32,1,1] = np.nan
    ds = xr.Dataset({'lwe_thickness':(('time','lat','lon'),data,{'Units':'cm'})},
                    coords={'time':time,'lat':[-1.,0.,1.],'lon':[0.125,0.375,0.625]})
    ds = ds.isel(time=np.delete(np.arange(293),190))
    source,destination = tmp_path/'csr.nc',tmp_path/'filled.nc'
    ds.to_netcdf(source,encoding={'time':{'units':'days since 2002-04-01 00:00:00',
                                        'calendar':'proleptic_gregorian','dtype':'int64'},
                                 'lwe_thickness':{'_FillValue':np.nan}})
    result = subprocess.run([sys.executable,str(Path(__file__).with_name('ssa_gapfill.py')),
        str(source),str(destination),'--variable','lwe_thickness','--regularize-months',
        '--isel','lat=1:2','--isel','lon=1:2','--cv'],capture_output=True,text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    with xr.open_dataset(destination) as out:
        assert dict(out.sizes) == {'time':293,'lat':1,'lon':1}
        assert out.lwe_thickness_ssa.attrs['Units'] == 'cm'
        assert out.lwe_thickness_ssa_cv_rmse.attrs['units'] == 'cm'
        assert np.isnan(out.lwe_thickness.values[190]).all()
        assert np.isfinite(out.lwe_thickness_ssa.values[190]).all()
        assert np.isfinite(out.lwe_thickness_ssa_cv_rmse.values).all()
