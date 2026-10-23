from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

gpd = pytest.importorskip('geopandas')
pytest.importorskip('rasterio')

from shapely.geometry import box
import xarray as xr

from extract_countries import country_mask,summarize


def test_area_weighting_missing_and_empty_country():
    countries=gpd.GeoDataFrame({'NAME':['A','B']},
        geometry=[box(-1,0,1,90),box(40,0,41,1)],crs='EPSG:4326')
    lat=np.array([67.5,22.5])
    lon=np.array([-0.5,0.5])
    names,labels,areas=country_mask(countries,'NAME',lat,lon)
    data=np.array([[10.,10.],[0.,0.]])
    mean,count,coverage=summarize(data,labels,areas,len(names))
    expected=10*(1-np.sin(np.pi/4))
    np.testing.assert_allclose(mean[0],expected)
    assert count[0]==4 and coverage[0]==1
    assert np.isnan(mean[1]) and np.isnan(coverage[1]) and count[1]==0
    data[0,:]=np.nan
    mean,count,coverage=summarize(data,labels,areas,len(names))
    assert mean[0]==0 and count[0]==2
    np.testing.assert_allclose(coverage[0],np.sin(np.pi/4))


def test_overlap_rejected():
    countries=gpd.GeoDataFrame({'NAME':['A','B']},
        geometry=[box(-1,0,1,90),box(-1,0,1,90)],crs='EPSG:4326')
    with pytest.raises(ValueError,match='overlap'):
        country_mask(countries,'NAME',np.array([67.5,22.5]),np.array([-0.5,0.5]))


def test_cli_longitude_wrap_and_country_files(tmp_path):
    # Ascending latitudes and 0..360 longitudes exercise both reorderings.
    data=np.array([[[10.,20.],[30.,40.]],[[np.nan,50.],[np.nan,70.]]])
    ds=xr.Dataset({'lwe_thickness_ssa':(('time','lat','lon'),data,{'Units':'cm'})},
        coords={'time':pd.date_range('2002-01-01',periods=2,freq='MS'),
                'lat':[-0.5,0.5], 'lon':[0.5,359.5]})
    nc=tmp_path/'in.nc'
    ds.to_netcdf(nc)
    shapes=gpd.GeoDataFrame({'NAME':['West','East','Tiny']},
        geometry=[box(-1,-1,0,1),box(0,-1,1,1),box(80,0,80.01,0.01)],crs='EPSG:4326')
    shp=tmp_path/'countries.shp'
    shapes.to_file(shp)
    out=tmp_path/'country_text'
    result=subprocess.run([sys.executable,str(Path(__file__).with_name('extract_countries.py')),
        str(nc),'--shapefile',str(shp),'--country-field','NAME','--output-dir',str(out)],
        capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    index=pd.read_csv(out/'country_index.txt',sep='\t')
    west=pd.read_csv(out/index.set_index('country').loc['West','file'],sep='\t')
    east=pd.read_csv(out/index.set_index('country').loc['East','file'],sep='\t')
    np.testing.assert_allclose(west['mean'],[30,60])
    assert east.loc[0,'mean']==20 and np.isnan(east.loc[1,'mean'])
    assert east.loc[1,'valid_area_fraction']==0
    assert list(west['month'])==['2002-01','2002-02']
    combined=pd.read_csv(out/'all_countries_monthly.txt',sep='\t')
    assert len(combined)==6
    assert '"units": "cm"' in (out/'metadata.json').read_text()


def test_projected_shapefile():
    countries=gpd.GeoDataFrame({'NAME':['A']},geometry=[box(-1,0,1,90)],crs='EPSG:4326')
    # Limit northern latitude for projection to avoid the Web Mercator pole.
    countries.geometry=[box(-1,0,1,80)]
    projected=countries.to_crs('EPSG:3857')
    args=('NAME',np.array([67.5,22.5]),np.array([-0.5,0.5]))
    a=country_mask(countries,*args)
    b=country_mask(projected,*args)
    np.testing.assert_array_equal(a[1],b[1])
