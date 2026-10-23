"""Monthly country means from a regular latitude/longitude GRACE NetCDF grid.

Country masks use pixel centres, not fractional border-cell coverage. Means
use spherical grid-cell areas and only finite data within each monthly mask.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.features
from rasterio.enums import MergeAlg
from rasterio.transform import from_origin
import xarray as xr


def country_mask(countries, field, latitude, longitude):
    """Return names, integer country labels and spherical cell areas (km²)."""
    if countries.crs is None:
        raise ValueError('Shapefile has no CRS. Supply its correct .prj file; do not guess the CRS.')
    if field not in countries.columns:
        raise ValueError(f'Unknown country field {field!r}; available: {list(countries.columns)}')
    if countries[field].isna().any():
        raise ValueError('Country name field has missing values; select a complete identifier field.')
    geometry_column = countries.geometry.name
    countries = countries[[field, geometry_column]].copy().to_crs('EPSG:4326')
    countries[field] = countries[field].astype(str)
    if countries.geometry.isna().any() or countries.geometry.is_empty.any():
        raise ValueError('Shapefile includes missing/empty geometries.')
    if not countries.geometry.is_valid.all():
        raise ValueError('Shapefile includes invalid geometries; repair them before extraction.')
    if not countries.geometry.geom_type.isin(['Polygon','MultiPolygon']).all():
        raise ValueError('Country boundaries must be Polygon or MultiPolygon geometries.')
    # Combine multiple records belonging to the same named country.
    countries = countries.dissolve(by=field, sort=True)
    names = countries.index.astype(str).to_numpy()
    lat, lon = np.asarray(latitude), np.asarray(longitude)
    if len(lat) < 2 or len(lon) < 2:
        raise ValueError('At least two latitude and longitude coordinates are required.')
    dy, dx = -np.diff(lat), np.diff(lon)
    if not (np.all(dy > 0) and np.all(dx > 0)
            and np.allclose(dy,dy[0],rtol=1e-4,atol=1e-7)
            and np.allclose(dx,dx[0],rtol=1e-4,atol=1e-7)):
        raise ValueError('Expected a regular grid with descending latitude and ascending longitude.')
    dy, dx = float(dy[0]), float(dx[0])
    if np.max(np.abs(lat)) > 90 or lon.min() < -180 or lon.max() >= 180:
        raise ValueError('Coordinates must be geographic degrees after longitude normalization.')
    transform = from_origin(float(lon[0]-dx/2),float(lat[0]+dy/2),dx,dy)
    shape = (len(lat),len(lon))
    shapes = [(geom,i+1) for i,geom in enumerate(countries.geometry)]
    labels = rasterio.features.rasterize(shapes,out_shape=shape,transform=transform,
                                         fill=0,all_touched=False,dtype='int32')
    # Do not silently double-count or arbitrarily assign overlapping countries.
    overlaps = rasterio.features.rasterize([(geom,1) for geom in countries.geometry],
        out_shape=shape,transform=transform,fill=0,all_touched=False,
        merge_alg=MergeAlg.add,dtype='int32')
    if np.any(overlaps > 1):
        raise ValueError('Country polygons overlap at grid-cell centres. Resolve overlapping country boundaries first.')
    north = np.minimum(lat+dy/2,90)
    south = np.maximum(lat-dy/2,-90)
    row_area = 6371.0088**2 * np.deg2rad(dx) * (
        np.sin(np.deg2rad(north))-np.sin(np.deg2rad(south)))
    areas = np.broadcast_to(row_area[:,None],shape)
    return names, labels, areas


def summarize(values, labels, areas, number_countries):
    """Reduce one spatial grid to weighted country means and coverage."""
    ids = labels.ravel()
    weight = areas.ravel()
    x = np.asarray(values,dtype=np.float64).ravel()
    valid = (ids > 0) & np.isfinite(x)
    bins = number_countries+1
    sums = np.bincount(ids[valid],weights=x[valid]*weight[valid],minlength=bins)[1:]
    valid_area = np.bincount(ids[valid],weights=weight[valid],minlength=bins)[1:]
    valid_count = np.bincount(ids[valid],minlength=bins)[1:]
    total_area = np.bincount(ids,weights=weight,minlength=bins)[1:]
    means = np.divide(sums,valid_area,out=np.full(number_countries,np.nan),where=valid_area>0)
    coverage = np.divide(valid_area,total_area,out=np.full(number_countries,np.nan),where=total_area>0)
    return means,valid_count,coverage


def extract(args):
    output = Path(args.output_dir)
    if output.exists():
        raise FileExistsError(f'Choose a new output directory; {output} already exists.')
    countries = gpd.read_file(args.shapefile)
    with xr.open_dataset(args.input) as ds:
        if args.variable not in ds:
            raise ValueError(f'Unknown variable {args.variable!r}; available: {list(ds.data_vars)}')
        data = ds[args.variable]
        if set(data.dims) != {args.time,args.lat,args.lon} or data.ndim != 3:
            raise ValueError('Variable must have exactly time, latitude and longitude dimensions.')
        # CSR uses 0..360 longitude. Convert to the -180..180 polygon convention
        # and reorder data together with its coordinates.
        data = data.assign_coords({args.lon:((data[args.lon]+180)%360)-180})
        data = data.sortby(args.lon).sortby(args.lat,ascending=False)
        data = data.transpose(args.time,args.lat,args.lon)
        dates = data[args.time]
        # dt access also supports nonstandard cftime calendars.
        years,months = dates.dt.year.values,dates.dt.month.values
        month_ids = years*12+months
        if len(np.unique(month_ids)) != len(month_ids) or np.any(np.diff(month_ids)<=0):
            raise ValueError('Input must have sorted, unique monthly records.')
        date_labels = [f'{int(y):04d}-{int(m):02d}' for y,m in zip(years,months)]
        names,labels,areas = country_mask(countries,args.country_field,
                                         data[args.lat].values,data[args.lon].values)
        selected = np.bincount(labels.ravel(),minlength=len(names)+1)[1:]
        selected_area = np.bincount(labels.ravel(),weights=areas.ravel(),minlength=len(names)+1)[1:]
        ntime,ncountry = len(date_labels),len(names)
        means = np.full((ntime,ncountry),np.nan)
        coverage = np.full_like(means,np.nan)
        counts = np.zeros((ntime,ncountry),dtype=np.int64)
        for start in range(0,ntime,args.time_chunk):
            stop = min(start+args.time_chunk,ntime)
            block = data.isel({args.time:slice(start,stop)}).values
            for offset,values in enumerate(block):
                i = start+offset
                means[i],counts[i],coverage[i] = summarize(values,labels,areas,ncountry)
            print(f'Processed {stop}/{ntime} months',flush=True)
        units = data.attrs.get('units',data.attrs.get('Units','unspecified'))
        source_attrs = dict(data.attrs)
    # Only create outputs after the spatial and temporal reductions succeed.
    output.mkdir(parents=True)
    manifest=[]
    all_rows=[]
    for j,name in enumerate(names):
        safe = re.sub(r'[^A-Za-z0-9_-]+','_',name).strip('_')[:100] or 'country'
        filename=f'{j+1:03d}_{safe}.txt'
        table=pd.DataFrame({'month':date_labels,'mean':means[:,j],
            'valid_cells':counts[:,j],'selected_cells':selected[j],
            'valid_area_fraction':coverage[:,j]})
        table.to_csv(output/filename,sep='\t',index=False,na_rep='NaN',float_format='%.8g')
        manifest.append({'country':name,'file':filename,'selected_cells':int(selected[j]),
                         'selected_grid_area_km2':float(selected_area[j])})
        all_rows.append(table.assign(country=name)[['country',*table.columns]])
    pd.DataFrame(manifest).to_csv(output/'country_index.txt',sep='\t',index=False,float_format='%.8g')
    pd.concat(all_rows,ignore_index=True).to_csv(output/'all_countries_monthly.txt',
        sep='\t',index=False,na_rep='NaN',float_format='%.8g')
    metadata=dict(input=str(args.input),shapefile=str(args.shapefile),variable=args.variable,
        country_field=args.country_field,units=units,months=ntime,countries=ncountry,
        spatial_method='Grid-cell centre within polygon; no fractional boundary-cell clipping',
        weights='Spherical grid-cell area using radius 6371.0088 km',
        missing_data='Mean uses finite values only; coverage denominator is all selected grid-cell area',
        source_variable_attributes=source_attrs)
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str),encoding='utf-8')
    print(f'Wrote {ncountry} country files to {output}; units: {units}')
    empty=names[selected==0]
    if len(empty):
        print('No grid-cell centres selected for: '+', '.join(empty))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',nargs='?',help='SSA output NetCDF')
    parser.add_argument('--shapefile',required=True,help='Country polygon shapefile (with .shx/.dbf/.prj files)')
    parser.add_argument('--country-field',help='Attribute used to identify countries, e.g. ADMIN or NAME')
    parser.add_argument('--list-fields',action='store_true',help='Show polygon attribute names and exit')
    parser.add_argument('--variable',default='lwe_thickness_ssa')
    parser.add_argument('--output-dir',default='country_text')
    parser.add_argument('--time',default='time')
    parser.add_argument('--lat',default='lat')
    parser.add_argument('--lon',default='lon')
    parser.add_argument('--time-chunk',type=int,default=12,help='Months loaded at once (default 12)')
    args=parser.parse_args()
    if args.list_fields:
        countries=gpd.read_file(args.shapefile)
        print('Attribute fields:',[c for c in countries.columns if c!=countries.geometry.name])
        print(countries.drop(columns=countries.geometry.name).head().to_string(index=False))
        return
    if not args.input or not args.country_field:
        parser.error('Provide input NetCDF and --country-field, or use --list-fields first')
    if args.time_chunk<1:
        parser.error('--time-chunk must be positive')
    extract(args)


if __name__=='__main__':
    main()
