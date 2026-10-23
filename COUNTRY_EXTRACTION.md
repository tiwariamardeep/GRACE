# Extract monthly country averages to text

`extract_countries.py` writes one tab-separated `.txt` file per country from
`lwe_thickness_ssa` in the completed NetCDF output. It produces monthly means,
not a list of every grid-cell value and not a sum of water-storage volume.
Physical units are unchanged: your CSR file uses centimetres.

## Install and update on Derecho

```bash
conda activate gdal_env
cd /glade/campaign/univ/umsu0014/tiwaria6/Other/Data/GRACE/GRACE_SSA_python
git pull --ff-only
conda install -c conda-forge "python=3.10" geopandas rasterio xarray netcdf4 numpy pandas
```

Review Conda's changes before accepting, especially changes to GDAL. This reuses
the Python 3.10 environment. pip users can install `requirements-countries.txt`;
on Derecho prefer the consistent conda-forge GDAL/geospatial stack.

## Identify the country attribute

Keep `.shp`, `.shx`, `.dbf`, and `.prj` files together. Supply your actual path:

```bash
python extract_countries.py --shapefile /path/to/world_countries.shp --list-fields
```

Find the field containing country names, such as `ADMIN`, `NAME`, or `NAME_EN`.
Multiple polygons with the same value are combined as one country. Distinct
countries must have distinct values. The code uses the boundaries and country
definitions provided in your shapefile; it does not impose its own definitions.

## Extract

For a country name field `ADMIN`:

```bash
python extract_countries.py \
  CSR_SSA_process_full.nc \
  --shapefile /path/to/world_countries.shp \
  --country-field ADMIN \
  --variable lwe_thickness_ssa \
  --output-dir country_text
```

Use the actual SSA output filename. Run within your site's compute allocation
for the full global file. This extractor uses one process; it does not need the
32-worker SSA allocation. Its default time block is 12 months; use
`--time-chunk 1` to reduce peak data-array memory if necessary.

## Outputs

- A country file such as `001_Afghanistan.txt` with columns `month`, `mean`,
  `valid_cells`, `selected_cells`, and `valid_area_fraction`.
- `all_countries_monthly.txt` containing all countries in long-table format.
- `country_index.txt` mapping the original country names to filenames and
  reporting selected grid-cell counts/areas.
- `metadata.json` recording units, source variable, paths, and method.

Country filenames have numeric prefixes to avoid collisions after sanitizing
names. Existing output directories are never overwritten. Date labels follow
the input, so a January 2002–December 2026 file produces 300 rows per country.
Input records must be sorted and have one record per calendar month. Absent
records are not invented by this extractor.

## Spatial method and limitations

- Supports a regular one-dimensional lat/lon grid, including your 0.25-degree
  CSR grid. It handles ascending/descending latitude and converts 0..360
  longitude to -180..180 while reordering values with coordinates.
- Reprojects country boundaries to EPSG:4326 using the shapefile's declared CRS.
  A missing CRS is an error; its identity is never guessed.
- Rasterizes countries once, selecting cells by their centres. Boundary cells
  are not clipped to fractional polygon areas. This is an approximation for
  national means, especially for small countries, islands, and coastal regions.
- Weights selected cells by their spherical grid-cell area. Equal weights on a
  latitude/longitude grid would overweight high-latitude cells.
- Mean = sum(value × selected valid cell area) / sum(selected valid cell area).
  NaNs and infinities are excluded. A zero-valued observation remains valid.
- `valid_area_fraction` is the valid selected-cell area divided by all selected
  grid-cell area for that country, not the fraction of its exact polygon area.
  A country with selected cells but no finite data has mean NaN and coverage 0.
- Countries with no grid-cell centres have NaN means, zero selected cells, and
  NaN coverage. They are reported rather than assigned invented estimates.
- Overlapping country polygons that select the same pixel cause an error; this
  avoids silent double counting. Invalid/empty geometries are rejected.
- Global polygons must be correctly split at the antimeridian in the -180..180
  geographic convention; the code does not repair dateline topology.
- The SSA variable includes reconstructed/extrapolated values. The extractor
  does not exclude them or estimate uncertainty in the country mean.
- A 0.25-degree output grid is not independent 0.25-degree GRACE information.
  Country averaging does not remove GRACE leakage or resolve native mascon
  resolution limitations.

Tests use synthetic polygons and NetCDF arrays to verify spherical weighting,
missing-data coverage, small-country handling, longitude wrapping, reprojection,
overlap rejection, and text-file output. Your shapefile and actual NetCDF data
have not been run in the development environment.
