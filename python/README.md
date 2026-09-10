# moc_py — Python side of moc

**Model vs Observation Glacier Velocity Comparison Module**, Python
implementation. The MATLAB implementation lives one folder up; see
[`../README.md`](../README.md) for the project overview, how `moc` relates to
the [green-ibis](../../green-ibis) workflows, and the shared path
configuration.

Python-native readers and plots for **observed** glacier velocities, for use in
Jupyter. This is **Phase 1**: the two observation formats are read natively and
compared/plotted in Python. The ISSM **model** side stays in MATLAB for now (see
"Phase 2" below) because a MATLAB-saved ISSM `model` object cannot be read
directly in Python.

## Why this split

| Data | Reader | Notes |
|---|---|---|
| Obs gridded velocity (`.nc`) | `load_obs_netcdf` | Files are netCDF4 (**HDF5**); read via `h5py` — **no `netCDF4`/`rioxarray` needed** |
| Obs flowline (`.mat`) | `load_obs_flowline` | MATLAB **v7** files, read via `scipy.io.loadmat` |
| ISSM model (`.mat`) | *(Phase 2)* | v7.3 MATLAB **class object** — export to netCDF from MATLAB, then read here |

Both readers return **`xarray.Dataset`** objects with decimal-year time coords
(converted from CF "days since" / MATLAB datenum), so you get `.sel`, `.isel`,
slicing, and plotting for free.

## Install

Required: `numpy scipy h5py xarray matplotlib pandas`. Optional extras —
`shapefile` (`geopandas shapely pyogrio`) for `moc_py.flowline_shp`,
`interactive` (`ipympl jupyterlab`) for notebooks and ROI drawing, `lazy`
(`dask rioxarray zarr ...`) for cloud reads. Everything else imports without them.

**pixi** (the project convention):

```bash
pixi add numpy scipy h5py xarray matplotlib pandas
pixi add geopandas shapely pyogrio ipympl jupyterlab   # optional extras
pixi run sync-env                 # keep environment.yml in step for conda users
pixi run pip install -e python/
```

**conda / mamba:**

```bash
conda create -n moc python=3.11 && conda activate moc
conda install -c conda-forge numpy scipy h5py xarray matplotlib pandas \
                             geopandas shapely pyogrio ipympl jupyterlab
pip install -e python/
```

**pip:**

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e "python/[shapefile]"      # core + the shapefile extra
```

Editable (`-e`) is recommended while the toolkit is moving. Check it:

```bash
python -c "import moc_py; print(moc_py.__file__)"
python -m moc_py               # show configured data paths
python -m moc_py.smoke_test    # load one of each obs type
```

## Data paths

Nothing is hard-coded to a machine. Each location resolves on first use:
explicit argument → environment variable (`MOC_OBS_NETCDF_DIR`, ...) → the
shared config file `~/.config/moc/paths.json` → an interactive prompt that
remembers the answer. The MATLAB side reads and writes the same file, so
configuring either one configures both.

```bash
python -m moc_py --set          # prompt for whatever is missing
```

```python
import moc_py as mp
mp.configure(obs_netcdf_dir="/data/velocity")   # this session
mp.save_config()                                # ... and remember it
mp.CONFIG.status()                              # where each setting came from
```

Unconfigured and unanswerable (a script, CI), a `MissingPathError` names the
setting and the three ways to set it. Full details in
[`../README.md`](../README.md#configuring-data-paths).

## Quick start

```python
import moc_py as mp

O = mp.load_obs_netcdf("Sentinel_Subset_Upernavik.nc", bands=["vv","vx","vy"],
                       trange=(2018, 2020))    # xarray.Dataset (y, x, time)
F = mp.load_obs_flowline("jakobshavn.mat")     # xarray.Dataset (point, time)

mp.plot_map(O, band="vv", tindex=0)
mp.plot_flowline(F, mode="hovmoller")
mp.plot_timeseries(F, dist=5000)

roi = mp.roi_from_obs(O, margin=2000, name="upernavik")   # saves .json + .exp
```

Interactive ROI in a notebook (needs `%matplotlib widget`):

```python
sel = mp.RoiSelector(O, band="vv", kind="rectangle", name="upernavik")
# drag on the figure, then:
roi = sel.save()
```

See `notebooks/moc_python_quickstart.ipynb` for the full walkthrough.

## Flowlines from the shapefiles

`moc_py.flowline_shp` turns the two archived shapefile collections into a
flowline the comparison functions accept. Give it an ROI and it loops the whole
Felikson collection (`CONFIG.flowline_shp_dir`), keeps the flowlines that reach
into the ROI, takes the centre one, splits it at the most recent Black &
Joughin terminus trace that crosses it (`CONFIG.termini_shp`), keeps the
upstream side, and clips that to the ROI:

```python
roi = mp.load_roi("sverdrups_trunk")        # or a Roi, polygon, or bbox tuple

FL = mp.flowline_from_roi(roi, spacing=200)  # xarray.Dataset (point,) x/y/d
FL.attrs["terminus_date"]                    # '2021-12-21' — what it was clipped at
FL.attrs["flowline"]                         # line used; cf. "center_flowline"

Fobs = mp.sample_obs_along_flowline(O, FL)   # -> (point, time) with `vel`
Fmod = mp.model_flowline(M, FL)              # the model along the same line
mp.plot_flowline(Fobs, mode="profile", time=2019.6, model=Fmod)
```

`FL` uses the same `(point,)` layout as `load_obs_flowline` — coords `x`, `y`
and `d` (along-line distance in metres, **zero at the terminus**, increasing
inland) — so it drops straight into `model_flowline`, `compare_timeseries` and
the plots, and `FL.to_netcdf(...)` saves it.

The steps are also usable on their own:

```python
F = mp.load_flowlines(roi)          # GeoDataFrame: glacier, flowline, geometry
T = mp.load_termini(roi)            # GeoDataFrame: date (parsed), Quality_Fl
c = mp.center_flowline(F, glacier="a045")   # mp.ordered_flowlines for the fan
up = mp.clip_upstream(c.geometry, mp.latest_terminus(T, line=c.geometry).geometry)
ds = mp.flowline_to_dataset(mp.clip_to_roi(up, roi), name="a045_06", spacing=200)
```

Notes on the defaults, all overridable:

- `load_flowlines` skips the `*_iterNN.shp` intermediate versions
  (`iterations=True` keeps them) and any shapefile without a `flowline` column.
- `load_termini` keeps only `Quality_Fl == 0` traces (`max_quality=None` keeps
  all) and takes a `trange` if you want the front position of a given epoch.
- "Centre" flowline means the one whose downstream end sits closest to the mean
  of all the downstream ends; pass `flowline="05"` to choose explicitly.
- If **no trace crosses the centre line**, the next line out is tried instead,
  working outwards and alternating sides (`ordered_flowlines`), and the first
  one a trace does cross is used — so a centre line that happens to fall in a
  gap between traces no longer sinks the call. `FL.attrs["center_flowline"]`
  records the centre of the fan, so a differing `flowline` attr means a swap
  happened. Passing `flowline=` explicitly disables the search.
- The split uses the **furthest-inland** crossing, so a wiggly terminus leaves
  no seaward remnant. Only when no trace crosses *any* line of the fan does the
  call raise, unless you pass `require_terminus=False` (which falls back to the
  un-split centre line).

The MATLAB mirror is `../flowline/` (`moc_flowline_from_roi` and the same
steps, `moc_`-prefixed) — same defaults, same results; it needs the Mapping
Toolbox.

Worked example:

```bash
python examples/example_flowline_from_shapefiles.py --save-figs
python examples/example_flowline_from_shapefiles.py --roi my_roi --glacier b045
```

## Guided tour on the real data

`examples/example_sverdrups_workflow.py` runs the whole package end to end on
files that are actually in the configured folders — the Python mirror of
`examples/example_sverdrups_workflow.m`. Sverdrups Glacier is the one place
where a gridded subset and a redrawn flowline cover the same ice:

```bash
python examples/example_sverdrups_workflow.py              # show the figures
python examples/example_sverdrups_workflow.py --save-figs  # write figs/python/*.png
python examples/example_sverdrups_workflow.py --no-model   # observations only
```

The same tour as a notebook — narrated section by section, with the figures and
tables already executed so it reads without running anything:
`notebooks/moc_python_sverdrups_workflow.ipynb` (runs end to end in ~11 s;
`jupyter lab notebooks/moc_python_sverdrups_workflow.ipynb`). It adds `python/`
to `sys.path` if `moc_py` is not pip-installed, so it works straight from a
clone.

It reads `SentinelMonthly_Sverdrups.nc` (72 monthly frames, 2015–2020) and
`sverdrups.mat` (31 points, 0–7.5 km, 2015–2024), builds an ROI from the
flowline footprint and round-trips it through `rois/` (`.json` + ISSM `.exp`),
maps the speed with ROI and flowline overlaid, runs `velocity_variance` three
ways (melt-season window, two individual slices, robust record-long + coverage),
prints a per-year table, and plots the flowline distance–time image beside a
point series cross-checked against the nearest gridded pixel.

**The model section runs by default here**, unlike the MATLAB tour: it reads the
exported `model_output/exported/Model_NW_fric1_Snapshot.nc` (232,756 vertices,
465,007 elements, 376 steps 2007–2022) in about a second, where MATLAB would
have to load a multi-GB ISSM `.mat`. It adds `spatial_diff` + `plot_diff`, a
modelled-vs-observed flowline profile, and a paired point series — currently
bias +48 m/yr and RMSE 157 m/yr over the ROI at t=2019.6, and bias +255 m/yr,
RMSE 338 m/yr, r=0.28 at d=3 km. The whole run takes ~7 s.

Its per-year table matches the MATLAB tour's to every printed digit (pixel
counts included), which is the practical parity check between the two
implementations.

## Normalised variability maps

`velocity_variance` maps how much each pixel's speed *changed* over a temporal
window, divided by that pixel's own reference speed. The absolute spread is
dominated by the fast trunk simply because it is fast; the normalised field is
dimensionless, so a pixel is comparable with its neighbours as a **fraction of
relative change**.

```python
V = mp.velocity_variance(O, roi, trange=(2019.0, 2019.5))   # roi optional
mp.plot_variance(V, "norm", roi=roi)      # spatial map of the fractional range
mp.plot_variance(V, "spread")             # the same window in m/yr
V.attrs["norm_mean"], V.attrs["norm_median"]   # ROI-mean fractional change

# two time slices: nearest frame to each date, plus a signed change field
V2 = mp.velocity_variance(O, roi, times=(2019.0, 2019.5))
mp.plot_variance(V2, "norm_change")       # (later - earlier) / speed, diverging
```

Fields: `spread` (m/yr), `reference` (m/yr), `norm = spread/reference`, `count`,
and `change`/`norm_change` when exactly two frames are used. Options:
`statistic` = `range` (default, peak-to-peak) / `std` / `iqr` / `mad`;
`normalize` = `mean` (default) / `median` / `first` / `max` / `none`;
`min_count` (frames a pixel needs) and `min_speed` (guard against huge ratios in
near-stagnant ice, off by default).

Mirrored on the MATLAB side by `compare/moc_velocity_variance.m` +
`plot/moc_plot_variance.m` (same options and numbers; `.normChange` there for
`norm_change` here).

## Lazy / cloud reads

`load_obs_netcdf` reads a whole local file into memory (via h5py). For cloud or
out-of-core work — remote COG stacks, zarr datacubes, HTTP netCDF — use
`open_obs_dataarray`, which returns a **lazy** (dask-backed) `xarray.DataArray`
and normalises it to moc_py conventions (dims `(y, x, time)`, `time` as decimal
year + `time_dt`, `epsg` attr). Nothing is read until you slice/compute, so the
downstream functions realise only the pixels/frames they touch.

It accepts an existing DataArray/Dataset **or** a path/URL/`.zarr`/glob. That
makes it the bridge for an already-lazy array such as the GrIMP
`myVelSeries.vv` (remote Sentinel-1 COGs) from your Flowlines notebook:

```python
# myVelSeries.vv : lazy xarray.DataArray, dims (time, y, x), datetime time coord
da = mp.open_obs_dataarray(myVelSeries.vv, band="vv")   # -> (y, x, time), lazy
O  = mp.to_obs_dataset(da)                               # obs Dataset (still lazy)

# now everything else just works, computing only what it needs:
roi = mp.roi_from_obs(O, margin=2000)
D   = mp.spatial_diff(M, O, roi, time=2019.6)            # realises 1 obs frame
C   = mp.compare_timeseries(M, O, xy=(x0, y0))           # realises 1 pixel column
mp.plot_map(O, band="vv", tindex=-1)
```

Or open a remote source directly:

```python
da = mp.open_obs_dataarray("https://.../cube.zarr", band="v", chunks="auto")
da = mp.open_obs_dataarray("s3://bucket/vv_*.tif", band="vv")   # COG stack
```

Lazy/cloud backends are optional extras: `pip install -e "python/[lazy]"`
(dask, rioxarray, zarr, h5netcdf, fsspec). `dask` is what keeps arrays chunked;
`rioxarray` supplies CRS (`.rio.crs` → epsg) and COG support.

## Layout

```
python/
  moc_py/
    __init__.py        public API
    config.py          path resolution + defaults (mirrors config_moc.m /
                       config/moc_path.m); mp.CONFIG, mp.configure
    __main__.py        `python -m moc_py` — show or set data paths
    util.py            datenum / CF-time -> datetime64 / decimal-year
    obs_netcdf.py      load_obs_netcdf  (h5py -> xarray, robust dim matching)
    obs_flowline.py    load_obs_flowline (scipy -> xarray, source normalising)
    roi.py             Roi, roi_from_obs, roi_mask, save/load, RoiSelector
    flowline_shp.py    flowline_from_roi + steps (Felikson lines x TBlack termini)
    variance.py        velocity_variance (normalised temporal variability maps)
    plot.py            plot_map, plot_flowline, plot_timeseries, plot_variance
    smoke_test.py      python -m moc_py.smoke_test  (no-plot checks)
  examples/
    example_sverdrups_workflow.py       end-to-end tour on the data that is on disk
    example_flowline_from_shapefiles.py ROI -> cropped, terminus-clipped flowline
  notebooks/
    moc_python_quickstart.ipynb
    moc_python_model_compare.ipynb
    moc_python_sverdrups_workflow.ipynb   the tour below, with outputs saved
  tests/
    test_model_synthetic.py             Phase 2 checks, synthetic mesh
    test_variance_synthetic.py          velocity_variance, synthetic grid
    test_flowline_shp_synthetic.py      flowline_shp, synthetic shapefiles
  requirements.txt
```

## Smoke test

```bash
python -m moc_py.smoke_test
```

Loads one of each obs type (auto-picking the smallest matching files), prints
shapes/ranges, and checks that netCDF bands are `(y, x, time)` with
`vv ≈ hypot(vx, vy)` (band order) and that flowline velocity is `(point, time)`
(catches the transposed `u2.v`). All checks pass against the current data.

The `tests/test_*_synthetic.py` scripts need no data at all — each fabricates
its own inputs and asserts against analytic answers:

```bash
python tests/test_flowline_shp_synthetic.py   # skips if geopandas is absent
```

## Conventions & parity with MATLAB

- Coordinates in metres, EPSG:3413. Time as **decimal year** (`time`) plus a
  `datetime64` (`time_dt`) coordinate.
- Flowline sources tried in order `filteredV → SentinelV → SentinelVmonthly → u2`;
  `u2.v` is auto-transposed to `(point, time)`.
- netCDF band order assumed `vv, vx, vy, ex, ey, dT` but read from the file's
  `component` variable when present. Two non-declared nodata sentinels are
  masked to NaN: `-2e9` on the signed bands (`vx`, `vy`) and `-1` on the
  physically non-negative ones (`vv`, `ex`, `ey`, `dT`).
- ROIs round-trip to both `.json` (Python) and `.exp` (ISSM/MATLAB) so an ROI
  drawn here can be reused by the MATLAB `moc_*` functions and vice versa.

## Phase 2 (model side) — built

The ISSM model stays in MATLAB only for a one-time **export**; everything after
is Python.

1. **MATLAB** `moc_export_model.m` writes a neutral netCDF4: `x(nv)`, `y(nv)`,
   `elements(ne,3)` (**1-based**, tagged `index_base=1`), `time(nt)` (decimal
   year), `vel/vx/vy(nv,nt)`, optional `surface`/`bed`, with `epsg` in the
   attrs. See `examples/example_export_model.m`. Run once per model.
2. **Python** `moc_py`:
   - `load_model(nc)` → :class:`Model` (elements auto-converted to 0-based;
     axes matched by size, so the MATLAB/HDF5 transpose quirk is handled).
   - `model_to_grid`, `model_at_points`, `model_flowline` — mesh interpolation
     via `matplotlib.tri.LinearTriInterpolator`, the same P1 linear
     interpolation as ISSM's `InterpFromMesh2d` / `InterpFromMeshToGrid`.
   - `spatial_diff(M, O, roi, ...)` → model/obs/diff on a common grid + stats;
     `compare_timeseries(M, obs, ...)` → paired series at a point/flowline.
   - plots: `plot_grid`, `plot_diff`, `plot_compare`; `plot_flowline` /
     `plot_timeseries` accept a `model=` overlay.

Walkthrough: `notebooks/moc_python_model_compare.ipynb`.

**Validation (no MATLAB needed):** `tests/test_model_synthetic.py` fabricates a
mesh with a known linear field, exercises the whole path, and asserts the
interpolation is exact (~1e-13) and a model-vs-itself difference is ~0 — for
both stored axis orders. All checks pass.

```bash
python tests/test_model_synthetic.py
```

`tests/test_variance_synthetic.py` does the same for `velocity_variance`: a
fabricated fast-ice/slow-ice grid where the absolute spread is 10x larger on the
fast side and the normalised spread 2x larger on the slow side, checking the
recovered fractions, ROI cropping and polygon masking, two-slice selection and
signed change, the `min_count`/`min_speed` guards, and the plots.

```bash
python tests/test_variance_synthetic.py
```
```
