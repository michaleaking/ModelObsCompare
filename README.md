# moc — Model vs Observation Glacier Velocity Comparison Module

`moc` (**Model vs Observation Glacier Velocity Comparison Module**) compares
**observed** surface velocities against **modelled** velocities from ISSM
over glacier basins of interest on the Greenland Ice Sheet.

This repo includes functionality for point-by-point, along-flow, and 2D 
comparisons. Tools for quantifying bias and rerunning bias-corrected differences are included, 
though in early development. 

- **Regions of interest** — draw a box or polygon on a velocity backdrop, or
  derive one from reading in boundaries from a glacier-specific netCDF file. Generated ROIs 
  can be formatted as JSON and or ISSM-compatible `.exp` contours, for easy pass-off to model
  runs if needed.
- **Time-series comparison** at a coordinate or along a flowline, with bias,
  RMSE, and correlation over common or interpolated-to-common timestamps.
- **Spatial-difference maps** (model − obs) over an ROI, regridding both fields
  onto a common grid.
- **Observed variability maps** — how much a pixel's speed swings over a melt
  season or a whole record, in absolute and speed-normalised terms.
- **Flowlines from archived shapefiles** — given an ROI, pull the centre
  flowline of the glacier through it and clip it at the most recently traced
  terminus, ready to compare along.

## Both Python and MATLAB implementations included in toolkit

The same capabilities exist in **MATLAB** and in **Python**, same conventions and 
data apply:

| | MATLAB | Python |
|---|---|---|
| location | this folder | [`python/`](python/) |
| entry point | `setup_moc` | `import moc_py as mp` |
| naming | `moc_load_obs_netcdf(...)` | `mp.load_obs_netcdf(...)` |
| model side | reads ISSM `.mat` directly (needs ISSM) | reads a netCDF exported by `moc_export_model.m` |
| detailed docs | this file | [`python/README.md`](python/README.md) |

Use MATLAB when you are already inside an ISSM session and want the model
object itself. Use Python for the observation side, for notebooks, and for
anything you want to run without a MATLAB licence. The model-side Python
functions read a netCDF that the MATLAB `moc_export_model.m` writes, so a
typical workflow crosses over once: export in MATLAB, analyse in Python.

Both read the **same configuration file**, so data locations are set once for
both (see [Configuring data paths](#configuring-data-paths)).

## Relationship to green-ibis

[`green-ibis`](../green-ibis) — *Greenland Ice Bed from Ice Surface* — infers
seasonal-timescale **basal** processes (basal shear stress, effective pressure,
friction coefficients) from ice-surface observations and ISSM model output.
Its `utils/` holds the physics: `basal_shear_stress_budd.m`,
`basal_shear_stress_schoof.m`, `driving_stress.m`,
`effective_pressure_full_connectivity.m`, and so on.

Those inversions are only as trustworthy as the velocity field underneath them.
That is where `moc` fits: it is the **validation and diagnostic step** on either
side of a green-ibis analysis.

```
   ISSM transient / inversion run
              |
              v
   +----------------------+        +--------------------------------+
   |  moc  (this repo)    |        |  green-ibis                    |
   |  does the model      | -----> |  infer basal conditions from   |
   |  reproduce observed  |  once  |  a velocity field you trust    |
   |  velocity?           |  it    |  (tau_b, N, C, driving stress) |
   +----------------------+  does  +--------------------------------+
        ^                                      |
        |______________________________________|
             where it does not, moc localises the
             misfit in space, time and along-flow
```

Concretely:

- **Before** a green-ibis inversion, use `moc` to check that the run being
  inverted actually matches observations over the glacier and period of
  interest — `moc_spatial_diff` for the map, `moc_compare_timeseries` for a
  point, `moc_plot_flowline` along the trunk.
- **After**, when an inferred basal field looks surprising, `moc` tells you
  whether the surprise is real or is inherited from a velocity misfit in that
  same place.
- **Shared conventions.** Both use EPSG:3413, metres, decimal years, and ISSM
  `.exp` contours; an ROI saved here can define a green-ibis analysis domain.
  Both also resolve machine-specific paths through environment variables —
  green-ibis uses `ISSM_DIR`, which `moc` reads too.
- **Complementary, not overlapping.** green-ibis `utils/extract_obs_from_model.m`
  pulls the observations *that were used to calibrate* a transient run; `moc`
  reads the observational products independently, so it can also judge periods
  and places the calibration never saw.

The two are separate repositories on purpose — `moc` has no green-ibis
dependency and vice versa — but they are designed to sit side by side in the
same parent folder.

## Installing the Python module (`moc_py`)

`moc_py` is a normal Python package; install it into whichever environment you
run notebooks and analyses from. From the repository root:

**pixi** (the project convention):

```bash
pixi add numpy scipy h5py xarray matplotlib pandas
pixi add geopandas shapely pyogrio        # optional: the flowline/shapefile tools
pixi add ipympl jupyterlab                # optional: notebooks + interactive ROIs
pixi run sync-env                         # keep environment.yml in step for conda users
pixi run pip install -e python/           # the package itself, editable
```

**conda / mamba:**

```bash
conda create -n moc python=3.11
conda activate moc
conda install -c conda-forge numpy scipy h5py xarray matplotlib pandas                              geopandas shapely pyogrio ipympl jupyterlab
pip install -e python/
```

**pip only** (geopandas wheels are fine on modern pip):

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e "python/[shapefile]"       # core + geopandas/shapely/pyogrio
# or: pip install -e python/              # core only
# or: pip install -r python/requirements.txt
```

Editable (`-e`) is recommended while the toolkit is still moving: your edits
take effect without reinstalling. Verify:

```bash
python -c "import moc_py; print(moc_py.__file__)"
python -m moc_py                # show configured data paths
python -m moc_py.smoke_test     # load one of each obs type and check orientation
```

Only `numpy scipy h5py xarray matplotlib pandas` are required. The
`shapefile` extra (`geopandas`, `shapely`, `pyogrio`) is needed **only** for
`moc_py.flowline_shp`; the `lazy` extra (`dask`, `rioxarray`, `zarr`, ...) only
for cloud/lazy reads. Everything else imports without them.

For **MATLAB**, there is nothing to install — run `setup_moc` from this folder,
which puts the toolkit (and ISSM, if configured) on the path.

## Configuring data paths

No path in either implementation is hard-coded to a machine. Each data location
is resolved the first time it is needed, in this order:

1. **an explicit argument** — `moc_load_flowlines(roi, 'dir', '/data/flowlines')`,
   or `mp.load_flowlines(roi, directory="/data/flowlines")`;
2. **an environment variable** — `MOC_MODEL_DIR`, `MOC_OBS_FLOWLINE_DIR`,
   `MOC_OBS_NETCDF_DIR`, `MOC_FLOWLINE_SHP_DIR`, `MOC_TERMINI_SHP`,
   `MOC_ROI_DIR`, and `ISSM_DIR` (the same variable green-ibis uses);
3. **the shared config file** — `~/.config/moc/paths.json` (override its
   location with `MOC_CONFIG`). **MATLAB and Python read and write the same
   file**, so answering once configures both;
4. **an interactive prompt** — when running in a terminal, `moc` asks for the
   path it needs, checks that it exists, and saves the answer to the config
   file so it is asked once.

If none of those apply (a script, a notebook with no console, CI), the error
names the missing setting and the three ways to supply it, rather than failing
on a stale path from someone else's machine.

Set them up front, or let the prompts happen as you go:

```matlab
moc_configure('-set')                              % prompt for anything missing
moc_configure('obs_netcdf_dir', '/data/velocity')  % this session
moc_configure('-save')                             % ... and remember it
moc_configure()                                    % show the table below
```

```bash
python -m moc_py --set          # prompt for anything missing
moc-config --set                # same, if the package is installed
```

```python
import moc_py as mp
mp.configure(obs_netcdf_dir="/data/velocity")
mp.save_config()                # persist to ~/.config/moc/paths.json
mp.CONFIG.status()              # {'obs_netcdf_dir': ('/data/velocity', 'argument'), ...}
```

Either `moc_configure()` or `python -m moc_py` prints where every setting came
from:

```
setting            source                 value
------------------ ---------------------- ----------------------------------------
issm_dir           ISSM_DIR               /opt/issm
model_dir          config file            /data/model_output
obs_netcdf_dir     session                /data/velocity
flowline_shp_dir   unset                  (not set)
roi_dir            built-in               .../ModelObsCompare/rois
```

| setting | what to point it at |
|---|---|
| `issm_dir` | root of the ISSM binary distribution (holds `bin/`, `lib/`) |
| `model_dir` | ISSM model output `.mat` / exported `.nc` files |
| `obs_flowline_dir` | per-glacier observed flowline `.mat` files |
| `obs_netcdf_dir` | gridded observed-velocity netCDF files |
| `flowline_shp_dir` | Felikson per-glacier flowline shapefiles |
| `termini_shp` | Black & Joughin traced-terminus shapefile |
| `roi_dir` | where selected ROIs are saved (defaults inside the repo) |

The datasets themselves are not in this repository — they are large and
externally archived. See [Data assumptions](#data-assumptions) for what each
one has to look like.

---

# MATLAB reference

The rest of this file documents the MATLAB implementation. For the Python
equivalent see [`python/README.md`](python/README.md).

## Requirements

- MATLAB (R2018b+ recommended for `drawrectangle`/`drawpolygon`; older releases
  fall back to `ginput`).
- ISSM **binary** distribution on disk — *optional*. It provides `loadmodel`,
  `InterpFromMeshToGrid`, `InterpFromMesh2d`, `SectionValues` and
  `expread/expwrite`, so without it the model-side functions are unavailable;
  the observation, ROI, variance and flowline tools all still work.
- **Mapping Toolbox** — only for `flowline/` (`shaperead`, `polyxpoly`); every
  other part of the toolkit runs without it.

## Setup

At the start of every session:

```matlab
cfg = setup_moc();   % adds this toolkit (and ISSM, if configured) to the path
```

`setup_moc` asks for nothing: each data location is resolved the first time it
is actually needed, so a session that only touches observations is never
interrupted about the model folder. It prints which locations are known and
which are not. To fill them in up front, see
[Configuring data paths](#configuring-data-paths) above:

```matlab
moc_configure('-set')     % prompt for anything missing, then remember it
```

`config_moc()` returns the conventions — `.epsg` (3413, Greenland Polar
Stereographic), `.grid_res` (200 m default regrid spacing), `.obs_bands` — plus
every path that is already configured. Because a path may legitimately be
unset, reach for one with `moc_path` rather than assuming the field exists:

```matlab
d   = moc_path('model_dir');    % resolves, prompting once if it must
cfg = config_moc();  cfg.epsg   % conventions are always present
```

## Data assumptions

**Model output** (`model_dir/Model_*.mat`) — saved ISSM `model` objects. The
loader supports both transient runs (`results.TransientSolution` with a `time`
array and `Vx/Vy/Vel/Surface` per step) and snapshot inversions
(`results.StressbalanceSolution`). Velocities are treated as **m/yr**.

**Observed gridded velocity** (`obs_netcdf_dir/*.nc`, e.g.
`Sentinel_Subset_Upernavik.nc`) — a 4-D `velocity(time, component, y, x)` array
with component bands `{vv, vx, vy, ex, ey, dT}`, coordinates in EPSG:3413, and
`time` as CF "days since …". The loader is robust to dimension order (it matches
axes by size) and converts time to decimal years. It also masks nodata to NaN:
these files declare a useless `_FillValue` **attribute** (NaN) and keep the real
per-band sentinels in a `_FillValue` **variable** — `-1` on the physically
non-negative bands (`vv`, `ex`, `ey`, `dT`) and `-2e9` on `vx`/`vy`. Without
that, 15–30% of pixels read as a `-1 m/yr` speed and quietly poison every map,
difference and variance downstream.

**Observed flowline files** (`obs_flowline_dir/<glacier>.mat`) — one glacier per
file: ordered `(x, y)` flowline points, along-line distance `d`, static geometry
(`zb`, `zthick`), and velocity time series in one of several layouts
(`filteredV` + `vti`, or structs `SentinelV`, `SentinelVmonthly`, `u2`). Times
are **MATLAB datenum** and are converted to decimal years. Point/time
orientation is auto-detected (`u2.v` is stored time × points and is transposed).

**Flowline shapefiles** (`flowline_shp_dir/glacier*.shp`) — Felikson et al.
(2020), one file per glacier holding several `Line` features identified by a
`flowline` attribute. Vertices run **from the terminus inland**, so point 1 is
the downstream end. `*_iterNN.shp` are intermediate versions and are skipped by
default.

**Terminus shapefiles** (`termini_shp`) — Black & Joughin traced terminus
positions, one `Line` per glacier per image with a `SourceDate` (`yyyy-mm-dd`)
and a `Quality_Fl` (non-zero flags lower-confidence traces; only `0` is kept by
default). Note its `Glacier_ID` is **unrelated** to the flowline collection's
glacier ids — the two are paired spatially, never by id.

Both are read as-is in `cfg.epsg` (EPSG:3413); nothing is reprojected.

## Directory layout

```
ModelObsCompare/
  config_moc.m            conventions + already-configured paths
  setup_moc.m             path setup (run first)
  config/                 machine-independent path resolution
    moc_path.m              resolve one location (env -> file -> prompt)
    moc_configure.m         show / set / save data paths
    moc_path_settings.m     the list of settings (shared with moc_py)
    moc_config_file.m       location + contents of ~/.config/moc/paths.json
    moc_clean_path.m        tidy a pasted path; moc_home_dir.m
  io/                     loaders
    moc_load_model.m        ISSM .mat -> lean mesh + velocity struct
    moc_load_obs_netcdf.m   gridded obs .nc -> (ny x nx x nt) bands
    moc_load_obs_flowline.m per-glacier flowline .mat -> normalised struct
    moc_list_models.m       table of available runs
  roi/                    region of interest
    moc_select_roi.m        interactive draw (rectangle/polygon)
    moc_roi_from_obs.m      ROI from a netCDF extent
    moc_save_roi.m          save .mat (+ ISSM .exp)
    moc_load_roi.m          load .mat or .exp
    moc_roi_mask.m          inpolygon test
  interp/                 mesh <-> grid/points
    moc_model_to_grid.m     regrid mesh field to regular grid (InterpFromMeshToGrid)
    moc_model_at_points.m   sample mesh at points over time (InterpFromMesh2d)
    moc_model_flowline.m    sample mesh along a flowline over time
  flowline/               flowlines from the shapefile archives (Mapping Toolbox)
    moc_flowline_from_roi.m ROI -> cropped, terminus-clipped flowline (the pipeline)
    moc_load_flowlines.m    loop every glacier*.shp, keep those in the ROI
    moc_load_termini.m      terminus traces in the ROI, dates parsed
    moc_center_flowline.m   middle line of a glacier's fan
    moc_ordered_flowlines.m rank a fan centre-outwards (the fallback order)
    moc_latest_terminus.m   newest trace, optionally one crossing a line
    moc_clip_upstream.m     split at the terminus, keep the upstream side
    moc_clip_to_roi.m       clip a line to the ROI polygon
    moc_roi_polygon.m       coerce struct/name/bbox/matrix -> polygon vertices
    moc_flowline_to_struct.m package (x,y) -> .x/.y/.d flowline struct
    moc_sample_obs_along_flowline.m  gridded obs -> (point x time) .vel
  compare/                pairing + metrics
    moc_spatial_diff.m      model-minus-obs on a common grid + stats
    moc_compare_timeseries.m paired series at a point/flowline + stats
    moc_velocity_variance.m normalised temporal variability of observed speed
  plot/
    moc_plot_map.m          gridded field map
    moc_plot_diff.m         3-panel model/obs/difference
    moc_plot_timeseries.m   paired time series
    moc_plot_flowline.m     profile or distance-time (Hovmoller) difference
    moc_plot_variance.m     normalised variability map
  util/
    moc_datenum2decyear.m / moc_decyear2datenum.m
  examples/
    example_select_roi.m
    example_spatial_diff.m
    example_flowline_timeseries.m
    example_velocity_variance.m
    example_sverdrups_workflow.m  end-to-end tour on the data that is on disk
    example_flowline_from_shapefiles.m  ROI -> cropped, terminus-clipped flowline
  tests/
    smoke_test_moc.m        no-plot loader + orientation checks
    test_variance_synthetic.m  synthetic check of moc_velocity_variance
    test_flowline_shp_synthetic.m  synthetic check of the flowline/ pipeline
  rois/                   saved ROIs (.mat + .exp)
```

## First run: smoke test

Before trusting any plots or metrics, run the smoke test. It loads one of each
data type (auto-picking the *smallest* model and netCDF files for speed), prints
shapes/ranges, and runs targeted orientation checks — no plotting:

```matlab
smoke_test_moc      % prints a PASS/WARN/FAIL summary
```

It verifies, in particular:
- obs netCDF bands come back as `(ny x nx x nt)` and `vv ≈ hypot(vx,vy)`
  (catches a wrong band order),
- flowline velocity is `(npts x nt)` (catches the transposed `u2.v`),
- `InterpFromMeshToGrid` output is `(numel(y) x numel(x))` and grid cell values
  match an independent `InterpFromMesh2d` sample at the same coordinates
  (catches a transposed/flipped regrid — the main thing to confirm).

Force specific files by editing the CONFIG block at the top of the script.

## Quick start

```matlab
cfg = setup_moc();

% --- spatial difference over an ROI ---
M   = moc_load_model('Model_NW_HindcastRun_TransientInversion.mat');
O   = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc','bands',{'vv'});
roi = moc_roi_from_obs(O,'margin',2000);
D   = moc_spatial_diff(M, O, roi, 'time',2019.0, 'obswindow',0.1);
moc_plot_diff(D, 'roi',roi);

% --- normalised variability over a temporal window ---
V = moc_velocity_variance(O, roi, 'trange',[2019.4 2019.8]);
moc_plot_variance(V, 'field','norm', 'roi',roi);   % change / local speed [-]
V2 = moc_velocity_variance(O, roi, 'times',[2019.45 2019.75]);  % two slices
moc_plot_variance(V2, 'field','normChange', 'roi',roi);

% --- flowline time series ---
Fobs = moc_load_obs_flowline('jakobshavn.mat');
Fmod = moc_model_flowline(M, Fobs);
C    = moc_compare_timeseries(M, Fobs, 'dist',5000);
moc_plot_timeseries(C,'datetime',true);
```

See `examples/` for the full walkthroughs.

## Guided tour on the real data

`examples/example_sverdrups_workflow.m` runs the toolkit end to end on files
that are actually in the configured folders. Sverdrups Glacier is the one place
where a gridded subset and a redrawn flowline cover the same ice:
`SentinelMonthly_Sverdrups.nc` (72 monthly frames, 2015–2020) and
`sverdrups.mat` (31 points, 0–7.5 km, 2015–2024).

```matlab
example_sverdrups_workflow                      % figures on screen
SAVE_FIGS = true;  example_sverdrups_workflow   % also writes figs/*.png
RUN_MODEL = true;  example_sverdrups_workflow   % adds the model panels
```

It loads both observation types, builds an ROI from the flowline footprint and
round-trips it through `rois/` (`.mat` + ISSM `.exp`), draws a speed map with
the ROI and flowline on top, then the variability maps (melt-season window, two
individual slices, robust record-long + coverage), a per-year table of relative
change, and the flowline distance–time image beside a point series compared
against the nearest gridded pixel. Section 7 (model difference, flowline
profile, paired series) is **off by default** because it loads a multi-GB ISSM
`.mat`; set `RUN_MODEL = true` to include it.

The Python mirror is `python/examples/example_sverdrups_workflow.py` (also as a
notebook, `python/notebooks/moc_python_sverdrups_workflow.ipynb`); it runs the
model section by default — the exported model netCDF loads in ~1 s — and
reproduces the table below to every digit.

Numbers from the current data, for reference — over the ROI, melt season
(Apr–Nov), `minspeed` 20 m/yr:

| year | frames | mean relative change | median | mean range |
|------|--------|----------------------|--------|------------|
| 2015 | 7 | 27.1% | 14.5% | 156 m/yr |
| 2016 | 7 | 39.6% | 15.7% | 170 m/yr |
| 2017 | 7 | 34.3% | 17.1% | 211 m/yr |
| 2018 | 7 | 34.8% | 15.2% | 204 m/yr |
| 2019 | 7 | 44.2% | 19.4% | 190 m/yr |
| 2020 | 7 | 38.8% | 18.5% | 183 m/yr |

## Normalised variability maps

`moc_velocity_variance` reduces the time axis of a gridded observation to a
per-pixel map of how much the speed *changed* over a window, divided by that
pixel's own reference speed. The absolute spread is dominated by the fast trunk
simply because it is fast; the normalised field is dimensionless, so a pixel is
comparable with its neighbours as a **fraction of relative change**.

```matlab
V = moc_velocity_variance(O, roi, 'trange',[2019.4 2019.8]);  % roi may be []
moc_plot_variance(V, 'field','norm',   'roi',roi);   % fractional range [-]
moc_plot_variance(V, 'field','spread', 'roi',roi);   % the same window in m/yr
V.stats.normMean, V.stats.normMedian                 % ROI-mean relative change

% two individual time slices: nearest frame to each date + a signed change field
V2 = moc_velocity_variance(O, roi, 'times',[2019.45 2019.75]);
moc_plot_variance(V2, 'field','normChange', 'roi',roi);   % diverging, centred on 0
```

Fields: `.spread` [m/yr], `.reference` [m/yr], `.norm = spread./reference`,
`.count`, and `.change`/`.normChange` when exactly two frames are used. Options:
`'statistic'` = `range` (default, peak-to-peak) / `std` / `iqr` / `mad`;
`'normalize'` = `mean` (default) / `median` / `first` / `max` / `none`;
`'mincount'` (frames a pixel needs, default 2) and `'minspeed'` (guard against
huge ratios in near-stagnant ice, off by default). Walkthrough:
`examples/example_velocity_variance.m`.

**Validation (no data, no ISSM):** `tests/test_variance_synthetic.m` fabricates a
fast-ice/slow-ice grid where the absolute spread is 10x larger on the fast side
and the normalised spread 2x larger on the slow side, then checks the recovered
fractions, ROI cropping and polygon masking, two-slice selection and signed
change, the `mincount`/`minspeed` guards, and the plots. All checks pass.

```matlab
test_variance_synthetic
```

This mirrors `moc_py.velocity_variance` / `mp.plot_variance` on the Python side
(same options, same numbers; MATLAB uses `.normChange` where Python uses
`norm_change`).

## Flowlines from the shapefiles

`flowline/` turns the two archived shapefile collections into a flowline the
comparison functions accept. Give it an ROI and it loops the whole Felikson
collection (`cfg.flowline_shp_dir`), keeps the flowlines that reach into the
ROI, takes the centre one, splits it at the most recent Black & Joughin terminus
trace that crosses it (`cfg.termini_shp`), keeps the upstream side, and clips
that to the ROI:

```matlab
roi = moc_load_roi('sverdrups_trunk');       % or a struct, bbox, or np x 2 matrix

FL = moc_flowline_from_roi(roi, 'spacing',200);
FL.terminus_date        % '2021-12-21' — the trace it was clipped at
FL.flowline             % the line used; compare with FL.center_flowline
FL.d(1)                 % 0, at the terminus; d increases inland

O    = moc_load_obs_netcdf('SentinelMonthly_Sverdrups.nc', 'bands',{'vv'});
Fobs = moc_sample_obs_along_flowline(O, FL);   % -> .vel (npoints x ntime)
M    = moc_load_model('Model_NW_HindcastRun_TransientInversion.mat', ...
                      'solution','transient', 'trange',[2018 2020]);
Fmod = moc_model_flowline(M, FL);
moc_plot_flowline(Fmod, Fobs, 'mode','profile', 'time',2019.6);
C    = moc_compare_timeseries(M, Fobs, 'dist',2500);
```

`FL` carries the same `.x`, `.y`, `.d` geometry fields as
`moc_load_obs_flowline`, so it goes straight into `moc_model_flowline`, and
`moc_sample_obs_along_flowline` adds the `.vel`/`.time`/`.source` fields the
plots and `moc_compare_timeseries` look for.

The steps are also usable on their own:

```matlab
F = moc_load_flowlines(roi);                     % struct array: glacier, flowline, x, y
T = moc_load_termini(roi);                       % struct array: date, quality, x, y
c = moc_center_flowline(F, 'glacier','a045');
t = moc_latest_terminus(T, 'line',[c.x c.y]);
[xu, yu] = moc_clip_upstream(c.x, c.y, t.x, t.y);
[xc, yc] = moc_clip_to_roi(xu, yu, roi);
FL = moc_flowline_to_struct(xc, yc, 'name','a045_06', 'spacing',200);
```

Notes on the defaults, all overridable:

- `moc_load_flowlines` skips the `*_iterNN.shp` intermediate versions
  (`'iterations',true` keeps them) and any shapefile without a `flowline`
  attribute (e.g. `sverdrup_masks.shp`).
- `moc_load_termini` keeps only `Quality_Fl == 0` traces (`'max_quality',[]`
  keeps all) and takes a `'trange'` if you want the front position of a given
  epoch.
- "Centre" flowline means the one whose downstream end sits closest to the mean
  of all the downstream ends; pass `'flowline','05'` to choose explicitly.
- If **no trace crosses the centre line**, the next line out is tried instead,
  working outwards and alternating sides (`moc_ordered_flowlines`), and the
  first one a trace does cross is used — so a centre line that happens to fall
  in a gap between traces no longer sinks the call. `FL.center_flowline` records
  the centre of the fan, so `FL.flowline` differing from it means a swap
  happened. Naming `'flowline'` explicitly disables the search.
- The split uses the **furthest-inland** crossing, so a wiggly terminus leaves
  no seaward remnant. Only when no trace crosses *any* line of the fan does the
  call error, unless you pass `'require_terminus',false` (which falls back to
  the un-split centre line).

Worked example: `examples/example_flowline_from_shapefiles.m` (set
`RUN_MODEL = true` for the model panel). The Python mirror is
`python/moc_py/flowline_shp.py` — same defaults, same results.

## Conventions & notes

- All coordinates are in metres, EPSG:3413. Model, obs, and flowlines must share
  this projection (they do for the configured datasets).
- Naming is `moc_*` for all functions to avoid clashing with ISSM's own
  functions already on the path.
- Model files are large (multi-GB). `moc_load_model` extracts only the mesh and
  velocity/geometry it needs; consider `'trange'` to subset transient runs.
- **Untested against the live data as shipped** — the code was written from the
  observed file structures and ISSM function signatures. Run the examples and
  verify the first plots (especially grid orientation from
  `InterpFromMeshToGrid` and obs band order) before trusting the metrics.
```
