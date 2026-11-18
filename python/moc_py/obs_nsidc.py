"""Adapter for GrIMP/NSIDC velocity subsets written by ``glacier_reader.py``.

The companion reader (``~/jupyter/glacierReader/glacier_reader.py``, lifted from
``matlab_glacier_reader.ipynb``) pulls velocity time series straight from NSIDC
with ``grimpfunc``/``nisardev`` and writes one netCDF per product:

    velocity(time, band, y, x)   bands {vx, vy, vv, ex, ey, ev, dT}, NaN no-data

That is the same physical content as the Sentinel/CSK subsets
:func:`moc_py.load_obs_netcdf` reads, in a slightly different dress: the band
axis is called ``band`` and carries its own string labels, where moc expects a
``component`` variable, and time is CF-encoded datetime rather than decimal
years. This module translates, so an NSIDC download drops straight into
:func:`moc_py.spatial_diff`, :func:`moc_py.compare_timeseries`,
:func:`moc_py.velocity_variance`, the ROI tools and the plotters.

The translation is not cosmetic. moc's eager reader looks for a ``component``
variable and, not finding one, falls back to the ``CONFIG.obs_bands`` order
``(vv, vx, vy, ex, ey, dT)`` — which is *not* the nisardev order
``(vx, vy, vv, ex, ey, ev, dT)``. Reading one of these files with
:func:`moc_py.load_obs_netcdf` directly would therefore return ``vx`` when you
asked for ``vv``, silently. Use :func:`load_nsidc_netcdf` instead.

Nothing here needs ``grimpfunc``, ``nisardev`` or GDAL — only an already
written file (plus ``h5netcdf``/``netCDF4`` to open it, and ``dask`` to stay
lazy). :func:`fetch_nsidc` is the one exception: it delegates to
``glacier_reader`` to download first, and so needs that environment.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from .config import CONFIG, resolve
from .obs_lazy import to_obs_dataset
from .util import datetime64_to_decyear

#: Band axis names used by the GrIMP/nisardev products, in preference order.
BAND_DIMS: tuple[str, ...] = ("band", "component")

#: GrIMP sentinel no-data values; anything at or below the sentinel is no data.
#: Mirrors ``glacier_reader.NO_DATA_SENTINELS``.
NO_DATA_SENTINELS: dict[str, float] = {
    "vx": -2.0e9, "vy": -2.0e9, "dT": -2.0e9,
    "vv": -1.0, "ex": -1.0, "ey": -1.0, "ev": -1.0,
}

#: Bands nisardev derives as ``sqrt(x**2 + y**2)``, mapped to their components.
#: Their no-data must be inherited, not thresholded: an unmasked sentinel pair
#: gives ``+2.83e9``, which survives any lower bound on speed.
DERIVED_BANDS: dict[str, tuple[str, str]] = {"vv": ("vx", "vy"), "ev": ("ex", "ey")}


def load_nsidc_netcdf(
    ncpath: str | Path,
    bands: list[str] | tuple[str, ...] | None = None,
    varname: str = "velocity",
    trange: tuple[float, float] | None = None,
    chunks: str | dict | None = "auto",
) -> xr.Dataset:
    """Load a ``glacier_reader`` NSIDC subset as a moc observation Dataset.

    Args:
        ncpath: Full path or bare filename (resolved against
            ``CONFIG.obs_netcdf_dir``, so downloads can live in the same place
            as the hand-made subsets).
        bands: Bands to keep, e.g. ``["vv"]``. ``None`` -> every band present.
        varname: Name of the stacked velocity variable (default ``"velocity"``).
        trange: ``(t0, t1)`` decimal-year window to keep. ``None`` -> all frames.
        chunks: Dask chunking, passed to :func:`xarray.open_dataset`. ``None``
            reads eagerly; ``"auto"`` (default) keeps the read lazy, which
            matters for a multi-year monthly mosaic stack.

    Returns:
        Dataset in the schema of :func:`moc_py.load_obs_netcdf`: one data var
        per band, each ``(y, x, time)``; coords ``x``, ``y``, ``time`` (decimal
        year) and ``time_dt`` (datetime64); attrs ``epsg``, ``bands``, ``file``,
        plus the NSIDC provenance attrs (``nsidc_product``, ``grimp_box``, ...)
        the writer recorded.

    Raises:
        FileNotFoundError: If the file cannot be found.
        KeyError: If ``varname`` is not a variable in the file.

    Example:
        >>> import moc_py as mp
        >>> from moc_py.obs_nsidc import load_nsidc_netcdf
        >>> O = load_nsidc_netcdf("Jakobshavn_NSIDC-0731.nc", bands=["vv"])
        >>> roi = mp.roi_from_obs(O, margin=2000)
    """
    path = resolve(ncpath, "obs_netcdf_dir")
    ds = xr.open_dataset(path, chunks=chunks)
    if varname not in ds.data_vars:
        candidates = [v for v in ds.data_vars if v != "spatial_ref"]
        if len(candidates) != 1:
            raise KeyError(
                f'Variable "{varname}" not in {path}; have {list(ds.data_vars)}'
            )
        varname = candidates[0]

    out = nsidc_to_obs(ds[varname], bands=bands, trange=trange,
                       epsg=_file_epsg(ds))
    out.attrs.update({k: v for k, v in ds.attrs.items() if k != "epsg"})
    out.attrs["file"] = str(path)
    return out


def nsidc_to_obs(
    source,
    bands: list[str] | tuple[str, ...] | None = None,
    trange: tuple[float, float] | None = None,
    epsg: int | None = None,
    clean: bool = True,
) -> xr.Dataset:
    """Convert a band-stacked NSIDC velocity array into a moc obs Dataset.

    Splits the ``band`` axis into one data var per band, puts the dims in moc's
    ``(y, x, time)`` order, and builds the decimal-year ``time`` coord that the
    model side is indexed on.

    Args:
        source: A band-stacked :class:`xarray.DataArray` — e.g.
            ``series.subset`` from a ``nisardev.nisarVelSeries``, or the
            ``velocity`` variable of a saved file — or a
            :class:`~xarray.Dataset` containing one.
        bands: Bands to keep. ``None`` -> every band present.
        trange: ``(t0, t1)`` decimal-year window to keep.
        epsg: Projection code. ``None`` -> from ``.rio.crs``/attrs, else
            ``CONFIG.epsg``.
        clean: Re-apply the GrIMP sentinel masking (see
            :func:`normalize_nodata`). Harmless on an already-clean file and
            the only protection for a file written by an older workflow.
            Default True.

    Returns:
        Dataset in the schema of :func:`moc_py.load_obs_netcdf`.

    Raises:
        ValueError: If no band axis (``band`` or ``component``) can be found.
    """
    da = source[_velocity_var(source)] if isinstance(source, xr.Dataset) else source

    band_dim = next((d for d in BAND_DIMS if d in da.dims), None)
    if band_dim is None:
        raise ValueError(
            f"No band axis among {BAND_DIMS} in dims {da.dims}; this does not "
            "look like a GrIMP/nisardev velocity stack."
        )
    if clean:
        da = normalize_nodata(da, band_dim=band_dim)

    labels = [str(b) for b in np.asarray(da[band_dim].values).ravel()]
    want = [b for b in (bands or labels) if b in labels]
    if bands and not want:
        raise KeyError(f"None of bands={list(bands)} present; file has {labels}")

    da = _with_decyear_time(da)
    if trange is not None and "time" in da.coords:
        tvals = np.asarray(da["time"].values, dtype="float64")
        da = da.isel(time=np.flatnonzero((tvals >= trange[0]) & (tvals <= trange[1])))

    order = ["y", "x"] + (["time"] if "time" in da.dims else [])
    arrays = {}
    for b in want:
        field = da.isel({band_dim: labels.index(b)}, drop=True).transpose(*order, ...)
        field.name = b
        field.attrs = {k: v for k, v in da.attrs.items() if k != "_FillValue"}
        field.attrs["band"] = b
        arrays[b] = field

    ds = to_obs_dataset(arrays)
    ds.attrs["epsg"] = int(epsg) if epsg is not None else _infer_epsg(da)
    return ds


def normalize_nodata(da: xr.DataArray, band_dim: str = "band") -> xr.DataArray:
    """Replace GrIMP no-data sentinels with NaN so NaN is the only no-data value.

    Mirrors ``glacier_reader.normalize_nodata`` so a file can be cleaned on the
    way out as well as on the way in.

    Args:
        da: Band-stacked velocity array whose ``band_dim`` labels are GrIMP band
            names (``vx``, ``vy``, ``vv``, ``ex``, ``ey``, ``ev``, ``dT``).
        band_dim: Name of the band axis. Default ``"band"``.

    Returns:
        ``float32`` copy in which every value at or below the band's sentinel,
        and every non-finite value, is NaN. The derived speed bands in
        :data:`DERIVED_BANDS` inherit their components' gaps instead of being
        thresholded, because ``sqrt((-2e9)**2 + (-2e9)**2)`` is a large
        *positive* number that no threshold on speed can catch.
    """
    labels = [str(b) for b in np.asarray(da[band_dim].values).ravel()]
    thresholds = xr.DataArray(
        [NO_DATA_SENTINELS.get(b, -np.inf) for b in labels],
        coords={band_dim: da[band_dim]}, dims=band_dim,
    )
    cleaned = da.where(np.isfinite(da) & (da > thresholds))

    for derived, (bx, by) in DERIVED_BANDS.items():
        if {derived, bx, by} <= set(labels):
            components_ok = cleaned.sel({band_dim: [bx, by]}).notnull().all(band_dim)
            cleaned = cleaned.where((cleaned[band_dim] != derived) | components_ok)

    cleaned = cleaned.astype("float32")
    stale = [c for c in cleaned.coords if str(c) == "_FillValue"]
    if stale:
        cleaned = cleaned.drop_vars(stale)
    cleaned.attrs = {k: v for k, v in da.attrs.items()
                     if k not in ("nodata", "_FillValue")}
    cleaned.attrs["nodata_convention"] = "NaN"
    return cleaned


def bbox_to_roi(bbox: dict[str, float], name: str = "", epsg: int | None = None):
    """Convert a picked bounding box into a moc :class:`~moc_py.roi.Roi`.

    The Greenland ROI picker in ``glacier_reader`` returns a plain
    ``{'minx', 'miny', 'maxx', 'maxy'}`` dict; this turns it into the rectangle
    ROI the rest of moc (and, via ``save_roi``, the MATLAB side) works with.

    Args:
        bbox: ``{'minx', 'miny', 'maxx', 'maxy'}`` in EPSG:3413 metres — e.g.
            ``selector.bbox`` or ``out["bbox"]`` from
            ``glacier_reader.glacier_velocity_workflow``.
        name: ROI name. If given, the ROI is also saved to ``CONFIG.roi_dir``
            (and an ISSM ``.exp`` contour alongside it).
        epsg: Projection code. ``None`` -> ``CONFIG.epsg``.

    Returns:
        A rectangular :class:`~moc_py.roi.Roi`.

    Example:
        >>> import glacier_reader as gr
        >>> from moc_py.obs_nsidc import bbox_to_roi
        >>> sel = gr.select_roi_interactive()      # drag a box
        >>> roi = bbox_to_roi(sel.bbox, name=sel.name)
    """
    from .roi import Roi, save_roi

    x0, x1 = float(bbox["minx"]), float(bbox["maxx"])
    y0, y1 = float(bbox["miny"]), float(bbox["maxy"])
    roi = Roi(name=name, kind="rectangle",
              x=np.array([x0, x1, x1, x0, x0]),
              y=np.array([y0, y0, y1, y1, y0]),
              epsg=int(epsg) if epsg is not None else int(CONFIG.epsg))
    return save_roi(roi) if name else roi


def fetch_nsidc(
    center: tuple[float, float] | None = None,
    dx_m: float | None = None,
    dy_m: float | None = None,
    first_date: str = "2015-01-01",
    last_date: str | None = None,
    products: list[str] | tuple[str, ...] = ("NSIDC-0725",),
    prefix: str | None = "glacier",
    out_dir: str | Path | None = None,
    bands: list[str] | tuple[str, ...] | None = None,
    reader_dir: str | Path | None = None,
    **kwargs,
) -> dict[str, xr.Dataset]:
    """Download an NSIDC velocity ROI and return it as moc obs Datasets.

    Thin delegation to ``glacier_reader.glacier_velocity_workflow`` (the NSIDC
    side) followed by :func:`load_nsidc_netcdf` on each file written. Local
    files are reused: a product whose netCDF already exists is not
    re-downloaded unless ``overwrite=True`` is passed through.

    Args:
        center: ROI centre ``(x, y)``, EPSG:3413 metres unless
            ``center_crs=`` is passed through. ``None`` opens the Greenland
            overview map (annual velocity mosaic over the standard outline) so
            the ROI can be drawn — use when the ROI is not known in advance.
        dx_m: Metres to add either side of the centre in x. ``None`` also opens
            the picker.
        dy_m: Same for y. ``None`` -> ``dx_m`` (square ROI).
        first_date: Inclusive NSIDC search start, ``"YYYY-MM-DD"``.
        last_date: Inclusive search end. ``None`` -> today.
        products: NSIDC product codes, e.g. ``["NSIDC-0725", "NSIDC-0481"]``.
        prefix: Leading part of each file name, typically the glacier name.
            ``None`` looks up the glacier the ROI covers in the named flux
            gates and asks for confirmation before using it.
        out_dir: Where the netCDFs are written. ``None`` ->
            ``CONFIG.obs_netcdf_dir``, so they land beside the other subsets.
        bands: Bands to return. ``None`` -> all.
        reader_dir: Directory holding ``glacier_reader.py``, if it is not
            importable already. ``None`` -> try the import as-is, then
            ``~/jupyter/glacierReader``.
        **kwargs: Passed to ``download_velocity_netcdfs``
            (``center_crs``, ``num_workers``, ``overwrite``, ``retries``, ...).

    Returns:
        ``{file stem: obs Dataset}``, one entry per product (per GrIMP box for
        the boxed products NSIDC-0481/0646), in the order written.

    Raises:
        ImportError: If ``glacier_reader`` (and its ``grimpfunc``/``nisardev``
            dependencies) cannot be imported. The download needs that
            environment; reading an already-downloaded file does not.

    Note:
        Needs NASA Earthdata credentials — call ``glacier_reader.login()`` once
        in the session first.
    """
    gr = _import_glacier_reader(reader_dir)
    target = Path(out_dir) if out_dir is not None else CONFIG.get("obs_netcdf_dir")
    out = gr.glacier_velocity_workflow(
        center=center, dx_m=dx_m, dy_m=dy_m, first_date=first_date,
        last_date=last_date, products=products, out_dir=target, prefix=prefix,
        **kwargs,
    )
    return {Path(f).stem: load_nsidc_netcdf(f, bands=bands) for f in out["files"]}


# ---------------------------------------------------------------------------
def _import_glacier_reader(reader_dir: str | Path | None):
    """Import the standalone ``glacier_reader`` module, searching ``reader_dir``."""
    candidates = [reader_dir] if reader_dir else [None, Path.home() / "jupyter" / "glacierReader"]
    last: Exception | None = None
    for cand in candidates:
        if cand is not None:
            path = str(Path(cand).expanduser())
            if path not in sys.path:
                sys.path.insert(0, path)
        try:
            return importlib.import_module("glacier_reader")
        except ImportError as err:  # module or its heavy deps absent
            last = err
    raise ImportError(
        "Could not import glacier_reader (the NSIDC download side). Pass "
        "reader_dir=... or run from an environment with grimpfunc + nisardev "
        f"installed; see glacierReader/pixi.toml. Original error: {last}"
    ) from last


def _velocity_var(ds: xr.Dataset) -> str:
    """Name of the band-stacked velocity variable in a Dataset."""
    if "velocity" in ds.data_vars:
        return "velocity"
    for name, var in ds.data_vars.items():
        if any(d in var.dims for d in BAND_DIMS):
            return str(name)
    raise KeyError(f"No band-stacked velocity variable in {list(ds.data_vars)}")


def _with_decyear_time(da: xr.DataArray) -> xr.DataArray:
    """Add decimal-year ``time`` and ``time_dt`` coords (nisardev stores datetimes)."""
    if "time" not in da.coords:
        return da
    tvals = np.asarray(da["time"].values)
    if np.issubdtype(tvals.dtype, np.datetime64):
        time_dt = tvals.astype("datetime64[ns]")
        da = da.assign_coords(time=("time", datetime64_to_decyear(time_dt)),
                              time_dt=("time", time_dt))
    return da


def _file_epsg(ds: xr.Dataset) -> int | None:
    """EPSG recorded by ``save_velocity_netcdf``, if present."""
    try:
        return int(ds.attrs["epsg"])
    except (KeyError, TypeError, ValueError):
        return None


def _infer_epsg(da: xr.DataArray) -> int:
    """Best-effort EPSG from a rioxarray CRS or an ``epsg`` attr."""
    try:
        crs = da.rio.crs
        if crs is not None and crs.to_epsg():
            return int(crs.to_epsg())
    except Exception:  # noqa: BLE001 - accessor missing / no CRS set
        pass
    try:
        return int(da.attrs["epsg"])
    except (KeyError, TypeError, ValueError):
        return int(CONFIG.epsg)
