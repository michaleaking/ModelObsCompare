"""Lazy / cloud-friendly observed-velocity reading as ``xarray.DataArray``.

The eager :func:`moc_py.load_obs_netcdf` pulls a whole local netCDF into memory
via h5py. This module is the complement: it produces (or wraps) a **lazy**
``xarray.DataArray`` — dask-backed, nothing read until you slice/compute — which
is what you want for cloud reads (remote COG stacks, zarr datacubes, HTTP
netCDF) such as the ``myVelSeries.vv`` array built in the Flowlines notebook.

Two entry points:

- :func:`open_obs_dataarray` — normalise ANY of {an existing DataArray, an
  existing Dataset, a path/URL/glob} into a lazy DataArray with moc_py
  conventions: dims ``(y, x, time)``, coords ``x``/``y``, ``time`` (decimal
  year) + ``time_dt`` (datetime64), and an ``epsg`` attr.
- :func:`to_obs_dataset` — wrap normalised DataArray(s) into a Dataset with the
  same schema as :func:`load_obs_netcdf`, so the ROI / spatial_diff /
  compare_timeseries / plotting functions consume it unchanged. Laziness is
  preserved: those functions realise only the slices they touch.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from .config import CONFIG
from .util import cf_time_to_datetime64, datetime64_to_decyear

_X_NAMES = ("x", "X", "easting", "xc")
_Y_NAMES = ("y", "Y", "northing", "yc")
_T_NAMES = ("time", "mid_date", "t", "date")


def open_obs_dataarray(source, band: str = "vv", varname: str = "velocity",
                       chunks="auto", x_dim: str | None = None,
                       y_dim: str | None = None, time_dim: str | None = None,
                       epsg: int | None = None, engine: str | None = None) -> xr.DataArray:
    """Open/normalise an observed-velocity source into a lazy ``DataArray``.

    Args:
        source: One of
            * ``xarray.DataArray`` — e.g. ``myVelSeries.vv`` (already a field), or
            * ``xarray.Dataset`` — a band is selected (see ``band``/``varname``), or
            * ``str``/``Path`` — a local/remote file, ``.zarr`` store, or glob;
              opened lazily with ``chunks``.
        band: Band to select when the source has a ``component`` dim or multiple
            data vars (default ``"vv"``; use ``"v"`` for itslive-style cubes).
        varname: Velocity variable name to look for in a Dataset (default
            ``"velocity"``); ignored if ``band`` is itself a data var.
        chunks: Dask chunking passed to ``xr.open_*`` / ``.chunk`` (default
            ``"auto"``; ``None`` disables chunking). Requires ``dask`` to stay lazy.
        x_dim, y_dim, time_dim: Override dimension names (else auto-detected).
        epsg: Projection code; if ``None`` inferred from ``.rio.crs`` / attrs /
            ``CONFIG.epsg``.
        engine: Optional xarray backend engine for remote/HTTP netCDF.

    Returns:
        A lazy ``xarray.DataArray`` with dims ``(y, x, time)``, coords ``x``,
        ``y``, ``time`` (decimal year) and ``time_dt`` (datetime64), plus attrs
        ``epsg`` and ``band``. **Not** loaded into memory.
    """
    obj = _open_source(source, chunks=chunks, engine=engine)
    da = _select_band(obj, band=band, varname=varname)

    xd = x_dim or _find_dim(da, _X_NAMES, "x")
    yd = y_dim or _find_dim(da, _Y_NAMES, "y")
    td = time_dim or _find_dim(da, _T_NAMES, "time", required=False)

    # canonical dim order (y, x, time); keep time optional (single-frame ok)
    order = [yd, xd] + ([td] if td else [])
    da = da.transpose(*order, ...)

    # standardise coord/dim names
    rename = {}
    if xd != "x":
        rename[xd] = "x"
    if yd != "y":
        rename[yd] = "y"
    if td and td != "time":
        rename[td] = "time"
    if rename:
        da = da.rename(rename)

    da = _attach_time(da)
    da.attrs["epsg"] = int(epsg) if epsg is not None else _infer_epsg(da, obj)
    da.attrs["band"] = band
    da.name = band
    return da


def to_obs_dataset(arrays, band: str | None = None) -> xr.Dataset:
    """Wrap normalised DataArray(s) into an obs Dataset (schema of load_obs_netcdf).

    Args:
        arrays: A single :class:`xarray.DataArray` from :func:`open_obs_dataarray`,
            or a list/tuple/dict of them (e.g. vv/vx/vy) sharing coords.
        band: Name to use for a single unnamed DataArray (default: its ``.name``
            or ``"vv"``).

    Returns:
        A lazy :class:`xarray.Dataset` with one data var per band, dims
        ``(y, x, time)``, coords ``x``/``y``/``time``/``time_dt``, and an
        ``epsg`` attr — directly usable by ``spatial_diff`` / ``plot_map`` /
        ``compare_timeseries``.
    """
    if isinstance(arrays, xr.DataArray):
        name = band or arrays.name or "vv"
        ds = arrays.to_dataset(name=name)
    elif isinstance(arrays, dict):
        ds = xr.Dataset({k: v for k, v in arrays.items()})
    else:  # list/tuple
        ds = xr.Dataset({(a.name or f"band{i}"): a for i, a in enumerate(arrays)})
    first = arrays if isinstance(arrays, xr.DataArray) else list(
        arrays.values() if isinstance(arrays, dict) else arrays)[0]
    ds.attrs["epsg"] = int(first.attrs.get("epsg", CONFIG.epsg))
    ds.attrs["bands"] = ",".join(ds.data_vars)
    return ds


# ---------------------------------------------------------------------------
def _open_source(source, chunks, engine):
    """Return an xarray object for a DataArray/Dataset/path, opened lazily."""
    if isinstance(source, (xr.DataArray, xr.Dataset)):
        obj = source
        if chunks is not None:
            try:
                obj = obj.chunk(chunks)
            except (ImportError, ValueError):
                pass  # dask absent or already chunked; leave as-is
        return obj

    src = str(source)
    if src.endswith(".zarr") or src.rstrip("/").endswith(".zarr"):
        return xr.open_zarr(src, chunks=chunks)
    # local glob of rasters/netcdf -> multi-file open along a new/So concat dim
    if any(ch in src for ch in "*?[") or isinstance(source, (list, tuple)):
        return xr.open_mfdataset(source, chunks=chunks, engine=engine,
                                 combine="by_coords")
    return xr.open_dataset(src, chunks=chunks, engine=engine)


def _select_band(obj, band, varname):
    """Pick a single velocity field (DataArray) from a DataArray/Dataset."""
    if isinstance(obj, xr.DataArray):
        da = obj
    elif band in obj.data_vars:
        da = obj[band]
    elif varname in obj.data_vars:
        da = obj[varname]
    elif len(obj.data_vars) == 1:
        da = obj[next(iter(obj.data_vars))]
    else:
        raise KeyError(
            f'Cannot find band "{band}" or var "{varname}" in Dataset with '
            f"vars {list(obj.data_vars)}; pass band=... explicitly."
        )
    # resolve a component dimension if present
    if "component" in da.dims:
        comp = _component_labels(da, obj)
        if band in comp:
            da = da.isel(component=comp.index(band))
        else:
            raise KeyError(f'Band "{band}" not in component labels {comp}')
    return da


def _component_labels(da, obj):
    """Component band names from a coord, else the configured default order."""
    for src in (da, obj):
        if "component" in getattr(src, "coords", {}):
            vals = np.asarray(src["component"].values).ravel()
            return [v.decode() if isinstance(v, (bytes, bytearray)) else str(v)
                    for v in vals]
    return list(CONFIG.obs_bands)


def _find_dim(da, candidates, kind, required=True):
    """Return the first matching dim name (by common aliases), else raise/None."""
    for c in candidates:
        if c in da.dims:
            return c
    # fall back: a coordinate that maps to exactly one dim
    for c in candidates:
        if c in da.coords and da[c].ndim == 1:
            return da[c].dims[0]
    if required:
        raise ValueError(
            f"Could not identify the {kind} dimension among {da.dims}; "
            f"pass {kind}_dim=... explicitly."
        )
    return None


def _attach_time(da):
    """Ensure a decimal-year ``time`` coord plus a ``time_dt`` datetime coord."""
    if "time" not in da.coords:
        return da
    tvals = da["time"].values
    if np.issubdtype(np.asarray(tvals).dtype, np.datetime64):
        time_dt = np.asarray(tvals, dtype="datetime64[ns]")
        time_yr = datetime64_to_decyear(time_dt)
    else:
        units = da["time"].attrs.get("units", "")
        if "since" in str(units):
            time_dt = cf_time_to_datetime64(tvals, units)
            time_yr = datetime64_to_decyear(time_dt)
        else:  # already decimal years
            time_yr = np.asarray(tvals, dtype="float64")
            time_dt = None
    da = da.assign_coords(time=("time", time_yr))
    if time_dt is not None:
        da = da.assign_coords(time_dt=("time", time_dt))
    return da


def _infer_epsg(da, obj):
    """Best-effort EPSG from a rioxarray CRS or common attrs."""
    for src in (da, obj):
        try:
            crs = src.rio.crs  # only if rioxarray is installed & CRS is set
            if crs is not None:
                code = crs.to_epsg()
                if code:
                    return int(code)
        except Exception:  # noqa: BLE001 - accessor missing / no CRS
            pass
        for key in ("epsg", "EPSG"):
            if key in getattr(src, "attrs", {}):
                try:
                    return int(src.attrs[key])
                except (TypeError, ValueError):
                    pass
    return CONFIG.epsg
