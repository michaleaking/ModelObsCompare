"""Reader for gridded observed-velocity netCDF files (Sentinel/CSK subsets).

These files store a 4-D ``velocity(time, component, y, x)`` array whose component
dimension holds bands like ``{vv, vx, vy, ex, ey, dT}``. They are netCDF4 (HDF5)
files, so this reader uses :mod:`h5py` directly and needs neither the
``netCDF4`` nor ``rioxarray`` packages. The reader is robust to dimension ORDER
(axes are located by size) and returns an :class:`xarray.Dataset` with each band
laid out as ``(y, x, time)``.

Mirrors ``io/moc_load_obs_netcdf.m``.
"""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np
import xarray as xr

from .config import CONFIG, resolve
from .util import cf_time_to_datetime64, datetime64_to_decyear


def load_obs_netcdf(
    ncpath: str | Path,
    bands: list[str] | tuple[str, ...] | None = None,
    varname: str = "velocity",
    trange: tuple[float, float] | None = None,
) -> xr.Dataset:
    """Load a gridded observed-velocity netCDF into an :class:`xarray.Dataset`.

    Args:
        ncpath: Full path or bare filename (resolved against
            ``CONFIG.obs_netcdf_dir``).
        bands: Subset of component bands to return. ``None`` -> all present.
        varname: Name of the velocity variable (default ``"velocity"``).
        trange: ``(t0, t1)`` decimal-year window to keep. ``None`` -> all frames.

    Returns:
        Dataset with dims ``(y, x, time)``:
            data vars: one per band (e.g. ``vv``, ``vx``, ``vy`` ...), each
            ``(y, x, time)`` in float32.
            coords: ``x`` (nx), ``y`` (ny), ``time`` (decimal year),
            ``time_dt`` (datetime64).
            attrs: ``epsg``, ``geotransform`` (if present), ``file``.

    Raises:
        FileNotFoundError: If the file cannot be found.
        KeyError: If ``varname`` is not in the file.
    """
    path = resolve(ncpath, "obs_netcdf_dir")

    with h5py.File(path, "r") as h:
        if varname not in h:
            raise KeyError(f'Variable "{varname}" not in {path}')

        x = np.asarray(h["x"][:], dtype="float64").ravel()
        y = np.asarray(h["y"][:], dtype="float64").ravel()
        t_raw = np.asarray(h["time"][:]).ravel()
        t_units = _attr_str(h["time"], "units")
        time_dt = cf_time_to_datetime64(t_raw, t_units)
        time_yr = datetime64_to_decyear(time_dt)

        band_names = _component_names(h)
        want = list(bands) if bands else list(band_names)

        vel = h[varname]
        shape = vel.shape
        role = _map_axes(shape, len(x), len(y), len(t_raw), len(band_names))
        fill = _fill_value(vel)  # explicit _FillValue/missing_value, if any

        # time subset (mask into the time axis)
        if trange is not None:
            tkeep = (time_yr >= trange[0]) & (time_yr <= trange[1])
        else:
            tkeep = np.ones(time_yr.size, dtype=bool)

        data_vars = {}
        for b in want:
            if b not in band_names:
                continue
            ci = band_names.index(b)
            sl = [slice(None)] * len(shape)
            sl[role["comp"]] = ci
            arr = np.asarray(vel[tuple(sl)], dtype="float32")  # component axis removed
            arr = _to_yxt(arr, shape, role)  # -> (y, x, time)
            arr = arr[:, :, tkeep]
            _mask_fill(arr, fill, band=b)  # nodata sentinels -> NaN (in place)
            data_vars[b] = (("y", "x", "time"), arr)

        geotransform = _attr_str(h.get("spatial_ref"), "GeoTransform") if "spatial_ref" in h else ""
        epsg = _parse_epsg(h)

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "x": ("x", x),
            "y": ("y", y),
            "time": ("time", time_yr[tkeep]),
            "time_dt": ("time", time_dt[tkeep]),
        },
        attrs={
            "epsg": int(epsg) if epsg else CONFIG.epsg,
            "geotransform": geotransform,
            "file": str(path),
            "bands": ",".join(k for k in data_vars),
        },
    )
    return ds


# ---------------------------------------------------------------------------
def _component_names(h: h5py.File) -> list[str]:
    """Read the 'component' band labels, falling back to the configured order."""
    if "component" in h:
        vals = h["component"][:]
        out = []
        for v in np.asarray(vals).ravel():
            out.append(v.decode() if isinstance(v, (bytes, bytearray)) else str(v))
        if out:
            return out
    return list(CONFIG.obs_bands)


def _map_axes(shape, nx, ny, nt, nc) -> dict:
    """Locate x / y / time / component axes of ``velocity`` by matching sizes."""
    role = {}
    for name, n in (("x", nx), ("y", ny), ("t", nt), ("comp", nc)):
        matches = [i for i, s in enumerate(shape) if s == n]
        if len(matches) != 1:
            raise ValueError(
                f"Cannot unambiguously match axis '{name}' (size {n}) in "
                f"velocity shape {shape}; sizes may be non-unique."
            )
        role[name] = matches[0]
    return role


def _to_yxt(arr: np.ndarray, shape, role) -> np.ndarray:
    """Transpose a single-band array (component axis removed) to ``(y, x, time)``."""
    # Axes remaining after removing the component axis, in original order.
    remaining = [ax for ax in range(len(shape)) if ax != role["comp"]]
    pos = {ax: i for i, ax in enumerate(remaining)}
    return np.transpose(arr, (pos[role["y"]], pos[role["x"]], pos[role["t"]]))


def _fill_value(dset):
    """Return an explicit fill/missing sentinel from a variable's attrs, or None."""
    for name in ("_FillValue", "missing_value"):
        if name in getattr(dset, "attrs", {}):
            v = np.asarray(dset.attrs[name]).ravel()
            if v.size:
                return float(v[0])
    return None


#: Bands that cannot physically be negative: ``vv`` is a speed magnitude,
#: ``ex``/``ey`` are error magnitudes, ``dT`` an image-pair separation. These
#: products flag their nodata with ``-1`` on these bands (and ``-2e9`` on the
#: signed ``vx``/``vy``), so any negative is a sentinel, not data.
NONNEGATIVE_BANDS = frozenset({"vv", "ex", "ey", "dT"})


def _mask_fill(arr: np.ndarray, fill, band: str = "") -> None:
    """Set fill/sentinel values to NaN in place.

    Catches an explicit ``_FillValue`` plus non-physical magnitudes (>=1e9),
    since some products use a large negative sentinel (e.g. -2e9) that is not
    declared as ``_FillValue`` and would otherwise leak through raw HDF5 reads.
    On a band in :data:`NONNEGATIVE_BANDS`, negatives are sentinels too (these
    files use ``-1``) and are masked as well.

    Args:
        arr: Band array, modified in place.
        fill: Declared fill value, or None.
        band: Band name, used to decide whether negatives are physical.
    """
    if fill is not None and np.isfinite(fill):
        arr[arr == np.float32(fill)] = np.nan
    arr[np.abs(arr) >= 1e9] = np.nan
    if band in NONNEGATIVE_BANDS:
        arr[arr < 0] = np.nan


def _attr_str(dset, name: str) -> str:
    """Return a string attribute from an h5py object, decoding bytes."""
    if dset is None or name not in getattr(dset, "attrs", {}):
        return ""
    v = dset.attrs[name]
    if isinstance(v, (bytes, bytearray)):
        return v.decode()
    return str(v)


def _parse_epsg(h: h5py.File):
    """Best-effort EPSG code from the spatial_ref WKT."""
    if "spatial_ref" not in h:
        return None
    for key in ("spatial_ref", "crs_wkt"):
        wkt = _attr_str(h["spatial_ref"], key)
        m = re.search(r'EPSG","(\d+)"\]\]\s*$', wkt)
        if m:
            return int(m.group(1))
    return None
