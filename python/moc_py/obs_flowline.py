"""Reader for per-glacier observed flowline ``.mat`` files (redrawn_Jan25).

These are MATLAB v7 files (readable by :func:`scipy.io.loadmat`). Each file
describes ONE glacier flowline: ordered ``(x, y)`` points, along-line distance
``d``, static geometry, and one or more observed velocity time series stored as
points-by-time (or time-by-points) matrices. This module normalises the common
variants to a single :class:`xarray.Dataset` with dims ``(point, time)``.

Mirrors ``io/moc_load_obs_flowline.m``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.io
import xarray as xr

from .config import CONFIG, resolve
from .util import datenum_to_datetime64, datetime64_to_decyear

#: Velocity sources tried in order when ``source="auto"``.
_SOURCE_ORDER = ("filteredV", "SentinelV", "SentinelVmonthly", "u2")


def load_obs_flowline(matpath: str | Path, source: str = "auto") -> xr.Dataset:
    """Load an observed flowline file into a tidy :class:`xarray.Dataset`.

    Args:
        matpath: Full path or bare filename (resolved against
            ``CONFIG.obs_flowline_dir``).
        source: Which velocity source to extract. One of ``"auto"`` (default),
            ``"filteredV"``, ``"SentinelV"``, ``"SentinelVmonthly"``, ``"u2"``.

    Returns:
        Dataset with dims ``(point, time)``:
            data vars: ``vel`` (m/yr), and where available ``err``, ``vx``, ``vy``.
            coords: ``x``, ``y``, ``d`` (all along ``point``); ``time`` (decimal
            year) and ``time_dt`` (datetime64) along ``time``.
            attrs: ``name``, ``source``, ``epsg``, ``file``.

    Raises:
        FileNotFoundError: If the file cannot be found.
        ValueError: If no usable velocity source is present.
    """
    path = resolve(matpath, "obs_flowline_dir")
    raw = scipy.io.loadmat(path, squeeze_me=False, struct_as_record=False)

    x = _col(raw, "x")
    y = _col(raw, "y")
    d = _col(raw, "d")
    npts = x.size
    if d.size != npts:  # some files store d transposed (1 x npts)
        d = d.ravel()[:npts]

    name = path.stem
    if "glrname" in raw:
        try:
            name = str(np.asarray(raw["glrname"]).ravel()[0]).strip()
        except Exception:  # pragma: no cover - cosmetic only
            pass

    order = _SOURCE_ORDER if source == "auto" else (source,)
    used, vel, err, vx, vy, tdn = "", None, None, None, None, None
    for src in order:
        if src not in raw:
            continue
        vel, err, vx, vy, tdn = _extract(raw, src, npts)
        if vel is not None and tdn is not None:
            used = src
            break

    if not used:
        raise ValueError(
            f"No usable velocity source in {path} (looked for {', '.join(order)})"
        )

    time_dt = datenum_to_datetime64(tdn.ravel())
    time_yr = datetime64_to_decyear(time_dt)

    data_vars = {"vel": (("point", "time"), vel)}
    if err is not None:
        data_vars["err"] = (("point", "time"), err)
    if vx is not None:
        data_vars["vx"] = (("point", "time"), vx)
    if vy is not None:
        data_vars["vy"] = (("point", "time"), vy)

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "x": ("point", x),
            "y": ("point", y),
            "d": ("point", d),
            "time": ("time", time_yr),
            "time_dt": ("time", time_dt),
        },
        attrs={"name": name, "source": used, "epsg": CONFIG.epsg, "file": str(path)},
    )
    return ds


# ---------------------------------------------------------------------------
def _col(raw: dict, key: str) -> np.ndarray:
    """Return ``raw[key]`` as a 1-D float array, or an empty array."""
    if key not in raw:
        return np.empty(0)
    return np.asarray(raw[key], dtype="float64").ravel()


def _getf(struct, field: str):
    """Fetch a field from a scipy mat_struct (struct_as_record=False)."""
    if struct is None:
        return None
    return getattr(struct, field, None)


def _extract(raw: dict, src: str, npts: int):
    """Pull (vel, err, vx, vy, datenum) as points-by-time arrays for one source."""
    vel = err = vx = vy = tdn = None

    if src == "filteredV":
        vel = raw.get("filteredV")
        err = raw.get("filteredVe")
        tdn = raw.get("vti")
    elif src in ("SentinelV", "SentinelVmonthly"):
        st = _first_struct(raw[src])
        vel = _getf(st, "velocity")
        err = _getf(st, "err")
        vx = _getf(st, "vx")
        vy = _getf(st, "vy")
        tdn = _getf(st, "t")
    elif src == "u2":
        st = _first_struct(raw["u2"])
        vel = _getf(st, "v")
        err = _getf(st, "ve")
        tdn = _getf(st, "t")

    if vel is None or tdn is None:
        return None, None, None, None, None

    nt = np.asarray(tdn).size
    vel = _orient(vel, npts, nt)
    err = _orient(err, npts, nt)
    vx = _orient(vx, npts, nt)
    vy = _orient(vy, npts, nt)
    return vel, err, vx, vy, np.asarray(tdn, dtype="float64")


def _first_struct(obj):
    """Unwrap a (1,1) object/struct array to its single mat_struct element."""
    arr = np.asarray(obj)
    flat = arr.ravel()
    return flat[0] if flat.size else None


def _orient(a, npts: int, nt: int):
    """Ensure ``a`` is (npts, nt); transpose (nt, npts) inputs (e.g. ``u2.v``)."""
    if a is None:
        return None
    a = np.asarray(a, dtype="float64")
    if a.ndim == 1:
        a = a.reshape(-1, 1) if a.size == npts else a.reshape(1, -1)
    if a.shape[0] == npts:
        return a
    if a.shape[1] == npts and a.shape[0] == nt:
        return a.T
    if a.shape[0] == nt:  # time on axis 0 -> points on rows
        return a.T
    return a
