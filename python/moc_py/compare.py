"""Model-vs-observation comparison on common grids and at points (Phase 2).

Builds on :mod:`moc_py.model` (mesh interpolation) and the observation readers
to produce spatial-difference maps and paired time series with error metrics.
Mirrors ``compare/moc_spatial_diff.m`` and ``compare/moc_compare_timeseries.m``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from . import model as _model


def spatial_diff(M, O, roi, field: str = "vel", obsband: str = "vv",
                 res: float | None = None, time: float | None = None,
                 obswindow: float = 0.0) -> xr.Dataset:
    """Model-minus-observation speed difference on a common ROI grid.

    Regrids the model mesh field and the gridded observations onto the SAME
    regular grid over an ROI, then differences them. The observed frame(s)
    nearest ``time`` (optionally averaged over ``+/- obswindow`` years) are used.

    Args:
        M: Model from :func:`moc_py.model.load_model`.
        O: Obs grid Dataset from :func:`moc_py.load_obs_netcdf`.
        roi: The region of interest.
        field: Model field to compare (``"vel"``).
        obsband: Obs band to compare (``"vv"`` speed).
        res: Common-grid spacing [m] (default ``CONFIG.grid_res``).
        time: Target decimal year (default: model's last step).
        obswindow: +/- years of obs frames to average (0 = nearest frame).

    Returns:
        Dataset with data vars ``model``, ``obs``, ``diff`` (dims ``y, x``),
        coords ``x``/``y``, and attrs including ``rmse``/``mae``/``bias``/``n``.
    """
    G = _model.model_to_grid(M, roi, res=res, field=field, time=time)
    Xq, Yq = np.meshgrid(G.x.values, G.y.values)

    ttarget = time
    if ttarget is None:
        ttarget = G.attrs["time"]
    if ttarget is None or np.isnan(ttarget):
        ttarget = float(np.median(O.time.values))

    obs_on_grid, obstime = _obs_on_grid(O, obsband, ttarget, obswindow, Xq, Yq)

    diff = G["z"].values - obs_on_grid
    valid = np.isfinite(diff)
    dv = diff[valid]
    stats = _stats(dv)

    ds = xr.Dataset(
        {
            "model": (("y", "x"), G["z"].values),
            "obs": (("y", "x"), obs_on_grid),
            "diff": (("y", "x"), diff),
        },
        coords={"x": ("x", G.x.values), "y": ("y", G.y.values)},
        attrs={
            "field": field, "obsband": obsband,
            "time": float(G.attrs["time"]) if G.attrs["time"] is not None else np.nan,
            "obstime_min": float(np.min(obstime)), "obstime_max": float(np.max(obstime)),
            "n": stats["n"], "bias": stats["mean"], "median": stats["median"],
            "rmse": stats["rmse"], "mae": stats["mae"], "epsg": M.epsg,
        },
    )
    return ds


@dataclass
class TimeSeriesComparison:
    """Paired modelled/observed series at one location (see :func:`compare_timeseries`)."""

    xy: tuple
    model_time: np.ndarray
    model_val: np.ndarray
    obs_time: np.ndarray
    obs_val: np.ndarray
    obs_err: np.ndarray | None
    model_at_obs: np.ndarray
    stats: dict


def compare_timeseries(M, obs, xy=None, dist: float | None = None,
                       field: str = "vel", obsband: str = "vv") -> TimeSeriesComparison:
    """Align a modelled and observed velocity series at one location.

    Args:
        M: Model from :func:`moc_py.model.load_model`.
        obs: An obs-flowline Dataset (dims ``point, time``) or a gridded obs
            Dataset (dims ``y, x, time``).
        xy: ``(x, y)`` location [m]. Required for gridded obs; optional for a
            flowline (nearest point).
        dist: Along-flowline distance [m] to sample (flowline obs only).
        field: Model field (``"vel"``).
        obsband: Obs band for gridded obs (``"vv"``).

    Returns:
        A :class:`TimeSeriesComparison`.
    """
    is_flow = "point" in getattr(obs, "dims", {})
    if is_flow:
        xy, obs_time, obs_val, obs_err = _flowline_pick(obs, xy, dist)
    else:
        xy, obs_time, obs_val, obs_err = _grid_pick(obs, xy, obsband)

    model_val = _model.model_at_points(M, xy[0], xy[1], field=field)[0]
    model_time = M.time

    if model_time.size >= 2 and not np.all(np.isnan(model_time)):
        model_at_obs = np.interp(obs_time, model_time, model_val,
                                 left=np.nan, right=np.nan)
    else:
        model_at_obs = np.full_like(obs_time, np.nan, dtype="float64")

    good = np.isfinite(model_at_obs) & np.isfinite(obs_val)
    resid = model_at_obs[good] - obs_val[good]
    st = _stats(resid)
    if good.sum() >= 2:
        st["r"] = float(np.corrcoef(model_at_obs[good], obs_val[good])[0, 1])
    else:
        st["r"] = np.nan

    return TimeSeriesComparison(
        xy=(float(xy[0]), float(xy[1])),
        model_time=model_time, model_val=model_val,
        obs_time=obs_time, obs_val=obs_val, obs_err=obs_err,
        model_at_obs=model_at_obs, stats=st,
    )


# ---------------------------------------------------------------------------
def _obs_on_grid(O, band, ttarget, window, Xq, Yq):
    """Select obs frame(s) near ttarget and interpolate onto (Xq, Yq)."""
    from scipy.interpolate import RegularGridInterpolator

    t = O.time.values
    if window > 0:
        sel = np.abs(t - ttarget) <= window
        if not sel.any():
            sel = np.zeros_like(t, dtype=bool)
            sel[np.argmin(np.abs(t - ttarget))] = True
    else:
        sel = np.zeros_like(t, dtype=bool)
        sel[np.argmin(np.abs(t - ttarget))] = True

    slc = O[band].isel(time=np.where(sel)[0]).mean("time").values  # (y, x)
    xo = O.x.values
    yo = O.y.values
    oy = np.argsort(yo)
    ox = np.argsort(xo)
    rgi = RegularGridInterpolator((yo[oy], xo[ox]), slc[oy][:, ox],
                                  bounds_error=False, fill_value=np.nan)
    pts = np.column_stack([Yq.ravel(), Xq.ravel()])
    return rgi(pts).reshape(Xq.shape), t[sel]


def _flowline_pick(F, xy, dist):
    """Choose a flowline point by distance or nearest coordinate."""
    if dist is not None:
        i = int(np.argmin(np.abs(F.d.values - dist)))
    elif xy is not None:
        i = int(np.argmin(np.hypot(F.x.values - xy[0], F.y.values - xy[1])))
    else:
        i = 0
    sub = F.isel(point=i)
    err = sub.err.values if "err" in sub else None
    return (float(sub.x), float(sub.y)), F.time.values, sub.vel.values, err


def _grid_pick(O, xy, band):
    """Nearest-pixel time series from a gridded obs Dataset."""
    if xy is None:
        raise ValueError("For gridded obs you must pass xy=(x, y).")
    ci = int(np.argmin(np.abs(O.x.values - xy[0])))
    ri = int(np.argmin(np.abs(O.y.values - xy[1])))
    v = O[band].isel(x=ci, y=ri).values
    err = None
    if "ex" in O and "ey" in O:
        err = np.hypot(O["ex"].isel(x=ci, y=ri).values, O["ey"].isel(x=ci, y=ri).values)
    return (float(O.x.values[ci]), float(O.y.values[ri])), O.time.values, v, err


def _stats(resid: np.ndarray) -> dict:
    """Bias/median/rmse/mae over a residual vector (NaNs already removed)."""
    resid = resid[np.isfinite(resid)]
    if resid.size == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan, "rmse": np.nan, "mae": np.nan}
    return {
        "n": int(resid.size),
        "mean": float(np.mean(resid)),
        "median": float(np.median(resid)),
        "rmse": float(np.sqrt(np.mean(resid ** 2))),
        "mae": float(np.mean(np.abs(resid))),
    }
