"""Plotting helpers for observed velocity grids and flowlines (Phase 1).

These operate on the :class:`xarray.Dataset` objects returned by the readers.
Model overlays are accepted where noted but are optional until Phase 2 adds the
ISSM model export/interpolation.
"""

from __future__ import annotations

import numpy as np


def _new_ax(ax):
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(figsize=(7, 6))
    return ax


def _pctl(v, pc):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, pc)) if v.size else 1.0


def plot_map(ds, band: str = "vv", tindex: int = 0, ax=None, clim=None,
             roi=None, title: str = "", cmap: str = "viridis"):
    """Plot one time slice of a gridded observed field as a map.

    Args:
        ds: Dataset from :func:`load_obs_netcdf`.
        band: Which band to show (default ``"vv"`` speed).
        tindex: Time index to display.
        ax: Existing axes (default: new figure).
        clim: ``(vmin, vmax)`` colour limits (default ``(0, 98th pct)``).
        roi: Optional :class:`~moc_py.roi.Roi` to overlay.
        title: Axis title (default auto).
        cmap: Colormap name.

    Returns:
        The matplotlib ``Axes``.
    """
    ax = _new_ax(ax)
    da = ds[band].isel(time=tindex)
    x = ds.x.values
    y = ds.y.values
    extent = [x.min(), x.max(), y.min(), y.max()]
    origin = "upper" if y[0] > y[-1] else "lower"
    Z = da.values
    vmin, vmax = (0, _pctl(Z, 98)) if clim is None else clim
    im = ax.imshow(Z, extent=extent, origin=origin, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="equal", interpolation="nearest")
    cb = ax.figure.colorbar(im, ax=ax)
    cb.set_label("speed [m/yr]")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if roi is not None:
        ax.plot(roi.x, roi.y, "r-", lw=1.5)
    t = float(ds.time.isel(time=tindex))
    ax.set_title(title or f"obs {band}  t={t:.3f}")
    return ax


def plot_timeseries(Fobs, dist: float | None = None, point: int | None = None,
                    ax=None, use_datetime: bool = True, model=None, title: str = ""):
    """Plot an observed velocity time series at one flowline point.

    Args:
        Fobs: Flowline Dataset from :func:`load_obs_flowline`.
        dist: Along-flowline distance [m] to sample (nearest point).
        point: Explicit point index (overrides ``dist``).
        ax: Existing axes.
        use_datetime: Plot on a datetime x-axis (else decimal year).
        model: Optional ``(t, v)`` tuple of a modelled series to overlay (Phase 2).
        title: Axis title.

    Returns:
        The matplotlib ``Axes``.
    """
    ax = _new_ax(ax)
    if point is None:
        if dist is not None:
            point = int(np.argmin(np.abs(Fobs.d.values - dist)))
        else:
            point = 0
    sub = Fobs.isel(point=point)
    t = sub.time_dt.values if use_datetime else sub.time.values
    v = sub.vel.values
    if "err" in sub:
        e = sub.err.values
        ax.errorbar(t, v, yerr=e, fmt="o", ms=3, color="0.2",
                    ecolor="0.7", capsize=0, label="observed")
    else:
        ax.plot(t, v, "o", ms=3, color="0.2", label="observed")
    if model is not None:
        mt, mv = model
        ax.plot(mt, mv, "-", color="#d81b1b", lw=1.8, label="model")
    ax.set_xlabel("date" if use_datetime else "decimal year")
    ax.set_ylabel("speed [m/yr]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    d_km = float(Fobs.d.isel(point=point)) / 1e3
    ax.set_title(title or f"{Fobs.attrs.get('name','')}  d={d_km:.1f} km")
    return ax


def plot_grid(G, ax=None, clim=None, roi=None, title: str = "", cmap: str = "viridis",
              var: str = "z"):
    """Plot a 2-D gridded field (e.g. model-on-grid from ``model_to_grid``).

    Args:
        G: Dataset with a 2-D data var (dims ``y, x``); default var ``"z"``.
        ax: Existing axes.
        clim: ``(vmin, vmax)`` colour limits (default ``(0, 98th pct)``).
        roi: Optional ROI polygon to overlay.
        title: Axis title.
        cmap: Colormap.
        var: Which data var to plot.

    Returns:
        The matplotlib ``Axes``.
    """
    ax = _new_ax(ax)
    Z = G[var].values
    x, y = G.x.values, G.y.values
    extent = [x.min(), x.max(), y.min(), y.max()]
    origin = "upper" if y[0] > y[-1] else "lower"
    vmin, vmax = (0, _pctl(Z, 98)) if clim is None else clim
    im = ax.imshow(Z, extent=extent, origin=origin, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="equal", interpolation="nearest")
    cb = ax.figure.colorbar(im, ax=ax)
    cb.set_label("speed [m/yr]")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if roi is not None:
        ax.plot(roi.x, roi.y, "r-", lw=1.5)
    t = G.attrs.get("time", None)
    ax.set_title(title or f"model {G.attrs.get('field','')}"
                 + (f"  t={t:.3f}" if isinstance(t, (int, float)) and np.isfinite(t) else ""))
    return ax


def plot_diff(D, clim=None, dlim=None, roi=None, cmap: str = "viridis"):
    """Three-panel model / obs / difference map from :func:`moc_py.compare.spatial_diff`.

    Args:
        D: Dataset from ``spatial_diff`` (vars ``model``, ``obs``, ``diff``).
        clim: Speed colour limits for the model/obs panels (default auto).
        dlim: Symmetric difference limit +/- [m/yr] (default ~95th pct of |diff|).
        roi: Optional ROI polygon to overlay on each panel.
        cmap: Colormap for the speed panels.

    Returns:
        The matplotlib ``Figure``.
    """
    import matplotlib.pyplot as plt

    if clim is None:
        both = np.concatenate([D["model"].values.ravel(), D["obs"].values.ravel()])
        clim = (0, _pctl(both, 98))
    if dlim is None:
        dlim = _pctl(np.abs(D["diff"].values), 95) or 1.0

    x, y = D.x.values, D.y.values
    extent = [x.min(), x.max(), y.min(), y.max()]
    origin = "upper" if y[0] > y[-1] else "lower"
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))

    for ax, key, ttl, cm, lims, lbl in [
        (axs[0], "model", f"model  t={D.attrs['time']:.3f}", cmap, clim, "speed [m/yr]"),
        (axs[1], "obs",
         f"obs  t=[{D.attrs['obstime_min']:.2f}..{D.attrs['obstime_max']:.2f}]",
         cmap, clim, "speed [m/yr]"),
        (axs[2], "diff",
         f"diff  RMSE={D.attrs['rmse']:.0f}  MAE={D.attrs['mae']:.0f}  bias={D.attrs['bias']:.0f}",
         "RdBu_r", (-dlim, dlim), "model - obs [m/yr]"),
    ]:
        im = ax.imshow(D[key].values, extent=extent, origin=origin,
                       vmin=lims[0], vmax=lims[1], cmap=cm, aspect="equal",
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax); cb.set_label(lbl)
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title(ttl)
        if roi is not None:
            ax.plot(roi.x, roi.y, "k-", lw=1.2)
    fig.tight_layout()
    return fig


#: Colorbar labels and default colour scaling per variable of a variance Dataset.
_VAR_STYLE = {
    "norm": ("{stat} / {ref} speed [-]", "magma", False),
    "spread": ("{stat} of speed [m/yr]", "magma", False),
    "reference": ("{ref} speed [m/yr]", "viridis", False),
    "count": ("valid frames [-]", "viridis", False),
    "change": ("later - earlier [m/yr]", "RdBu_r", True),
    "norm_change": ("(later - earlier) / {ref} speed [-]", "RdBu_r", True),
}


def plot_variance(V, var: str = "norm", ax=None, clim=None, roi=None,
                  title: str = "", cmap: str | None = None):
    """Map the per-pixel temporal variability from :func:`moc_py.velocity_variance`.

    Args:
        V: Dataset from ``velocity_variance``.
        var: Which field to map — ``"norm"`` (spread as a fraction of the
            pixel's own speed, the comparable-between-pixels one), ``"spread"``
            (absolute [m/yr]), ``"reference"``, ``"count"``, or, for a two-slice
            window, the signed ``"change"`` / ``"norm_change"``.
        ax: Existing axes (default: new figure).
        clim: ``(vmin, vmax)`` colour limits. Default ``(0, 98th pct)``, or
            symmetric about zero for the signed fields.
        roi: Optional ROI polygon to overlay.
        title: Axis title (default auto, includes the ROI-mean).
        cmap: Colormap (default: sequential, or diverging for signed fields).

    Returns:
        The matplotlib ``Axes``.

    Raises:
        ValueError: If ``var`` is not in ``V``.
    """
    if var not in V:
        raise ValueError(f"{var!r} not in variance Dataset (has {list(V.data_vars)})")

    ax = _new_ax(ax)
    stat = V.attrs.get("statistic", "spread")
    ref = V.attrs.get("normalize", "mean")
    lbl, default_cmap, signed = _VAR_STYLE.get(var, ("", "viridis", False))
    Z = V[var].values.astype("float64")
    if var == "count":
        Z[Z == 0] = np.nan

    if clim is None:
        if signed:
            lim = _pctl(np.abs(Z), 98) or 1.0
            clim = (-lim, lim)
        else:
            clim = (0, _pctl(Z, 98))

    x, y = V.x.values, V.y.values
    extent = [x.min(), x.max(), y.min(), y.max()]
    origin = "upper" if y[0] > y[-1] else "lower"
    im = ax.imshow(Z, extent=extent, origin=origin, vmin=clim[0], vmax=clim[1],
                   cmap=cmap or default_cmap, aspect="equal",
                   interpolation="nearest")
    cb = ax.figure.colorbar(im, ax=ax)
    cb.set_label(lbl.format(stat=stat, ref=ref))
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if roi is not None:
        ax.plot(roi.x, roi.y, "c-", lw=1.5)

    mean_key = f"{var}_mean" if f"{var}_mean" in V.attrs else None
    summary = f"  mean={V.attrs[mean_key]:.3g}" if mean_key else ""
    ax.set_title(title or
                 f"{V.attrs.get('roi', '') or V.attrs.get('band', '')} {var}"
                 f" ({stat}/{ref})\n"
                 f"t=[{V.attrs['time_min']:.2f}..{V.attrs['time_max']:.2f}]"
                 f"  n={V.attrs.get('nframes', 0)}{summary}", fontsize=10)
    return ax


def plot_compare(C, ax=None, use_datetime: bool = True, title: str = ""):
    """Plot a paired model/obs series from :func:`moc_py.compare.compare_timeseries`.

    Args:
        C: A ``TimeSeriesComparison``.
        ax: Existing axes.
        use_datetime: Plot on a datetime x-axis (else decimal year).
        title: Axis title.

    Returns:
        The matplotlib ``Axes``.
    """
    from .util import decyear_to_datetime64

    ax = _new_ax(ax)
    mt = decyear_to_datetime64(C.model_time) if use_datetime else C.model_time
    ot = decyear_to_datetime64(C.obs_time) if use_datetime else C.obs_time
    if C.obs_err is not None and np.any(np.isfinite(C.obs_err)):
        ax.errorbar(ot, C.obs_val, yerr=C.obs_err, fmt="o", ms=3, color="0.2",
                    ecolor="0.7", capsize=0, label="observed")
    else:
        ax.plot(ot, C.obs_val, "o", ms=3, color="0.2", label="observed")
    ax.plot(mt, C.model_val, "-", color="#d81b1b", lw=1.8, label="model")
    ax.set_xlabel("date" if use_datetime else "decimal year")
    ax.set_ylabel("speed [m/yr]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    s = C.stats
    ax.set_title(title or f"({C.xy[0]:.0f}, {C.xy[1]:.0f})  "
                 f"RMSE={s['rmse']:.0f}  bias={s['mean']:.0f}  r={s.get('r', float('nan')):.2f}")
    return ax


def plot_flowline(Fobs, mode: str = "profile", time: float | None = None,
                  ax=None, model=None, dlim=None):
    """Plot velocity along a flowline: a profile at one time or a distance-time image.

    Args:
        Fobs: Flowline Dataset from :func:`load_obs_flowline`.
        mode: ``"profile"`` (speed vs distance at one time) or ``"hovmoller"``
            (distance-time image of observed speed, or model-obs if ``model`` given).
        time: Decimal year for profile mode (default: median obs time).
        ax: Existing axes.
        model: Optional modelled flowline Dataset with matching ``d``/``time``/``vel``
            (Phase 2). If given in hovmoller mode, plots model - obs.
        dlim: Symmetric colour limit for a difference image.

    Returns:
        The matplotlib ``Axes``.
    """
    ax = _new_ax(ax)
    d_km = Fobs.d.values / 1e3

    if mode == "profile":
        tt = time if time is not None else float(np.median(Fobs.time.values))
        oi = int(np.argmin(np.abs(Fobs.time.values - tt)))
        ax.plot(d_km, Fobs.vel.isel(time=oi).values, "o-", ms=3, color="0.3",
                label=f"obs t={float(Fobs.time.isel(time=oi)):.2f}")
        if model is not None:
            mi = int(np.argmin(np.abs(model.time.values - tt)))
            ax.plot(model.d.values / 1e3, model.vel.isel(time=mi).values, "-",
                    color="#d81b1b", lw=1.8, label=f"model t={float(model.time.isel(time=mi)):.2f}")
        ax.set_xlabel("along-flowline distance [km]")
        ax.set_ylabel("speed [m/yr]")
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_title(f"{Fobs.attrs.get('name','')} profile near t={tt:.2f}")

    elif mode == "hovmoller":
        if model is not None:
            Z = model.vel.values - Fobs.vel.values  # (point, time)
            label = "model - obs [m/yr]"
            cmap = "RdBu_r"
            lim = dlim if dlim is not None else _pctl(np.abs(Z), 95)
            vmin, vmax = -lim, lim
        else:
            Z = Fobs.vel.values
            label = "speed [m/yr]"
            cmap = "viridis"
            vmin, vmax = 0, _pctl(Z, 98)
        t = Fobs.time.values
        im = ax.pcolormesh(t, d_km, Z, shading="nearest", cmap=cmap,
                           vmin=vmin, vmax=vmax)
        cb = ax.figure.colorbar(im, ax=ax)
        cb.set_label(label)
        ax.set_xlabel("decimal year")
        ax.set_ylabel("along-flowline distance [km]")
        ax.set_title(f"{Fobs.attrs.get('name','')} flowline")
    else:
        raise ValueError(f"unknown mode: {mode}")
    return ax
