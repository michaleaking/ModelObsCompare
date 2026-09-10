"""Per-pixel temporal variability of observed velocity over a time window.

Answers "where, inside this ROI, did the ice speed change the most *relative to
how fast it normally goes*". The absolute spread (e.g. peak-to-peak speed over a
window) is dominated by the fast trunk of a glacier simply because it is fast;
dividing each pixel's spread by its own reference speed gives a dimensionless
map in which a pixel can be compared with its neighbours as a *fraction of
relative change* rather than an absolute magnitude change.

Works on any Dataset with dims ``(y, x, time)`` — the observed grids from
:func:`moc_py.load_obs_netcdf` in particular. There is no MATLAB counterpart.
"""

from __future__ import annotations

import warnings

import numpy as np
import xarray as xr

from .roi import Roi, roi_mask

#: Spread statistics accepted by :func:`velocity_variance`.
STATISTICS = ("range", "std", "iqr", "mad")

#: Per-pixel normalisers accepted by :func:`velocity_variance`.
NORMALIZERS = ("mean", "median", "first", "max", "none")


def velocity_variance(
    ds: xr.Dataset,
    roi: Roi | None = None,
    band: str = "vv",
    trange: tuple[float, float] | None = None,
    times: tuple[float, float] | list[float] | None = None,
    statistic: str = "range",
    normalize: str = "mean",
    min_count: int = 2,
    min_speed: float = 0.0,
) -> xr.Dataset:
    """Map the normalised temporal spread of speed over a window, per pixel.

    For every pixel the frames inside the temporal window are reduced to a
    spread (``statistic``) and a reference speed (``normalize``); the returned
    ``norm`` field is ``spread / reference``, i.e. the change expressed as a
    fraction of that pixel's own speed. Pixels are therefore comparable with
    their neighbours regardless of how fast the ice is there.

    The window is chosen in one of three ways:
      * ``times=(t0, t1)`` — the two frames nearest those decimal years (the
        "two time slices" case; a signed ``change`` field is added),
      * ``trange=(t0, t1)`` — every frame inside that decimal-year window,
      * neither — every frame in ``ds``.

    Args:
        ds: Gridded Dataset with dims ``(y, x, time)``, e.g. from
            :func:`moc_py.load_obs_netcdf`.
        roi: Region of interest to restrict to. ``None`` -> the full grid.
        band: Which data var to use (default ``"vv"`` speed).
        trange: ``(t0, t1)`` decimal-year window, inclusive.
        times: Target decimal years; the nearest frame to each is used. Pass two
            values to compare two time slices.
        statistic: Spread measure over the window — ``"range"`` (max - min),
            ``"std"``, ``"iqr"`` (75th - 25th percentile) or ``"mad"`` (median
            absolute deviation from the median).
        normalize: Per-pixel reference speed — ``"mean"``, ``"median"``,
            ``"first"`` (earliest frame), ``"max"``, or ``"none"`` (then
            ``norm`` is the raw spread).
        min_count: Minimum number of finite frames a pixel needs; pixels with
            fewer are NaN everywhere in the output.
        min_speed: Reference speeds below this [m/yr] are treated as invalid
            (0 = keep all). Useful to stop near-stagnant ice from producing huge
            ratios that swamp the colour scale.

    Returns:
        Dataset over the ROI with dims ``(y, x)``:
            ``spread`` — absolute spread [m/yr],
            ``reference`` — per-pixel normaliser [m/yr],
            ``norm`` — ``spread / reference`` (dimensionless fraction),
            ``count`` — finite frames per pixel,
            ``change`` / ``norm_change`` — signed later-minus-earlier difference,
            present only when exactly two frames are used.
        Coords ``x``, ``y``, and ``time`` (the frames used). Attrs carry the
        ROI-mean/median/p90 and finite-pixel count of ``norm`` and ``spread``
        (``norm_mean``, ``norm_median``, ``norm_p90``, ``norm_npix``, and the
        same for ``spread``), ``npix`` (= ``norm_npix``, the pixels the
        normalised map reports on — fewer than ``spread_npix`` when
        ``min_speed`` bites), ``time_min``/``time_max``, ``statistic``,
        ``normalize``, ``band``, ``epsg`` and ``roi``.

    Raises:
        ValueError: On an unknown ``statistic``/``normalize``, a missing band,
            or a window that selects no frames.

    Example:
        >>> V = velocity_variance(O, roi, trange=(2019.0, 2019.5))
        >>> V.attrs["norm_mean"]          # mean fractional change over the ROI
        >>> plot_variance(V, roi=roi)     # spatial map of it
    """
    if statistic not in STATISTICS:
        raise ValueError(f"statistic must be one of {STATISTICS}, got {statistic!r}")
    if normalize not in NORMALIZERS:
        raise ValueError(f"normalize must be one of {NORMALIZERS}, got {normalize!r}")
    if band not in ds:
        raise ValueError(f"band {band!r} not in Dataset (has {list(ds.data_vars)})")

    sub = _roi_subset(ds, roi)
    tidx = _time_indices(sub.time.values, trange, times)
    da = sub[band].isel(time=tidx).transpose("y", "x", "time")

    A = np.asarray(da.values, dtype="float64")  # (ny, nx, nt)
    count = np.isfinite(A).sum(axis=2)
    keep = count >= max(int(min_count), 1)

    spread = _spread(A, statistic)
    reference = _reference(A, normalize)

    bad_ref = ~np.isfinite(reference) | (reference <= 0) | (reference < min_speed)
    if normalize == "none":
        norm = spread.copy()
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            norm = spread / reference
        norm[bad_ref] = np.nan

    if roi is not None:
        inside = roi_mask(roi, *np.meshgrid(sub.x.values, sub.y.values))
        keep &= inside

    for arr in (spread, reference, norm):
        arr[~keep] = np.nan

    data_vars = {
        "spread": (("y", "x"), spread),
        "reference": (("y", "x"), reference),
        "norm": (("y", "x"), norm),
        "count": (("y", "x"), count.astype("int32")),
    }

    t_used = np.asarray(da.time.values, dtype="float64")
    if t_used.size == 2:
        change = A[:, :, 1] - A[:, :, 0]  # later minus earlier
        change[~keep] = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            norm_change = change / reference
        norm_change[bad_ref | ~keep] = np.nan
        data_vars["change"] = (("y", "x"), change)
        data_vars["norm_change"] = (("y", "x"), norm_change)

    V = xr.Dataset(
        data_vars=data_vars,
        coords={
            "x": ("x", sub.x.values),
            "y": ("y", sub.y.values),
            "time": ("time", t_used),
        },
        attrs={
            "band": band,
            "statistic": statistic,
            "normalize": normalize,
            "min_count": int(min_count),
            "min_speed": float(min_speed),
            "nframes": int(t_used.size),
            "time_min": float(np.min(t_used)),
            "time_max": float(np.max(t_used)),
            "roi": roi.name if roi is not None else "",
            "epsg": int(ds.attrs.get("epsg", 3413)),
        },
    )
    V.attrs.update(_summary(norm, "norm"))
    V.attrs.update(_summary(spread, "spread"))
    V.attrs["npix"] = V.attrs["norm_npix"]   # pixels the normalised map reports on
    return V


# ---------------------------------------------------------------------------
def _roi_subset(ds: xr.Dataset, roi: Roi | None) -> xr.Dataset:
    """Crop a Dataset to an ROI's bounding box (no-op if ``roi`` is None).

    Args:
        ds: Dataset with ``x``/``y`` coords (either axis may be descending).
        roi: The region of interest, or None.

    Returns:
        The cropped Dataset.

    Raises:
        ValueError: If the ROI does not overlap the grid.
    """
    if roi is None:
        return ds
    xmin, xmax, ymin, ymax = roi.bbox
    x = ds.x.values
    y = ds.y.values
    xi = np.where((x >= xmin) & (x <= xmax))[0]
    yi = np.where((y >= ymin) & (y <= ymax))[0]
    if xi.size == 0 or yi.size == 0:
        raise ValueError(
            f"ROI bbox ({xmin:.0f}..{xmax:.0f}, {ymin:.0f}..{ymax:.0f}) "
            "does not overlap the grid"
        )
    return ds.isel(x=xi, y=yi)


def _time_indices(t: np.ndarray, trange, times) -> np.ndarray:
    """Resolve the time-window arguments to frame indices (sorted, unique).

    Args:
        t: Decimal-year coordinate of the Dataset.
        trange: ``(t0, t1)`` inclusive window, or None.
        times: Target decimal years matched to the nearest frame, or None.

    Returns:
        Integer index array into the time axis.

    Raises:
        ValueError: If both are given, or the selection is empty/too small.
    """
    if times is not None and trange is not None:
        raise ValueError("Pass either 'times' or 'trange', not both.")

    if times is not None:
        targets = np.atleast_1d(np.asarray(times, dtype="float64"))
        idx = np.unique([int(np.argmin(np.abs(t - tt))) for tt in targets])
        if idx.size < targets.size:
            warnings.warn(
                f"{targets.size} target times collapsed onto {idx.size} distinct "
                "frame(s); the requested slices are not resolved by this series.",
                stacklevel=3,
            )
    elif trange is not None:
        idx = np.where((t >= trange[0]) & (t <= trange[1]))[0]
        if idx.size == 0:
            raise ValueError(
                f"No frames in trange {trange}; series spans "
                f"{t.min():.3f}..{t.max():.3f}"
            )
    else:
        idx = np.arange(t.size)

    if idx.size < 2:
        raise ValueError("Need at least 2 time slices to measure variability.")
    return idx[np.argsort(t[idx])]  # time-ascending, so 'first'/'change' are ordered


def _spread(A: np.ndarray, statistic: str) -> np.ndarray:
    """Reduce a ``(y, x, time)`` stack to a per-pixel spread, ignoring NaNs.

    Args:
        A: Speed stack [m/yr].
        statistic: One of :data:`STATISTICS`.

    Returns:
        ``(y, x)`` spread [m/yr]; all-NaN pixels come back NaN.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN pixels
        if statistic == "range":
            return np.nanmax(A, axis=2) - np.nanmin(A, axis=2)
        if statistic == "std":
            return np.nanstd(A, axis=2)
        if statistic == "iqr":
            q75, q25 = np.nanpercentile(A, [75, 25], axis=2)
            return q75 - q25
        med = np.nanmedian(A, axis=2)
        return np.nanmedian(np.abs(A - med[:, :, None]), axis=2)


def _reference(A: np.ndarray, normalize: str) -> np.ndarray:
    """Per-pixel reference speed used to normalise the spread.

    Args:
        A: ``(y, x, time)`` speed stack [m/yr], time ascending.
        normalize: One of :data:`NORMALIZERS`.

    Returns:
        ``(y, x)`` reference speed [m/yr]. For ``"none"`` this is the mean,
        reported for context but not divided by.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if normalize == "median":
            return np.nanmedian(A, axis=2)
        if normalize == "max":
            return np.nanmax(A, axis=2)
        if normalize == "first":
            # Earliest finite frame per pixel (frames may be individually NaN).
            first = np.argmax(np.isfinite(A), axis=2)
            ref = np.take_along_axis(A, first[:, :, None], axis=2)[:, :, 0]
            return ref
        return np.nanmean(A, axis=2)


def _summary(field: np.ndarray, prefix: str) -> dict:
    """Mean/median/percentile summary of a 2-D field over its finite pixels.

    Args:
        field: The ``(y, x)`` field to summarise.
        prefix: Attribute-name prefix, e.g. ``"norm"``.

    Returns:
        Dict of ``{prefix}_mean``, ``_median``, ``_p90``, ``_npix``. The keys are
        prefixed so summarising a second field cannot overwrite the first.
    """
    v = field[np.isfinite(field)]
    if v.size == 0:
        return {f"{prefix}_mean": np.nan, f"{prefix}_median": np.nan,
                f"{prefix}_p90": np.nan, f"{prefix}_npix": 0}
    return {
        f"{prefix}_mean": float(np.mean(v)),
        f"{prefix}_median": float(np.median(v)),
        f"{prefix}_p90": float(np.percentile(v, 90)),
        f"{prefix}_npix": int(v.size),
    }
