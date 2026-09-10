"""Build flowline inputs from shapefiles: crop to an ROI, clip at the terminus.

Two archived shapefile collections are combined here:

* **Flowlines** — Felikson et al. (2020), one shapefile per glacier
  (``glacier0001.shp`` ... ``glacierd216.shp``), each holding several
  ``LineString`` flowlines identified by a ``flowline`` attribute. Vertices run
  from the terminus inland, so index 0 is the downstream end.
* **Termini** — Black & Joughin traced terminus positions
  (``glacier_termini_v01.0.shp``), one ``LineString`` per glacier per image
  with a ``SourceDate`` and a ``Quality_Fl``.

The pipeline is :func:`flowline_from_roi`:

    ROI -> flowlines intersecting it -> pick the centre one -> split it at the
    most recent overlapping terminus trace -> keep the upstream part -> clip to
    the ROI -> resample -> :class:`xarray.Dataset`

The Dataset uses the same ``(point,)`` layout as :func:`moc_py.load_obs_flowline`
(coords ``x``, ``y``, ``d``), so it drops straight into
:func:`moc_py.model_flowline` and :func:`moc_py.compare_timeseries`. Pair it
with a gridded observation subset through :func:`sample_obs_along_flowline` to
get the full ``(point, time)`` observed flowline the plotting functions expect.

Requires the optional ``shapefile`` extra (``geopandas``, ``shapely``).
"""

from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import xarray as xr

from .config import CONFIG

#: Filename pattern for a Felikson per-glacier flowline shapefile.
_GLACIER_RE = re.compile(r"^glacier(?P<glacier>[a-z]?\d{3,4})(?:_iter(?P<iter>\d+))?$")


# ---------------------------------------------------------------------------
# ROI handling
# ---------------------------------------------------------------------------
def roi_polygon(roi, epsg: int = CONFIG.epsg):
    """Coerce anything ROI-shaped into a shapely polygon in ``epsg``.

    Args:
        roi: A :class:`~moc_py.roi.Roi`, a shapely polygon, a
            ``(xmin, xmax, ymin, ymax)`` bbox as returned by ``Roi.bbox``, a
            ``GeoDataFrame``/``GeoSeries`` (its unary union is used), or an ROI
            name / path resolvable by :func:`moc_py.load_roi`.
        epsg: Projection the returned polygon must be in (default
            ``CONFIG.epsg``, 3413).

    Returns:
        A shapely ``Polygon``/``MultiPolygon`` in ``epsg``.

    Raises:
        TypeError: If ``roi`` is not one of the supported forms.
    """
    import geopandas as gpd
    from shapely.geometry import Polygon, box
    from shapely.geometry.base import BaseGeometry

    src_epsg = epsg

    if isinstance(roi, BaseGeometry):
        poly = roi
    elif isinstance(roi, (gpd.GeoDataFrame, gpd.GeoSeries)):
        poly = roi.union_all() if hasattr(roi, "union_all") else roi.unary_union
        if roi.crs is not None:
            src_epsg = int(roi.crs.to_epsg())
    elif hasattr(roi, "x") and hasattr(roi, "y"):  # moc_py.Roi
        poly = Polygon(np.column_stack([np.asarray(roi.x), np.asarray(roi.y)]))
        src_epsg = int(getattr(roi, "epsg", epsg))
    elif isinstance(roi, (str, Path)):
        from .roi import load_roi

        return roi_polygon(load_roi(roi), epsg=epsg)
    elif np.ndim(roi) == 1 and len(roi) == 4:
        xmin, xmax, ymin, ymax = map(float, roi)
        poly = box(xmin, ymin, xmax, ymax)
    else:
        raise TypeError(f"Cannot interpret {type(roi)!r} as an ROI")

    if src_epsg != epsg:
        poly = (gpd.GeoSeries([poly], crs=f"EPSG:{src_epsg}")
                .to_crs(epsg=epsg).iloc[0])
    return poly


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def load_flowlines(roi=None, directory: str | Path | None = None,
                   pattern: str = "glacier*.shp", iterations: bool = False,
                   clip: bool = False, predicate: str = "intersects"):
    """Read every per-glacier flowline shapefile and keep those inside an ROI.

    Loops over the whole Felikson collection (the ``glob`` in the *Flowlines*
    notebook), concatenating the files into one table rather than a dict of
    per-glacier objects.

    Args:
        roi: Region of interest, in any form accepted by :func:`roi_polygon`.
            ``None`` returns the full collection.
        directory: Folder of flowline shapefiles (default
            ``CONFIG.flowline_shp_dir``).
        pattern: Glob for the shapefiles to read.
        iterations: Keep the ``*_iterNN.shp`` intermediate versions too
            (default ``False``: only the final flowline per glacier).
        clip: Clip geometries to the ROI boundary instead of returning whole
            flowlines that merely reach into it.
        predicate: Spatial test against the ROI — ``"intersects"`` (any
            overlap, default) or ``"within"`` (entirely inside).

    Returns:
        ``GeoDataFrame`` with columns ``glacier``, ``flowline``, ``iteration``,
        ``file`` and ``geometry``, in EPSG:3413.

    Raises:
        FileNotFoundError: If the directory holds no matching shapefiles.
    """
    import geopandas as gpd
    import pandas as pd

    directory = Path(directory) if directory else CONFIG.flowline_shp_dir
    files = sorted(glob.glob(str(directory / pattern)))
    if not files:
        raise FileNotFoundError(f"No shapefiles matching {pattern} in {directory}")

    poly = roi_polygon(roi) if roi is not None else None
    frames = []
    for f in files:
        m = _GLACIER_RE.match(Path(f).stem)
        if m is None:  # e.g. sverdrup_masks.shp — not a flowline file
            continue
        if m.group("iter") and not iterations:
            continue
        g = gpd.read_file(f)
        if "flowline" not in g.columns:
            continue
        if poly is not None:
            g = g[getattr(g.geometry, predicate)(poly)]
            if g.empty:
                continue
        g = g.assign(glacier=m.group("glacier"),
                     iteration=int(m.group("iter") or 0),
                     file=f)
        frames.append(g)

    if not frames:
        cols = ["glacier", "flowline", "iteration", "file", "geometry"]
        return gpd.GeoDataFrame(pd.DataFrame(columns=cols), geometry="geometry",
                                crs=f"EPSG:{CONFIG.epsg}")

    out = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True),
                           crs=frames[0].crs)
    out = out.to_crs(epsg=CONFIG.epsg)
    if clip and poly is not None:
        out = gpd.clip(out, poly)
        out = out[~out.geometry.is_empty]
    cols = ["glacier", "flowline", "iteration", "file", "geometry"]
    return out[cols + [c for c in out.columns if c not in cols]].reset_index(drop=True)


def load_termini(roi=None, path: str | Path | None = None,
                 max_quality: int | None = 0, trange=None):
    """Read traced terminus positions and keep those inside an ROI.

    Args:
        roi: Region of interest, in any form accepted by :func:`roi_polygon`.
            ``None`` returns every trace in the file.
        path: Terminus shapefile (default ``CONFIG.termini_shp``).
        max_quality: Keep traces whose ``Quality_Fl`` is at most this. The
            product flags lower-confidence traces with non-zero values, so the
            default ``0`` keeps only the cleanest; pass ``None`` to keep all.
        trange: Optional ``(start, end)`` date bound, anything
            :func:`pandas.to_datetime` accepts; either end may be ``None``.

    Returns:
        ``GeoDataFrame`` of terminus ``LineString``s in EPSG:3413, with an added
        ``date`` column (``datetime64``), sorted oldest to newest.
    """
    import geopandas as gpd
    import pandas as pd

    path = Path(path) if path else CONFIG.termini_shp
    t = gpd.read_file(path).to_crs(epsg=CONFIG.epsg)
    t["date"] = pd.to_datetime(t["SourceDate"], errors="coerce")

    if max_quality is not None and "Quality_Fl" in t.columns:
        t = t[t["Quality_Fl"] <= max_quality]
    if trange is not None:
        lo, hi = trange
        if lo is not None:
            t = t[t["date"] >= pd.to_datetime(lo)]
        if hi is not None:
            t = t[t["date"] <= pd.to_datetime(hi)]
    if roi is not None:
        t = t[t.intersects(roi_polygon(roi))]

    return (t.dropna(subset=["date"]).sort_values("date")
            .reset_index(drop=True))


# ---------------------------------------------------------------------------
# Geometry operations
# ---------------------------------------------------------------------------
def ordered_flowlines(flowlines, glacier: str | None = None):
    """Order a glacier's fan of flowlines from the centre outwards.

    Felikson's flowlines are seeded side by side across the trunk, so ranking
    them by how far their downstream end sits from the mean of all the
    downstream ends walks the fan centre first, then alternately out to either
    side. That order is what :func:`flowline_from_roi` follows when the centre
    line has no terminus trace crossing it.

    Args:
        flowlines: ``GeoDataFrame`` from :func:`load_flowlines`.
        glacier: Restrict to this glacier id (e.g. ``"a045"``).

    Returns:
        The rows, centre first (a ``GeoDataFrame``).

    Raises:
        ValueError: If ``flowlines`` is empty, or ``glacier`` matches nothing.
    """
    g = flowlines if glacier is None else flowlines[flowlines["glacier"] == glacier]
    if len(g) == 0:
        raise ValueError(f"No flowlines to choose from (glacier={glacier!r})")
    if len(g) == 1:
        return g

    ends = np.array([line.coords[0] for line in g.geometry])
    d = np.hypot(*(ends - ends.mean(axis=0)).T)
    return g.iloc[np.argsort(d, kind="stable")]


def center_flowline(flowlines, glacier: str | None = None):
    """Pick the middle flowline of a glacier's fan of flowlines.

    The "centre" line is the one whose downstream end sits closest to the mean
    of all the downstream ends — the first of :func:`ordered_flowlines`.

    Args:
        flowlines: ``GeoDataFrame`` from :func:`load_flowlines`.
        glacier: Restrict to this glacier id (e.g. ``"a045"``). Required when
            more than one glacier is present, unless you want the centre of the
            pooled set.

    Returns:
        The selected row (a ``GeoSeries``).

    Raises:
        ValueError: If ``flowlines`` is empty, or ``glacier`` matches nothing.
    """
    return ordered_flowlines(flowlines, glacier).iloc[0]


def latest_terminus(termini, line=None, max_distance: float | None = None):
    """Return the most recent terminus trace, optionally one crossing a line.

    Args:
        termini: ``GeoDataFrame`` from :func:`load_termini`.
        line: If given, only traces that intersect this flowline are considered
            (falling back to the whole set if none do, so a flowline that has
            already retreated past every trace still gets an answer).
        max_distance: When falling back, ignore traces further than this [m]
            from ``line``.

    Returns:
        The newest matching row, or ``None`` if there is nothing to return.
    """
    if len(termini) == 0:
        return None

    cand = termini
    if line is not None:
        crossing = termini[termini.intersects(line)]
        if len(crossing):
            cand = crossing
        elif max_distance is not None:
            cand = termini[termini.distance(line) <= max_distance]
            if len(cand) == 0:
                return None
    return cand.loc[cand["date"].idxmax()]


def clip_upstream(line, terminus, downstream_end: str = "auto"):
    """Cut a flowline at a terminus trace and keep the upstream part.

    The terminus acts as a splitting feature: everything seaward of the
    furthest-inland crossing is discarded.

    Args:
        line: Flowline ``LineString``.
        terminus: Terminus ``LineString`` (or any geometry that crosses it).
        downstream_end: Which end of ``line`` is the seaward one — ``"first"``,
            ``"last"``, or ``"auto"`` (default: whichever endpoint lies closer
            to ``terminus``; Felikson flowlines start at the terminus, so this
            resolves to ``"first"``).

    Returns:
        A ``LineString`` running from the terminus inland, or ``None`` if the
        two geometries never cross.
    """
    import shapely
    from shapely import get_coordinates
    from shapely.geometry import Point
    from shapely.ops import substring

    if downstream_end == "auto":
        first, last = Point(line.coords[0]), Point(line.coords[-1])
        downstream_end = ("first" if terminus.distance(first) <= terminus.distance(last)
                          else "last")
    if downstream_end == "last":
        line = type(line)(list(line.coords)[::-1])
    elif downstream_end != "first":
        raise ValueError("downstream_end must be 'first', 'last' or 'auto'")

    crossings = get_coordinates(line.intersection(terminus))
    if crossings.size == 0:
        return None

    # Furthest-inland crossing, so no seaward remnant survives a wiggly terminus.
    # errstate: GEOS raises a spurious invalid-value warning locating a point on
    # a line this long (the flowlines run hundreds of km); the result is exact.
    with np.errstate(invalid="ignore"):
        s_cut = float(np.max(shapely.line_locate_point(
            line, shapely.points(crossings))))
    if s_cut >= line.length:
        return None
    return substring(line, s_cut, line.length)


def clip_to_roi(line, roi, keep: str = "longest"):
    """Clip a line to an ROI, resolving the pieces a concave ROI can produce.

    Args:
        line: A ``LineString``.
        roi: Region of interest, in any form accepted by :func:`roi_polygon`.
        keep: ``"longest"`` (default) keeps the longest connected piece;
            ``"first"`` keeps the piece nearest the line's start.

    Returns:
        The clipped ``LineString``, or ``None`` if nothing falls inside.
    """
    from shapely.ops import linemerge

    clipped = line.intersection(roi_polygon(roi))
    if clipped.is_empty:
        return None
    if clipped.geom_type == "MultiLineString":
        clipped = linemerge(clipped)
    if clipped.geom_type == "LineString":
        return clipped

    parts = [g for g in getattr(clipped, "geoms", []) if g.geom_type == "LineString"]
    if not parts:
        return None
    if keep == "first":
        return min(parts, key=lambda g: line.project(g.interpolate(0.0)))
    return max(parts, key=lambda g: g.length)


# ---------------------------------------------------------------------------
# Packaging for the comparison functions
# ---------------------------------------------------------------------------
def flowline_to_dataset(line, name: str = "flowline", spacing: float | None = None,
                        **attrs) -> xr.Dataset:
    """Turn a ``LineString`` into the ``(point,)`` layout the toolkit expects.

    Args:
        line: The flowline, ordered from its downstream end inland.
        name: Value for the ``name`` attribute.
        spacing: Resample to this even along-line spacing [m]. ``None``
            (default) keeps the shapefile's own vertices.
        **attrs: Extra attributes to record on the Dataset (e.g.
            ``terminus_date``).

    Returns:
        Dataset with dims ``(point,)``, coords ``x``, ``y`` and ``d``
        (along-line distance [m], zero at the downstream end), and attrs
        ``name``, ``epsg``, ``length``. It is accepted anywhere an obs-flowline
        Dataset is — :func:`moc_py.model_flowline`,
        :func:`moc_py.compare_timeseries` — and gains ``vel``/``time`` from
        :func:`sample_obs_along_flowline`.

    Raises:
        ValueError: If ``line`` is empty or has fewer than two vertices.
    """
    if line is None or line.is_empty:
        raise ValueError("Cannot build a Dataset from an empty line")

    if spacing:
        n = max(int(np.floor(line.length / spacing)) + 1, 2)
        d = np.arange(n, dtype="float64") * float(spacing)
        xy = np.array([line.interpolate(float(s)).coords[0] for s in d])
        x, y = xy[:, 0], xy[:, 1]
    else:
        xy = np.asarray(line.coords, dtype="float64")
        x, y = xy[:, 0], xy[:, 1]
        step = np.hypot(np.diff(x), np.diff(y))
        d = np.concatenate([[0.0], np.cumsum(step)])

    if x.size < 2:
        raise ValueError("Flowline has fewer than two points after processing")

    return xr.Dataset(
        coords={"x": ("point", x), "y": ("point", y), "d": ("point", d)},
        attrs={"name": name, "epsg": CONFIG.epsg,
               "length": float(d[-1]), **attrs},
    )


def flowline_from_roi(roi, glacier: str | None = None, flowline: str | None = None,
                      spacing: float | None = None, name: str | None = None,
                      flowline_dir: str | Path | None = None,
                      termini_path: str | Path | None = None,
                      max_quality: int | None = 0, trange=None,
                      require_terminus: bool = True) -> xr.Dataset:
    """Build one ROI-cropped, terminus-clipped flowline ready for comparison.

    Runs the whole pipeline: read every flowline shapefile, keep those reaching
    into the ROI, take the centre flowline of the glacier, split it at the most
    recent terminus trace that crosses it, keep the upstream side, clip that to
    the ROI, and package it as a Dataset.

    If no terminus trace crosses the centre flowline, the next line out is tried
    instead — working outwards from the centre, alternating sides — and the
    first one a trace does cross is used. Compare the ``flowline`` and
    ``center_flowline`` attributes to see whether that happened. Naming a
    ``flowline`` explicitly disables the search: that line is used or the call
    fails.

    Args:
        roi: Region of interest, in any form accepted by :func:`roi_polygon`.
        glacier: Glacier id (e.g. ``"a045"``). Default: the glacier with the
            most flowlines in the ROI.
        flowline: Flowline id (e.g. ``"05"``). Default: the centre one, per
            :func:`center_flowline`, falling outwards to its neighbours if no
            terminus trace crosses it.
        spacing: Even resampling interval [m]; ``None`` keeps native vertices.
        name: Dataset ``name`` attribute (default ``"<glacier>_<flowline>"``).
        flowline_dir: Folder of flowline shapefiles (default
            ``CONFIG.flowline_shp_dir``).
        termini_path: Terminus shapefile (default ``CONFIG.termini_shp``).
        max_quality: Terminus quality cut-off, see :func:`load_termini`.
        trange: Optional ``(start, end)`` date bound on the terminus traces —
            use it to clip at the front position of a particular epoch.
        require_terminus: Raise if no terminus trace crosses any candidate
            flowline. Set ``False`` to fall back to the un-split centre
            flowline (the returned ``terminus_date`` attribute is then empty).

    Returns:
        Dataset from :func:`flowline_to_dataset`, with extra attrs ``glacier``,
        ``flowline`` (the line actually used), ``center_flowline`` (the centre
        of the fan, which differs when the search fell outwards),
        ``terminus_date``, ``terminus_file`` and ``file``.

    Raises:
        ValueError: If the ROI contains no flowlines, no terminus trace crosses
            any candidate, or the chosen flowline is clipped away entirely.
    """
    F = load_flowlines(roi, directory=flowline_dir)
    if len(F) == 0:
        raise ValueError("No flowlines intersect this ROI")

    if glacier is None:
        glacier = F["glacier"].value_counts().idxmax()
    sub = F[F["glacier"] == glacier]
    if len(sub) == 0:
        raise ValueError(f"Glacier {glacier!r} has no flowlines in this ROI "
                         f"(available: {sorted(F['glacier'].unique())})")

    ranked = ordered_flowlines(sub)          # centre first, then outwards
    center_id = str(ranked.iloc[0]["flowline"])
    if flowline is None:
        candidates = ranked
    else:
        match = sub[sub["flowline"] == flowline]
        if len(match) == 0:
            raise ValueError(f"Flowline {flowline!r} not in glacier {glacier!r} "
                             f"(available: {sorted(sub['flowline'])})")
        candidates = match               # named explicitly: no outward search

    T = load_termini(roi, path=termini_path, max_quality=max_quality, trange=trange)

    # Walk out from the centre until a terminus trace actually crosses a line.
    row, line, tdate, tfile = None, None, "", ""
    for _, cand in candidates.iterrows():
        term = latest_terminus(T, line=cand.geometry)
        if term is None:
            continue
        upstream = clip_upstream(cand.geometry, term.geometry)
        if upstream is None:
            continue
        row, line = cand, upstream
        tdate = str(np.datetime64(term["date"], "D"))
        tfile = str(termini_path or CONFIG.termini_shp)
        break

    if row is None:
        if require_terminus:
            tried = ", ".join(str(f) for f in candidates["flowline"])
            raise ValueError(
                f"No terminus trace in the ROI crosses any flowline of glacier "
                f"{glacier!r} (tried {tried}); pass require_terminus=False to "
                f"skip the clip"
            )
        row = candidates.iloc[0]         # un-split, from the centre outwards
        line = row.geometry

    line = clip_to_roi(line, roi)
    if line is None or line.length == 0:
        raise ValueError("Nothing of the flowline is left inside the ROI once "
                         "the part downstream of the terminus is removed")

    return flowline_to_dataset(
        line, name=name or f"{glacier}_{row['flowline']}", spacing=spacing,
        glacier=str(glacier), flowline=str(row["flowline"]),
        center_flowline=center_id,
        terminus_date=tdate, terminus_file=tfile, file=str(row["file"]),
    )


def sample_obs_along_flowline(O: xr.Dataset, F: xr.Dataset, band: str = "vv",
                              method: str = "linear") -> xr.Dataset:
    """Sample a gridded observation Dataset along a flowline.

    Turns the geometry-only Dataset from :func:`flowline_from_roi` into the
    ``(point, time)`` form :func:`moc_py.load_obs_flowline` produces, so it can
    be handed to :func:`moc_py.plot_flowline` and :func:`moc_py.compare_timeseries`.

    Args:
        O: Gridded obs Dataset from :func:`moc_py.load_obs_netcdf` (dims
            ``y, x, time``).
        F: Flowline Dataset with ``x``/``y``/``d`` coords.
        band: Band to sample as ``vel`` (default ``"vv"``). ``vx``/``vy`` are
            carried across too when present in ``O``.
        method: Interpolation passed to :meth:`xarray.Dataset.interp`
            (``"linear"`` or ``"nearest"``).

    Returns:
        Dataset with ``vel`` (and where available ``vx``, ``vy``) on dims
        ``(point, time)``, coords ``x``, ``y``, ``d``, ``time``, and the
        flowline's attrs plus ``source``.

    Raises:
        KeyError: If ``band`` is not a variable of ``O``.
    """
    if band not in O:
        raise KeyError(f"Band {band!r} not in the obs Dataset ({list(O.data_vars)})")

    px = xr.DataArray(F.x.values, dims="point")
    py = xr.DataArray(F.y.values, dims="point")
    wanted = [band] + [b for b in ("vx", "vy", "ex", "ey") if b in O and b != band]
    S = O[wanted].interp(x=px, y=py, method=method)

    ds = S.rename({band: "vel"}).transpose("point", "time")
    ds = ds.assign_coords(x=("point", F.x.values), y=("point", F.y.values),
                          d=("point", F.d.values))
    ds.attrs = {**F.attrs, "source": band,
                "obs_file": O.attrs.get("file", ""), "epsg": CONFIG.epsg}
    return ds
