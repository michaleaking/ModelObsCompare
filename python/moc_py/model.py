"""Load an exported ISSM model netCDF and interpolate it (Phase 2).

Reads the neutral netCDF written by ``moc_export_model.m`` (schema
``moc_py-model-1.0``) with :mod:`h5py`, builds a
:class:`matplotlib.tri.Triangulation`, and samples the model velocity onto
regular grids, points, or flowlines using ``LinearTriInterpolator`` — the same
P1 linear interpolation as ISSM's ``InterpFromMesh2d`` /
``InterpFromMeshToGrid``.

Axes are matched by size, so the reader is robust to the MATLAB/HDF5 transpose
quirk. Element indices are converted to 0-based using the file's ``index_base``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import xarray as xr

from .config import CONFIG, resolve


@dataclass
class Model:
    """A loaded model mesh + velocity time series.

    Attributes:
        x, y: Vertex coordinates ``(nv,)`` [m].
        elements: Triangulation ``(ne, 3)``, **0-based**.
        time: ``(nt,)`` decimal years (may be ``[nan]`` for a snapshot).
        fields: Dict of ``name -> (nv, nt)`` arrays (e.g. ``vel``, ``vx``, ``vy``).
        bed: Optional ``(nv,)`` bed elevation [m].
        epsg: Projection code.
        file: Source netCDF path.
    """

    x: np.ndarray
    y: np.ndarray
    elements: np.ndarray
    time: np.ndarray
    fields: dict
    bed: np.ndarray | None = None
    epsg: int = CONFIG.epsg
    file: str = ""
    _tri: object = field(default=None, repr=False)

    @property
    def triangulation(self):
        """A cached :class:`matplotlib.tri.Triangulation` for the mesh."""
        if self._tri is None:
            from matplotlib.tri import Triangulation
            self._tri = Triangulation(self.x, self.y, self.elements)
        return self._tri

    def field(self, name: str) -> np.ndarray:
        """Return the ``(nv, nt)`` array for a field, raising if absent."""
        if name not in self.fields:
            raise KeyError(f'Model has no field "{name}"; have {list(self.fields)}')
        return self.fields[name]

    def time_index(self, time: float | None, tindex: int | None) -> int:
        """Resolve a time selector to a column index (nearest ``time`` wins)."""
        nt = self.time.size
        if time is not None and not np.all(np.isnan(self.time)):
            return int(np.argmin(np.abs(self.time - time)))
        if tindex is not None:
            return int(np.clip(tindex, 0, nt - 1))
        return nt - 1  # default: last (most recent) step


def load_model(ncpath: str | Path) -> Model:
    """Load an exported model netCDF into a :class:`Model`.

    Args:
        ncpath: Full path or bare filename (resolved against ``CONFIG.model_dir``).

    Returns:
        A :class:`Model`.
    """
    path = resolve(ncpath, "model_dir")
    with h5py.File(path, "r") as h:
        x = np.asarray(h["x"][:], dtype="float64").ravel()
        y = np.asarray(h["y"][:], dtype="float64").ravel()
        nv = x.size
        time = np.asarray(h["time"][:], dtype="float64").ravel()
        nt = time.size

        elements = np.asarray(h["elements"][:])
        elements = _orient_elements(elements)
        base = _attr(h["elements"], "index_base", 1)
        elements = (elements - int(base)).astype("int64")

        skip = {"x", "y", "elements", "time", "surface"}
        fields = {}
        for name in h.keys():
            if name in skip:
                continue
            arr = np.asarray(h[name][:], dtype="float64")
            if arr.ndim == 2:
                fields[name] = _orient_field(arr, nv, nt)
        if "surface" in h:
            fields["surface"] = _orient_field(np.asarray(h["surface"][:], "float64"), nv, nt)

        bed = None
        if "bed" in h:
            bed = np.asarray(h["bed"][:], dtype="float64").ravel()
        epsg = int(_attr(h, "epsg", CONFIG.epsg, is_group=True))

    return Model(x=x, y=y, elements=elements, time=time, fields=fields,
                 bed=bed, epsg=epsg, file=str(path))


def model_to_grid(M: Model, roi, res: float | None = None, field: str = "vel",
                  time: float | None = None, tindex: int | None = None,
                  mask: bool = True) -> xr.Dataset:
    """Interpolate a model field onto a regular grid over an ROI.

    Args:
        M: The model.
        roi: :class:`~moc_py.roi.Roi` defining the bbox (and polygon mask).
        res: Grid spacing [m] (default ``CONFIG.grid_res``).
        field: Field name (``"vel"`` default).
        time: Target decimal year (nearest step). Overrides ``tindex``.
        tindex: Column index (default: last).
        mask: NaN-out cells outside the ROI polygon (default True).

    Returns:
        Dataset with data var ``z`` (dims ``y, x``), coords ``x``/``y``, and
        attrs ``field``, ``time``, ``res``.
    """
    from matplotlib.tri import LinearTriInterpolator
    from .roi import roi_mask

    res = res if res is not None else CONFIG.grid_res
    ti = M.time_index(time, tindex)
    z = M.field(field)[:, ti]

    xmin, xmax, ymin, ymax = roi.bbox
    xg = np.arange(xmin, xmax + res, res)
    yg = np.arange(ymin, ymax + res, res)
    X, Y = np.meshgrid(xg, yg)

    interp = LinearTriInterpolator(M.triangulation, z)
    grid = np.asarray(interp(X, Y).filled(np.nan))

    if mask:
        grid[~roi_mask(roi, X, Y)] = np.nan

    tval = float(M.time[ti]) if not np.all(np.isnan(M.time)) else np.nan
    return xr.Dataset(
        {"z": (("y", "x"), grid)},
        coords={"x": ("x", xg), "y": ("y", yg)},
        attrs={"field": field, "time": tval, "res": float(res), "epsg": M.epsg},
    )


def model_at_points(M: Model, xq, yq, field: str = "vel") -> np.ndarray:
    """Interpolate a model field at points for every time step.

    Args:
        M: The model.
        xq, yq: Query coordinates ``(np,)`` [m].
        field: Field name.

    Returns:
        ``(np, nt)`` array; NaN where a point is outside the mesh.
    """
    from matplotlib.tri import LinearTriInterpolator

    xq = np.atleast_1d(np.asarray(xq, dtype="float64"))
    yq = np.atleast_1d(np.asarray(yq, dtype="float64"))
    Z = M.field(field)
    out = np.full((xq.size, M.time.size), np.nan)
    for k in range(M.time.size):
        interp = LinearTriInterpolator(M.triangulation, Z[:, k])
        out[:, k] = np.asarray(interp(xq, yq).filled(np.nan))
    return out


def model_flowline(M: Model, flowline, field: str = "vel") -> xr.Dataset:
    """Sample a model field along a flowline for every time step.

    Args:
        M: The model.
        flowline: An obs-flowline Dataset (uses ``x``, ``y``, ``d``), a
            :class:`~moc_py.roi.Roi`, or an ``(np, 2)`` array of coordinates.
        field: Field name.

    Returns:
        Dataset with dims ``(point, time)`` matching the observed-flowline
        layout: ``vel`` data var, coords ``x``/``y``/``d`` and ``time``. This
        drops straight into :func:`moc_py.plot.plot_flowline` as ``model=``.
    """
    x, y, d = _line_coords(flowline)
    vals = model_at_points(M, x, y, field=field)
    return xr.Dataset(
        {"vel": (("point", "time"), vals)},
        coords={"x": ("point", x), "y": ("point", y), "d": ("point", d),
                "time": ("time", M.time)},
        attrs={"name": "model", "field": field, "epsg": M.epsg},
    )


# ---------------------------------------------------------------------------
def _orient_elements(elements: np.ndarray) -> np.ndarray:
    """Return elements as ``(ne, 3)`` regardless of stored axis order."""
    elements = np.asarray(elements)
    if elements.ndim != 2:
        raise ValueError(f"elements must be 2-D, got shape {elements.shape}")
    if elements.shape[1] == 3:
        return elements
    if elements.shape[0] == 3:
        return elements.T
    raise ValueError(f"elements has no axis of length 3: {elements.shape}")


def _orient_field(arr: np.ndarray, nv: int, nt: int) -> np.ndarray:
    """Return a 2-D field as ``(nv, nt)`` regardless of stored axis order."""
    if arr.shape == (nv, nt):
        return arr
    if arr.shape == (nt, nv):
        return arr.T
    if arr.shape[0] == nv:
        return arr
    if arr.shape[1] == nv:
        return arr.T
    return arr


def _line_coords(fl):
    """Extract (x, y, d) from a flowline Dataset, ROI, or (np,2) array."""
    if isinstance(fl, np.ndarray) and fl.ndim == 2 and fl.shape[1] == 2:
        x, y = fl[:, 0].astype("float64"), fl[:, 1].astype("float64")
        d = None
    elif isinstance(fl, xr.Dataset) and "x" in fl.coords and "y" in fl.coords:
        x = np.asarray(fl.x.values, dtype="float64")
        y = np.asarray(fl.y.values, dtype="float64")
        d = np.asarray(fl.d.values, dtype="float64") if "d" in fl.coords else None
    elif hasattr(fl, "x") and hasattr(fl, "y"):  # Roi
        x = np.asarray(fl.x, dtype="float64")
        y = np.asarray(fl.y, dtype="float64")
        d = None
    else:
        raise TypeError("flowline must be an obs Dataset, a Roi, or an (np,2) array")
    if d is None:
        d = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    return x, y, d


def _attr(obj, name, default, is_group: bool = False):
    """Fetch an HDF5 attribute (from a dataset or the root group), with default."""
    attrs = obj.attrs if not is_group else obj["/"].attrs if "/" in obj else obj.attrs
    if name in attrs:
        v = np.asarray(attrs[name]).ravel()
        if v.size:
            return v[0]
    return default
