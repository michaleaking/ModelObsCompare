"""Region-of-interest helpers: build, mask, save/load, and interactively draw.

An ROI is a small, serialisable polygon in projected (EPSG:3413) coordinates.
It can be derived from an observed grid's extent, drawn interactively in a
notebook, tested against points, and round-tripped to JSON or an ISSM ``.exp``
contour for interoperability with the MATLAB toolkit.

Mirrors the ``roi/`` functions on the MATLAB side.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from .config import CONFIG


@dataclass
class Roi:
    """A region of interest.

    Attributes:
        name: Short identifier (used as the saved filename stem).
        kind: ``"rectangle"`` or ``"polygon"``.
        x: Closed polygon vertex x-coordinates [m].
        y: Closed polygon vertex y-coordinates [m].
        epsg: Projection code of the coordinates.
        file: Path the ROI was saved to (empty if unsaved).
    """

    name: str
    kind: str
    x: np.ndarray
    y: np.ndarray
    epsg: int = CONFIG.epsg
    file: str = ""

    @property
    def bbox(self) -> tuple[float, float, float, float]:
        """Bounding box ``(xmin, xmax, ymin, ymax)``."""
        return (float(np.min(self.x)), float(np.max(self.x)),
                float(np.min(self.y)), float(np.max(self.y)))


def roi_from_obs(ds, margin: float = 0.0, name: str = "") -> Roi:
    """Build a rectangular ROI from the extent of an observed grid Dataset.

    Args:
        ds: Dataset from :func:`load_obs_netcdf` (uses ``ds.x``, ``ds.y``, epsg).
        margin: Inset from each edge [m] (avoids no-data border pixels).
        name: If given, the ROI is also saved to ``CONFIG.roi_dir``.

    Returns:
        A rectangular :class:`Roi`.
    """
    x = np.asarray(ds.x.values, dtype="float64")
    y = np.asarray(ds.y.values, dtype="float64")
    xmin, xmax = x.min() + margin, x.max() - margin
    ymin, ymax = y.min() + margin, y.max() - margin
    roi = Roi(
        name=name,
        kind="rectangle",
        x=np.array([xmin, xmax, xmax, xmin, xmin]),
        y=np.array([ymin, ymin, ymax, ymax, ymin]),
        epsg=int(ds.attrs.get("epsg", CONFIG.epsg)),
    )
    if name:
        save_roi(roi)
    return roi


def roi_mask(roi: Roi, x, y) -> np.ndarray:
    """Test which points fall inside an ROI polygon.

    Args:
        roi: The region of interest.
        x, y: Coordinate arrays (matching shape; grids from ``meshgrid`` are fine).

    Returns:
        Boolean array (shape of ``x``) that is True inside the ROI.
    """
    from matplotlib.path import Path as MplPath

    x = np.asarray(x)
    y = np.asarray(y)
    pts = np.column_stack([x.ravel(), y.ravel()])
    poly = MplPath(np.column_stack([roi.x, roi.y]))
    return poly.contains_points(pts).reshape(x.shape)


def save_roi(roi: Roi, directory: str | Path | None = None) -> Roi:
    """Save an ROI to ``<dir>/<name>.json`` and an ISSM ``<name>.exp`` contour.

    Args:
        roi: The ROI to save (``roi.name`` must be set).
        directory: Output folder (default ``CONFIG.roi_dir``).

    Returns:
        The ROI with ``file`` set to the JSON path.
    """
    directory = Path(directory) if directory else CONFIG.roi_dir
    directory.mkdir(parents=True, exist_ok=True)
    if not roi.name:
        raise ValueError("roi.name must be set before saving")

    jpath = directory / f"{roi.name}.json"
    payload = {**asdict(roi), "x": list(map(float, roi.x)), "y": list(map(float, roi.y))}
    jpath.write_text(json.dumps(payload, indent=2))
    roi.file = str(jpath)

    _write_exp(roi, directory / f"{roi.name}.exp")
    return roi


def load_roi(name: str | Path, directory: str | Path | None = None) -> Roi:
    """Load an ROI from a ``.json`` (or ``.exp``) file.

    Args:
        name: ROI name or a path to a ``.json``/``.exp`` file.
        directory: Folder to search (default ``CONFIG.roi_dir``).

    Returns:
        The loaded :class:`Roi`.
    """
    directory = Path(directory) if directory else CONFIG.roi_dir
    p = Path(name)
    for cand in (p, directory / p, directory / f"{p}.json", directory / f"{p}.exp"):
        if cand.is_file():
            if cand.suffix == ".json":
                d = json.loads(cand.read_text())
                return Roi(name=d["name"], kind=d["kind"],
                          x=np.array(d["x"]), y=np.array(d["y"]),
                          epsg=int(d.get("epsg", CONFIG.epsg)), file=str(cand))
            if cand.suffix == ".exp":
                x, y = _read_exp(cand)
                return Roi(name=cand.stem, kind="polygon", x=x, y=y, file=str(cand))
    raise FileNotFoundError(f"No ROI found for '{name}' in {directory}")


class RoiSelector:
    """Interactive ROI drawing for Jupyter (needs ``%matplotlib widget``).

    Draw a rectangle (or polygon) on an observed-velocity backdrop; the result
    is available as :attr:`roi` once you finish the selection.

    Example:
        >>> sel = RoiSelector(O, band="vv", tindex=0, kind="rectangle")
        >>> # drag on the figure ...
        >>> roi = sel.roi                     # then use / save it
    """

    def __init__(self, ds, band: str = "vv", tindex: int = 0,
                 kind: str = "rectangle", name: str = "", clim=None):
        from matplotlib.widgets import RectangleSelector, PolygonSelector
        from . import plot as _plot

        self.name = name
        self.kind = kind
        self.epsg = int(ds.attrs.get("epsg", CONFIG.epsg))
        self.roi: Roi | None = None
        self.ax = _plot.plot_map(ds, band=band, tindex=tindex, clim=clim,
                                 title=f"Draw {kind} ROI")
        self._fig = self.ax.figure
        if kind == "rectangle":
            self._sel = RectangleSelector(self.ax, self._on_rect,
                                          useblit=True, interactive=True)
        else:
            self._sel = PolygonSelector(self.ax, self._on_poly, useblit=True)

    def _on_rect(self, eclick, erelease):
        x1, x2 = sorted((eclick.xdata, erelease.xdata))
        y1, y2 = sorted((eclick.ydata, erelease.ydata))
        self.roi = Roi(self.name, "rectangle",
                       np.array([x1, x2, x2, x1, x1]),
                       np.array([y1, y1, y2, y2, y1]), self.epsg)

    def _on_poly(self, verts):
        v = np.asarray(verts)
        vx = np.r_[v[:, 0], v[0, 0]]
        vy = np.r_[v[:, 1], v[0, 1]]
        self.roi = Roi(self.name, "polygon", vx, vy, self.epsg)

    def save(self) -> Roi:
        """Save the current selection (raises if nothing drawn yet)."""
        if self.roi is None:
            raise RuntimeError("No ROI drawn yet.")
        if not self.roi.name:
            raise ValueError("Set a name before saving (RoiSelector(..., name=...)).")
        return save_roi(self.roi)


# ---------------------------------------------------------------------------
def _write_exp(roi: Roi, path: Path) -> None:
    """Write an ISSM Argus ``.exp`` contour (interop with the MATLAB toolkit)."""
    x = np.asarray(roi.x, dtype="float64")
    y = np.asarray(roi.y, dtype="float64")
    lines = [f"## Name:{roi.name or 'roi'}", "## Icon:0",
             "# Points Count Value", f"{x.size} 1.", "# X pos Y pos"]
    lines += [f"{xi:.6f} {yi:.6f}" for xi, yi in zip(x, y)]
    path.write_text("\n".join(lines) + "\n")


def _read_exp(path: Path):
    """Read x, y vertex arrays from an ISSM ``.exp`` contour."""
    xs, ys = [], []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(xs), np.array(ys)
