"""moc_py — Python-native readers, model interpolation, and plots for
comparing observed vs ISSM-modelled glacier velocities.

Phase 1 (observations, pure Python):
  - :func:`load_obs_netcdf`   gridded Sentinel/CSK velocity subsets -> xarray
  - :func:`load_obs_flowline` per-glacier flowline .mat files       -> xarray
  - ROI tools (:class:`Roi`, :func:`roi_from_obs`, :class:`RoiSelector`, ...)
  - :func:`flowline_from_roi`  Felikson flowline shapefiles cropped to an ROI
    and clipped at the most recent Black & Joughin terminus trace

Phase 2 (model side, needs a netCDF from ``moc_export_model.m``):
  - :func:`load_model`        exported ISSM mesh + velocity          -> Model
  - :func:`model_to_grid` / :func:`model_at_points` / :func:`model_flowline`
    (mesh interpolation via matplotlib.tri, == ISSM InterpFromMesh2d)
  - :func:`spatial_diff` / :func:`compare_timeseries` and their plots

Data locations are not hard-coded: set them with :func:`configure`, the
``MOC_*`` environment variables, or ``python -m moc_py.config --set`` (which
prompts once and remembers). See :mod:`moc_py.config`.

Typical use::

    import moc_py as mp
    mp.configure(obs_netcdf_dir="~/data/velocity")   # once per machine
    O = mp.load_obs_netcdf("Sentinel_Subset_Upernavik.nc", bands=["vv"])
    M = mp.load_model("nw_fric1.nc")            # exported from MATLAB
    roi = mp.roi_from_obs(O, margin=2000)
    D = mp.spatial_diff(M, O, roi, time=2019.6)
    mp.plot_diff(D, roi=roi)
"""

from __future__ import annotations

from .config import CONFIG, Config, configure, save_config, MissingPathError
from .obs_flowline import load_obs_flowline
from .obs_netcdf import load_obs_netcdf
from .obs_lazy import open_obs_dataarray, to_obs_dataset
from .roi import Roi, RoiSelector, load_roi, roi_from_obs, roi_mask, save_roi
from .flowline_shp import (
    roi_polygon, load_flowlines, load_termini, center_flowline, ordered_flowlines,
    latest_terminus,
    clip_upstream, clip_to_roi, flowline_to_dataset, flowline_from_roi,
    sample_obs_along_flowline,
)
from .model import (
    Model, load_model, model_to_grid, model_at_points, model_flowline,
)
from .compare import spatial_diff, compare_timeseries, TimeSeriesComparison
from .variance import velocity_variance
from .plot import (
    plot_flowline, plot_map, plot_timeseries,
    plot_grid, plot_diff, plot_compare, plot_variance,
)
from . import util

__all__ = [
    "CONFIG", "Config", "configure", "save_config", "MissingPathError",
    # obs
    "load_obs_netcdf", "load_obs_flowline",
    "open_obs_dataarray", "to_obs_dataset",
    # roi
    "Roi", "RoiSelector", "roi_from_obs", "roi_mask", "save_roi", "load_roi",
    # flowlines from shapefiles (needs the "shapefile" extra)
    "roi_polygon", "load_flowlines", "load_termini", "center_flowline",
    "ordered_flowlines", "latest_terminus", "clip_upstream", "clip_to_roi", "flowline_to_dataset",
    "flowline_from_roi", "sample_obs_along_flowline",
    # model (Phase 2)
    "Model", "load_model", "model_to_grid", "model_at_points", "model_flowline",
    # compare (Phase 2)
    "spatial_diff", "compare_timeseries", "TimeSeriesComparison",
    # temporal variability of the observations
    "velocity_variance",
    # plotting
    "plot_map", "plot_flowline", "plot_timeseries",
    "plot_grid", "plot_diff", "plot_compare", "plot_variance",
    "util",
]
