"""End-to-end tour of moc_py on the data that is actually on disk.

Python mirror of ``examples/example_sverdrups_workflow.m`` on the MATLAB side.
Sverdrups Glacier (NW Greenland) is the one place in the configured folders
where a gridded observation subset and a redrawn flowline cover the same ice:

    obs grid : SentinelMonthly_Sverdrups.nc  (72 monthly frames, 2015-2020)
    flowline : sverdrups.mat                 (31 points, 0-7.5 km, 2015-2024)
    model    : exported/Model_NW_fric1_Snapshot.nc  (376 steps, 2007-2022)

Demonstrates, in order:
    1. the two observation readers,
    2. an ROI built from the flowline footprint, saved (.json + ISSM .exp) and
       read back,
    3. a speed map with the ROI and flowline overlaid,
    4. :func:`moc_py.velocity_variance` — melt-season window, two individual
       time slices, and a robust record-long variant + coverage map,
    5. a per-year table of relative change over the ROI,
    6. the flowline: distance-time image and a point series checked against the
       nearest gridded pixel,
    7. the model: spatial difference, flowline profile, and a paired series.

Unlike the MATLAB tour, section 7 runs by default: the exported model netCDF
loads in ~1 s, where the MATLAB path would read a multi-GB ISSM ``.mat``.

Run::

    python examples/example_sverdrups_workflow.py                  # show figures
    python examples/example_sverdrups_workflow.py --save-figs      # write PNGs
    python examples/example_sverdrups_workflow.py --no-model       # obs only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402

OBS_NC = "SentinelMonthly_Sverdrups.nc"
OBS_FL = "sverdrups.mat"
MODEL_NC = "exported/Model_NW_fric1_Snapshot.nc"

MELT_WINDOW = (2019.35, 2019.85)   # May-Sep 2019
TWO_SLICES = (2019.15, 2019.6)     # late winter -> late summer
POINT_DIST = 3000.0                # [m] along the flowline
MIN_SPEED = 20.0                   # [m/yr] below this a ratio is meaningless


def build_roi(F: xr.Dataset, pad: float = 3000.0,
              name: str = "sverdrups_trunk") -> mp.Roi:
    """Build, save and re-load a rectangular ROI around a flowline's footprint.

    Args:
        F: Flowline Dataset from :func:`moc_py.load_obs_flowline`.
        pad: Margin around the flowline bounding box [m].
        name: ROI name; also the saved filename stem in ``CONFIG.roi_dir``.

    Returns:
        The ROI, round-tripped through the saved ``.json`` (an ``.exp`` contour
        is written alongside it for the MATLAB/ISSM side).
    """
    x0, x1 = F.x.values.min() - pad, F.x.values.max() + pad
    y0, y1 = F.y.values.min() - pad, F.y.values.max() + pad
    roi = mp.Roi(name=name, kind="rectangle",
                 x=np.array([x0, x1, x1, x0, x0]),
                 y=np.array([y0, y0, y1, y1, y0]))
    mp.save_roi(roi)
    return mp.load_roi(name)


def yearly_table(O: xr.Dataset, roi: mp.Roi, years: range) -> None:
    """Print per-year relative change over the ROI for the melt season.

    Args:
        O: Gridded obs Dataset.
        roi: Region of interest.
        years: Years to tabulate (Apr-Nov of each is used).
    """
    print("\n  melt season (Apr-Nov) over the ROI")
    print("  year | frames | px    | mean rel. | median rel. | mean range [m/yr]")
    print("  -----+--------+-------+-----------+-------------+------------------")
    for yr in years:
        try:
            V = mp.velocity_variance(O, roi, trange=(yr + 0.25, yr + 0.85),
                                     min_count=3, min_speed=MIN_SPEED)
        except ValueError as e:
            print(f"  {yr:4d} |  skipped ({e})")
            continue
        a = V.attrs
        print(f"  {yr:4d} |  {a['nframes']:5d} | {a['npix']:5d} |"
              f"   {100 * a['norm_mean']:6.1f}% |     {100 * a['norm_median']:6.1f}% |"
              f"   {a['spread_mean']:8.0f}")


def flowline_figure(F: xr.Dataset, O: xr.Dataset, plt) -> "plt.Figure":
    """Draw the flowline distance-time image beside one point's time series.

    The gridded product is overlaid on the point series as an independent check
    that the two readers agree where they overlap.

    Args:
        F: Flowline Dataset.
        O: Gridded obs Dataset.
        plt: The pyplot module (passed in so the backend is chosen by the caller).

    Returns:
        The matplotlib ``Figure``.
    """
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    mp.plot_flowline(F, mode="hovmoller", ax=axs[0])

    k = int(np.argmin(np.abs(F.d.values - POINT_DIST)))
    mp.plot_timeseries(F, point=k, ax=axs[1], use_datetime=False)
    ci = int(np.argmin(np.abs(O.x.values - float(F.x.values[k]))))
    ri = int(np.argmin(np.abs(O.y.values - float(F.y.values[k]))))
    axs[1].plot(O.time.values, O["vv"].isel(x=ci, y=ri).values, "s-", ms=3,
                color="#d81b1b", lw=1, label="gridded obs (nearest px)")
    axs[1].legend()
    fig.tight_layout()
    return fig


def model_figures(M, O: xr.Dataset, F: xr.Dataset, roi: mp.Roi, plt) -> list:
    """Model-vs-observation panels: spatial difference, profile, paired series.

    Args:
        M: Model from :func:`moc_py.load_model`.
        O: Gridded obs Dataset.
        F: Flowline Dataset.
        roi: Region of interest.
        plt: The pyplot module.

    Returns:
        List of ``(name, Figure)`` pairs.
    """
    D = mp.spatial_diff(M, O, roi, time=2019.6, obswindow=0.05)
    print(f"\nmodel - obs at t={D.attrs['time']:.2f} over the ROI: "
          f"bias={D.attrs['bias']:.0f}  RMSE={D.attrs['rmse']:.0f}  "
          f"MAE={D.attrs['mae']:.0f} m/yr over {D.attrs['n']} px")
    f6 = mp.plot_diff(D, roi=roi)

    Fmod = mp.model_flowline(M, F)
    f7, ax = plt.subplots(figsize=(7.5, 5))
    mp.plot_flowline(F, mode="profile", time=2019.6, model=Fmod, ax=ax)
    f7.tight_layout()

    C = mp.compare_timeseries(M, F, dist=POINT_DIST)
    print(f"at d={POINT_DIST / 1e3:.1f} km: bias={C.stats['mean']:.0f}  "
          f"RMSE={C.stats['rmse']:.0f} m/yr  r={C.stats['r']:.2f} "
          f"over {C.stats['n']} paired times")
    f8, ax = plt.subplots(figsize=(8, 4.8))
    mp.plot_compare(C, ax=ax, use_datetime=False)
    f8.tight_layout()

    return [("06_model_obs_diff", f6), ("07_flowline_profile", f7),
            ("08_point_timeseries", f8)]


def main(argv: list[str] | None = None) -> int:
    """Run the tour.

    Args:
        argv: Command-line arguments (default ``sys.argv[1:]``).

    Returns:
        Process exit code.
    """
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--obs", default=OBS_NC, help="gridded obs netCDF")
    ap.add_argument("--flowline", default=OBS_FL, help="flowline .mat")
    ap.add_argument("--model", default=MODEL_NC, help="exported model netCDF")
    ap.add_argument("--no-model", action="store_true", help="skip section 7")
    ap.add_argument("--save-figs", action="store_true",
                    help="write PNGs to <repo>/figs/python instead of showing them")
    args = ap.parse_args(argv)

    import matplotlib
    if args.save_figs:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figs = []

    # --- 1. observations ----------------------------------------------------
    O = mp.load_obs_netcdf(args.obs, bands=["vv"])
    F = mp.load_obs_flowline(args.flowline)
    t = O.time.values
    print(f"\nobs grid  : {O.sizes['x']}x{O.sizes['y']} px, {t.size} frames "
          f"{t.min():.2f}-{t.max():.2f}, "
          f"{100 * np.isfinite(O['vv'].values).mean():.0f}% of pixel-frames valid")
    ft = F.time.values
    print(f"flowline  : {F.attrs.get('name', args.flowline)}, {F.sizes['point']} points "
          f"(0-{F.d.values.max() / 1e3:.1f} km), {ft.size} times "
          f"{ft.min():.2f}-{ft.max():.2f} (source {F.attrs.get('source', '?')})")

    # --- 2. ROI from the flowline footprint --------------------------------
    roi = build_roi(F)
    Xo, Yo = np.meshgrid(O.x.values, O.y.values)
    xmin, xmax, ymin, ymax = roi.bbox
    print(f"ROI       : {(xmax - xmin) / 1e3:.1f} x {(ymax - ymin) / 1e3:.1f} km, "
          f"{int(mp.roi_mask(roi, Xo, Yo).sum())} of {Xo.size} grid pixels inside")

    # --- 3. speed map with the ROI and flowline on top ---------------------
    ti = int(np.argmin(np.abs(t - 2019.6)))
    f1, ax = plt.subplots(figsize=(7.5, 6))
    mp.plot_map(O, band="vv", tindex=ti, ax=ax, clim=(0, 3000), roi=roi)
    ax.plot(F.x.values, F.y.values, "w-", lw=2.5)
    ax.plot(F.x.values, F.y.values, "k--", lw=1)
    f1.tight_layout()
    figs.append(("01_speed_map", f1))

    # --- 4a. melt-season variability, normalised ---------------------------
    # min_speed matters here: without it the ROI mean jumps from ~0.4 to ~1.2 as
    # near-stagnant margin pixels divide a few m/yr of noise by a few m/yr of speed.
    V = mp.velocity_variance(O, roi, trange=MELT_WINDOW,
                             min_count=3, min_speed=MIN_SPEED)
    f2, axs = plt.subplots(1, 2, figsize=(13, 5))
    mp.plot_variance(V, "norm", ax=axs[0], roi=roi,
                     title="relative: range / mean speed [-]")
    mp.plot_variance(V, "spread", ax=axs[1], roi=roi,
                     title="absolute: range of speed [m/yr]")
    f2.tight_layout()
    figs.append(("02_melt_season_variability", f2))

    # --- 4b. between two individual time slices ----------------------------
    V2 = mp.velocity_variance(O, roi, times=TWO_SLICES)
    t0, t1 = V2.time.values
    f3, axs = plt.subplots(1, 2, figsize=(13, 5))
    mp.plot_variance(V2, "norm_change", ax=axs[0], roi=roi,
                     title=f"relative change {t0:.2f} -> {t1:.2f} [-]")
    mp.plot_variance(V2, "change", ax=axs[1], roi=roi,
                     title=f"absolute change {t0:.2f} -> {t1:.2f} [m/yr]")
    f3.tight_layout()
    figs.append(("03_two_slice_change", f3))
    print(f"\nlate winter ({t0:.2f}) -> late summer ({t1:.2f}): "
          f"mean |change| = {100 * V2.attrs['norm_mean']:.0f}% of local speed")

    # --- 4c. robust, record-long variant + coverage ------------------------
    Vr = mp.velocity_variance(O, roi, statistic="iqr", normalize="median",
                              min_speed=MIN_SPEED, min_count=24)
    f4, axs = plt.subplots(1, 2, figsize=(13, 5))
    mp.plot_variance(Vr, "norm", ax=axs[0], roi=roi,
                     title=f"{t.min():.0f}-{t.max():.0f} IQR / median speed [-]"
                           f"  (mean {Vr.attrs['norm_mean']:.2f})")
    mp.plot_variance(Vr, "count", ax=axs[1], roi=roi,
                     title=f"valid frames per pixel (of {t.size})")
    f4.tight_layout()
    figs.append(("04_record_long_robust", f4))

    # --- 5. per-year table --------------------------------------------------
    yearly_table(O, roi, range(int(np.ceil(t.min())), int(np.floor(t.max())) + 1))

    # --- 6. the flowline ----------------------------------------------------
    figs.append(("05_flowline", flowline_figure(F, O, plt)))

    # --- 7. the model -------------------------------------------------------
    if args.no_model:
        print(f"\n(model section skipped — drop --no-model to load {args.model})")
    else:
        M = mp.load_model(args.model)
        print(f"\nmodel     : {M.x.size} vertices, {M.elements.shape[0]} elements, "
              f"{M.time.size} steps {M.time.min():.2f}-{M.time.max():.2f}")
        figs.extend(model_figures(M, O, F, roi, plt))

    # --- output -------------------------------------------------------------
    if args.save_figs:
        outdir = Path(__file__).resolve().parents[2] / "figs" / "python"
        outdir.mkdir(parents=True, exist_ok=True)
        for name, fig in figs:
            fig.savefig(outdir / f"{name}.png", dpi=120, bbox_inches="tight")
        print(f"\n{len(figs)} figures written to {outdir}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
