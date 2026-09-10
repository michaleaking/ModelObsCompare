"""Build a comparison-ready flowline from the archived shapefiles.

Takes a region of interest and produces the single centre-glacier flowline that
runs through it, clipped at the most recently traced terminus position:

    ROI  ->  Felikson flowlines that reach into it
         ->  centre flowline of the busiest glacier
         ->  split at the newest Black & Joughin terminus trace, upstream side
         ->  clipped to the ROI, resampled
         ->  xarray Dataset  ->  model_flowline / plot_flowline / compare_timeseries

The default ROI is ``sverdrups_trunk``, the one saved ROI that overlaps a
gridded observation subset and an exported model, so sections 3-4 have data to
work with. Any saved ROI name, ROI ``.json``/``.exp`` path, or bbox works::

    python examples/example_flowline_from_shapefiles.py
    python examples/example_flowline_from_shapefiles.py --roi my_roi --spacing 100
    python examples/example_flowline_from_shapefiles.py --glacier b045 --flowline 05
    python examples/example_flowline_from_shapefiles.py --save-figs
    python examples/example_flowline_from_shapefiles.py --write flowline.nc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402

DEFAULT_ROI = "sverdrups_trunk"
OBS_NC = "SentinelMonthly_Sverdrups.nc"
MODEL_NC = "exported/Model_NW_fric1_Snapshot.nc"
PROFILE_TIME = 2019.6


def survey(roi) -> None:
    """Print what the ROI catches before any clipping happens.

    Args:
        roi: Region of interest, in any form :func:`moc_py.roi_polygon` accepts.
    """
    F = mp.load_flowlines(roi)
    T = mp.load_termini(roi)
    print(f"\nflowlines : {len(F)} in the ROI, "
          f"glaciers {sorted(F['glacier'].unique())}")
    for g, sub in F.groupby("glacier"):
        print(f"    {g}: flowlines {sorted(sub['flowline'])}")
    if len(T):
        print(f"termini   : {len(T)} traces, {T['date'].min().date()} to "
              f"{T['date'].max().date()}")
    else:
        print("termini   : none in the ROI")


def overview_figure(FL, roi, plt):
    """Map the ROI, every flowline in it, and the terminus used for the clip.

    Args:
        FL: Packaged flowline Dataset from :func:`moc_py.flowline_from_roi`.
        roi: The region of interest.
        plt: The pyplot module.

    Returns:
        The matplotlib ``Figure``.
    """
    F = mp.load_flowlines(roi)
    T = mp.load_termini(roi)
    poly = mp.roi_polygon(roi)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(*poly.exterior.xy, "k-", lw=1.5, label="ROI")
    for _, r in F.iterrows():
        ax.plot(*r.geometry.xy, "-", color="0.75", lw=0.8)
    for _, r in T.iterrows():
        ax.plot(*r.geometry.xy, "-", color="#7fb2e5", lw=0.5)
    if len(T):
        newest = T.loc[T["date"].idxmax()]
        ax.plot(*newest.geometry.xy, "-", color="#1f6fb4", lw=2.5,
                label=f"terminus {FL.attrs['terminus_date'] or newest['date'].date()}")
    ax.plot(FL.x.values, FL.y.values, "-", color="#d81b1b", lw=2.5,
            label=f"clipped flowline ({FL.attrs['name']})")
    ax.plot(FL.x.values[0], FL.y.values[0], "o", color="#d81b1b", ms=7,
            label="d = 0 (terminus)")

    b = poly.bounds
    pad = 0.15 * max(b[2] - b[0], b[3] - b[1])
    ax.set_xlim(b[0] - pad, b[2] + pad)
    ax.set_ylim(b[1] - pad, b[3] + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m, EPSG:3413]")
    ax.set_ylabel("y [m, EPSG:3413]")
    ax.set_title("flowlines and termini in the ROI")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def main(argv: list[str] | None = None) -> int:
    """Run the example.

    Args:
        argv: Command-line arguments (default ``sys.argv[1:]``).

    Returns:
        Process exit code.
    """
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--roi", default=DEFAULT_ROI,
                    help="saved ROI name, .json/.exp path, or 'xmin,xmax,ymin,ymax'")
    ap.add_argument("--glacier", default=None, help="glacier id, e.g. a045")
    ap.add_argument("--flowline", default=None,
                    help="flowline id, e.g. 05 (default: the centre one)")
    ap.add_argument("--spacing", type=float, default=200.0,
                    help="even resampling interval [m]; 0 keeps native vertices")
    ap.add_argument("--obs", default=OBS_NC, help="gridded obs netCDF")
    ap.add_argument("--model", default=MODEL_NC, help="exported model netCDF")
    ap.add_argument("--write", default=None, help="save the flowline to this netCDF")
    ap.add_argument("--save-figs", action="store_true",
                    help="write PNGs to <repo>/figs/python instead of showing them")
    args = ap.parse_args(argv)

    import matplotlib
    if args.save_figs:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    roi = ([float(v) for v in args.roi.split(",")] if "," in args.roi
           else mp.load_roi(args.roi))

    # --- 1. what is in the ROI ---------------------------------------------
    survey(roi)

    # --- 2. the pipeline ----------------------------------------------------
    FL = mp.flowline_from_roi(roi, glacier=args.glacier, flowline=args.flowline,
                              spacing=args.spacing or None)
    a = FL.attrs
    print(f"\nflowline  : glacier {a['glacier']} line {a['flowline']}, "
          f"{FL.sizes['point']} points over {a['length'] / 1e3:.2f} km")
    print(f"            clipped at the {a['terminus_date']} terminus trace")
    print(f"            d = 0 at ({FL.x.values[0]:.0f}, {FL.y.values[0]:.0f})")

    figs = [("10_flowline_from_shapefiles", overview_figure(FL, roi, plt))]

    if args.write:
        FL.to_netcdf(args.write)
        print(f"            written to {args.write}")

    # --- 3. hand it the observations ---------------------------------------
    try:
        O = mp.load_obs_netcdf(args.obs, bands=["vv"])
    except FileNotFoundError:
        print(f"\n(no gridded obs at {args.obs} — stopping after the geometry)")
        O = None
    if O is not None:
        Fobs = mp.sample_obs_along_flowline(O, FL)
        good = np.isfinite(Fobs.vel.values)
        print(f"\nobs       : {Fobs.sizes['time']} frames sampled onto the line, "
              f"{100 * good.mean():.0f}% of point-frames valid")
        f, axs = plt.subplots(1, 2, figsize=(13, 4.8))
        mp.plot_flowline(Fobs, mode="hovmoller", ax=axs[0])
        mp.plot_flowline(Fobs, mode="profile", time=PROFILE_TIME, ax=axs[1])
        f.tight_layout()
        figs.append(("11_flowline_obs", f))

        # --- 4. and the model ----------------------------------------------
        try:
            M = mp.load_model(args.model)
        except FileNotFoundError:
            print(f"(no exported model at {args.model} — skipping the comparison)")
            M = None
        if M is not None:
            Fmod = mp.model_flowline(M, FL)
            C = mp.compare_timeseries(M, Fobs, dist=0.5 * a["length"])
            print(f"model-obs : at d={0.5 * a['length'] / 1e3:.1f} km  "
                  f"bias={C.stats['mean']:.0f}  RMSE={C.stats['rmse']:.0f} m/yr  "
                  f"r={C.stats['r']:.2f} over {C.stats['n']} paired times")
            f, ax = plt.subplots(figsize=(7.5, 5))
            mp.plot_flowline(Fobs, mode="profile", time=PROFILE_TIME,
                             model=Fmod, ax=ax)
            f.tight_layout()
            figs.append(("12_flowline_model_obs", f))

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
