"""Self-contained validation of :func:`moc_py.velocity_variance` (no data needed).

Fabricates a two-domain grid — fast ice (2000 m/yr, +/-100) beside slow ice
(100 m/yr, +/-10) — so the ABSOLUTE spread is 10x larger on the fast side while
the RELATIVE (normalised) spread is 2x larger on the slow side, and asserts:

  * the normalised map recovers the known fractions (0.10 vs 0.20) while the
    absolute spread recovers 200 vs 20 m/yr,
  * an ROI restricts the output grid, a polygon ROI masks pixels inside its
    bbox but outside the polygon,
  * two-slice mode picks the nearest frames and the signed change is
    later-minus-earlier regardless of the ROI/time argument order,
  * pixels with too few finite frames, and references below `min_speed`, are
    NaN,
  * every statistic x normalize combination runs, and the plots draw.

Run:  python tests/test_variance_synthetic.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import xarray as xr

matplotlib.use("Agg")  # headless: the plot check must not need a display

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402

TIMES = np.array([2019.0, 2019.25, 2019.5, 2019.75])
FAST, SLOW = 2000.0, 100.0      # mean speeds [m/yr] of the two halves
FAST_AMP, SLOW_AMP = 100.0, 10.0  # +/- seasonal amplitude [m/yr]


def _make_obs(nx: int = 50, ny: int = 50) -> xr.Dataset:
    """Build a synthetic obs grid: fast ice on the left, slow ice on the right.

    The y coordinate DESCENDS, as in the real Sentinel subsets. One pixel is
    all-NaN to exercise the ``min_count`` path.

    Args:
        nx: Number of columns.
        ny: Number of rows.

    Returns:
        Dataset with ``vv`` on dims ``(y, x, time)``.
    """
    x = -300000 + 200.0 * np.arange(nx)
    y = -1810000 - 200.0 * np.arange(ny)
    fast = np.arange(nx) < nx // 2
    base = np.where(fast, FAST, SLOW)[None, :, None]
    amp = np.where(fast, FAST_AMP, SLOW_AMP)[None, :, None]
    season = np.array([-1.0, 0.0, 1.0, 0.0])[None, None, :]  # peak-to-peak = 2*amp
    vv = np.broadcast_to(base + amp * season, (ny, nx, TIMES.size)).astype("float32").copy()
    vv[0, 0, :] = np.nan  # dead pixel
    return xr.Dataset(
        {"vv": (("y", "x", "time"), vv)},
        coords={"x": ("x", x), "y": ("y", y), "time": ("time", TIMES)},
        attrs={"epsg": 3413},
    )


def _check(msg, cond):
    print(f"  [{'PASS' if cond else 'FAIL':4s}] {msg}")
    return bool(cond)


def _raises(fn, *args, **kwargs) -> bool:
    """True if ``fn`` raises ValueError."""
    try:
        fn(*args, **kwargs)
    except ValueError:
        return True
    return False


def run() -> bool:
    ok = True
    O = _make_obs()
    nx = O.sizes["x"]

    # --- full grid, whole series -------------------------------------------
    V = mp.velocity_variance(O)
    f_norm = float(V["norm"].values[5, 5])            # fast half
    s_norm = float(V["norm"].values[5, nx - 5])       # slow half
    f_spread = float(V["spread"].values[5, 5])
    s_spread = float(V["spread"].values[5, nx - 5])
    ok &= _check(f"absolute spread is fast-biased ({f_spread:.0f} vs {s_spread:.0f} m/yr)",
                 np.isclose(f_spread, 2 * FAST_AMP) and np.isclose(s_spread, 2 * SLOW_AMP))
    ok &= _check(f"normalised spread inverts it ({f_norm:.3f} vs {s_norm:.3f})",
                 np.isclose(f_norm, 2 * FAST_AMP / FAST)
                 and np.isclose(s_norm, 2 * SLOW_AMP / SLOW))
    ok &= _check("all-NaN pixel is NaN with count 0",
                 np.isnan(V["norm"].values[0, 0]) and V["count"].values[0, 0] == 0)
    ok &= _check(f"ROI-mean summary in attrs (norm_mean={V.attrs['norm_mean']:.3f})",
                 np.isclose(V.attrs["norm_mean"], 0.5 * (2 * FAST_AMP / FAST + 2 * SLOW_AMP / SLOW),
                            atol=1e-3))
    ok &= _check("full grid used when roi is None",
                 (V.sizes["y"], V.sizes["x"]) == (O.sizes["y"], nx))

    # --- rectangular ROI + two time slices ---------------------------------
    roi = mp.Roi("rightbox", "rectangle",
                 np.array([-295000.0, -290500, -290500, -295000, -295000]),
                 np.array([-1819000.0, -1819000, -1812000, -1812000, -1819000]))
    V2 = mp.velocity_variance(O, roi, times=(2019.02, 2019.48))
    ok &= _check(f"ROI crops the grid ({V2.sizes['y']}x{V2.sizes['x']} of "
                 f"{O.sizes['y']}x{nx})",
                 V2.sizes["x"] < nx and V2.sizes["y"] < O.sizes["y"])
    ok &= _check(f"two nearest frames selected (t={list(np.round(V2.time.values, 2))})",
                 V2.attrs["nframes"] == 2
                 and np.allclose(V2.time.values, [2019.0, 2019.5]))
    ok &= _check(f"mean normalised change over the ROI = {V2.attrs['norm_mean']:.3f}",
                 np.isclose(V2.attrs["norm_mean"], 2 * SLOW_AMP / SLOW))
    ok &= _check("signed change is later minus earlier (+2*amp here)",
                 np.isclose(np.nanmean(V2["change"].values), 2 * SLOW_AMP)
                 and np.isclose(np.nanmean(V2["norm_change"].values), 2 * SLOW_AMP / SLOW))
    ok &= _check("no signed change field for a >2-frame window",
                 "change" not in mp.velocity_variance(O, roi))

    # --- polygon ROI masks inside its own bbox -----------------------------
    tri = mp.Roi("tri", "polygon",
                 np.array([-300000.0, -290200, -300000, -300000]),
                 np.array([-1819800.0, -1819800, -1810000, -1819800]))
    V3 = mp.velocity_variance(O, tri, trange=(2019.0, 2019.6),
                              statistic="std", normalize="median")
    inside = np.isfinite(V3["norm"].values).mean()
    ok &= _check(f"polygon masks ~half its bbox ({inside * 100:.0f}% kept)",
                 0.35 < inside < 0.65)

    # --- guards -------------------------------------------------------------
    V4 = mp.velocity_variance(O, min_speed=500.0)
    ok &= _check("min_speed drops slow-ice references",
                 np.isnan(V4["norm"].values[5, nx - 5])
                 and np.isfinite(V4["norm"].values[5, 5]))
    ok &= _check("min_count drops thin pixels",
                 np.isnan(mp.velocity_variance(O, min_count=5)["norm"].values).all())
    ok &= _check("bad arguments raise ValueError",
                 _raises(mp.velocity_variance, O, statistic="nope")
                 and _raises(mp.velocity_variance, O, normalize="nope")
                 and _raises(mp.velocity_variance, O, band="zz")
                 and _raises(mp.velocity_variance, O, trange=(2050, 2060))
                 and _raises(mp.velocity_variance, O, times=(2019.0,))
                 and _raises(mp.velocity_variance, O, times=(2019.0, 2019.5),
                             trange=(2019.0, 2019.5)))

    # --- every statistic x normalize combination ----------------------------
    combos_ok = True
    for stat in mp.variance.STATISTICS:
        for nrm in mp.variance.NORMALIZERS:
            v = mp.velocity_variance(O, statistic=stat, normalize=nrm)
            combos_ok &= np.isfinite(v.attrs["norm_mean"]) and v.attrs["npix"] > 0
    ok &= _check(f"{len(mp.variance.STATISTICS)}x{len(mp.variance.NORMALIZERS)} "
                 "statistic/normalize combinations run", combos_ok)

    # --- plotting -----------------------------------------------------------
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    mp.plot_variance(V, "norm", ax=axs[0], roi=roi)
    mp.plot_variance(V, "spread", ax=axs[1])
    mp.plot_variance(V2, "norm_change", ax=axs[2])
    plt.close(fig)
    ok &= _check("plot_variance draws norm / spread / norm_change", True)
    ok &= _check("plot_variance rejects an absent field",
                 _raises(mp.plot_variance, V, "change"))
    return ok


if __name__ == "__main__":
    print("=== velocity_variance (synthetic) ===")
    passed = run()
    print("\n" + ("ALL VARIANCE CHECKS PASSED" if passed else "SOME CHECKS FAILED"))
    sys.exit(0 if passed else 1)
