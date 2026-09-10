"""No-plot smoke test for the moc_py readers (Phase 1).

Loads one of each observation type, prints shapes/ranges, and runs orientation
and band-mapping checks. Run as a module::

    python -m moc_py.smoke_test
    python -m moc_py.smoke_test --nc Sentinel_Subset_Upernavik.nc --fl jakobshavn.mat

Exits non-zero if any check FAILs.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from .config import CONFIG, MissingPathError
from .obs_flowline import load_obs_flowline
from .obs_netcdf import load_obs_netcdf


def _pick_first(directory, patterns):
    for pat in patterns:
        hits = sorted(directory.glob(pat), key=lambda p: p.stat().st_size)
        if hits:
            return hits[0]
    return None


def _pick_first_in(setting: str, patterns):
    """Smallest matching file in a configured folder, or None if unconfigured.

    Args:
        setting: Path setting name, e.g. ``"obs_netcdf_dir"``.
        patterns: Glob patterns, tried in order.

    Returns:
        A :class:`~pathlib.Path`, or ``None`` if the setting is not configured
        (the smoke test then reports that rather than prompting).
    """
    try:
        directory = CONFIG.get(setting, prompt=False)
    except MissingPathError:
        return None
    return _pick_first(directory, patterns)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="moc_py reader smoke test")
    ap.add_argument("--nc", default=None, help="obs netCDF file (bare name ok)")
    ap.add_argument("--fl", default=None, help="obs flowline .mat (bare name ok)")
    args = ap.parse_args(argv)

    results: list[tuple[str, str, str]] = []

    def rec(name, status, msg):
        results.append((name, status, msg))
        print(f"  [{status:4s}] {name:22s} {msg}")

    # --- netCDF ---------------------------------------------------------
    nc = args.nc or _pick_first_in("obs_netcdf_dir", ["*Subset*.nc", "*.nc"])
    print(f"\n--- OBS netCDF: {nc} ---")
    try:
        O = load_obs_netcdf(nc, bands=["vv", "vx", "vy"], trange=(2018, 2020))
        ny, nx, nt = O.sizes["y"], O.sizes["x"], O.sizes["time"]
        print(f"  dims y={ny} x={nx} time={nt} | epsg={O.attrs['epsg']}")
        print(f"  time {float(O.time.min()):.3f}..{float(O.time.max()):.3f}")
        for b in ("vv", "vx", "vy"):
            v = O[b].values
            print(f"  {b}: shape {v.shape} range {np.nanmin(v):.0f}..{np.nanmax(v):.0f} "
                  f"nNaN={100*np.mean(np.isnan(v)):.0f}%")
        if O["vv"].dims == ("y", "x", "time") and O["vv"].shape == (ny, nx, nt):
            rec("obsnc.orient", "PASS", f"vv is (y,x,time) {(ny,nx,nt)}")
        else:
            rec("obsnc.orient", "FAIL", f"vv dims/shape {O['vv'].dims}/{O['vv'].shape}")
        A = O["vv"].isel(time=0).values
        B = np.hypot(O["vx"].isel(time=0).values, O["vy"].isel(time=0).values)
        g = np.isfinite(A) & np.isfinite(B) & (A > 1)
        rel = float(np.median(np.abs(A[g] - B[g]) / A[g])) if g.sum() > 50 else np.nan
        if rel < 0.05:
            rec("obsnc.bands", "PASS", f"vv~hypot(vx,vy) rel.err {100*rel:.2f}%")
        else:
            rec("obsnc.bands", "WARN", f"vv vs hypot rel.err {100*rel:.2f}% - check band order")
    except Exception as e:  # noqa: BLE001
        rec("obsnc.load", "FAIL", str(e))

    # --- flowline -------------------------------------------------------
    fl = args.fl
    if fl is None:
        fl = _pick_first_in("obs_flowline_dir", ["jakobshavn.mat", "*.mat"])
    print(f"\n--- OBS flowline: {fl} ---")
    try:
        F = load_obs_flowline(fl)
        npts, ntf = F.sizes["point"], F.sizes["time"]
        print(f"  name={F.attrs['name']} source={F.attrs['source']} "
              f"pts={npts} time={ntf}")
        print(f"  d {float(F.d.min()):.0f}..{float(F.d.max()):.0f} m "
              f"monotonic={bool(np.all(np.diff(F.d.values) >= 0))}")
        print(f"  time {float(F.time.min()):.2f}..{float(F.time.max()):.2f}")
        print(f"  vel shape {F.vel.shape} range "
              f"{np.nanmin(F.vel.values):.0f}..{np.nanmax(F.vel.values):.0f}")
        if F.vel.dims == ("point", "time") and F.vel.shape == (npts, ntf):
            rec("flowline.orient", "PASS", f"vel is (point,time) {(npts, ntf)}")
        else:
            rec("flowline.orient", "FAIL", f"vel dims/shape {F.vel.dims}/{F.vel.shape}")
    except Exception as e:  # noqa: BLE001
        rec("flowline.load", "FAIL", str(e))

    # --- summary --------------------------------------------------------
    nP = sum(s == "PASS" for _, s, _ in results)
    nW = sum(s == "WARN" for _, s, _ in results)
    nF = sum(s == "FAIL" for _, s, _ in results)
    print(f"\n{'='*50}\n  {nP} PASS   {nW} WARN   {nF} FAIL")
    return 1 if nF else 0


if __name__ == "__main__":
    sys.exit(main())
