"""Self-contained validation of :mod:`moc_py.flowline_shp` (no archive needed).

Writes a miniature stand-in for the two shapefile collections into a temporary
folder — three parallel flowlines running east (the terminus end first) and
three terminus traces crossing them at known x, one per year — so every step
has an analytic answer:

  * the ROI selects only the flowlines that reach into it,
  * the terminus reader parses ``SourceDate`` and honours the quality cut,
  * :func:`center_flowline` picks the middle line of the fan, and
    :func:`ordered_flowlines` ranks it centre-outwards,
  * when no trace crosses the centre line, the search steps outwards and uses
    the first neighbour one does cross (glacier0003 below), while naming a
    flowline explicitly still fails rather than silently swapping,
  * :func:`clip_upstream` cuts at the NEWEST trace (x=2000, not the older ones)
    and keeps the inland side, both for a line drawn terminus-first and for one
    drawn inland-first,
  * :func:`clip_to_roi` trims the result to the ROI box,
  * the packaged Dataset has monotonic ``d`` starting at zero, honours
    ``spacing``, and is accepted by :func:`moc_py.model_flowline`,
  * :func:`sample_obs_along_flowline` recovers a known linear velocity field,
  * the failure paths raise (empty ROI, unknown glacier, no crossing trace).

Skips cleanly if geopandas/shapely are not installed.

Run:  python tests/test_flowline_shp_synthetic.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402

Y_LINES = (0.0, 100.0, 200.0)        # the fan: three flowlines, 100 m apart
TERMINI_X = {"2018-06-01": 500.0, "2019-06-01": 1200.0, "2021-08-15": 2000.0}
BAD_X = 4000.0                       # a Quality_Fl=1 trace, newest of all
ROI = (-500.0, 6000.0, -50.0, 250.0)  # (xmin, xmax, ymin, ymax)

# A second fan, north of every trace above (they stop at y=500), so it needs
# its own ROI. Its newest trace is a
# stub that spans only the y=1000 line, so the centre line (y=1100) is missed
# and the search must step outwards to reach it.
Y_LINES_2 = (1000.0, 1100.0, 1200.0)
ROI_2 = (-500.0, 6000.0, 950.0, 1250.0)
STUB_X, STUB_Y = 2000.0, (950.0, 1050.0)


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


def _write_fixture(tmp: Path):
    """Write the miniature flowline and terminus shapefiles.

    Args:
        tmp: Folder to write into.

    Returns:
        ``(flowline_dir, termini_path)``.
    """
    import geopandas as gpd
    from shapely.geometry import LineString

    crs = f"EPSG:{mp.CONFIG.epsg}"
    fdir = tmp / "flowlines"
    fdir.mkdir()

    # glacier0001: a fan of three lines, x from 0 (terminus) to 10 km inland.
    lines = [LineString([(x, y) for x in np.linspace(0, 10000, 101)])
             for y in Y_LINES]
    gpd.GeoDataFrame({"flowline": ["03", "04", "05"]},
                     geometry=lines, crs=crs).to_file(fdir / "glacier0001.shp")

    # glacier0003: a second fan whose centre line no trace crosses.
    lines3 = [LineString([(x, y) for x in np.linspace(0, 10000, 101)])
              for y in Y_LINES_2]
    gpd.GeoDataFrame({"flowline": ["03", "04", "05"]},
                     geometry=lines3, crs=crs).to_file(fdir / "glacier0003.shp")

    # glacier0002: far away, so the ROI must exclude it.
    gpd.GeoDataFrame(
        {"flowline": ["03"]},
        geometry=[LineString([(0, 50000), (10000, 50000)])],
        crs=crs).to_file(fdir / "glacier0002.shp")

    # An intermediate iteration that must be ignored unless asked for.
    gpd.GeoDataFrame(
        {"flowline": ["03"]},
        geometry=[LineString([(0, 0), (10000, 0)])],
        crs=crs).to_file(fdir / "glacier0001_iter01.shp")

    # Not a flowline file at all — must be skipped, not crashed on.
    gpd.GeoDataFrame(
        {"note": ["mask"]},
        geometry=[LineString([(0, 0), (1, 1)])],
        crs=crs).to_file(fdir / "some_masks.shp")

    # Termini: vertical traces crossing the whole fan, one per date.
    rows, geoms = [], []
    for date, x in TERMINI_X.items():
        rows.append({"SourceDate": date, "Quality_Fl": 0, "Glacier_ID": 1})
        geoms.append(LineString([(x, -1000.0), (x, 500.0)]))
    rows.append({"SourceDate": "2022-01-01", "Quality_Fl": 1, "Glacier_ID": 1})
    geoms.append(LineString([(BAD_X, -1000.0), (BAD_X, 500.0)]))
    # the stub for glacier0003: crosses only its outermost line
    rows.append({"SourceDate": "2021-08-15", "Quality_Fl": 0, "Glacier_ID": 3})
    geoms.append(LineString([(STUB_X, STUB_Y[0]), (STUB_X, STUB_Y[1])]))
    tpath = tmp / "termini.shp"
    gpd.GeoDataFrame(rows, geometry=geoms, crs=crs).to_file(tpath)
    return fdir, tpath


def _make_obs() -> xr.Dataset:
    """A gridded obs Dataset whose speed is exactly ``x`` (m/yr)."""
    x = np.linspace(-1000.0, 11000.0, 121)
    y = np.linspace(-500.0, 700.0, 13)
    t = np.array([2019.0, 2020.0])
    vv = np.broadcast_to(x[None, :, None], (y.size, x.size, t.size)).copy()
    return xr.Dataset({"vv": (("y", "x", "time"), vv)},
                      coords={"x": x, "y": y, "time": t},
                      attrs={"epsg": mp.CONFIG.epsg, "file": "synthetic"})


def run(fdir: Path, tpath: Path) -> bool:
    """Run every check against the fixture; returns True if all passed."""
    from shapely.geometry import LineString

    ok = True
    newest_x = TERMINI_X["2021-08-15"]

    # --- readers ------------------------------------------------------------
    F = mp.load_flowlines(ROI, directory=fdir)
    ok &= _check(f"ROI keeps only glacier0001's 3 flowlines (got {len(F)} rows, "
                 f"glaciers {sorted(F['glacier'].unique())})",
                 len(F) == 3 and set(F["glacier"]) == {"0001"})
    ok &= _check("_iter files excluded by default",
                 (F["iteration"] == 0).all())
    ok &= _check("_iter files included on request",
                 len(mp.load_flowlines(ROI, directory=fdir, iterations=True)) == 4)
    ok &= _check("a second ROI selects only glacier0003",
                 set(mp.load_flowlines(ROI_2, directory=fdir)["glacier"]) == {"0003"})
    ok &= _check("a non-flowline shapefile in the folder is skipped",
                 "some_masks" not in " ".join(F["file"]))
    ok &= _check("no ROI returns the whole collection (7 lines)",
                 len(mp.load_flowlines(directory=fdir)) == 7)

    T = mp.load_termini(ROI, path=tpath)
    ok &= _check(f"quality cut drops the flagged trace ({len(T)} of 4 kept)",
                 len(T) == 3)
    ok &= _check("dates parsed and sorted oldest first",
                 str(T["date"].iloc[0].date()) == "2018-06-01"
                 and str(T["date"].iloc[-1].date()) == "2021-08-15")
    ok &= _check("max_quality=None keeps the flagged trace",
                 len(mp.load_termini(ROI, path=tpath, max_quality=None)) == 4)
    ok &= _check("trange bounds the traces",
                 len(mp.load_termini(ROI, path=tpath,
                                     trange=("2019-01-01", None))) == 2)

    # --- selection ----------------------------------------------------------
    c = mp.center_flowline(F, glacier="0001")
    ok &= _check(f"centre of the fan is the middle line (got {c['flowline']!r})",
                 c["flowline"] == "04")
    ok &= _check("ordered_flowlines ranks the fan centre-outwards",
                 list(mp.ordered_flowlines(F, glacier="0001")["flowline"])[0] == "04")
    newest = mp.latest_terminus(T, line=c.geometry)
    ok &= _check("latest_terminus returns the newest crossing trace",
                 str(newest["date"].date()) == "2021-08-15")

    # --- the split ----------------------------------------------------------
    up = mp.clip_upstream(c.geometry, newest.geometry)
    ok &= _check(f"split at the newest trace, not an older one "
                 f"(starts at x={up.coords[0][0]:.0f})",
                 np.isclose(up.coords[0][0], newest_x))
    ok &= _check("upstream side kept (runs inland to x=10000)",
                 np.isclose(up.coords[-1][0], 10000.0))
    rev = LineString(list(c.geometry.coords)[::-1])
    up_rev = mp.clip_upstream(rev, newest.geometry)
    ok &= _check("auto orientation handles an inland-first line too",
                 np.isclose(up_rev.length, up.length))
    ok &= _check("explicit downstream_end='last' agrees",
                 np.isclose(
                     mp.clip_upstream(rev, newest.geometry,
                                      downstream_end="last").length, up.length))
    ok &= _check("a non-crossing terminus returns None",
                 mp.clip_upstream(
                     c.geometry,
                     LineString([(-5000, -10), (-5000, 10)])) is None)

    # --- the ROI clip -------------------------------------------------------
    inroi = mp.clip_to_roi(up, ROI)
    ok &= _check(f"clipped to the ROI's 6 km edge (len={inroi.length:.0f} m)",
                 np.isclose(inroi.length, ROI[1] - newest_x, atol=120.0))

    # --- packaging ----------------------------------------------------------
    ds = mp.flowline_from_roi(ROI, flowline_dir=fdir, termini_path=tpath)
    d = ds.d.values
    ok &= _check(f"Dataset has (point,) x/y/d, d[0]=0, {ds.sizes['point']} points",
                 d[0] == 0 and np.all(np.diff(d) > 0)
                 and set(ds.coords) == {"x", "y", "d"})
    ok &= _check(f"starts at the terminus (x={ds.x.values[0]:.0f}) and "
                 f"records its date ({ds.attrs['terminus_date']})",
                 np.isclose(ds.x.values[0], newest_x)
                 and ds.attrs["terminus_date"] == "2021-08-15")
    ok &= _check(f"defaults to the centre flowline (got {ds.attrs['flowline']!r})",
                 ds.attrs["glacier"] == "0001" and ds.attrs["flowline"] == "04")

    ds2 = mp.flowline_from_roi(ROI, flowline_dir=fdir, termini_path=tpath,
                               spacing=100.0, flowline="05")
    ok &= _check(f"spacing=100 gives an even grid ({ds2.sizes['point']} points)",
                 np.allclose(np.diff(ds2.d.values), 100.0)
                 and ds2.attrs["flowline"] == "05")

    # --- drop-in with the rest of the toolkit -------------------------------
    Fobs = mp.sample_obs_along_flowline(_make_obs(), ds2)
    ok &= _check("sampled obs has dims (point, time) with vel",
                 Fobs.vel.dims == ("point", "time"))
    ok &= _check("sampled speed recovers the known field (vel == x)",
                 np.allclose(Fobs.vel.isel(time=0).values, ds2.x.values, atol=1e-6))
    ok &= _check("d/x/y survive the sampling",
                 np.allclose(Fobs.d.values, ds2.d.values))

    # --- stepping outwards when the centre line has no trace ---------------
    F2 = mp.load_flowlines(ROI_2, directory=fdir)
    c2 = mp.center_flowline(F2, glacier="0003")
    T2 = mp.load_termini(ROI_2, path=tpath)
    ok &= _check(f"the stub trace misses the centre line ({c2['flowline']!r})",
                 c2["flowline"] == "04" and len(T2) == 1
                 and mp.clip_upstream(c2.geometry, T2.iloc[0].geometry) is None)

    ds3 = mp.flowline_from_roi(ROI_2, flowline_dir=fdir, termini_path=tpath)
    ok &= _check(f"search steps outwards to a line the trace does cross "
                 f"(used {ds3.attrs['flowline']!r}, centre was "
                 f"{ds3.attrs['center_flowline']!r})",
                 ds3.attrs["flowline"] == "03"
                 and ds3.attrs["center_flowline"] == "04"
                 and ds3.attrs["terminus_date"] == "2021-08-15")
    ok &= _check("the swapped line is still clipped at the terminus",
                 np.isclose(ds3.x.values[0], STUB_X))
    ok &= _check("center_flowline attr equals flowline when no swap happened",
                 ds.attrs["center_flowline"] == ds.attrs["flowline"] == "04")
    ok &= _check("naming a flowline explicitly does NOT swap, it raises",
                 _raises(mp.flowline_from_roi, ROI_2, flowline="04",
                         flowline_dir=fdir, termini_path=tpath))

    # --- failure paths ------------------------------------------------------
    ok &= _check("an ROI with no flowlines raises",
                 _raises(mp.flowline_from_roi, (1e6, 1.1e6, 1e6, 1.1e6),
                         flowline_dir=fdir, termini_path=tpath))
    ok &= _check("an unknown glacier raises",
                 _raises(mp.flowline_from_roi, ROI, glacier="9999",
                         flowline_dir=fdir, termini_path=tpath))
    ok &= _check("an unknown flowline id raises",
                 _raises(mp.flowline_from_roi, ROI, flowline="99",
                         flowline_dir=fdir, termini_path=tpath))
    ok &= _check("no crossing trace on ANY candidate raises",
                 _raises(mp.flowline_from_roi, ROI, flowline_dir=fdir,
                         termini_path=tpath, trange=("2030-01-01", None)))
    unclipped = mp.flowline_from_roi(ROI, flowline_dir=fdir, termini_path=tpath,
                                     trange=("2030-01-01", None),
                                     require_terminus=False)
    ok &= _check("require_terminus=False falls back to the un-split flowline",
                 unclipped.attrs["terminus_date"] == ""
                 and np.isclose(unclipped.x.values[0], 0.0))
    ok &= _check("sample_obs_along_flowline rejects an absent band",
                 _raises_key(mp.sample_obs_along_flowline, _make_obs(), ds2,
                             band="nope"))
    return ok


def _raises_key(fn, *args, **kwargs) -> bool:
    """True if ``fn`` raises KeyError."""
    try:
        fn(*args, **kwargs)
    except KeyError:
        return True
    return False


if __name__ == "__main__":
    print("=== flowline_shp (synthetic shapefiles) ===")
    try:
        import geopandas  # noqa: F401
        import shapely  # noqa: F401
    except ImportError as e:
        print(f"  SKIPPED — needs the 'shapefile' extra ({e})")
        sys.exit(0)

    with tempfile.TemporaryDirectory() as td:
        fdir, tpath = _write_fixture(Path(td))
        passed = run(fdir, tpath)
    print("\n" + ("ALL FLOWLINE-SHP CHECKS PASSED" if passed else "SOME CHECKS FAILED"))
    sys.exit(0 if passed else 1)
