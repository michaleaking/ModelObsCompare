"""Self-contained validation of :mod:`moc_py.obs_nsidc` (no NSIDC access needed).

Fabricates the netCDF that ``glacier_reader.save_velocity_netcdf`` writes — a
``velocity(time, band, y, x)`` stack in nisardev band order
``(vx, vy, vv, ex, ey, ev, dT)``, datetime64 times, GrIMP sentinels left in on
purpose — and asserts that:

  * :func:`moc_py.load_nsidc_netcdf` labels bands by the file's own ``band``
    coord, so ``vv`` is the speed and not ``vx`` (the whole point of the
    adapter: moc's ``CONFIG.obs_bands`` order is different),
  * the output matches the schema of :func:`moc_py.load_obs_netcdf` — dims
    ``(y, x, time)``, decimal-year ``time`` plus ``time_dt``, an ``epsg`` attr —
    and feeds the downstream ROI / variance / spatial tools unchanged,
  * GrIMP sentinels (``-2e9`` on ``vx``/``vy``/``dT``, ``-1`` on the
    non-negative bands) become NaN, and the derived ``vv``/``ev`` inherit their
    components' gaps rather than surviving as ``+2.83e9``,
  * band subsetting, a decimal-year ``trange``, and provenance attrs survive,
  * :func:`moc_py.load_obs_netcdf` now WARNS instead of silently guessing when
    a file carries no band labels at all.

Run:  python tests/test_nsidc_synthetic.py
"""

from __future__ import annotations

import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402
from moc_py.obs_nsidc import nsidc_to_obs  # noqa: E402

#: nisardev band order — deliberately NOT moc's CONFIG.obs_bands order.
BANDS = ("vx", "vy", "vv", "ex", "ey", "ev", "dT")
VX, VY = 300.0, 400.0            # so vv = hypot(vx, vy) = 500 exactly
TIMES = np.array(["2019-01-15", "2019-04-15", "2019-07-15", "2019-10-15"],
                 dtype="datetime64[ns]")


def _make_file(path: Path, nx: int = 12, ny: int = 10) -> Path:
    """Write a synthetic glacier_reader-style NSIDC subset, sentinels included.

    The bottom row of every frame is no-data, flagged the way the GrIMP COGs do
    it: ``-2e9`` on the signed/elapsed bands and ``-1`` on the non-negative
    ones. ``vv``/``ev`` are left as ``hypot`` of the raw sentinels, which is how
    an unmasked read really behaves (+2.83e9), so the inheritance logic is
    exercised rather than assumed.

    Args:
        path: Output ``.nc`` path.
        nx: Grid columns. Default 12.
        ny: Grid rows. Default 10.

    Returns:
        The file written.
    """
    nt = TIMES.size
    stack = np.zeros((nt, len(BANDS), ny, nx), dtype="float32")
    base = {"vx": VX, "vy": VY, "vv": np.hypot(VX, VY),
            "ex": 10.0, "ey": 10.0, "ev": np.hypot(10.0, 10.0), "dT": 6.0}
    for bi, band in enumerate(BANDS):
        stack[:, bi] = base[band]
        stack[:, bi, -1, :] = -2.0e9 if band in ("vx", "vy", "dT") else -1.0
    # vv/ev as nisardev would compute them from UNMASKED sentinels
    stack[:, BANDS.index("vv"), -1, :] = np.hypot(-2.0e9, -2.0e9)
    stack[:, BANDS.index("ev"), -1, :] = np.hypot(-1.0, -1.0)

    da = xr.DataArray(
        stack, dims=("time", "band", "y", "x"),
        coords={"time": TIMES, "band": list(BANDS),
                "y": np.linspace(-2.28e6, -2.26e6, ny),
                "x": np.linspace(-1.93e5, -1.73e5, nx)},
        name="velocity",
    )
    ds = da.to_dataset()
    ds.attrs.update({"epsg": 3413, "nsidc_product": "NSIDC-0731",
                     "grimp_box": "", "glacier": "Synthetic"})
    ds.to_netcdf(path)
    return path


def _make_unlabelled(path: Path, nx: int = 12, ny: int = 10) -> Path:
    """Write a band-stacked velocity file with NO band labels at all.

    Built with :mod:`h5py` rather than ``xarray.to_netcdf`` so the fixture is
    real HDF5 regardless of which netCDF engine happens to be installed —
    :func:`moc_py.load_obs_netcdf` reads via h5py and cannot open a netCDF3
    file.

    Args:
        path: Output path.
        nx: Grid columns. Default 12.
        ny: Grid rows. Default 10.

    Returns:
        The file written.
    """
    import h5py

    with h5py.File(path, "w") as h:
        h.create_dataset("x", data=np.linspace(-1.93e5, -1.73e5, nx))
        h.create_dataset("y", data=np.linspace(-2.28e6, -2.26e6, ny))
        t = h.create_dataset("time", data=np.array([14.0, 104.0, 195.0, 287.0]))
        t.attrs["units"] = "days since 2019-01-01"
        h.create_dataset(
            "velocity",
            data=np.zeros((TIMES.size, len(BANDS), ny, nx), dtype="float32"),
        )
    return path


def _check(label: str, ok: bool) -> bool:
    """Print a pass/fail line and return the result unchanged."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    return ok


def _raises(fn, *args, **kwargs) -> bool:
    """True if calling ``fn`` raises anything."""
    try:
        fn(*args, **kwargs)
    except Exception:  # noqa: BLE001 - any failure is the expected outcome
        return True
    return False


def run() -> bool:
    """Run every check and return True only if all of them passed."""
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        path = _make_file(Path(tmp) / "Synthetic_NSIDC-0731.nc")

        O = mp.load_nsidc_netcdf(path)
        ok &= _check("all 7 nisardev bands come through",
                     set(O.data_vars) == set(BANDS))
        ok &= _check("dims are (y, x, time)",
                     O["vv"].dims == ("y", "x", "time"))

        # The adapter's reason to exist: band labels read from the file, not guessed.
        good = O["vv"].isel(y=0, x=0, time=0)
        ok &= _check("vv is the speed, not vx (band labels honoured)",
                     np.isclose(float(good), np.hypot(VX, VY))
                     and np.isclose(float(O["vx"].isel(y=0, x=0, time=0)), VX))

        ok &= _check("time is decimal years with a datetime64 companion",
                     2019.0 < float(O.time[0]) < 2019.2
                     and np.issubdtype(O["time_dt"].dtype, np.datetime64))
        ok &= _check("epsg and provenance attrs survive",
                     O.attrs["epsg"] == 3413
                     and O.attrs["nsidc_product"] == "NSIDC-0731")

        # No-data: sentinels gone, derived bands inheriting their components.
        bottom = {b: O[b].isel(y=-1).values for b in BANDS}
        ok &= _check("GrIMP sentinels are NaN on every band",
                     all(np.all(np.isnan(v)) for v in bottom.values()))
        ok &= _check("derived vv/ev did not survive as +2.83e9",
                     float(np.nanmax(O["vv"].values)) < 1e4)
        ok &= _check("valid pixels are untouched",
                     np.isclose(float(np.nanmean(O["vv"].values)), np.hypot(VX, VY)))

        # Subsetting.
        one = mp.load_nsidc_netcdf(path, bands=["vv"])
        ok &= _check("bands= subsets", list(one.data_vars) == ["vv"])
        win = mp.load_nsidc_netcdf(path, bands=["vv"], trange=(2019.4, 2019.8))
        ok &= _check("trange= subsets in decimal years", win.sizes["time"] == 2)
        ok &= _check("an absent band raises",
                     _raises(mp.load_nsidc_netcdf, path, bands=["nope"]))

        # Downstream moc tools consume it unchanged.
        roi = mp.roi_from_obs(O, margin=0.0)
        V = mp.velocity_variance(O, roi=roi)
        ok &= _check("roi_from_obs + velocity_variance accept the Dataset",
                     V is not None and roi is not None)

        # In-memory conversion (what a live nisarVelSeries.subset goes through).
        with xr.open_dataset(path) as raw:
            mem = nsidc_to_obs(raw["velocity"], bands=["vv", "vx"])
        ok &= _check("nsidc_to_obs works on an in-memory stack",
                     set(mem.data_vars) == {"vv", "vx"}
                     and mem["vv"].dims == ("y", "x", "time")
                     and np.isclose(float(mem["vv"].isel(y=0, x=0, time=0)),
                                    np.hypot(VX, VY)))

        # And the eager reader no longer guesses in silence.
        unlabelled = _make_unlabelled(Path(tmp) / "unlabelled.nc")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                mp.load_obs_netcdf(unlabelled)
            except Exception:  # noqa: BLE001 - the warning is what is under test
                pass
        ok &= _check("load_obs_netcdf warns when band labels are missing",
                     any("component" in str(w.message) for w in caught))
    return ok


if __name__ == "__main__":
    print("=== obs_nsidc (synthetic) ===")
    passed = run()
    print("\n" + ("ALL NSIDC CHECKS PASSED" if passed else "SOME CHECKS FAILED"))
    sys.exit(0 if passed else 1)
