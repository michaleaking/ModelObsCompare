"""Self-contained validation of the Phase-2 model path (no MATLAB needed).

Fabricates a triangulated mesh with a known LINEAR velocity field, writes it in
the exported-model netCDF schema, then checks that:

  * load_model returns 0-based elements and a working triangulation,
  * the mesh interpolation reproduces the linear field exactly (P1 interp is
    exact for linear fields) at points, on a grid, and along a flowline,
  * spatial_diff of the model against an observation equal to the same field is
    ~0 (verifies obs regridding + differencing sign),
  * the loader is robust to either stored axis order (the MATLAB/HDF5 transpose
    quirk) for both `elements` and the field variables.

Run:  python tests/test_model_synthetic.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import xarray as xr
from scipy.spatial import Delaunay

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import moc_py as mp  # noqa: E402

# A known linear field:  z(x, y, t) = A*x + B*y + C + t_offset
A, B, C = 1.0e-3, 2.0e-3, 500.0
TIMES = np.array([2018.0, 2019.0])


def linear_field(x, y, t):
    return A * x + B * y + C + (t - TIMES[0]) * 100.0


def _make_mesh(n=1500, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-316000, -267000, n)
    y = rng.uniform(-1860000, -1809000, n)
    tri = Delaunay(np.column_stack([x, y]))
    return x, y, tri.simplices  # elements 0-based (ne, 3)


def _write_model_nc(path, x, y, elements0, transpose=False):
    """Write the exported-model schema with 1-based elements (like MATLAB)."""
    nv, nt = x.size, TIMES.size
    vel = np.stack([linear_field(x, y, t) for t in TIMES], axis=1)  # (nv, nt)
    elem1 = (elements0 + 1).astype("int32")  # 1-based, like MATLAB export
    with h5py.File(path, "w") as h:
        h.create_dataset("x", data=x)
        h.create_dataset("y", data=y)
        h.create_dataset("time", data=TIMES)
        e = h.create_dataset("elements", data=(elem1.T if transpose else elem1))
        e.attrs["index_base"] = np.int32(1)
        for name in ("vel", "vx", "vy"):
            h.create_dataset(name, data=(vel.T if transpose else vel))
        h.attrs["epsg"] = np.int32(3413)
        h.attrs["index_base"] = np.int32(1)


def _check(msg, cond):
    print(f"  [{'PASS' if cond else 'FAIL':4s}] {msg}")
    return bool(cond)


def run(transpose: bool) -> bool:
    ok = True
    x, y, elements0 = _make_mesh()
    with tempfile.TemporaryDirectory() as td:
        nc = Path(td) / "synthetic_model.nc"
        _write_model_nc(nc, x, y, elements0, transpose=transpose)
        M = mp.load_model(nc)

        ok &= _check("elements 0-based (min index == 0)", M.elements.min() == 0)
        ok &= _check("elements shape (ne, 3)", M.elements.shape[1] == 3)
        ok &= _check("field shape (nv, nt)", M.field("vel").shape == (x.size, TIMES.size))

        # interior sample points (away from the convex hull to stay inside)
        xs = np.array([-3.0e5, -2.9e5, -2.85e5])
        ys = np.array([-1.85e6, -1.84e6, -1.83e6])
        got = mp.model_at_points(M, xs, ys, field="vel")  # (npts, nt)
        want = np.stack([linear_field(xs, ys, t) for t in TIMES], axis=1)
        finite = np.isfinite(got)
        err = np.nanmax(np.abs(got[finite] - want[finite]))
        ok &= _check(f"model_at_points exact on linear field (max err {err:.2e})", err < 1e-6)

        # grid
        roi = mp.Roi("t", "rectangle",
                     np.array([-3.05e5, -2.75e5, -2.75e5, -3.05e5, -3.05e5]),
                     np.array([-1.855e6, -1.855e6, -1.815e6, -1.815e6, -1.855e6]))
        G = mp.model_to_grid(M, roi, res=500, field="vel", time=2019.0)
        X, Y = np.meshgrid(G.x.values, G.y.values)
        gwant = linear_field(X, Y, 2019.0)
        v = np.isfinite(G["z"].values)
        gerr = np.nanmax(np.abs(G["z"].values[v] - gwant[v]))
        ok &= _check(f"model_to_grid exact inside ROI (max err {gerr:.2e}, {v.mean()*100:.0f}% valid)",
                     gerr < 1e-6 and v.mean() > 0.5)

        # spatial_diff vs an observation equal to the same field  -> ~0
        ox = np.linspace(-3.05e5, -2.75e5, 80)
        oy = np.linspace(-1.855e6, -1.815e6, 90)
        OX, OY = np.meshgrid(ox, oy)
        vv = linear_field(OX, OY, 2019.0)[:, :, None]  # (y, x, 1)
        O = xr.Dataset({"vv": (("y", "x", "time"), vv)},
                       coords={"x": ox, "y": oy, "time": [2019.0]},
                       attrs={"epsg": 3413})
        D = mp.spatial_diff(M, O, roi, res=500, time=2019.0)
        ok &= _check(f"spatial_diff(model vs same field) ~ 0 (RMSE {D.attrs['rmse']:.2e})",
                     D.attrs["rmse"] < 1e-3)

        # flowline sampling matches the field
        fl = xr.Dataset(
            {"vel": (("point", "time"), np.zeros((len(xs), len(TIMES))))},
            coords={"x": ("point", xs), "y": ("point", ys),
                    "d": ("point", np.array([0.0, 1e4, 2e4])), "time": ("time", TIMES)},
            attrs={"name": "obs"})
        Fm = mp.model_flowline(M, fl, field="vel")
        ferr = np.nanmax(np.abs(Fm.vel.values - want))
        ok &= _check(f"model_flowline exact (max err {ferr:.2e})", ferr < 1e-6)
    return ok


if __name__ == "__main__":
    all_ok = True
    for tp in (False, True):
        print(f"\n=== stored axis order: {'TRANSPOSED' if tp else 'canonical'} ===")
        all_ok &= run(tp)
    print("\n" + ("ALL PHASE-2 CHECKS PASSED" if all_ok else "SOME CHECKS FAILED"))
    sys.exit(0 if all_ok else 1)
