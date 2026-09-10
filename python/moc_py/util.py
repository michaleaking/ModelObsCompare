"""Time-conversion helpers shared by the moc_py readers.

Observed flowline ``.mat`` files store time as MATLAB *datenum*; observed
netCDF files store CF ``"<step> since <ref>"`` times. Both are converted to a
common ``numpy.datetime64`` representation and to decimal years (to align with
ISSM model output, which uses decimal years).
"""

from __future__ import annotations

import re

import numpy as np

# MATLAB datenum of the Unix epoch (1970-01-01). datenum(1970,1,1) == 719529.
_MATLAB_UNIX_EPOCH_DAYS = 719529


def datenum_to_datetime64(dn: np.ndarray | float) -> np.ndarray:
    """Convert MATLAB datenum(s) to ``datetime64[ns]``.

    Args:
        dn: MATLAB serial date number(s).

    Returns:
        ``numpy.datetime64[ns]`` array (scalar in -> 0-d array).
    """
    dn = np.asarray(dn, dtype="float64")
    seconds = np.round((dn - _MATLAB_UNIX_EPOCH_DAYS) * 86400.0)
    return np.datetime64("1970-01-01T00:00:00") + seconds.astype("timedelta64[s]")


def datetime64_to_decyear(dt: np.ndarray) -> np.ndarray:
    """Convert ``datetime64`` value(s) to decimal years (e.g. 2018.5).

    Args:
        dt: ``numpy.datetime64`` array (any precision).

    Returns:
        ``float64`` array of decimal years.
    """
    dt = np.asarray(dt, dtype="datetime64[ns]")
    year_start = dt.astype("datetime64[Y]")  # floor to Jan 1 of each year
    year = year_start.astype(int) + 1970
    next_start = (year_start + np.timedelta64(1, "Y")).astype("datetime64[ns]")
    ys = year_start.astype("datetime64[ns]")
    frac = (dt - ys).astype("timedelta64[s]").astype("float64") / (
        (next_start - ys).astype("timedelta64[s]").astype("float64")
    )
    return year.astype("float64") + frac


def datenum_to_decyear(dn: np.ndarray | float) -> np.ndarray:
    """Convert MATLAB datenum(s) directly to decimal years."""
    return datetime64_to_decyear(datenum_to_datetime64(dn))


def cf_time_to_datetime64(values: np.ndarray, units: str | None) -> np.ndarray:
    """Convert CF ``"<step> since <ref>"`` times to ``datetime64[ns]``.

    Args:
        values: Numeric time offsets.
        units: CF units string, e.g. ``"days since 2015-01-06 12:00:00"``.
            If ``None`` or unparseable, values are returned as days since Unix
            epoch (best-effort).

    Returns:
        ``numpy.datetime64[ns]`` array.
    """
    values = np.asarray(values, dtype="float64")
    if not units:
        base = np.datetime64("1970-01-01T00:00:00")
        step = "days"
    else:
        m = re.match(r"\s*(\w+)\s+since\s+(.+)", units)
        if not m:
            base = np.datetime64("1970-01-01T00:00:00")
            step = "days"
        else:
            step = m.group(1).lower()
            ref = m.group(2).strip().replace("T", " ")
            ref = re.sub(r"\s*(UTC|Z)\s*$", "", ref)
            try:
                base = np.datetime64(ref.replace(" ", "T"))
            except ValueError:
                # fall back to date-only
                base = np.datetime64(ref.split()[0])
    per = {
        "days": 86400.0, "day": 86400.0,
        "hours": 3600.0, "hour": 3600.0,
        "minutes": 60.0, "minute": 60.0,
        "seconds": 1.0, "second": 1.0,
    }.get(step, 86400.0)
    seconds = np.round(values * per)
    return base + seconds.astype("timedelta64[s]")


def decyear_to_datetime64(yr: np.ndarray | float) -> np.ndarray:
    """Convert decimal year(s) to ``datetime64[ns]`` (inverse of decyear)."""
    yr = np.asarray(yr, dtype="float64")
    year = np.floor(yr).astype(int)
    frac = yr - year
    start = np.array([np.datetime64(f"{y:04d}-01-01") for y in year], dtype="datetime64[ns]")
    nxt = np.array([np.datetime64(f"{y + 1:04d}-01-01") for y in year], dtype="datetime64[ns]")
    span = (nxt - start).astype("timedelta64[s]").astype("float64")
    return start + np.round(frac * span).astype("timedelta64[s]")
