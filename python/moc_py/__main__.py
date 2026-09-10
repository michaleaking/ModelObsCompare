"""Command-line entry point: show or set the moc_py data paths.

    python -m moc_py            # show every setting and where it came from
    python -m moc_py --set      # prompt for the ones still missing

Installing the package also provides the same thing as ``moc-config``.
"""

from __future__ import annotations

from .config import main

if __name__ == "__main__":
    raise SystemExit(main())
