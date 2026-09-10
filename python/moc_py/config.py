"""Central configuration (paths and defaults) for the moc_py toolkit.

Nothing here is hard-coded to one machine. Every data path is resolved on
first use, in this order:

1. an explicit argument — the ``directory=`` / ``path=`` keyword most readers
   take, or :func:`configure` to set one globally;
2. an environment variable (``MOC_MODEL_DIR``, ``MOC_OBS_NETCDF_DIR``, ...);
3. the shared config file, ``~/.config/moc/paths.json`` (override its location
   with ``MOC_CONFIG``) — the MATLAB side reads the same file, so configuring
   once covers both;
4. an interactive prompt, when running in a terminal — the answer is saved to
   the config file, so each path is asked for once.

If none of those produce a path (a script, a notebook without a console, CI),
a :class:`MissingPathError` explains exactly which setting is missing and the
three ways to supply it.

    python -m moc_py                 # show every setting and where it came from
    python -m moc_py --set           # prompt for the ones still missing
    moc-config --set                 # same, once the package is installed

Mirrors ``config_moc.m`` on the MATLAB side.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

#: Where a resolved path is remembered across sessions (shared with MATLAB).
#: ``MOC_CONFIG`` overrides it.
DEFAULT_CONFIG_FILE = (
    Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "moc" / "paths.json"
)

#: The configurable data locations: setting -> (kind, env var, description).
PATH_SETTINGS: dict[str, tuple[str, str, str]] = {
    "model_dir": (
        "dir", "MOC_MODEL_DIR",
        "folder of ISSM model output .mat / exported .nc files",
    ),
    "obs_flowline_dir": (
        "dir", "MOC_OBS_FLOWLINE_DIR",
        "folder of per-glacier observed flowline .mat files (e.g. redrawn_Jan25)",
    ),
    "obs_netcdf_dir": (
        "dir", "MOC_OBS_NETCDF_DIR",
        "folder of gridded observed-velocity netCDF files",
    ),
    "flowline_shp_dir": (
        "dir", "MOC_FLOWLINE_SHP_DIR",
        "folder of Felikson per-glacier flowline shapefiles (glacier*.shp)",
    ),
    "termini_shp": (
        "file", "MOC_TERMINI_SHP",
        "Black & Joughin traced-terminus shapefile (glacier_termini_v01.0.shp)",
    ),
    "roi_dir": (
        "dir", "MOC_ROI_DIR",
        "folder to save selected ROIs in (created if missing)",
    ),
}


class MissingPathError(RuntimeError):
    """Raised when a data path is needed but has not been configured."""


def config_file() -> Path:
    """Return the config file location (``MOC_CONFIG``, else the default)."""
    override = os.environ.get("MOC_CONFIG")
    return Path(override).expanduser() if override else DEFAULT_CONFIG_FILE


def _read_config_file() -> dict:
    """Load the saved paths, or an empty dict if there is no readable file."""
    p = config_file()
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return {}


def _clean(text: str) -> str:
    """Tidy a pasted or drag-and-dropped path (quotes, whitespace, ``~``)."""
    return str(Path(text.strip().strip('"').strip("'")).expanduser())


def _is_interactive() -> bool:
    """True when there is a console a prompt could actually reach."""
    try:
        return sys.stdin is not None and sys.stdin.isatty()
    except (AttributeError, ValueError):  # detached stdin
        return False


class Config:
    """Toolkit configuration: data locations plus projection/grid defaults.

    Path attributes (``model_dir``, ``obs_flowline_dir``, ``obs_netcdf_dir``,
    ``flowline_shp_dir``, ``termini_shp``, ``roi_dir``) resolve on first access
    through the chain described in the module docstring, and are cached
    afterwards. The non-path settings are plain values:

    Attributes:
        epsg: Projection EPSG code the data are stored in (3413 = Greenland PS).
        grid_res: Default regular-grid spacing in metres.
        obs_bands: Component order in the Sentinel-type netCDF ``velocity`` var.
        prompt: Whether an unset path may be asked for interactively.

    Example:
        >>> import moc_py as mp
        >>> mp.configure(obs_netcdf_dir="~/data/velocity")   # or set MOC_OBS_NETCDF_DIR
        >>> mp.CONFIG.obs_netcdf_dir
        PosixPath('/home/you/data/velocity')
    """

    def __init__(self, *, epsg: int = 3413, grid_res: float = 200.0,
                 obs_bands: tuple[str, ...] = ("vv", "vx", "vy", "ex", "ey", "dT"),
                 prompt: bool = True, **paths):
        self.epsg = epsg
        self.grid_res = grid_res
        self.obs_bands = obs_bands
        self.prompt = prompt
        self._explicit: dict[str, Path] = {}
        self._cache: dict[str, Path] = {}
        if paths:
            self.set(**paths)

    # -- resolution ---------------------------------------------------------
    def __getattr__(self, name: str) -> Path:
        """Resolve a path setting on first access (called only if not an attr)."""
        if name in PATH_SETTINGS:
            return self.get(name)
        raise AttributeError(f"{type(self).__name__!r} has no attribute {name!r}")

    def get(self, name: str, prompt: bool | None = None) -> Path:
        """Resolve one path setting.

        Args:
            name: A key of :data:`PATH_SETTINGS`.
            prompt: Override whether an interactive prompt is allowed.

        Returns:
            The configured :class:`~pathlib.Path`.

        Raises:
            KeyError: If ``name`` is not a known setting.
            MissingPathError: If it is not configured and cannot be asked for.
        """
        if name not in PATH_SETTINGS:
            raise KeyError(f"Unknown path setting {name!r}; "
                           f"known: {', '.join(PATH_SETTINGS)}")
        if name in self._cache:
            return self._cache[name]

        kind, env, desc = PATH_SETTINGS[name]
        value = self._explicit.get(name)
        if value is None and os.environ.get(env):
            value = Path(_clean(os.environ[env]))
        if value is None:
            saved = _read_config_file().get(name)
            if saved:
                value = Path(_clean(str(saved)))
        if value is None:
            value = self._default(name)

        if value is None:
            allow = self.prompt if prompt is None else prompt
            if allow and _is_interactive():
                value = self._ask(name, kind, env, desc)
        if value is None:
            raise MissingPathError(
                f"The path setting '{name}' ({desc}) is not configured.\n"
                f"Set it in any of these ways:\n"
                f"  export {env}=/path/to/data\n"
                f"  python -c \"import moc_py; moc_py.configure({name}='/path/to/data')\"\n"
                f"  python -m moc_py --set            (prompts, then remembers)"
            )

        if name == "roi_dir":
            value.mkdir(parents=True, exist_ok=True)
        self._cache[name] = value
        return value

    def _default(self, name: str):
        """Built-in fallback for settings that ship with the repo."""
        if name == "roi_dir":
            return Path(__file__).resolve().parent.parent / "rois"
        return None

    def _ask(self, name: str, kind: str, env: str, desc: str):
        """Prompt for one path, validate it, and remember the answer."""
        what = "folder" if kind == "dir" else "file"
        print(f"\nmoc_py needs the {what} for '{name}':\n    {desc}", file=sys.stderr)
        print(f"(leave blank to cancel; set {env} to skip this next time)",
              file=sys.stderr)
        for _ in range(3):
            try:
                reply = input(f"  path to {name}: ")
            except (EOFError, KeyboardInterrupt):
                print(file=sys.stderr)
                return None
            if not reply.strip():
                return None
            p = Path(_clean(reply))
            if (kind == "dir" and p.is_dir()) or (kind == "file" and p.is_file()):
                self._explicit[name] = p
                try:
                    save_config(self)
                    print(f"  saved to {config_file()}", file=sys.stderr)
                except OSError as e:                      # read-only home, etc.
                    print(f"  (could not save: {e})", file=sys.stderr)
                return p
            print(f"  no such {what}: {p}", file=sys.stderr)
        return None

    # -- mutation / inspection ---------------------------------------------
    def set(self, **paths) -> "Config":
        """Set one or more path settings explicitly (highest precedence).

        Args:
            **paths: ``setting=path`` pairs; a value of ``None`` clears it.

        Returns:
            This config, so calls can be chained.

        Raises:
            KeyError: If a name is not a known path setting.
        """
        for name, value in paths.items():
            if name not in PATH_SETTINGS:
                raise KeyError(f"Unknown path setting {name!r}; "
                               f"known: {', '.join(PATH_SETTINGS)}")
            self._cache.pop(name, None)
            if value is None:
                self._explicit.pop(name, None)
            else:
                self._explicit[name] = Path(_clean(str(value)))
        return self

    def status(self) -> dict[str, tuple[str, str]]:
        """Report each path setting as ``(value, source)`` without prompting.

        Returns:
            Mapping of setting name to its resolved value (``""`` if unset) and
            where it came from (``"argument"``, the env var name, ``"config
            file"``, ``"built-in"``, or ``"unset"``).
        """
        saved = _read_config_file()
        out: dict[str, tuple[str, str]] = {}
        for name, (_, env, _desc) in PATH_SETTINGS.items():
            if name in self._explicit:
                out[name] = (str(self._explicit[name]), "argument")
            elif os.environ.get(env):
                out[name] = (_clean(os.environ[env]), env)
            elif saved.get(name):
                out[name] = (_clean(str(saved[name])), "config file")
            elif self._default(name) is not None:
                out[name] = (str(self._default(name)), "built-in")
            else:
                out[name] = ("", "unset")
        return out


def save_config(cfg: "Config | None" = None, path: str | Path | None = None) -> Path:
    """Write the currently set paths to the shared config file.

    Existing entries for settings that are not set are left untouched, so
    saving never silently drops another session's answer.

    Args:
        cfg: Config to take values from (default the module-level ``CONFIG``).
        path: Destination (default :func:`config_file`).

    Returns:
        The file written.
    """
    cfg = cfg if cfg is not None else CONFIG
    dest = Path(path).expanduser() if path else config_file()
    data = _read_config_file() if dest == config_file() else {}
    for name, (value, source) in cfg.status().items():
        if value and source in ("argument", "config file"):
            data[name] = value
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return dest


def configure(**paths) -> "Config":
    """Set data paths on the module-level :data:`CONFIG`.

    Args:
        **paths: ``setting=path`` pairs, e.g.
            ``configure(obs_netcdf_dir="~/data/vel", termini_shp="~/data/t.shp")``.

    Returns:
        The updated :data:`CONFIG`.

    Example:
        >>> import moc_py as mp
        >>> mp.configure(flowline_shp_dir="~/supp_data/Felikson_flowlines")
    """
    return CONFIG.set(**paths)


#: Module-level default configuration. Import and mutate, or build your own.
CONFIG = Config()


def resolve(path_like: str | Path, setting: str) -> Path:
    """Resolve a bare filename against a configured directory.

    A path that already points at a file is returned untouched — the directory
    setting is only consulted (and only then possibly prompted for) when the
    lookup actually needs it.

    Args:
        path_like: A full path or a bare filename.
        setting: Name of the directory setting to search, e.g.
            ``"obs_netcdf_dir"``.

    Returns:
        An existing :class:`~pathlib.Path`.

    Raises:
        FileNotFoundError: If neither the path nor ``<setting>/path_like`` exists.
        MissingPathError: If the directory setting is needed but unconfigured.
    """
    p = Path(path_like).expanduser()
    if p.is_file():
        return p
    cand = CONFIG.get(setting) / p
    if cand.is_file():
        return cand
    raise FileNotFoundError(f"Not found: {path_like} (also tried {cand})")


def main(argv: list[str] | None = None) -> int:
    """Show the configured paths, or prompt for the missing ones.

    Args:
        argv: Command-line arguments (default ``sys.argv[1:]``).

    Returns:
        Process exit code (0 if every setting is configured).
    """
    import argparse

    ap = argparse.ArgumentParser(
        prog="python -m moc_py",
        description="Show or set the moc_py data paths.")
    ap.add_argument("--set", action="store_true",
                    help="prompt for every setting that is still missing")
    ap.add_argument("--all", action="store_true",
                    help="with --set, re-prompt for settings that are already set")
    args = ap.parse_args(argv)

    if args.set:
        for name, (_, source) in list(CONFIG.status().items()):
            if source == "unset" or (args.all and source != "built-in"):
                kind, env, desc = PATH_SETTINGS[name]
                CONFIG._ask(name, kind, env, desc)

    print(f"\nmoc_py paths   (config file: {config_file()})")
    print(f"{'setting':<18} {'source':<22} value")
    print(f"{'-' * 18} {'-' * 22} {'-' * 40}")
    missing = 0
    for name, (value, source) in CONFIG.status().items():
        print(f"{name:<18} {source:<22} {value or '(not set)'}")
        missing += source == "unset"
    if missing:
        print(f"\n{missing} setting(s) not configured — "
              f"run with --set, or export the MOC_* variables.")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
