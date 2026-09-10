function S = moc_path_settings()
%MOC_PATH_SETTINGS  The configurable data locations (single source of truth).
%
%   Every path the toolkit needs is described here once: its kind, the
%   environment variable that can supply it, and a human-readable description
%   used by the interactive prompt. The Python side keeps the same list in
%   moc_py/config.py and reads/writes the same config file, so configuring
%   either language configures both.
%
%   Usage:
%      S = moc_path_settings();
%
%   Output (1 x n struct array S):
%      .name (char) setting name, e.g. 'model_dir'
%      .kind (char) 'dir' | 'file'
%      .env  (char) environment variable that supplies it
%      .desc (char) what the user should point it at

S = struct('name',{},'kind',{},'env',{},'desc',{});

S(end+1) = struct('name','issm_dir', 'kind','dir', 'env','ISSM_DIR', ...
    'desc','root of the ISSM binary distribution (the folder holding bin/ and lib/)');
S(end+1) = struct('name','model_dir', 'kind','dir', 'env','MOC_MODEL_DIR', ...
    'desc','folder of ISSM model output .mat / exported .nc files');
S(end+1) = struct('name','obs_flowline_dir', 'kind','dir', 'env','MOC_OBS_FLOWLINE_DIR', ...
    'desc','folder of per-glacier observed flowline .mat files (e.g. redrawn_Jan25)');
S(end+1) = struct('name','obs_netcdf_dir', 'kind','dir', 'env','MOC_OBS_NETCDF_DIR', ...
    'desc','folder of gridded observed-velocity netCDF files');
S(end+1) = struct('name','flowline_shp_dir', 'kind','dir', 'env','MOC_FLOWLINE_SHP_DIR', ...
    'desc','folder of Felikson per-glacier flowline shapefiles (glacier*.shp)');
S(end+1) = struct('name','termini_shp', 'kind','file', 'env','MOC_TERMINI_SHP', ...
    'desc','Black & Joughin traced-terminus shapefile (glacier_termini_v01.0.shp)');
S(end+1) = struct('name','roi_dir', 'kind','dir', 'env','MOC_ROI_DIR', ...
    'desc','folder to save selected ROIs in (created if missing)');

end
