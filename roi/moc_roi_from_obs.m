function roi = moc_roi_from_obs(O, varargin)
%MOC_ROI_FROM_OBS  Build a rectangular ROI from the extent of an observed grid.
%
%   Convenient when a glacier-specific netCDF (e.g. Sentinel_Subset_Upernavik.nc)
%   already defines the region you care about: use its bounding box directly,
%   optionally shrunk by a margin to avoid edge/no-data pixels.
%
%   Usage:
%      roi = moc_roi_from_obs(O)                          % O from moc_load_obs_netcdf
%      roi = moc_roi_from_obs(O, 'margin',2000, 'name','upernavik')
%
%   Input:
%      O : obs grid struct with fields .x (1xnx), .y (1xny), and .epsg
%   Name/value options:
%      'margin' : (double) inset from each edge [m] (default 0)
%      'name'   : (char)   ROI name; if given, the ROI is saved
%
%   Output:
%      roi : rectangular ROI struct (see moc_select_roi)

p = inputParser;
p.addParameter('margin',0, @isscalar);
p.addParameter('name','', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;
m = opt.margin;

xmin = min(O.x)+m; xmax = max(O.x)-m;
ymin = min(O.y)+m; ymax = max(O.y)-m;

roi = struct();
roi.name = char(opt.name);
roi.type = 'rectangle';
roi.x = [xmin xmax xmax xmin xmin]';
roi.y = [ymin ymin ymax ymax ymin]';
roi.bbox = [xmin xmax ymin ymax];
if isfield(O,'epsg'), roi.epsg = O.epsg; else, cfg = config_moc(); roi.epsg = cfg.epsg; end
roi.file = '';

if ~isempty(roi.name)
    roi = moc_save_roi(roi);
end

end
