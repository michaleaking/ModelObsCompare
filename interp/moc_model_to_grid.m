function G = moc_model_to_grid(M, roi, varargin)
%MOC_MODEL_TO_GRID  Regrid an ISSM model field from its mesh onto a regular grid.
%
%   Wraps ISSM's InterpFromMeshToGrid to sample an unstructured model field on
%   a regular (x,y) grid covering an ROI, masking points outside the ROI polygon.
%   This produces the gridded model field needed for spatial-difference maps.
%
%   Usage:
%      G = moc_model_to_grid(M, roi)
%      G = moc_model_to_grid(M, roi, 'field','vel', 'res',200, 'tindex',[])
%      G = moc_model_to_grid(M, roi, 'time',2019.0)     % nearest transient step
%
%   Input:
%      M   : model struct from moc_load_model
%      roi : ROI struct (defines the bounding box + polygon mask)
%   Name/value options:
%      'field'  : (char)   'vel'(default) | 'vx' | 'vy' | 'surface'
%      'res'    : (double) grid spacing [m] (default config_moc().grid_res)
%      'tindex' : (double) column index into the time series (default: last)
%      'time'   : (double) decimal year; overrides tindex by nearest match
%      'mask'   : (logical) mask outside the ROI polygon to NaN (default true)
%
%   Output (struct G):
%      .x, .y   (1xnx / 1xny)   grid coordinates [m]
%      .X, .Y   (ny x nx)       meshgrid matrices
%      .grid    (ny x nx)       interpolated field (NaN outside mesh/ROI)
%      .field   (char)          field name
%      .time    (double)        decimal year of the slice used
%      .res     (double)        grid spacing used

p = inputParser;
p.addParameter('field','vel', @(s) ischar(s)||isstring(s));
p.addParameter('res',[], @(v) isempty(v)||isscalar(v));
p.addParameter('tindex',[], @(v) isempty(v)||isscalar(v));
p.addParameter('time',[], @(v) isempty(v)||isscalar(v));
p.addParameter('mask',true, @islogical);
p.parse(varargin{:});
opt = p.Results;
field = char(opt.field);

if ~isfield(M, field)
    error('moc_model_to_grid:badField', 'Model struct has no field "%s".', field);
end
res = opt.res;
if isempty(res), cfg = config_moc(); res = cfg.grid_res; end

% Pick the time slice
nt = size(M.(field), 2);
if ~isempty(opt.time) && ~all(isnan(M.time))
    [~, ti] = min(abs(M.time - opt.time));
elseif ~isempty(opt.tindex)
    ti = opt.tindex;
else
    ti = nt;   % default to last (most recent) step
end
ti = max(1, min(ti, nt));
data = M.(field)(:, ti);

% Build the regular grid over the ROI bounding box
bb = roi.bbox;
xg = bb(1):res:bb(2);
yg = bb(3):res:bb(4);

% InterpFromMeshToGrid returns a [numel(yg) x numel(xg)] grid (rows=y, cols=x)
grid = InterpFromMeshToGrid(M.elements, M.x, M.y, data, xg, yg, NaN);

[X, Y] = meshgrid(xg, yg);

% Mask outside the ROI polygon
if opt.mask
    inside = moc_roi_mask(roi, X, Y);
    grid(~inside) = NaN;
end

G = struct();
G.x = xg; G.y = yg;
G.X = X;  G.Y = Y;
G.grid = grid;
G.field = field;
if all(isnan(M.time)), G.time = NaN; else, G.time = M.time(ti); end
G.res = res;

end
