function roi = moc_select_roi(base, varargin)
%MOC_SELECT_ROI  Interactively draw a region of interest on a background field.
%
%   Displays a background velocity field and lets you draw a rectangle or
%   polygon; returns an ROI struct and (optionally) saves it to disk as both a
%   .mat and an ISSM Argus .exp contour for reuse in mesh/interp tools.
%
%   Usage:
%      roi = moc_select_roi(O)                       % O from moc_load_obs_netcdf
%      roi = moc_select_roi(G)                        % G from moc_model_to_grid
%      roi = moc_select_roi(base, 'shape','polygon', 'name','upernavik')
%
%   Input:
%      base : a struct describing the backdrop. Either
%             - an obs grid struct with .x (1xnx), .y (1xny) and a 2-D field, or
%             - {x, y, field2d} passed as three separate args (see options), or
%             - a model struct M (from moc_load_model): shown as scattered speed.
%   Name/value options:
%      'shape'   : (char) 'rectangle'(default) | 'polygon'
%      'field'   : (char) which field name in an obs struct to display (default 'vv')
%      'tindex'  : (double) time slice to display for 3-D fields (default 1)
%      'name'    : (char) ROI name; if given, saves rois/<name>.mat and .exp
%      'save'    : (logical) force save even without a name (default false)
%      'clim'    : (1x2) color limits for the backdrop
%
%   Output (struct roi):
%      .name  (char)            ROI name
%      .type  (char)            'rectangle' | 'polygon'
%      .x,.y  (Nx1 double)      closed polygon vertices [m]
%      .bbox  (1x4 double)      [xmin xmax ymin ymax]
%      .epsg  (double)          projection code inherited from the backdrop
%      .file  (char)            saved .mat path ('' if not saved)

p = inputParser;
p.addParameter('shape','rectangle', @(s) any(strcmpi(s,{'rectangle','polygon'})));
p.addParameter('field','vv', @(s) ischar(s)||isstring(s));
p.addParameter('tindex',1, @isscalar);
p.addParameter('name','', @(s) ischar(s)||isstring(s));
p.addParameter('save',false, @islogical);
p.addParameter('clim',[], @(v) isempty(v)||numel(v)==2);
p.parse(varargin{:});
opt = p.Results;
shape = lower(char(opt.shape));

% --- draw the backdrop ---
fig = figure('Name','Select ROI','Color','w');
ax  = axes(fig); hold(ax,'on');
epsg = NaN;
[xb, yb, Zb, epsg] = local_backdrop(base, char(opt.field), opt.tindex, epsg);
if ~isempty(Zb)
    imagesc(ax, xb, yb, Zb, 'AlphaData', ~isnan(Zb));
    set(ax,'YDir','normal');
    if ~isempty(opt.clim), caxis(ax, opt.clim); else
        caxis(ax, [0, local_prctile(Zb(:), 98)]);
    end
    colormap(ax, parula); cb = colorbar(ax); ylabel(cb,'speed [m/yr]');
end
axis(ax,'equal'); axis(ax,'tight');
xlabel(ax,'x [m]'); ylabel(ax,'y [m]');
title(ax, sprintf('Draw %s ROI, then double-click / press Enter to accept', shape));

% --- interactive drawing (modern ROI tools with ginput fallback) ---
[px, py, usedShape] = local_draw(ax, shape);

% Close the polygon
if px(1)~=px(end) || py(1)~=py(end)
    px(end+1) = px(1); py(end+1) = py(1);
end

roi = struct();
roi.name = char(opt.name);
roi.type = usedShape;
roi.x = px(:);
roi.y = py(:);
roi.bbox = [min(px) max(px) min(py) max(py)];
roi.epsg = epsg;
roi.file = '';

% Overlay the final ROI
plot(ax, roi.x, roi.y, 'r-', 'LineWidth', 2);

% --- optional save ---
if opt.save || ~isempty(roi.name)
    if isempty(roi.name)
        roi.name = ['roi_' datestr(now,'yyyymmdd_HHMMSS')];
    end
    roi = moc_save_roi(roi);
end

end

% ===========================================================================
function [xb, yb, Zb, epsg] = local_backdrop(base, field, tindex, epsg)
%LOCAL_BACKDROP  Extract (x, y, 2-D field) to display behind the ROI.
xb = []; yb = []; Zb = [];
if isstruct(base) && isfield(base,'x') && isfield(base,'y')
    if isfield(base,'epsg') && ~isnan(base.epsg), epsg = base.epsg; end
    xb = base.x(:)'; yb = base.y(:)';
    if isfield(base, field)
        Z = base.(field);
    elseif isfield(base,'grid')      % model-to-grid struct
        Z = base.grid;
    elseif isfield(base,'vv')
        Z = base.vv;
    else
        Z = [];
    end
    if ~isempty(Z)
        if ndims(Z)==3, Z = Z(:,:,min(tindex,size(Z,3))); end
        Zb = Z;
    end
elseif isstruct(base) && isfield(base,'vel')   % model struct: scatter speed
    % Build a coarse background grid from the mesh vertices for context.
    xb = linspace(min(base.x), max(base.x), 400);
    yb = linspace(min(base.y), max(base.y), 400);
    v  = base.vel(:, min(tindex, size(base.vel,2)));
    Zb = griddata(base.x, base.y, v, xb, yb');
    if isfield(base,'epsg'), epsg = base.epsg; end
end
end

% ---------------------------------------------------------------------------
function [px, py, usedShape] = local_draw(ax, shape)
%LOCAL_DRAW  Use drawrectangle/drawpolygon if available, else ginput.
usedShape = shape;
hasDraw = exist('drawrectangle','file') || exist('drawpolygon','file');
if hasDraw
    if strcmp(shape,'rectangle')
        h = drawrectangle(ax);
        wait(h);
        pos = h.Position;   % [x y w h]
        px = [pos(1), pos(1)+pos(3), pos(1)+pos(3), pos(1)];
        py = [pos(2), pos(2), pos(2)+pos(4), pos(2)+pos(4)];
    else
        h = drawpolygon(ax);
        wait(h);
        px = h.Position(:,1)'; py = h.Position(:,2)';
    end
else
    % Fallback: manual clicking
    title(ax, 'Click ROI vertices; press Enter when done');
    [px, py] = ginput;
    px = px'; py = py';
    if strcmp(shape,'rectangle') && numel(px)>=2
        % interpret first two clicks as opposite corners
        x1 = px(1); x2 = px(2); y1 = py(1); y2 = py(2);
        px = [x1 x2 x2 x1]; py = [y1 y1 y2 y2];
    end
end
end

% ---------------------------------------------------------------------------
function q = local_prctile(v, pc)
%LOCAL_PRCTILE  Percentile without the Statistics Toolbox dependency.
v = v(~isnan(v));
if isempty(v), q = 1; return; end
v = sort(v);
q = v(max(1, round(pc/100*numel(v))));
end
