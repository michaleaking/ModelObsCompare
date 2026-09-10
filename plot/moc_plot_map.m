function ax = moc_plot_map(G, varargin)
%MOC_PLOT_MAP  Plot a gridded field (model-on-grid or an observed slice) as a map.
%
%   Usage:
%      ax = moc_plot_map(G)                          % G from moc_model_to_grid
%      ax = moc_plot_map(O, 'band','vv', 'tindex',1) % O from moc_load_obs_netcdf
%      ax = moc_plot_map(G, 'ax',ax, 'clim',[0 3000], 'roi',roi)
%
%   Input:
%      G : gridded struct. Either
%          - moc_model_to_grid output (.x,.y,.grid), or
%          - moc_load_obs_netcdf output (.x,.y,.<band> as ny x nx x nt).
%   Name/value options:
%      'band'   : (char)    band name for an obs struct (default 'vv')
%      'tindex' : (double)  time slice for a 3-D obs field (default 1)
%      'ax'     : (axes)    target axes (default: new figure)
%      'clim'   : (1x2)     color limits (default: [0, 98th pct])
%      'roi'    : (struct)  overlay this ROI polygon
%      'title'  : (char)    axis title
%      'cmap'   : (char/Nx3) colormap (default parula)
%
%   Output:
%      ax : handle to the axes drawn on

p = inputParser;
p.addParameter('band','vv', @(s) ischar(s)||isstring(s));
p.addParameter('tindex',1, @isscalar);
p.addParameter('ax',[], @(h) isempty(h)||isgraphics(h));
p.addParameter('clim',[], @(v) isempty(v)||numel(v)==2);
p.addParameter('roi',[], @(v) isempty(v)||isstruct(v));
p.addParameter('title','', @(s) ischar(s)||isstring(s));
p.addParameter('cmap','parula', @(v) ischar(v)||isnumeric(v));
p.parse(varargin{:});
opt = p.Results;

[x, y, Z, deflabel] = local_get_field(G, char(opt.band), opt.tindex);

ax = opt.ax;
if isempty(ax), figure('Color','w'); ax = axes(); end
imagesc(ax, x, y, Z, 'AlphaData', ~isnan(Z));
set(ax,'YDir','normal'); axis(ax,'equal'); axis(ax,'tight'); hold(ax,'on');

if isempty(opt.clim)
    caxis(ax, [0, local_prctile(Z(:), 98)]);
else
    caxis(ax, opt.clim);
end
colormap(ax, opt.cmap);
cb = colorbar(ax); ylabel(cb, 'speed [m/yr]');
xlabel(ax,'x [m]'); ylabel(ax,'y [m]');

if ~isempty(opt.roi)
    plot(ax, opt.roi.x, opt.roi.y, 'r-', 'LineWidth', 1.5);
end
if ~isempty(char(opt.title))
    title(ax, char(opt.title), 'Interpreter','none');
else
    title(ax, deflabel, 'Interpreter','none');
end

end

% ===========================================================================
function [x, y, Z, lbl] = local_get_field(G, band, tindex)
x = G.x(:)'; y = G.y(:)'; lbl = '';
if isfield(G,'grid')
    Z = G.grid;
    lbl = sprintf('model %s', local_str(G,'field','vel'));
    if isfield(G,'time') && ~isnan(G.time), lbl = sprintf('%s  t=%.3f', lbl, G.time); end
elseif isfield(G, band)
    Z = G.(band);
    if ndims(Z)==3, Z = Z(:,:,min(tindex,size(Z,3))); end
    lbl = sprintf('obs %s', band);
    if isfield(G,'time') && numel(G.time)>=tindex, lbl = sprintf('%s  t=%.3f', lbl, G.time(tindex)); end
else
    error('moc_plot_map:noField', 'Struct has neither .grid nor band "%s".', band);
end
end

function s = local_str(S, f, def)
if isfield(S,f) && (ischar(S.(f))||isstring(S.(f))), s = char(S.(f)); else, s = def; end
end

function q = local_prctile(v, pc)
v = v(isfinite(v));
if isempty(v), q = 1; return; end
v = sort(v); q = v(max(1, round(pc/100*numel(v))));
end
