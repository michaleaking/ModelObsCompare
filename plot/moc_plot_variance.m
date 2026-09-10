function ax = moc_plot_variance(V, varargin)
%MOC_PLOT_VARIANCE  Map the per-pixel temporal variability from moc_velocity_variance.
%
%   Draws one field of a variance struct as a map: the normalised spread (the
%   one that is comparable between neighbouring pixels), the absolute spread in
%   m/yr, the reference speed, the frame count, or — for a two-slice window —
%   the signed change, which is drawn with a diverging colormap centred on zero.
%
%   Usage:
%      ax = moc_plot_variance(V)                            % V.norm
%      ax = moc_plot_variance(V, 'field','spread')
%      ax = moc_plot_variance(V, 'field','normChange', 'roi',roi)
%      ax = moc_plot_variance(V, 'ax',ax, 'clim',[0 0.5])
%
%   Input:
%      V : struct from moc_velocity_variance
%   Name/value options:
%      'field' : (char)     which map to draw (default 'norm'):
%                           'norm' | 'spread' | 'reference' | 'count' |
%                           'change' | 'normChange'
%      'ax'    : (axes)     target axes (default: new figure)
%      'clim'  : (1x2)      color limits (default [0, 98th pct], or symmetric
%                           about zero for the signed fields)
%      'roi'   : (struct)   overlay this ROI polygon
%      'title' : (char)     axis title (default auto, includes the ROI mean)
%      'cmap'  : (char/Nx3) colormap (default: parula, or blue-white-red for
%                           the signed fields)
%
%   Output:
%      ax : handle to the axes drawn on
%
%   See also MOC_VELOCITY_VARIANCE, MOC_PLOT_MAP, MOC_PLOT_DIFF.

p = inputParser;
p.addParameter('field','norm', @(s) ischar(s)||isstring(s));
p.addParameter('ax',[], @(h) isempty(h)||isgraphics(h));
p.addParameter('clim',[], @(v) isempty(v)||numel(v)==2);
p.addParameter('roi',[], @(v) isempty(v)||isstruct(v));
p.addParameter('title','', @(s) ischar(s)||isstring(s));
p.addParameter('cmap',[], @(v) isempty(v)||ischar(v)||isnumeric(v));
p.parse(varargin{:});
opt = p.Results;

fld = char(opt.field);
if ~isfield(V, fld)
    error('moc_plot_variance:noField', ...
          'Variance struct has no field "%s" (a signed change field exists only for a 2-slice window).', fld);
end

Z = double(V.(fld));
if strcmp(fld,'count'), Z(Z == 0) = NaN; end
[lbl, signed] = local_style(fld, V);

ax = opt.ax;
if isempty(ax), figure('Color','w'); ax = axes(); end
imagesc(ax, V.x, V.y, Z, 'AlphaData', ~isnan(Z));
set(ax,'YDir','normal'); axis(ax,'equal'); axis(ax,'tight'); hold(ax,'on');

clim = opt.clim;
if isempty(clim)
    if signed
        lim = local_prctile(abs(Z(:)), 98);
        if ~isfinite(lim) || lim == 0, lim = 1; end
        clim = [-lim lim];
    else
        clim = [0, local_prctile(Z(:), 98)];
    end
end
caxis(ax, clim);

cmap = opt.cmap;
if isempty(cmap)
    if signed, cmap = local_diverging(); else, cmap = 'parula'; end
end
colormap(ax, cmap);
cb = colorbar(ax); ylabel(cb, lbl);
xlabel(ax,'x [m]'); ylabel(ax,'y [m]');

if ~isempty(opt.roi)
    plot(ax, opt.roi.x, opt.roi.y, 'c-', 'LineWidth', 1.5);
end

if ~isempty(char(opt.title))
    title(ax, char(opt.title), 'Interpreter','none');
else
    who = V.roi;
    if isempty(who), who = V.band; end
    % One line: a second line would collide with the y-axis 1e6 exponent label.
    % The statistic/normaliser are already spelled out on the colorbar.
    ttl = sprintf('%s %s  t=[%.2f..%.2f]  n=%d', who, fld, ...
                  min(V.time), max(V.time), numel(V.time));
    switch fld
        case 'norm',   ttl = sprintf('%s  mean=%.3g', ttl, V.stats.normMean);
        case 'spread', ttl = sprintf('%s  mean=%.3g m/yr', ttl, V.stats.spreadMean);
    end
    title(ax, ttl, 'Interpreter','none');
end

end

% ===========================================================================
function [lbl, signed] = local_style(fld, V)
%LOCAL_STYLE  Colorbar label and whether the field is signed (diverging).
signed = false;
switch fld
    case 'norm'
        lbl = sprintf('%s / %s speed [-]', V.statistic, V.normalize);
    case 'spread'
        lbl = sprintf('%s of speed [m/yr]', V.statistic);
    case 'reference'
        lbl = sprintf('%s speed [m/yr]', V.normalize);
    case 'count'
        lbl = 'valid frames [-]';
    case 'change'
        lbl = 'later - earlier [m/yr]'; signed = true;
    case 'normChange'
        lbl = sprintf('(later - earlier) / %s speed [-]', V.normalize); signed = true;
    otherwise
        lbl = fld;
end
end

function cmap = local_diverging()
%LOCAL_DIVERGING  Blue-white-red diverging colormap (no toolbox needed).
n = 256; h = (n-1)/2;
r = [linspace(0.23,1,h+1), linspace(1,0.71,h)]';
g = [linspace(0.30,1,h+1), linspace(1,0.09,h)]';
b = [linspace(0.75,1,h+1), linspace(1,0.17,h)]';
cmap = [r g b];
end

function q = local_prctile(v, pc)
v = v(isfinite(v));
if isempty(v), q = 1; return; end
v = sort(v); q = v(max(1, round(pc/100*numel(v))));
end
