function fig = moc_plot_diff(D, varargin)
%MOC_PLOT_DIFF  Three-panel model / obs / difference map from moc_spatial_diff.
%
%   Draws the gridded model field, the observed field, and their difference
%   (model - obs) on a shared spatial grid, with a diverging colormap for the
%   difference panel centred on zero.
%
%   Usage:
%      fig = moc_plot_diff(D)                       % D from moc_spatial_diff
%      fig = moc_plot_diff(D, 'clim',[0 3000], 'dlim',500, 'roi',roi)
%
%   Input:
%      D : struct from moc_spatial_diff (.x,.y,.model,.obs,.diff,.stats)
%   Name/value options:
%      'clim' : (1x2)    speed color limits for model/obs panels (default auto)
%      'dlim' : (double) symmetric difference limit +/- [m/yr] (default auto ~95th pct)
%      'roi'  : (struct) overlay ROI polygon on each panel
%      'cmap' : (char)   colormap for the speed panels (default parula)
%
%   Output:
%      fig : figure handle

p = inputParser;
p.addParameter('clim',[], @(v) isempty(v)||numel(v)==2);
p.addParameter('dlim',[], @(v) isempty(v)||isscalar(v));
p.addParameter('roi',[], @(v) isempty(v)||isstruct(v));
p.addParameter('cmap','parula', @(v) ischar(v)||isnumeric(v));
p.parse(varargin{:});
opt = p.Results;

clim = opt.clim;
if isempty(clim)
    both = [D.model(:); D.obs(:)];
    clim = [0, local_prctile(both, 98)];
end
dlim = opt.dlim;
if isempty(dlim)
    dlim = local_prctile(abs(D.diff(:)), 95);
    if ~isfinite(dlim) || dlim==0, dlim = 1; end
end

fig = figure('Color','w','Position',[100 100 1500 460]);

ax1 = subplot(1,3,1);
moc_plot_map(struct('x',D.x,'y',D.y,'grid',D.model,'field','model','time',D.time), ...
    'ax',ax1,'clim',clim,'cmap',opt.cmap,'roi',opt.roi,'title', ...
    sprintf('model  t=%.3f', D.time));

ax2 = subplot(1,3,2);
moc_plot_map(struct('x',D.x,'y',D.y,'grid',D.obs,'field','obs'), ...
    'ax',ax2,'clim',clim,'cmap',opt.cmap,'roi',opt.roi,'title', ...
    sprintf('obs  t=[%.2f..%.2f]', min(D.obstime), max(D.obstime)));

ax3 = subplot(1,3,3);
imagesc(ax3, D.x, D.y, D.diff, 'AlphaData', ~isnan(D.diff));
set(ax3,'YDir','normal'); axis(ax3,'equal'); axis(ax3,'tight'); hold(ax3,'on');
caxis(ax3, [-dlim dlim]);
colormap(ax3, local_diverging());
cb = colorbar(ax3); ylabel(cb,'model - obs [m/yr]');
xlabel(ax3,'x [m]'); ylabel(ax3,'y [m]');
if ~isempty(opt.roi), plot(ax3, opt.roi.x, opt.roi.y, 'k-','LineWidth',1.5); end
title(ax3, sprintf('diff  RMSE=%.0f  MAE=%.0f  bias=%.0f m/yr', ...
    D.stats.rmse, D.stats.mae, D.stats.meanDiff), 'Interpreter','none');

end

% ===========================================================================
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
