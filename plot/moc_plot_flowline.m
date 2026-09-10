function ax = moc_plot_flowline(Fmod, Fobs, varargin)
%MOC_PLOT_FLOWLINE  Compare modelled and observed velocity along a flowline.
%
%   Two modes:
%     'profile' : speed vs along-flowline distance at one time (model & obs).
%     'hovmoller': model-minus-obs difference as a distance-time image.
%
%   Usage:
%      ax = moc_plot_flowline(Fmod, Fobs)                     % profile at matched time
%      ax = moc_plot_flowline(Fmod, Fobs, 'mode','hovmoller')
%      ax = moc_plot_flowline(Fmod, Fobs, 'mode','profile', 'time',2019.0)
%
%   Input:
%      Fmod : modelled flowline struct (moc_model_flowline): .d,.time,.val
%      Fobs : observed flowline struct (moc_load_obs_flowline): .d,.time,.vel
%   Name/value options:
%      'mode'  : (char) 'profile'(default) | 'hovmoller'
%      'time'  : (double) decimal year for profile mode (default: median obs time)
%      'ax'    : (axes) target axes
%      'dlim'  : (double) symmetric diff limit for hovmoller [m/yr]
%
%   Output:
%      ax : axes handle

p = inputParser;
p.addParameter('mode','profile', @(s) any(strcmpi(s,{'profile','hovmoller'})));
p.addParameter('time',[], @(v) isempty(v)||isscalar(v));
p.addParameter('ax',[], @(h) isempty(h)||isgraphics(h));
p.addParameter('dlim',[], @(v) isempty(v)||isscalar(v));
p.parse(varargin{:});
opt = p.Results;
mode = lower(char(opt.mode));

ax = opt.ax;
if isempty(ax), figure('Color','w'); ax = axes(); end

switch mode
    case 'profile'
        tt = opt.time;
        if isempty(tt), tt = median(Fobs.time); end
        % nearest model and obs time columns
        [~, mi] = min(abs(Fmod.time - tt));
        [~, oi] = min(abs(Fobs.time - tt));
        hold(ax,'on'); grid(ax,'on');
        plot(ax, Fobs.d/1e3, Fobs.vel(:,oi), 'o-', 'Color',[0.3 0.3 0.3], ...
            'MarkerSize',3, 'DisplayName',sprintf('obs t=%.2f',Fobs.time(oi)));
        plot(ax, Fmod.d/1e3, Fmod.val(:,mi), '-', 'Color',[0.85 0.1 0.1], ...
            'LineWidth',1.8, 'DisplayName',sprintf('model t=%.2f',Fmod.time(mi)));
        xlabel(ax,'along-flowline distance [km]');
        ylabel(ax,'speed [m/yr]');
        legend(ax,'Location','best');
        title(ax, sprintf('flowline profile near t=%.2f', tt), 'Interpreter','none');

    case 'hovmoller'
        % interpolate model onto obs times, difference on the obs distance axis
        modI = local_interp_time(Fmod.d, Fmod.time, Fmod.val, Fobs.d, Fobs.time);
        Dfl = modI - Fobs.vel;             % (npts x nobs)
        dlim = opt.dlim;
        if isempty(dlim)
            dlim = local_prctile(abs(Dfl(:)), 95);
            if ~isfinite(dlim)||dlim==0, dlim = 1; end
        end
        imagesc(ax, Fobs.time, Fobs.d/1e3, Dfl, 'AlphaData', ~isnan(Dfl));
        set(ax,'YDir','normal'); caxis(ax,[-dlim dlim]);
        colormap(ax, local_diverging());
        cb = colorbar(ax); ylabel(cb,'model - obs [m/yr]');
        xlabel(ax,'decimal year'); ylabel(ax,'along-flowline distance [km]');
        title(ax,'flowline difference (model - obs)','Interpreter','none');
end

end

% ===========================================================================
function modI = local_interp_time(dmod, tmod, vmod, dobs, tobs)
%LOCAL_INTERP_TIME  Put model (dmod x tmod) onto the obs (dobs x tobs) grid.
% 1) interpolate in distance to obs points, 2) interpolate in time to obs times.
tmp = interp1(dmod, vmod, dobs, 'linear', NaN);   % (nobsPts x tmod)
modI = nan(numel(dobs), numel(tobs));
for i = 1:numel(dobs)
    if all(isnan(tmod)) || numel(tmod)<2
        modI(i,:) = NaN;
    else
        modI(i,:) = interp1(tmod, tmp(i,:), tobs, 'linear', NaN);
    end
end
end

function cmap = local_diverging()
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
