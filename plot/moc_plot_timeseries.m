function ax = moc_plot_timeseries(C, varargin)
%MOC_PLOT_TIMESERIES  Plot paired modelled vs observed velocity time series.
%
%   Usage:
%      ax = moc_plot_timeseries(C)                  % C from moc_compare_timeseries
%      ax = moc_plot_timeseries(C, 'ax',ax, 'datetime',true)
%
%   Input:
%      C : struct from moc_compare_timeseries
%   Name/value options:
%      'ax'       : (axes)    target axes (default: new figure)
%      'datetime' : (logical) plot on a datetime x-axis (default false = decimal yr)
%      'title'    : (char)    axis title
%
%   Output:
%      ax : axes handle

p = inputParser;
p.addParameter('ax',[], @(h) isempty(h)||isgraphics(h));
p.addParameter('datetime',false, @islogical);
p.addParameter('title','', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

ax = opt.ax;
if isempty(ax), figure('Color','w'); ax = axes(); end
hold(ax,'on'); grid(ax,'on');

% x-axis conversion
mt = C.model_time; ot = C.obs_time;
if opt.datetime
    mt = datetime(moc_decyear2datenum(mt),'ConvertFrom','datenum');
    ot = datetime(moc_decyear2datenum(ot),'ConvertFrom','datenum');
end

% observed with error bars where available
if ~isempty(C.obs_err) && any(isfinite(C.obs_err))
    errorbar(ax, ot, C.obs_val, C.obs_err, 'o', 'Color',[0.2 0.2 0.2], ...
        'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerSize',4, 'CapSize',0, ...
        'DisplayName','observed');
else
    plot(ax, ot, C.obs_val, 'o', 'Color',[0.2 0.2 0.2], ...
        'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerSize',4, 'DisplayName','observed');
end

% modelled line
plot(ax, mt, C.model_val, '-', 'Color',[0.85 0.1 0.1], 'LineWidth',1.8, ...
    'DisplayName','model');

xlabel(ax, ternary(opt.datetime,'date','decimal year'));
ylabel(ax,'speed [m/yr]');
legend(ax,'Location','best');

if ~isempty(char(opt.title))
    ttl = char(opt.title);
else
    ttl = sprintf('(%.0f, %.0f)  RMSE=%.0f  bias=%.0f  r=%.2f', ...
        C.xy(1), C.xy(2), C.stats.rmse, C.stats.bias, C.stats.r);
end
title(ax, ttl, 'Interpreter','none');

end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
