%% example_velocity_variance.m
% Map how much the observed speed varies over a temporal window, normalised so
% neighbouring pixels are comparable as a FRACTION of relative change rather
% than an absolute magnitude change.
% Pipeline: load gridded obs -> pick/define ROI -> reduce the time axis to a
% per-pixel spread + reference speed -> plot the normalised map -> repeat for
% two individual time slices.

cfg = setup_moc();

%% 1. Load gridded observations --------------------------------------------
O = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc', ...
                        'bands',{'vv'}, 'trange',[2018 2020]);

%% 2. Define the ROI (or pass [] for the full grid) -------------------------
if isfile(fullfile(moc_path('roi_dir'),'upernavik_box.mat'))
    roi = moc_load_roi('upernavik_box');
else
    roi = moc_roi_from_obs(O, 'margin',2000);
end

%% 3. Variability over a window ---------------------------------------------
% Peak-to-peak speed over the 2019 melt season, divided by each pixel's own
% mean speed over the same window.
V = moc_velocity_variance(O, roi, 'trange',[2019.4 2019.8]);

figure('Color','w','Position',[100 100 1100 460]);
moc_plot_variance(V, 'field','norm',   'ax',subplot(1,2,1), 'roi',roi);
moc_plot_variance(V, 'field','spread', 'ax',subplot(1,2,2), 'roi',roi);
% Left: comparable between pixels. Right: the same window in m/yr, which the
% fast trunk dominates simply because it is fast.

fprintf('mean relative change over the ROI: %.1f%% (median %.1f%%, p90 %.1f%%)\n', ...
    100*V.stats.normMean, 100*V.stats.normMedian, 100*V.stats.normP90);

%% 4. Between two individual time slices ------------------------------------
% The frames nearest each date are used; this also gives a SIGNED change field.
V2 = moc_velocity_variance(O, roi, 'times',[2019.45 2019.75]);

figure('Color','w','Position',[100 100 1100 460]);
moc_plot_variance(V2, 'field','normChange', 'ax',subplot(1,2,1), 'roi',roi);
moc_plot_variance(V2, 'field','change',     'ax',subplot(1,2,2), 'roi',roi);
fprintf('t=[%.3f %.3f]: mean |change| = %.1f%% of local speed\n', ...
    V2.time(1), V2.time(2), 100*V2.stats.normMean);

%% 5. Variants ---------------------------------------------------------------
% Robust spread (IQR) about the median, ignoring near-stagnant ice where a
% ratio would explode:
Vr = moc_velocity_variance(O, roi, 'trange',[2019 2020], ...
        'statistic','iqr', 'normalize','median', 'minspeed',20, 'mincount',5);
moc_plot_variance(Vr, 'field','norm', 'roi',roi);

% Coverage check: how many frames actually contributed to each pixel.
moc_plot_variance(Vr, 'field','count', 'roi',roi);
