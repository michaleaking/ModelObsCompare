%% test_variance_synthetic.m
% Self-contained validation of MOC_VELOCITY_VARIANCE / MOC_PLOT_VARIANCE.
%
% Needs no data and no ISSM: it fabricates a two-domain obs grid — fast ice
% (2000 m/yr, +/-100) beside slow ice (100 m/yr, +/-10) — so the ABSOLUTE spread
% is 10x larger on the fast side while the RELATIVE (normalised) spread is 2x
% larger on the slow side, then checks that
%   - the normalised map recovers the known fractions (0.10 vs 0.20) while the
%     absolute spread recovers 200 vs 20 m/yr,
%   - an ROI crops the output grid and a polygon ROI masks pixels that are
%     inside its bounding box but outside the polygon,
%   - two-slice mode picks the nearest frames and the signed change is
%     later-minus-earlier,
%   - pixels with too few finite frames, and references below 'minspeed', are NaN,
%   - every statistic x normalize combination runs, and the plots draw.
%
% Usage:
%   >> test_variance_synthetic          % prints a PASS/FAIL summary

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(here)));      % toolkit only; ISSM is not needed here

fprintf('==========  moc_velocity_variance (synthetic)  ==========\n');

TIMES    = [2019.0 2019.25 2019.5 2019.75];
FAST     = 2000;  SLOW     = 100;       % mean speeds [m/yr] of the two halves
FAST_AMP = 100;   SLOW_AMP = 10;        % +/- seasonal amplitude [m/yr]

R = struct('name',{},'ok',{},'msg',{});

%% Synthetic obs grid (y DESCENDS, as in the real Sentinel subsets) ---------
nx = 50; ny = 50;
O = struct();
O.x = -300000 + 200*(0:nx-1);
O.y = -1810000 - 200*(0:ny-1);
O.time = TIMES;
O.epsg = 3413;
fast   = (1:nx) <= floor(nx/2);
base   = reshape(FAST*fast + SLOW*~fast, [1 nx 1]);
amp    = reshape(FAST_AMP*fast + SLOW_AMP*~fast, [1 nx 1]);
season = reshape([-1 0 1 0], [1 1 numel(TIMES)]);          % peak-to-peak = 2*amp
O.vv   = repmat(base, [ny 1 numel(TIMES)]) + repmat(amp, [ny 1 1]) .* repmat(season, [ny nx 1]);
O.vv(1,1,:) = NaN;                                          % one dead pixel
O.bands = {'vv'};

%% 1. Full grid, whole series ----------------------------------------------
V = moc_velocity_variance(O);
fNorm = V.norm(5,5);         sNorm = V.norm(5,nx-4);
fSpread = V.spread(5,5);     sSpread = V.spread(5,nx-4);

R = rec(R, 'spread is fast-biased', ...
    near(fSpread, 2*FAST_AMP) && near(sSpread, 2*SLOW_AMP), ...
    sprintf('%.0f vs %.0f m/yr', fSpread, sSpread));
R = rec(R, 'normalised spread inverts it', ...
    near(fNorm, 2*FAST_AMP/FAST) && near(sNorm, 2*SLOW_AMP/SLOW), ...
    sprintf('%.3f vs %.3f', fNorm, sNorm));
R = rec(R, 'all-NaN pixel is NaN with count 0', ...
    isnan(V.norm(1,1)) && V.count(1,1)==0, '');
R = rec(R, 'ROI-mean summary in .stats', ...
    near(V.stats.normMean, 0.5*(2*FAST_AMP/FAST + 2*SLOW_AMP/SLOW), 1e-3), ...
    sprintf('normMean=%.3f', V.stats.normMean));
R = rec(R, 'full grid used when roi is empty', ...
    isequal(size(V.norm), [ny nx]), mat2str(size(V.norm)));

%% 2. Rectangular ROI + two time slices ------------------------------------
roi = struct('name','rightbox', 'kind','rectangle', ...
    'x',[-295000 -290500 -290500 -295000 -295000], ...
    'y',[-1819000 -1819000 -1812000 -1812000 -1819000]);
V2 = moc_velocity_variance(O, roi, 'times',[2019.02 2019.48]);

R = rec(R, 'ROI crops the grid', ...
    size(V2.norm,1) < ny && size(V2.norm,2) < nx, ...
    sprintf('%s of %s', mat2str(size(V2.norm)), mat2str([ny nx])));
R = rec(R, 'two nearest frames selected', ...
    numel(V2.time)==2 && near(V2.time(1),2019.0) && near(V2.time(2),2019.5), ...
    mat2str(V2.time));
R = rec(R, 'mean normalised change over the ROI', ...
    near(V2.stats.normMean, 2*SLOW_AMP/SLOW), sprintf('%.3f', V2.stats.normMean));
R = rec(R, 'signed change is later minus earlier', ...
    near(mean(V2.change(isfinite(V2.change))), 2*SLOW_AMP) && ...
    near(mean(V2.normChange(isfinite(V2.normChange))), 2*SLOW_AMP/SLOW), '');
R = rec(R, 'no signed change field for a >2-frame window', ...
    ~isfield(moc_velocity_variance(O, roi), 'change'), '');

%% 3. Polygon ROI masks inside its own bounding box ------------------------
tri = struct('name','tri', 'kind','polygon', ...
    'x',[-300000 -290200 -300000 -300000], ...
    'y',[-1819800 -1819800 -1810000 -1819800]);
V3 = moc_velocity_variance(O, tri, 'trange',[2019.0 2019.6], ...
                           'statistic','std', 'normalize','median');
inside = mean(isfinite(V3.norm(:)));
R = rec(R, 'polygon masks ~half its bbox', inside > 0.35 && inside < 0.65, ...
    sprintf('%.0f%% kept', 100*inside));

%% 4. Guards ----------------------------------------------------------------
V4 = moc_velocity_variance(O, [], 'minspeed',500);
R = rec(R, 'minspeed drops slow-ice references', ...
    isnan(V4.norm(5,nx-4)) && isfinite(V4.norm(5,5)), '');
V5 = moc_velocity_variance(O, [], 'mincount',5);
R = rec(R, 'mincount drops thin pixels', all(isnan(V5.norm(:))), '');
R = rec(R, 'bad arguments error', ...
    raises(@() moc_velocity_variance(O,[],'statistic','nope')) && ...
    raises(@() moc_velocity_variance(O,[],'normalize','nope')) && ...
    raises(@() moc_velocity_variance(O,[],'band','zz')) && ...
    raises(@() moc_velocity_variance(O,[],'trange',[2050 2060])) && ...
    raises(@() moc_velocity_variance(O,[],'times',2019.0)) && ...
    raises(@() moc_velocity_variance(O,[],'times',[2019 2019.5],'trange',[2019 2019.5])), '');

%% 5. Every statistic x normalize combination -------------------------------
STATS = {'range','std','iqr','mad'};
NORMS = {'mean','median','first','max','none'};
combosOk = true;
for i = 1:numel(STATS)
    for j = 1:numel(NORMS)
        v = moc_velocity_variance(O, [], 'statistic',STATS{i}, 'normalize',NORMS{j});
        combosOk = combosOk && isfinite(v.stats.normMean) && v.stats.npix > 0;
    end
end
R = rec(R, 'statistic x normalize combinations run', combosOk, ...
    sprintf('%dx%d', numel(STATS), numel(NORMS)));

%% 6. Plots ------------------------------------------------------------------
plotOk = true; plotMsg = '';
try
    fig = figure('Color','w','Visible','off','Position',[100 100 1500 460]);
    moc_plot_variance(V,  'field','norm',       'ax',subplot(1,3,1), 'roi',roi);
    moc_plot_variance(V,  'field','spread',     'ax',subplot(1,3,2));
    moc_plot_variance(V2, 'field','normChange', 'ax',subplot(1,3,3));
    close(fig);
catch ME
    plotOk = false; plotMsg = ME.message;
end
R = rec(R, 'moc_plot_variance draws norm / spread / normChange', plotOk, plotMsg);
R = rec(R, 'moc_plot_variance rejects an absent field', ...
    raises(@() moc_plot_variance(V, 'field','change')), '');

%% Summary -------------------------------------------------------------------
fprintf('\n----------------------------  summary  ----------------------------\n');
for i = 1:numel(R)
    if R(i).ok, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  [%-4s] %-45s %s\n', tag, R(i).name, R(i).msg);
end
if all([R.ok])
    fprintf('\nALL VARIANCE CHECKS PASSED\n');
else
    fprintf(2, '\n%d CHECK(S) FAILED\n', sum(~[R.ok]));
end

%% ===========================================================================
function R = rec(R, name, ok, msg)
%REC  Append one check result to the accumulator.
R(end+1) = struct('name',name, 'ok',logical(ok), 'msg',msg);
end

function tf = near(a, b, tol)
%NEAR  Absolute-tolerance scalar comparison (default 1e-9).
if nargin < 3, tol = 1e-9; end
tf = isfinite(a) && abs(a - b) <= tol;
end

function tf = raises(fh)
%RAISES  True if calling FH throws.
try
    fh();
    tf = false;
catch
    tf = true;
end
end
