function passed = test_flowline_shp_synthetic()
%TEST_FLOWLINE_SHP_SYNTHETIC  Self-contained validation of the flowline/ pipeline.
%
%   Writes a miniature stand-in for the two shapefile collections into a temp
%   folder — three parallel flowlines running east (the terminus end first) and
%   three terminus traces crossing them at known x, one per year — so every step
%   has an analytic answer:
%
%     * the ROI selects only the flowlines that reach into it,
%     * the terminus reader parses SourceDate and honours the quality cut,
%     * moc_center_flowline picks the middle line of the fan, and
%       moc_ordered_flowlines ranks it centre-outwards,
%     * when no trace crosses the centre line, the search steps outwards and
%       uses the first neighbour one does cross (glacier0003 below), while
%       naming a flowline explicitly still fails rather than silently swapping,
%     * moc_clip_upstream cuts at the NEWEST trace (x=2000, not the older ones)
%       and keeps the inland side, for a line drawn either way round,
%     * moc_clip_to_roi trims the result to the ROI box,
%     * the packaged struct has monotonic d starting at zero and honours spacing,
%     * moc_sample_obs_along_flowline recovers a known linear velocity field,
%     * the failure paths error (empty ROI, unknown glacier, no crossing trace).
%
%   MATLAB mirror of python/tests/test_flowline_shp_synthetic.py.
%
%   Usage:
%      >> setup_moc
%      >> test_flowline_shp_synthetic

Y_LINES   = [0 100 200];                       % the fan: three lines, 100 m apart
TERM_DATE = {'2018-06-01','2019-06-01','2021-08-15'};
TERM_X    = [500 1200 2000];
BAD_X     = 4000;                              % a Quality_Fl=1 trace, newest of all
ROI       = [-500 6000 -50 250];               % [xmin xmax ymin ymax]
NEWEST_X  = TERM_X(end);

% A second fan, north of every trace above (they stop at y=500), so it needs its
% own ROI. Its newest trace is a stub spanning only the y=1000 line, so the
% centre line (y=1100) is missed and the search must step outwards to reach it.
Y_LINES_2 = [1000 1100 1200];
ROI_2     = [-500 6000 950 1250];
STUB_X    = 2000;
STUB_Y    = [950 1050];

fprintf('=== flowline_shp (synthetic shapefiles) ===\n');
if isempty(which('shaperead'))
    fprintf('  SKIPPED — needs the Mapping Toolbox\n');
    passed = true; return;
end

tmp = tempname(); mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
[fdir, tpath] = local_write_fixture(tmp, Y_LINES, TERM_DATE, TERM_X, BAD_X, ...
    Y_LINES_2, STUB_X, STUB_Y);

ok = true;

%% readers
F = moc_load_flowlines(ROI, 'dir', fdir);
ok = ok & check(sprintf('ROI keeps only glacier0001''s 3 flowlines (got %d rows, glaciers %s)', ...
    numel(F), strjoin(unique({F.glacier}),',')), ...
    numel(F)==3 && isequal(unique({F.glacier}), {'0001'}));
ok = ok & check('_iter files excluded by default', all([F.iteration]==0));
ok = ok & check('_iter files included on request', ...
    numel(moc_load_flowlines(ROI, 'dir', fdir, 'iterations', true))==4);
ok = ok & check('a non-flowline shapefile in the folder is skipped', ...
    ~any(contains({F.file}, 'some_masks')));
ok = ok & check('no ROI returns the whole collection (7 lines)', ...
    numel(moc_load_flowlines([], 'dir', fdir))==7);
ok = ok & check('a second ROI selects only glacier0003', ...
    isequal(unique({moc_load_flowlines(ROI_2, 'dir', fdir).glacier}), {'0003'}));

T = moc_load_termini(ROI, 'path', tpath);
ok = ok & check(sprintf('quality cut drops the flagged trace (%d of 4 kept)', numel(T)), ...
    numel(T)==3);
ok = ok & check('dates parsed and sorted oldest first', ...
    strcmp(datestr(T(1).date,'yyyy-mm-dd'), '2018-06-01') && ...
    strcmp(datestr(T(end).date,'yyyy-mm-dd'), '2021-08-15')); %#ok<DATST>
ok = ok & check('max_quality=[] keeps the flagged trace', ...
    numel(moc_load_termini(ROI, 'path', tpath, 'max_quality', []))==4);
ok = ok & check('trange bounds the traces', ...
    numel(moc_load_termini(ROI, 'path', tpath, 'trange', {'2019-01-01',''}))==2);

%% selection
[c, ~] = moc_center_flowline(F, 'glacier', '0001');
ok = ok & check(sprintf('centre of the fan is the middle line (got %s)', c.flowline), ...
    strcmp(c.flowline, '04'));
ordIdx = moc_ordered_flowlines(F, 'glacier', '0001');
ok = ok & check('ordered_flowlines ranks the fan centre-outwards', ...
    strcmp(F(ordIdx(1)).flowline, '04'));
newest = moc_latest_terminus(T, 'line', [c.x c.y]);
ok = ok & check('latest_terminus returns the newest crossing trace', ...
    strcmp(datestr(newest.date,'yyyy-mm-dd'), '2021-08-15')); %#ok<DATST>

%% the split
[xu, yu] = moc_clip_upstream(c.x, c.y, newest.x, newest.y);
ok = ok & check(sprintf('split at the newest trace, not an older one (starts at x=%.0f)', xu(1)), ...
    abs(xu(1) - NEWEST_X) < 1e-6);
ok = ok & check('upstream side kept (runs inland to x=10000)', abs(xu(end) - 10000) < 1e-6);
lenUp = sum(hypot(diff(xu), diff(yu)));
[xr, yr] = moc_clip_upstream(flipud(c.x), flipud(c.y), newest.x, newest.y);
ok = ok & check('auto orientation handles an inland-first line too', ...
    abs(sum(hypot(diff(xr), diff(yr))) - lenUp) < 1e-6);
[xe, ye] = moc_clip_upstream(flipud(c.x), flipud(c.y), newest.x, newest.y, ...
    'downstream_end', 'last');
ok = ok & check('explicit downstream_end=''last'' agrees', ...
    abs(sum(hypot(diff(xe), diff(ye))) - lenUp) < 1e-6);
[xn, ~] = moc_clip_upstream(c.x, c.y, [-5000; -5000], [-10; 10]);
ok = ok & check('a non-crossing terminus returns empty', isempty(xn));

%% the ROI clip
[xc, yc] = moc_clip_to_roi(xu, yu, ROI);
lenIn = sum(hypot(diff(xc), diff(yc)));
ok = ok & check(sprintf('clipped to the ROI''s 6 km edge (len=%.0f m)', lenIn), ...
    abs(lenIn - (ROI(2) - NEWEST_X)) < 120);

%% packaging
FL = moc_flowline_from_roi(ROI, 'dir', fdir, 'termini_path', tpath);
ok = ok & check(sprintf('struct has x/y/d, d(1)=0, %d points', numel(FL.x)), ...
    FL.d(1)==0 && all(diff(FL.d) > 0) && numel(FL.x)==numel(FL.d));
ok = ok & check(sprintf('starts at the terminus (x=%.0f) and records its date (%s)', ...
    FL.x(1), FL.terminus_date), ...
    abs(FL.x(1) - NEWEST_X) < 1e-6 && strcmp(FL.terminus_date, '2021-08-15'));
ok = ok & check(sprintf('defaults to the centre flowline (got %s)', FL.flowline), ...
    strcmp(FL.glacier,'0001') && strcmp(FL.flowline,'04'));

FL2 = moc_flowline_from_roi(ROI, 'dir', fdir, 'termini_path', tpath, ...
    'spacing', 100, 'flowline', '05');
ok = ok & check(sprintf('spacing=100 gives an even grid (%d points)', numel(FL2.d)), ...
    all(abs(diff(FL2.d) - 100) < 1e-9) && strcmp(FL2.flowline,'05'));

%% drop-in with the rest of the toolkit
O = local_make_obs();
Fobs = moc_sample_obs_along_flowline(O, FL2);
ok = ok & check('sampled obs is (npoints x ntime) with .vel', ...
    isequal(size(Fobs.vel), [numel(FL2.x) numel(O.time)]));
ok = ok & check('sampled speed recovers the known field (vel == x)', ...
    max(abs(Fobs.vel(:,1) - FL2.x)) < 1e-6);
ok = ok & check('d/x/y survive the sampling', isequal(Fobs.d, FL2.d));
ok = ok & check('output satisfies moc_compare_timeseries''s flowline test', ...
    isfield(Fobs,'d') && isfield(Fobs,'vel') && isfield(Fobs,'source'));

%% stepping outwards when the centre line has no trace
F2 = moc_load_flowlines(ROI_2, 'dir', fdir);
c2 = moc_center_flowline(F2, 'glacier', '0003');
T2 = moc_load_termini(ROI_2, 'path', tpath);
[xMiss, ~] = moc_clip_upstream(c2.x, c2.y, T2(1).x, T2(1).y);
ok = ok & check(sprintf('the stub trace misses the centre line (%s)', c2.flowline), ...
    strcmp(c2.flowline,'04') && numel(T2)==1 && isempty(xMiss));

FL4 = moc_flowline_from_roi(ROI_2, 'dir', fdir, 'termini_path', tpath);
ok = ok & check(sprintf(['search steps outwards to a line the trace does cross ' ...
    '(used %s, centre was %s)'], FL4.flowline, FL4.center_flowline), ...
    strcmp(FL4.flowline,'03') && strcmp(FL4.center_flowline,'04') && ...
    strcmp(FL4.terminus_date,'2021-08-15'));
ok = ok & check('the swapped line is still clipped at the terminus', ...
    abs(FL4.x(1) - STUB_X) < 1e-6);
ok = ok & check('center_flowline equals flowline when no swap happened', ...
    strcmp(FL.center_flowline, FL.flowline) && strcmp(FL.flowline,'04'));
ok = ok & check('naming a flowline explicitly does NOT swap, it errors', ...
    raises(@() moc_flowline_from_roi(ROI_2, 'flowline','04', 'dir', fdir, ...
        'termini_path', tpath)));

%% failure paths
ok = ok & check('an ROI with no flowlines errors', ...
    raises(@() moc_flowline_from_roi([1e6 1.1e6 1e6 1.1e6], 'dir', fdir, 'termini_path', tpath)));
ok = ok & check('an unknown glacier errors', ...
    raises(@() moc_flowline_from_roi(ROI, 'glacier','9999', 'dir', fdir, 'termini_path', tpath)));
ok = ok & check('an unknown flowline id errors', ...
    raises(@() moc_flowline_from_roi(ROI, 'flowline','99', 'dir', fdir, 'termini_path', tpath)));
ok = ok & check('no crossing trace on ANY candidate errors', ...
    raises(@() moc_flowline_from_roi(ROI, 'dir', fdir, 'termini_path', tpath, ...
        'trange', {'2030-01-01',''})));
FL3 = moc_flowline_from_roi(ROI, 'dir', fdir, 'termini_path', tpath, ...
    'trange', {'2030-01-01',''}, 'require_terminus', false);
ok = ok & check('require_terminus=false falls back to the un-split flowline', ...
    isempty(FL3.terminus_date) && abs(FL3.x(1)) < 1e-6);
ok = ok & check('moc_sample_obs_along_flowline rejects an absent band', ...
    raises(@() moc_sample_obs_along_flowline(O, FL2, 'band','nope')));

fprintf('\n%s\n', ternary(ok, 'ALL FLOWLINE-SHP CHECKS PASSED', 'SOME CHECKS FAILED'));
if nargout > 0, passed = ok; end

end

% ===========================================================================
function [fdir, tpath] = local_write_fixture(tmp, yLines, termDate, termX, badX, ...
    yLines2, stubX, stubY)
%LOCAL_WRITE_FIXTURE  Write the miniature flowline and terminus shapefiles.
fdir = fullfile(tmp, 'flowlines');
mkdir(fdir);
xs = linspace(0, 10000, 101);

% glacier0001: a fan of three lines, x from 0 (terminus) to 10 km inland
S = struct('Geometry',{},'X',{},'Y',{},'flowline',{});
ids = {'03','04','05'};
for k = 1:numel(yLines)
    S(k).Geometry = 'Line';
    S(k).X = xs;
    S(k).Y = yLines(k) * ones(size(xs));
    S(k).flowline = ids{k};
end
shapewrite(S, fullfile(fdir, 'glacier0001.shp'));

% glacier0003: a second fan whose centre line no trace crosses
S3b = struct('Geometry',{},'X',{},'Y',{},'flowline',{});
for k = 1:numel(yLines2)
    S3b(k).Geometry = 'Line';
    S3b(k).X = xs;
    S3b(k).Y = yLines2(k) * ones(size(xs));
    S3b(k).flowline = ids{k};
end
shapewrite(S3b, fullfile(fdir, 'glacier0003.shp'));

% glacier0002: far away, so the ROI must exclude it
S2 = struct('Geometry','Line','X',[0 10000],'Y',[50000 50000],'flowline','03');
shapewrite(S2, fullfile(fdir, 'glacier0002.shp'));

% an intermediate iteration that must be ignored unless asked for
S3 = struct('Geometry','Line','X',[0 10000],'Y',[0 0],'flowline','03');
shapewrite(S3, fullfile(fdir, 'glacier0001_iter01.shp'));

% not a flowline file at all — must be skipped, not crashed on
S4 = struct('Geometry','Line','X',[0 1],'Y',[0 1],'note','mask');
shapewrite(S4, fullfile(fdir, 'some_masks.shp'));

% termini: vertical traces crossing the whole fan, one per date
T = struct('Geometry',{},'X',{},'Y',{},'SourceDate',{},'Quality_Fl',{},'Glacier_ID',{});
for k = 1:numel(termX)
    T(k).Geometry = 'Line';
    T(k).X = [termX(k) termX(k)];
    T(k).Y = [-1000 500];
    T(k).SourceDate = termDate{k};
    T(k).Quality_Fl = 0;
    T(k).Glacier_ID = 1;
end
T(end+1) = struct('Geometry','Line','X',[badX badX],'Y',[-1000 500], ...
    'SourceDate','2022-01-01','Quality_Fl',1,'Glacier_ID',1);
% the stub for glacier0003: crosses only its outermost line
T(end+1) = struct('Geometry','Line','X',[stubX stubX],'Y',stubY, ...
    'SourceDate','2021-08-15','Quality_Fl',0,'Glacier_ID',3);
tpath = fullfile(tmp, 'termini.shp');
shapewrite(T, tpath);
end

% ---------------------------------------------------------------------------
function O = local_make_obs()
%LOCAL_MAKE_OBS  A gridded obs struct whose speed is exactly x (m/yr).
O = struct();
O.x = linspace(-1000, 11000, 121);
O.y = linspace(-500, 700, 13);
O.time = [2019 2020];
O.bands = {'vv'};
O.file = 'synthetic';
cfg = config_moc();
O.epsg = cfg.epsg;
O.vv = repmat(reshape(O.x, 1, [], 1), numel(O.y), 1, numel(O.time));
end

% ---------------------------------------------------------------------------
function ok = check(msg, cond)
%CHECK  Print a PASS/FAIL line and return the condition.
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%-4s] %s\n', tag, msg);
end

% ---------------------------------------------------------------------------
function tf = raises(fn)
%RAISES  True if calling fn throws.
try
    fn();
    tf = false;
catch
    tf = true;
end
end

% ---------------------------------------------------------------------------
function out = ternary(cond, a, b)
%TERNARY  Inline conditional for the summary line.
if cond, out = a; else, out = b; end
end
