%% smoke_test_moc.m
% Fast, no-plot sanity check for the ModelObsCompare toolkit.
%
% Loads one of each data type (ISSM model, gridded obs netCDF, observed
% flowline) and prints shapes / ranges, then runs targeted ORIENTATION checks:
%   - obs netCDF band field is (ny x nx x nt) and vv ~= hypot(vx,vy)
%   - flowline velocity is (npts x nt)
%   - InterpFromMeshToGrid output orientation matches (numel(y) x numel(x)),
%     verified independently by sampling the mesh at grid-cell centres
%
% Nothing is plotted. Each section is wrapped in try/catch so one failure does
% not stop the rest. A PASS/WARN/FAIL summary is printed at the end.
%
% Usage:
%   >> smoke_test_moc                 % auto-pick smallest model + first obs files
% Optionally set the variables in the CONFIG block below first.

clear; clc;
fprintf('==================  ModelObsCompare smoke test  ==================\n');

% ---- results accumulator (nested-friendly via a struct handle) -----------
R = struct('name',{},'status',{},'msg',{});
rec = @(n,s,m) struct('name',n,'status',s,'msg',m);

%% CONFIG -------------------------------------------------------------------
% Leave empty to auto-discover. Set to a filename (bare or full path) to force.
MODEL_FILE    = '';    % e.g. 'Model_NW_fric1_SnapshotInversion_2007.mat'
MODEL_SOLUTION= 'auto';% 'auto' | 'transient' | 'stressbalance'
OBS_NC_FILE   = '';    % e.g. 'Sentinel_Subset_Upernavik.nc'
OBS_FL_FILE   = '';    % e.g. 'jakobshavn.mat'

%% 0. Path setup ------------------------------------------------------------
try
    cfg = setup_moc();
    R(end+1) = rec('setup_moc','PASS','ISSM + toolbox on path');
catch ME
    R(end+1) = rec('setup_moc','FAIL',ME.message);
    fprintf(2,'setup_moc failed: %s\nAborting.\n',ME.message);
    return
end

%% 1. Model loader ----------------------------------------------------------
M = [];
try
    if isempty(MODEL_FILE)
        T = moc_list_models();
        assert(~isempty(T),'no Model_* files found in %s',moc_path('model_dir'));
        [~,ord] = sort(T.bytes,'ascend');      % smallest first = fastest to load
        MODEL_FILE = T.file{ord(1)};
    end
    fprintf('\n--- MODEL: %s ---\n', MODEL_FILE);
    M = moc_load_model(MODEL_FILE,'solution',MODEL_SOLUTION);

    nv = numel(M.x); ne = size(M.elements,1); nt = numel(M.time);
    prn('mesh.x',        M.x);
    prn('mesh.y',        M.y);
    fprintf('  %-14s size=%s  (expect ne x 3)\n','elements', mat2str(size(M.elements)));
    fprintf('  %-14s [%s]  (%d step(s))\n','time', rngstr(M.time), nt);
    prn('vel',           M.vel);
    fprintf('  %-14s size=%s  (expect nv x nt = %d x %d)\n','vel', mat2str(size(M.vel)), nv, nt);
    if ~isempty(M.bed), prn('bed', M.bed); end

    % orientation: velocity rows must equal #vertices
    if size(M.vel,1)==nv && size(M.elements,2)==3
        R(end+1)=rec('model.orient','PASS',sprintf('vel is nv x nt (%dx%d), elements nx3',nv,nt));
    else
        R(end+1)=rec('model.orient','FAIL',sprintf('vel size %s vs nv=%d, elements %s', ...
            mat2str(size(M.vel)),nv,mat2str(size(M.elements))));
    end
    % plausibility of speeds (m/yr)
    vmax = max(M.vel(:,end));
    if vmax>10 && vmax<1e5
        R(end+1)=rec('model.units','PASS',sprintf('max speed %.0f m/yr looks like m/yr',vmax));
    else
        R(end+1)=rec('model.units','WARN',sprintf('max speed %.3g — check units (m/yr expected)',vmax));
    end
catch ME
    R(end+1)=rec('model.load','FAIL',ME.message);
    fprintf(2,'  model load failed: %s\n',ME.message);
end

%% 2. Observed gridded netCDF ----------------------------------------------
O = [];
try
    if isempty(OBS_NC_FILE)
        ncdir = moc_path('obs_netcdf_dir');
        d = dir(fullfile(ncdir,'*Subset*.nc'));
        if isempty(d), d = dir(fullfile(ncdir,'*.nc')); end
        assert(~isempty(d),'no .nc files found in %s',ncdir);
        [~, order]=sort([d.bytes],'ascend');   % smallest netCDF first
        OBS_NC_FILE = fullfile(d(order(1)).folder,d(order(1)).name);
    end
    fprintf('\n--- OBS netCDF: %s ---\n', OBS_NC_FILE);
    O = moc_load_obs_netcdf(OBS_NC_FILE,'bands',{'vv','vx','vy'});

    nx=numel(O.x); ny=numel(O.y); nto=numel(O.time);
    fprintf('  %-14s n=%d  [%s] m\n','x', nx, rngstr(O.x));
    fprintf('  %-14s n=%d  [%s] m\n','y', ny, rngstr(O.y));
    fprintf('  %-14s n=%d  [%s] (decimal yr)\n','time', nto, rngstr(O.time));
    fprintf('  %-14s EPSG:%g\n','crs', O.epsg);
    for b = {'vv','vx','vy'}
        if isfield(O,b{1})
            fprintf('  %-14s size=%s  [%s]\n', ['band ' b{1}], ...
                mat2str(size(O.(b{1}))), rngstr(O.(b{1})(:)));
        end
    end

    % orientation: each band must be (ny x nx x nt)
    okdim = isfield(O,'vv') && isequal(size(O.vv),[ny nx nto]);
    if okdim
        R(end+1)=rec('obsnc.orient','PASS',sprintf('vv is ny x nx x nt (%dx%dx%d)',ny,nx,nto));
    else
        R(end+1)=rec('obsnc.orient','FAIL',sprintf('vv size %s != [%d %d %d]', ...
            mat2str(size(O.vv)),ny,nx,nto));
    end
    % band-mapping check: speed vv should ~ hypot(vx,vy)
    if all(isfield(O,{'vv','vx','vy'}))
        f = min(1,size(O.vv,3));   % first frame
        A = O.vv(:,:,f); B = hypot(O.vx(:,:,f),O.vy(:,:,f));
        good = isfinite(A)&isfinite(B)&(A>1);
        if nnz(good)>50
            relerr = median(abs(A(good)-B(good))./A(good));
            if relerr < 0.05
                R(end+1)=rec('obsnc.bands','PASS',sprintf('vv ~= hypot(vx,vy), median rel.err %.1f%%',100*relerr));
            else
                R(end+1)=rec('obsnc.bands','WARN',sprintf('vv vs hypot(vx,vy) median rel.err %.1f%% — check band order %s', ...
                    100*relerr, strjoin(cfg.obs_bands,',')));
            end
        else
            R(end+1)=rec('obsnc.bands','WARN','too few valid pixels to check band mapping');
        end
    end
catch ME
    R(end+1)=rec('obsnc.load','FAIL',ME.message);
    fprintf(2,'  obs netCDF load failed: %s\n',ME.message);
end

%% 3. Observed flowline -----------------------------------------------------
F = [];
try
    if isempty(OBS_FL_FILE)
        fldir = moc_path('obs_flowline_dir');
        pref = fullfile(fldir,'jakobshavn.mat');
        if isfile(pref)
            OBS_FL_FILE = pref;
        else
            d = dir(fullfile(fldir,'*.mat'));
            assert(~isempty(d),'no flowline .mat files in %s',fldir);
            OBS_FL_FILE = fullfile(d(1).folder,d(1).name);
        end
    end
    fprintf('\n--- OBS flowline: %s ---\n', OBS_FL_FILE);
    F = moc_load_obs_flowline(OBS_FL_FILE);

    npts=numel(F.x); ntf=numel(F.time);
    fprintf('  %-14s ''%s''  source=%s\n','name', F.name, F.source);
    prn('x', F.x); prn('y', F.y);
    fprintf('  %-14s [%s] m  (monotonic: %d)\n','d', rngstr(F.d), all(diff(F.d)>=0));
    fprintf('  %-14s n=%d  [%s] (decimal yr)\n','time', ntf, rngstr(F.time));
    fprintf('  %-14s size=%s  [%s]  (expect npts x nt = %d x %d)\n','vel', ...
        mat2str(size(F.vel)), rngstr(F.vel(:)), npts, ntf);

    if isequal(size(F.vel),[npts ntf])
        R(end+1)=rec('flowline.orient','PASS',sprintf('vel is npts x nt (%dx%d)',npts,ntf));
    else
        R(end+1)=rec('flowline.orient','FAIL',sprintf('vel size %s != [%d %d]', ...
            mat2str(size(F.vel)),npts,ntf));
    end
catch ME
    R(end+1)=rec('flowline.load','FAIL',ME.message);
    fprintf(2,'  flowline load failed: %s\n',ME.message);
end

%% 4. ROI + regrid orientation (the important one) --------------------------
try
    assert(~isempty(M),'need a model for the regrid check');
    if ~isempty(O)
        roi = moc_roi_from_obs(O,'margin',1000);
    else
        % fall back to a small box around the model's busiest area
        roi = struct('name','smoke','type','rectangle', ...
            'x',[],'y',[],'bbox',[min(M.x) min(M.x)+2e4 min(M.y) min(M.y)+2e4],'epsg',cfg.epsg);
        roi.x=[roi.bbox(1) roi.bbox(2) roi.bbox(2) roi.bbox(1) roi.bbox(1)]';
        roi.y=[roi.bbox(3) roi.bbox(3) roi.bbox(4) roi.bbox(4) roi.bbox(3)]';
    end
    fprintf('\n--- REGRID over ROI bbox [%.0f %.0f %.0f %.0f] @ %g m ---\n', ...
        roi.bbox, cfg.grid_res);

    G = moc_model_to_grid(M, roi, 'field','vel', 'res',cfg.grid_res);
    nxg=numel(G.x); nyg=numel(G.y);
    fprintf('  grid.x n=%d, grid.y n=%d\n', nxg, nyg);
    fprintf('  grid size=%s  (expect ny x nx = %d x %d)\n', mat2str(size(G.grid)), nyg, nxg);
    fprintf('  valid (non-NaN) fraction = %.1f%%\n', 100*mean(isfinite(G.grid(:))));

    % Shape check
    if isequal(size(G.grid),[nyg nxg])
        R(end+1)=rec('regrid.shape','PASS',sprintf('grid is ny x nx (%dx%d)',nyg,nxg));
    elseif isequal(size(G.grid),[nxg nyg])
        R(end+1)=rec('regrid.shape','FAIL','grid is TRANSPOSED (nx x ny) — fix InterpFromMeshToGrid usage in moc_model_to_grid');
    else
        R(end+1)=rec('regrid.shape','FAIL',sprintf('unexpected grid size %s',mat2str(size(G.grid))));
    end

    % Independent orientation check: sample the mesh at a few grid-cell centres
    % via InterpFromMesh2d and compare to the grid value at the same (r,c).
    fin = find(isfinite(G.grid));
    if numel(fin) >= 5
        pick = fin(round(linspace(1,numel(fin),min(8,numel(fin)))));
        xs = G.X(pick); ys = G.Y(pick);   % coordinates of those grid cells
        ref = InterpFromMesh2d(M.elements, M.x, M.y, M.vel(:,end), xs, ys);
        gv  = G.grid(pick);
        d = abs(gv(:)-ref(:)); denom = max(abs(ref(:)),1);
        relerr = median(d./denom,'omitnan');
        if relerr < 0.02
            R(end+1)=rec('regrid.orient','PASS',sprintf('grid(r,c) matches mesh sample at (X,Y), median rel.err %.2f%%',100*relerr));
        else
            R(end+1)=rec('regrid.orient','FAIL',sprintf('grid(r,c) vs mesh sample median rel.err %.1f%% — X/Y or row/col mismatch',100*relerr));
        end
    else
        R(end+1)=rec('regrid.orient','WARN','too few valid grid cells inside ROI to check orientation (ROI over ice?)');
    end
catch ME
    R(end+1)=rec('regrid','FAIL',ME.message);
    fprintf(2,'  regrid check failed: %s\n',ME.message);
end

%% 5. Point interpolation over time ----------------------------------------
try
    assert(~isempty(M),'need a model');
    if ~isempty(F)
        px = F.x(min(20,numel(F.x))); py = F.y(min(20,numel(F.y)));
        where = 'flowline pt 20';
    else
        px = mean(M.x); py = mean(M.y); where = 'mesh centroid';
    end
    P = moc_model_at_points(M, px, py, 'field','vel');
    fprintf('\n--- POINT interp (%s) ---\n', where);
    fprintf('  val size=%s  [%s]  (expect 1 x nt)\n', mat2str(size(P.val)), rngstr(P.val(:)));
    if size(P.val,1)==1
        R(end+1)=rec('points.orient','PASS','val is 1 x nt');
    else
        R(end+1)=rec('points.orient','WARN',sprintf('val size %s',mat2str(size(P.val))));
    end
catch ME
    R(end+1)=rec('points','FAIL',ME.message);
    fprintf(2,'  point interp failed: %s\n',ME.message);
end

%% Summary ------------------------------------------------------------------
fprintf('\n==============================  SUMMARY  ==========================\n');
nP=0;nW=0;nF=0;
for i=1:numel(R)
    switch R(i).status
        case 'PASS', nP=nP+1;
        case 'WARN', nW=nW+1;
        case 'FAIL', nF=nF+1;
    end
    fprintf('  [%-4s] %-18s %s\n', R(i).status, R(i).name, R(i).msg);
end
fprintf('-------------------------------------------------------------------\n');
fprintf('  %d PASS   %d WARN   %d FAIL\n', nP, nW, nF);
if nF==0
    fprintf('  Smoke test OK%s.\n', ternary(nW>0,' (with warnings to review)',''));
else
    fprintf(2,'  %d FAILURE(S) — review orientation notes above.\n', nF);
end

%% local helpers ------------------------------------------------------------
function prn(name, v)
%PRN  Print min/max/#NaN for a numeric array field.
v = v(:);
fprintf('  %-14s [%s]  nNaN=%d/%d\n', name, rngstr(v), nnz(isnan(v)), numel(v));
end

function s = rngstr(v)
%RNGSTR  Compact "min .. max" string, NaN-safe.
v = v(:);
if isempty(v), s='empty'; return; end
s = sprintf('%.4g .. %.4g', min(v,[],'omitnan'), max(v,[],'omitnan'));
end

function out = ternary(c,a,b)
if c, out=a; else, out=b; end
end
