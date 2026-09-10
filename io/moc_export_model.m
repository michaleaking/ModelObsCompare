function ncpath = moc_export_model(M, ncpath, varargin)
%MOC_EXPORT_MODEL  Write an ISSM model's mesh + velocity to a neutral netCDF.
%
%   Bridges the MATLAB and Python sides: an ISSM `model` object cannot be read
%   in Python, but its mesh and velocity time series are just arrays. This
%   writes them to a language-neutral netCDF4 (HDF5) file that `moc_py.model`
%   reads with h5py — no ISSM or MATLAB needed on the Python side.
%
%   Usage:
%      moc_export_model(M, 'nw_fric1.nc')                 % M from moc_load_model
%      moc_export_model('Model_NW_..._2007.mat', out.nc)  % path -> loads first
%      moc_export_model(M, out, 'fields',{'vel','vx','vy'}, 'trange',[2018 2020])
%
%   Input:
%      M      : struct from MOC_LOAD_MODEL, or a path/filename to load.
%      ncpath : output netCDF path.
%   Name/value options:
%      'fields' : (cellstr) velocity fields to write (default {'vel','vx','vy'})
%      'trange' : (1x2)     [t0 t1] decimal-year window (transient only)
%      'overwrite' : (logical) delete an existing output first (default true)
%
%   Output:
%      ncpath : the written file path.
%
%   Schema (dims: vertex, elem, corner=3, time):
%      x(vertex), y(vertex)              [m]        mesh vertex coordinates
%      elements(elem, corner)  int32     1-based    triangulation (attr index_base=1)
%      time(time)              double    dec. year  (NaN for a snapshot)
%      vel/vx/vy(vertex, time) double    m/yr       requested velocity fields
%      surface(vertex, time)   double    m          if present in M
%      bed(vertex)             double    m          if present in M
%   Global attrs: epsg, solution, source_file, index_base, Conventions.

p = inputParser;
p.addParameter('fields', {'vel','vx','vy'}, @iscell);
p.addParameter('trange', [], @(v) isempty(v)||numel(v)==2);
p.addParameter('overwrite', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

% Accept a path/filename instead of a loaded struct
if ischar(M) || isstring(M)
    M = moc_load_model(char(M), 'trange', opt.trange);
elseif ~isempty(opt.trange) && ~all(isnan(M.time))
    keep = M.time >= opt.trange(1) & M.time <= opt.trange(2);
    M = local_subset_time(M, keep);
end

if opt.overwrite && isfile(ncpath)
    delete(ncpath);
end

nv = numel(M.x);
ne = size(M.elements,1);
nt = numel(M.time);

% --- create dimensions/variables (netcdf4 so h5py/xarray can read it) ------
nccreate(ncpath,'x','Dimensions',{'vertex',nv},'Datatype','double','Format','netcdf4');
ncwrite(ncpath,'x',M.x(:));
nccreate(ncpath,'y','Dimensions',{'vertex',nv},'Datatype','double');
ncwrite(ncpath,'y',M.y(:));

nccreate(ncpath,'elements','Dimensions',{'elem',ne,'corner',3},'Datatype','int32');
ncwrite(ncpath,'elements',int32(M.elements));

nccreate(ncpath,'time','Dimensions',{'time',nt},'Datatype','double');
ncwrite(ncpath,'time',M.time(:));

for f = opt.fields(:)'
    name = f{1};
    if ~isfield(M, name) || isempty(M.(name))
        warning('moc_export_model:missingField','Field "%s" not in M; skipping.', name);
        continue
    end
    nccreate(ncpath,name,'Dimensions',{'vertex',nv,'time',nt},'Datatype','double');
    ncwrite(ncpath,name,M.(name));
    ncwriteatt(ncpath,name,'units','m/yr');
end

if isfield(M,'surface') && ~isempty(M.surface) && size(M.surface,2)==nt
    nccreate(ncpath,'surface','Dimensions',{'vertex',nv,'time',nt},'Datatype','double');
    ncwrite(ncpath,'surface',M.surface);
    ncwriteatt(ncpath,'surface','units','m');
end
if isfield(M,'bed') && ~isempty(M.bed)
    nccreate(ncpath,'bed','Dimensions',{'vertex',nv},'Datatype','double');
    ncwrite(ncpath,'bed',M.bed(:));
    ncwriteatt(ncpath,'bed','units','m');
end

% --- coordinate/units attrs ------------------------------------------------
ncwriteatt(ncpath,'x','units','m'); ncwriteatt(ncpath,'y','units','m');
ncwriteatt(ncpath,'time','units','decimal_year');
ncwriteatt(ncpath,'elements','index_base',int32(1));
ncwriteatt(ncpath,'elements','long_name','triangulation vertex indices (1-based)');

% --- global attributes -----------------------------------------------------
cfg = config_moc();
ncwriteatt(ncpath,'/','epsg',int32(cfg.epsg));
ncwriteatt(ncpath,'/','index_base',int32(1));
sol = 'unknown'; if isfield(M,'solution'), sol = M.solution; end
ncwriteatt(ncpath,'/','solution',sol);
src = ''; if isfield(M,'file'), src = M.file; end
ncwriteatt(ncpath,'/','source_file',src);
ncwriteatt(ncpath,'/','Conventions','moc_py-model-1.0');
ncwriteatt(ncpath,'/','created_by','moc_export_model.m');

fprintf('Exported model -> %s\n  %d vertices, %d elements, %d time step(s), fields: %s\n', ...
    ncpath, nv, ne, nt, strjoin(opt.fields,','));

end

% ---------------------------------------------------------------------------
function M = local_subset_time(M, keep)
%LOCAL_SUBSET_TIME  Keep a subset of transient time steps in the model struct.
M.time = M.time(keep);
for f = {'vx','vy','vel','surface'}
    if isfield(M,f{1}) && ~isempty(M.(f{1})) && size(M.(f{1}),2)>1
        M.(f{1}) = M.(f{1})(:, keep);
    end
end
end
