function out = moc_configure(varargin)
%MOC_CONFIGURE  Show, set, or interactively fill in the toolkit's data paths.
%
%   The MATLAB counterpart of `python -m moc_py` on the Python side, writing
%   the same shared config file, so a path set here is also found there.
%
%   Usage:
%      moc_configure()                                % show every setting
%      moc_configure('-set')                          % prompt for missing ones
%      moc_configure('-set','-all')                   % re-prompt for all of them
%      moc_configure('model_dir','/data/models', ...) % set one or more
%      moc_configure('-save')                         % write current paths to disk
%      moc_configure('-clear')                        % drop this session's overrides
%      s = moc_configure('-status')                   % struct array, no printing
%
%   Setting a path with a name/value pair applies for this MATLAB session and
%   takes precedence over the environment and the config file. Follow it with
%   moc_configure('-save') to make it permanent.
%
%   Output:
%      out : for '-save', the config file written (char); for '-status', a
%            struct array with .name/.value/.source; otherwise omitted.
%
%   See also MOC_PATH, MOC_PATH_SETTINGS, MOC_CONFIG_FILE, CONFIG_MOC.

S = moc_path_settings();
args = varargin;

doSet = false; doAll = false; doShow = true;

% --- flags ---------------------------------------------------------------
keep = true(1, numel(args));
for i = 1:numel(args)
    if ~(ischar(args{i}) || isstring(args{i})), continue; end
    switch lower(char(args{i}))
        case '-set',    doSet = true;  keep(i) = false;
        case '-all',    doAll = true;  keep(i) = false;
        case '-clear',  moc_path_state('clear'); keep(i) = false;
        case '-status'
            out = local_status(S); return;
        case '-save'
            out = local_save(S); keep(i) = false; doShow = false;
    end
end
args = args(keep);

% --- name/value pairs ----------------------------------------------------
if mod(numel(args), 2) ~= 0
    error('moc_configure:badArgs', ...
        'Path settings must be given as name/value pairs.');
end
for i = 1:2:numel(args)
    name = char(args{i});
    if ~any(strcmp({S.name}, name))
        error('moc_configure:unknown', ...
            'Unknown path setting "%s"; known: %s.', name, strjoin({S.name}, ', '));
    end
    moc_path_state('set', name, moc_clean_path(args{i+1}));
    doShow = doShow && isempty(args);   % setting explicitly: stay quiet
end
if ~isempty(args), doShow = false; end

% --- interactive fill-in -------------------------------------------------
if doSet
    st = local_status(S);
    for k = 1:numel(S)
        needs = strcmp(st(k).source, 'unset') || ...
                (doAll && ~strcmp(st(k).source, 'built-in'));
        if needs
            try
                moc_path(S(k).name);         % prompts and saves
            catch
                % cancelled or not answerable — leave it unset
            end
        end
    end
    doShow = true;
end

% --- report --------------------------------------------------------------
if doShow
    st = local_status(S);
    fprintf('\nModelObsCompare paths   (config file: %s)\n', moc_config_file());
    fprintf('%-18s %-22s %s\n', 'setting', 'source', 'value');
    fprintf('%-18s %-22s %s\n', repmat('-',1,18), repmat('-',1,22), repmat('-',1,40));
    missing = 0;
    for k = 1:numel(st)
        v = st(k).value;
        if isempty(v), v = '(not set)'; missing = missing + 1; end
        fprintf('%-18s %-22s %s\n', st(k).name, st(k).source, v);
    end
    if missing > 0
        fprintf(['\n%d setting(s) not configured — run moc_configure(''-set''), ' ...
                 'or set the MOC_* environment variables.\n'], missing);
    end
end

end

% ===========================================================================
function st = local_status(S)
%LOCAL_STATUS  Where each setting currently comes from, without prompting.
ov = moc_path_state('get');
[~, data] = moc_config_file();
st = struct('name',{},'value',{},'source',{});
for k = 1:numel(S)
    name = S(k).name;
    value = ''; source = 'unset';
    if isfield(ov, name) && ~isempty(ov.(name))
        value = ov.(name); source = 'session';
    elseif ~isempty(getenv(S(k).env))
        value = moc_clean_path(getenv(S(k).env)); source = S(k).env;
    elseif isfield(data, name) && ~isempty(data.(name))
        value = moc_clean_path(char(data.(name))); source = 'config file';
    else
        try
            value = moc_path(name, 'prompt', false); source = 'built-in';
        catch
            value = ''; source = 'unset';
        end
    end
    st(end+1) = struct('name', name, 'value', value, 'source', source); %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------
function f = local_save(S)
%LOCAL_SAVE  Write the currently set paths to the shared config file.
%   Entries for settings that are not set are left as they were, so saving
%   never drops a path another session (or the Python side) put there.
[f, data] = moc_config_file();
st = local_status(S);
for k = 1:numel(st)
    if ~isempty(st(k).value) && any(strcmp(st(k).source, {'session','config file'}))
        data.(st(k).name) = st(k).value;
    end
end
d = fileparts(f);
if ~isfolder(d), mkdir(d); end
fid = fopen(f, 'w');
if fid < 0, error('moc_configure:saveFailed', 'Could not write %s', f); end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(data, 'PrettyPrint', true));
end
