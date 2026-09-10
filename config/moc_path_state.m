function out = moc_path_state(action, name, value)
%MOC_PATH_STATE  Session-level path overrides shared by moc_path/moc_configure.
%
%   Internal helper: holds the paths set with MOC_CONFIGURE for the life of the
%   MATLAB session (they take precedence over the environment and the config
%   file, and are not written to disk unless you ask MOC_CONFIGURE to save).
%
%   Usage:
%      s = moc_path_state('get');
%      s = moc_path_state('set', 'model_dir', '/path/to/models');
%      s = moc_path_state('clear');
%
%   Input:
%      action : (char) 'get' | 'set' | 'clear'
%      name   : (char) setting name (for 'set')
%      value  : (char) path        (for 'set')
%   Output:
%      out : (struct) the current overrides

persistent OVERRIDES
if isempty(OVERRIDES), OVERRIDES = struct(); end

switch lower(action)
    case 'get'
        % nothing to do
    case 'set'
        OVERRIDES.(name) = char(value);
    case 'clear'
        if nargin > 1 && ~isempty(name)
            if isfield(OVERRIDES, name), OVERRIDES = rmfield(OVERRIDES, name); end
        else
            OVERRIDES = struct();
        end
    otherwise
        error('moc_path_state:badAction', 'action must be get, set or clear.');
end
out = OVERRIDES;

end
