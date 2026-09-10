function [p, data] = moc_config_file()
%MOC_CONFIG_FILE  Location and contents of the shared path config file.
%
%   The same file the Python side reads and writes (moc_py/config.py), so a
%   path answered once is known to both languages. Location:
%      $MOC_CONFIG                        if set, else
%      $XDG_CONFIG_HOME/moc/paths.json    if set, else
%      ~/.config/moc/paths.json
%
%   Usage:
%      p = moc_config_file();
%      [p, data] = moc_config_file();
%
%   Output:
%      p    : (char)   full path to the config file (may not exist yet)
%      data : (struct) its decoded contents (empty struct if absent/unreadable)

override = getenv('MOC_CONFIG');
if ~isempty(override)
    p = char(override);
else
    xdg = getenv('XDG_CONFIG_HOME');
    if isempty(xdg)
        xdg = fullfile(moc_home_dir(), '.config');
    end
    p = fullfile(char(xdg), 'moc', 'paths.json');
end

data = struct();
if nargout > 1 && isfile(p)
    try
        data = jsondecode(fileread(p));
    catch
        data = struct();   % unreadable or malformed: treat as empty
    end
end

end
