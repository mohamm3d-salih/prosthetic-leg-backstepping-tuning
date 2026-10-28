function [root, res, logs] = p1_root()
%P1_ROOT  Repository root and its results and logs folders (created if missing).
root = fileparts(fileparts(mfilename('fullpath')));
res = fullfile(root, 'results');   if ~exist(res, 'dir'), mkdir(res); end
logs = fullfile(root, 'logs');     if ~exist(logs, 'dir'), mkdir(logs); end
end
