function logfile = p1_log(action, name)
%P1_LOG  Start or stop a log file (diary) in the logs folder.
persistent current
logfile = '';
switch lower(action)
    case 'start'
        [~, ~, logdir] = p1_root();
        logfile = fullfile(logdir, sprintf('%s_%s.log', name, datestr(now, 'yyyy-mm-dd_HHMMSS')));
        diary off
        diary(logfile);
        current = logfile;
        fprintf('==== %s ====\n', name);
        fprintf('started : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf('MATLAB  : %s\n', version);
        fprintf('computer: %s\n\n', computer);
    case 'stop'
        fprintf('\nfinished: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        diary off
        logfile = current;
        current = '';
end
end
