function run_paper1(what, algs)
%RUN_PAPER1  Reproduce the results of the paper.
%
%   Backstepping control of the shank angle of a belt-driven underactuated
%   prosthetic leg (one pelvis-mounted motor, passive hip). The two gains are
%   tuned by six metaheuristics on a noisy 1 kHz digital loop and compared with
%   trial-and-error gains, simple alternatives and a set of robustness tests.
%
%   run_paper1                    all steps 0-5 (several hours on 8 workers)
%   run_paper1('check')           step 0: quick checks and a run-time estimate (~5 min)
%   run_paper1('landscape')       step 1: cost map over (c3, c4)               (~12 min)
%   run_paper1('tune')            step 2: six optimisers x 15 runs             (~2-3 h)
%   run_paper1('tune', {'DE'})    step 2 for the listed optimisers only
%   run_paper1('stats')           step 3: statistics of the optimisers         (~1 min)
%   run_paper1('evaluate')        step 4: tuned gains vs trial and error       (~15 min)
%   run_paper1('extra')           step 5: alternatives and practical limits    (~25 min)
%
%   The results folder of the repository already holds every finished run, so
%   'stats', 'evaluate' and 'extra' can be repeated without tuning again; the
%   tuning skips runs whose result file exists. Each step writes a log file to
%   the logs folder. The Parallel Computing Toolbox is used when available
%   (parfor falls back to a normal loop otherwise).
if nargin < 1, what = 'all'; end
if nargin < 2, algs = {}; end

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(genpath(fullfile(root, 'lib')));
pool = [];
if exist('gcp', 'file'), pool = gcp('nocreate'); end   % workers need the same path
if ~isempty(pool)
    pctRunOnAll(sprintf('addpath(''%s''); addpath(genpath(''%s''));', ...
        fullfile(root, 'src'), fullfile(root, 'lib')));
end

steps = {'p1_step0_check', 'p1_step1_landscape', 'p1_step2_tuning', ...
         'p1_step3_stats', 'p1_step4_evaluate', 'p1_step5_revision'};
switch lower(what)
    case 'all',       todo = steps;
    case 'check',     todo = steps(1);
    case 'landscape', todo = steps(2);
    case 'tune',      todo = steps(3);
    case 'stats',     todo = steps(4);
    case 'evaluate',  todo = steps(5);
    case 'extra',     todo = steps(6);
    otherwise
        error('unknown option "%s" (all, check, landscape, tune, stats, evaluate, extra)', what);
end

status = cell(numel(todo), 3);
for i = 1:numel(todo)
    f = p1_log('start', todo{i});
    t0 = tic;
    try
        if strcmp(todo{i}, 'p1_step2_tuning')
            feval(todo{i}, algs);
        else
            feval(todo{i});
        end
        status(i, :) = {todo{i}, 'OK', ''};
    catch err
        fprintf('\n*** step failed: %s\n', err.message);
        for s = 1:numel(err.stack)
            fprintf('    in %s (line %d)\n', err.stack(s).name, err.stack(s).line);
        end
        status(i, :) = {todo{i}, 'FAILED', err.message};
        if strcmp(todo{i}, 'p1_step0_check')
            fprintf('\n*** the check failed, so the other steps were skipped\n');
            p1_log('stop');
            status = status(1:i, :);
            break
        end
    end
    fprintf('\nstep time: %.1f min\n', toc(t0)/60);
    p1_log('stop');
    fprintf('log saved: %s\n', f);
end

fprintf('\nsummary\n');
for i = 1:size(status, 1)
    fprintf('  %-20s %-7s %s\n', status{i, :});
end
end
