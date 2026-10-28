function p1_step2_tuning(algs)
%P1_STEP2_TUNING  Tune [c3 c4] with the six optimisers (runs and seeds from
%p1_settings). Same population, budget, box, cost and seeds for every optimiser;
%run k uses rng(k).
%   Each finished run is saved to results/p1_run_<alg>_<seed>.mat (x, J, Jhist,
%   evals, minutes), so an interrupted tuning continues with the next run.
%   p1_step2_tuning({'DE', 'PSO'}) tunes only the listed optimisers.
S = p1_settings();  [~, res] = p1_root();
if nargin < 1 || isempty(algs), algs = S.algs; end
scen = S.scen;
f = @(X) p1_popcost(X, scen);
t_all = tic;
for a = 1:numel(algs)
    ia = find(strcmp(S.algs, algs{a}));
    if isempty(ia), error('unknown optimiser "%s"', algs{a}); end
    fprintf('\n==== %s ====\n', S.algs{ia});
    for s = S.seeds
        fr = fullfile(res, sprintf('p1_run_%s_%d.mat', tag(S.algs{ia}), s));
        if exist(fr, 'file')
            Q = load(fr);
            fprintf('run %2d already done: J = %.6f  c3 = %.3f  c4 = %.3f\n', s, Q.J, Q.x(1), Q.x(2));
            continue
        end
        t0 = tic;
        R = S.funs{ia}(f, S.lb, S.ub, S.n_pop, S.n_iter, s);
        minutes = toc(t0)/60;
        x = R.x;  J = R.J;  Jhist = R.Jhist(:);  evals = R.evals(:);  alg = S.algs{ia};  seed = s; %#ok<NASGU>
        save(fr, 'x', 'J', 'Jhist', 'evals', 'minutes', 'alg', 'seed');
        fprintf('run %2d: J = %.6f  c3 = %.3f  c4 = %.3f  (%d evaluations, %.1f min)\n', ...
                s, J, x(1), x(2), evals(end), minutes);
    end
end
fprintf('\ntotal tuning time of this call: %.1f h\n', toc(t_all)/3600);
end

function t = tag(name)
t = strrep(strrep(name, '+', '_'), '-', '');
end
