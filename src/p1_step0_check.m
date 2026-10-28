function p1_step0_check()
%P1_STEP0_CHECK  Quick checks before the long runs:
%  (a) the simulator reproduces the reference cost of an earlier study
%      (J = 0.860736 for the gains [91.99 259.13] with T = 40 s);
%  (b) cost of the trial-and-error gains c3 = c4 = 9;
%  (c) every optimiser uses exactly the budget on a cheap test function;
%  (d) time of one cost evaluation and of one population -> run-time estimate.
S = p1_settings();
fprintf('Box: c3 in [%g, %g], c4 in [%g, %g]; population %d, %d generations, budget %d evaluations\n', ...
        S.lb(1), S.ub(1), S.lb(2), S.ub(2), S.n_pop, S.n_iter, S.budget);
fprintf('Scenarios of the robust cost (shank mass scale, noise factor, delay ms):\n');
for k = 1:numel(S.scen)
    fprintf('   %d: %-24s %.1f  %g  %g\n', k, S.scen_names{k}, S.scen(k).mass, S.scen(k).noise, S.scen(k).delay);
end
% (a)
pref = [91.99051419868184 259.1313409162851];
[Jref, Jsref] = p1_cost(pref, S.scen, 40);
fprintf('\n(a) reference gains of the earlier study (T = 40 s): J = %.6f (expected 0.860736), scenarios %s\n', Jref, mat2str(Jsref, 6));
if abs(Jref - 0.860736061962433) > 1e-6
    error('The cost differs from the reference value (%.9f instead of 0.860736062) - check the MATLAB path.', Jref);
end
t0 = tic;  [J20, Js20] = p1_cost(pref, S.scen);  t1 = toc(t0);
fprintf('    the same gains with T = %g s (cost used in this paper): J = %.6f, scenarios %s\n', S.T_tune, J20, mat2str(Js20, 6));
% (b)
[Jh, Jsh] = p1_cost(S.hand, S.scen);
fprintf('(b) trial-and-error gains c3 = c4 = 9 (T = %g s): J = %.6f, scenarios %s\n', S.T_tune, Jh, mat2str(Jsh, 6));
% (c)
fsph = @(X) sum(((X - [250 750])./[500 1500]).^2, 2);
for a = 1:numel(S.algs)
    R = S.funs{a}(fsph, S.lb, S.ub, S.n_pop, S.n_iter, 1);
    fprintf('(c) %-7s test function: best x = %s, J = %.2e, evaluations = %d\n', S.algs{a}, mat2str(R.x, 4), R.J, R.evals(end));
    if R.evals(end) ~= S.budget
        error('%s used %d evaluations instead of %d', S.algs{a}, R.evals(end), S.budget);
    end
end
% (d)
X = S.lb + rand(S.n_pop, 2).*(S.ub - S.lb);
t0 = tic;  Jp = p1_popcost(X, S.scen);  tp = toc(t0);
fprintf('\n(d) one cost evaluation (serial): %.1f s; one population of %d (parallel): %.1f s (%d diverged)\n', ...
        t1, S.n_pop, tp, sum(Jp >= 1e6));
nrun = numel(S.algs)*numel(S.seeds);
fprintf('    estimated time: one optimiser run %.0f min, all %d runs %.1f h, cost map %.0f min\n', ...
        tp*(S.n_iter + 1)/60, nrun, nrun*tp*(S.n_iter + 1)/3600, tp*S.grid_n^2/S.n_pop/60);
p1_save_json('p1_step0_check', struct('J_ref', Jref, 'J_ref_T', J20, 'T_tune', S.T_tune, 'J_hand', Jh, 'Js_hand', Jsh, ...
    't_eval_s', t1, 't_pop_s', tp, 'budget', S.budget, 'n_pop', S.n_pop, 'n_iter', S.n_iter));
end
