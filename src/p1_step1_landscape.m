function p1_step1_landscape()
%P1_STEP1_LANDSCAPE  Cost map over the search box: 20 x 20 log-spaced (c3, c4)
%points with the three scenario costs. Shows where the optimum is, how flat the
%valley is and where the loop diverges.
%   results/p1_landscape.csv (c3, c4, J, J1, J2, J3; 1e6 = diverged), p1_landscape.json
S = p1_settings();  [~, res] = p1_root();
c3 = logspace(log10(S.lb(1)), log10(S.ub(1)), S.grid_n);
c4 = logspace(log10(S.lb(2)), log10(S.ub(2)), S.grid_n);
[A, B] = ndgrid(c3, c4);  X = [A(:) B(:)];  n = size(X, 1);
J = zeros(n, 1);  Js = NaN(n, numel(S.scen));  scen = S.scen;
t0 = tic;
parfor i = 1:n
    [Ji, js] = p1_cost(X(i, :), scen);
    J(i) = Ji;  row = NaN(1, numel(scen));  row(1:numel(js)) = js;  Js(i, :) = row;
end
fprintf('%d points in %.1f min\n', n, toc(t0)/60);
fid = fopen(fullfile(res, 'p1_landscape.csv'), 'w');
fprintf(fid, 'c3,c4,J,J1,J2,J3\n');
for i = 1:n, fprintf(fid, '%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', X(i, :), J(i), Js(i, :)); end
fclose(fid);
ok = J < 1e6;
[Jmin, k] = min(J);
fprintf('Grid minimum: c3 = %.2f, c4 = %.2f, J = %.6f\n', X(k, 1), X(k, 2), Jmin);
fprintf('Diverged points: %d of %d (%.1f %%)\n', sum(~ok), n, 100*mean(~ok));
for tol = [0.01 0.05 0.10 0.5]
    fprintf('Points within %4.0f %% of the grid minimum: %d (%.1f %%)\n', 100*tol, sum(J <= (1 + tol)*Jmin), 100*mean(J <= (1 + tol)*Jmin));
end
[Jh, Jsh] = p1_cost(S.hand, S.scen);
fprintf('Trial-and-error gains: J = %.4f (%.0f times the grid minimum)\n', Jh, Jh/Jmin);
p1_save_json('p1_landscape', struct('grid_c3', c3, 'grid_c4', c4, 'J_min', Jmin, 'x_min', X(k, :), ...
    'n_points', n, 'n_diverged', sum(~ok), 'within_1pct', sum(J <= 1.01*Jmin), 'within_5pct', sum(J <= 1.05*Jmin), ...
    'within_10pct', sum(J <= 1.10*Jmin), 'J_hand', Jh, 'Js_hand', Jsh));
end
