function p1_step3_stats()
%P1_STEP3_STATS  Statistical comparison of the optimisers.
%No toolbox needed; every test is written out below.
%   * descriptive statistics, bootstrap 95 % CI of the median (10000 resamples)
%   * Kruskal-Wallis test (Friedman with the seed as block is printed for reference only)
%   * 15 pairwise exact Wilcoxon rank-sum tests, Holm correction, Vargha-Delaney A12
%   * equal spread: Brown-Forsythe (all), Fligner-Killeen (DE vs each), Holm
%   * runs within 0.1 % and 1 % of the best cost, Fisher exact tests vs DE
%   * evaluations needed to come within 1 % and 0.1 % of the best cost
%     (a run that never gets there counts as budget + 1)
%   * median best-so-far cost at 25/50/75/100 % of the budget
%   * spread of the tuned gains between runs
%   results/p1_stats.json, p1_convergence.csv
S = p1_settings();  [~, res] = p1_root();
k = numel(S.algs);  n = numel(S.seeds);
J = cell(1, k);  X = cell(1, k);  H = cell(1, k);  E = cell(1, k);  T = cell(1, k);
fid = fopen(fullfile(res, 'p1_convergence.csv'), 'w');
fprintf(fid, 'alg,seed,evals,Jbest\n');
for a = 1:k
    J{a} = zeros(n, 1);  X{a} = zeros(n, 2);  T{a} = zeros(n, 1);  H{a} = cell(n, 1);  E{a} = cell(n, 1);
    for i = 1:n
        fr = fullfile(res, sprintf('p1_run_%s_%d.mat', tag(S.algs{a}), S.seeds(i)));
        if ~exist(fr, 'file'), error('missing run %s seed %d - run the tuning first', S.algs{a}, S.seeds(i)); end
        Q = load(fr);
        J{a}(i) = Q.J;  X{a}(i, :) = Q.x;  T{a}(i) = Q.minutes;
        H{a}{i} = Q.Jhist(:);  E{a}{i} = Q.evals(:);
        for q = 1:numel(Q.evals), fprintf(fid, '%s,%d,%d,%.10g\n', S.algs{a}, S.seeds(i), Q.evals(q), Q.Jhist(q)); end
    end
end
fclose(fid);
Jbest = min(cellfun(@min, J));
[~, ab] = min(cellfun(@min, J));  [~, ib] = min(J{ab});  xbest = X{ab}(ib, :);
fprintf('Best cost found by any optimiser: J = %.6f (%s, run %d): c3 = %.4f, c4 = %.4f\n', ...
        Jbest, S.algs{ab}, S.seeds(ib), xbest(1), xbest(2));

% ---- descriptive statistics
fprintf('\n=== Final cost of every run (%d runs per optimiser, %d evaluations per run) ===\n', n, S.budget);
fprintf('%-7s %10s %10s %10s %10s %10s %10s  %-22s %8s\n', 'alg', 'best', 'median', 'mean', 'SD', 'worst', 'gap med %', '95% CI of median', 'min/run');
D = struct();
for a = 1:k
    ci = boot_ci(J{a});
    D(a).alg = S.algs{a};  D(a).best = min(J{a});  D(a).median = median(J{a});  D(a).mean = mean(J{a});
    D(a).sd = std(J{a});  D(a).worst = max(J{a});  D(a).ci = ci;  D(a).minutes = mean(T{a});
    D(a).gap_median_pct = 100*(median(J{a})/Jbest - 1);  D(a).gap_worst_pct = 100*(max(J{a})/Jbest - 1);
    D(a).x = X{a};  D(a).J = J{a};
    D(a).c3_sd = std(X{a}(:, 1));  D(a).c4_sd = std(X{a}(:, 2));
    D(a).c3_range = [min(X{a}(:, 1)) max(X{a}(:, 1))];  D(a).c4_range = [min(X{a}(:, 2)) max(X{a}(:, 2))];
    fprintf('%-7s %10.6f %10.6f %10.6f %10.2e %10.6f %10.3f  [%.6f, %.6f] %8.1f\n', S.algs{a}, D(a).best, D(a).median, ...
            D(a).mean, D(a).sd, D(a).worst, D(a).gap_median_pct, ci(1), ci(2), D(a).minutes);
end
fprintf('\nTuned gains (range over the %d runs):\n', n);
for a = 1:k
    fprintf('%-7s c3 %8.3f - %8.3f (SD %.3f)   c4 %8.3f - %8.3f (SD %.3f)\n', S.algs{a}, D(a).c3_range, D(a).c3_sd, D(a).c4_range, D(a).c4_sd);
end

% ---- omnibus tests
[Hkw, pkw] = kruskal_wallis(J);
[Qf, pf, mr] = friedman_test(cell2mat(J));
fprintf('\nKruskal-Wallis: H = %.3f, p = %.4g\n', Hkw, pkw);
fprintf('Friedman (blocks = seeds): chi2 = %.3f, p = %.4g\n', Qf, pf);
fprintf('Mean ranks over the seeds (1 = best):');
[~, o2] = sort(mr);
for a = o2, fprintf('  %s %.2f', S.algs{a}, mr(a)); end
fprintf('\n');

% ---- pairwise Wilcoxon + Holm + A12
pairs = nchoosek(1:k, 2);  np = size(pairs, 1);
p = zeros(np, 1);  A = p;
for q = 1:np
    p(q) = ranksum_exact(J{pairs(q, 1)}, J{pairs(q, 2)});
    A(q) = a12(J{pairs(q, 1)}, J{pairs(q, 2)});
end
ph = holm(p);
fprintf('\nPairwise Wilcoxon rank-sum (exact), Holm-corrected; A12 = P(first gives lower J than second)\n');
[~, o3] = sort(p);
PW = struct('a', {}, 'b', {}, 'p', {}, 'p_holm', {}, 'A12', {});
for q = o3'
    fprintf('  %-7s vs %-7s p = %.4f  p_Holm = %.4f  A12 = %.2f (%s)\n', S.algs{pairs(q, 1)}, ...
            S.algs{pairs(q, 2)}, p(q), ph(q), A(q), a12_label(A(q)));
    PW(end+1) = struct('a', S.algs{pairs(q, 1)}, 'b', S.algs{pairs(q, 2)}, 'p', p(q), 'p_holm', ph(q), 'A12', A(q)); %#ok<AGROW>
end

% ---- dispersion
iref = find(strcmp(S.algs, 'DE'));  others = setdiff(1:k, iref);
[Wbf, pbf] = brown_forsythe(J);
fprintf('\nBrown-Forsythe equal-spread test: W = %.3f, p = %.4g\n', Wbf, pbf);
pf2 = zeros(numel(others), 1);
for q = 1:numel(others), [~, pf2(q)] = fligner_killeen({J{iref}, J{others(q)}}); end
pfh = holm(pf2);
FK = struct('alg', {}, 'p', {}, 'p_holm', {}, 'sd_ratio', {});
for q = 1:numel(others)
    sr = std(J{others(q)})/max(std(J{iref}), eps);
    fprintf('  DE vs %-7s Fligner-Killeen p = %.4f  Holm %.4f  SD ratio %.1f\n', S.algs{others(q)}, pf2(q), pfh(q), sr);
    FK(end+1) = struct('alg', S.algs{others(q)}, 'p', pf2(q), 'p_holm', pfh(q), 'sd_ratio', sr); %#ok<AGROW>
end

% ---- success rates
SR = struct();
for tol = [0.001 0.01]
    thr = Jbest*(1 + tol);  s = cellfun(@(x) sum(x <= thr), J);
    pfi = zeros(numel(others), 1);
    for q = 1:numel(others), pfi(q) = fisher_exact([s(iref) n - s(iref); s(others(q)) n - s(others(q))]); end
    pfih = holm(pfi);
    fprintf('\nsuccess (J <= %.6f, within %.1f %%): ', thr, 100*tol);
    for a = 1:k, fprintf('%s %d/%d  ', S.algs{a}, s(a), n); end
    fprintf('\n   Fisher vs DE (p / Holm):');
    for q = 1:numel(others), fprintf('  %s %.3f/%.3f', S.algs{others(q)}, pfi(q), pfih(q)); end
    fprintf('\n');
    SR.(sprintf('tol_%s', strrep(sprintf('%g', 100*tol), '.', 'p'))) = struct('count', s, 'fisher_p', pfi', 'fisher_holm', pfih');
end

% ---- efficiency: evaluations to reach within 1 % and 0.1 % of the best cost
EF = struct();
for tol = [0.01 0.001]
    ev = cell(1, k);
    for a = 1:k
        ev{a} = zeros(n, 1);
        for i = 1:n
            r = find(H{a}{i} <= Jbest*(1 + tol), 1, 'first');
            if isempty(r), ev{a}(i) = S.budget + 1; else, ev{a}(i) = E{a}{i}(r); end
        end
    end
    [He, pe] = kruskal_wallis(ev);
    nm = sprintf('tol_%s', strrep(sprintf('%g', 100*tol), '.', 'p'));
    fprintf('\nEvaluations to come within %.1f %% of the best cost (median [min, max]; %d = never):\n', 100*tol, S.budget + 1);
    for a = 1:k
        fprintf('  %-7s %6.0f [%d, %d]  reached in %d/%d runs\n', S.algs{a}, median(ev{a}), min(ev{a}), max(ev{a}), sum(ev{a} <= S.budget), n);
    end
    fprintf('  Kruskal-Wallis: H = %.3f, p = %.4g\n', He, pe);
    EF.(nm) = struct('evals', {ev}, 'median', cellfun(@median, ev), 'reached', cellfun(@(v) sum(v <= S.budget), ev), 'kw_p', pe, 'kw_H', He);
end

% ---- budget sensitivity
fprintf('\nMedian best-so-far cost at 25/50/75/100 %% of the budget (%d evaluations):\n', S.budget);
BS = zeros(k, 4);  fr = [0.25 0.5 0.75 1];
for a = 1:k
    fprintf('  %-7s', S.algs{a});
    for j = 1:4
        v = zeros(n, 1);
        for i = 1:n
            r = find(E{a}{i} <= fr(j)*S.budget + 1e-9, 1, 'last');
            if isempty(r), v(i) = NaN; else, v(i) = H{a}{i}(r); end
        end
        BS(a, j) = median(v);
        fprintf('  %3.0f%%: %.6f', 100*fr(j), BS(a, j));
    end
    fprintf('\n');
end

p1_save_json('p1_stats', struct('algs', {S.algs}, 'budget', S.budget, 'n_runs', n, 'J_best', Jbest, ...
    'x_best', xbest, 'alg_best', S.algs{ab}, 'seed_best', S.seeds(ib), 'desc', D, ...
    'kruskal_H', Hkw, 'kruskal_p', pkw, 'friedman_chi2', Qf, 'friedman_p', pf, 'mean_ranks', mr, ...
    'pairwise', PW, 'brown_forsythe_W', Wbf, 'brown_forsythe_p', pbf, 'fligner', FK, ...
    'success', SR, 'efficiency', EF, 'budget_median', BS, 'budget_fracs', fr));
end

% ======================================================================
function t = tag(name)
t = strrep(strrep(name, '+', '_'), '-', '');
end

function ci = boot_ci(x)
rs = RandStream('mt19937ar', 'Seed', 0);  nb = 10000;  m = zeros(nb, 1);  nx = numel(x);
for b = 1:nb, m(b) = median(x(ceil(nx*rand(rs, nx, 1)))); end
m = sort(m);  ci = [m(round(0.025*nb)) m(round(0.975*nb))];
end

function [r, tc] = tied_ranks(v)
[s, i] = sort(v(:));  n = numel(v);  r = zeros(n, 1);  tc = 0;  a = 1;
while a <= n
    b = a;
    while b < n && s(b + 1) == s(a), b = b + 1; end
    r(i(a:b)) = (a + b)/2;  t = b - a + 1;  tc = tc + t^3 - t;  a = b + 1;
end
end

function p = chi2_upper(x, df)
p = gammainc(x/2, df/2, 'upper');
end

function [Hs, p] = kruskal_wallis(G)
x = cell2mat(G(:));  N = numel(x);  [r, tc] = tied_ranks(x);
Hs = 0;  a = 0;
for i = 1:numel(G)
    ni = numel(G{i});  Hs = Hs + sum(r(a + 1:a + ni))^2/ni;  a = a + ni;
end
Hs = 12/(N*(N + 1))*Hs - 3*(N + 1);
if tc < N^3 - N, Hs = Hs/(1 - tc/(N^3 - N)); else, Hs = 0; end
p = chi2_upper(Hs, numel(G) - 1);
end

function [Q, p, meanrank] = friedman_test(X)
[n, k] = size(X);  R = zeros(n, k);  tcs = 0;
for i = 1:n, [R(i, :), tc] = tied_ranks(X(i, :)); tcs = tcs + tc; end
Rj = sum(R, 1);
Q = 12/(n*k*(k + 1))*sum(Rj.^2) - 3*n*(k + 1);
if tcs < n*(k^3 - k), Q = Q/(1 - tcs/(n*(k^3 - k))); else, Q = 0; end
p = chi2_upper(Q, k - 1);  meanrank = Rj/n;
end

function p = ranksum_exact(x, y)
% two-sided exact Wilcoxon rank-sum test: full enumeration for small samples (5 runs),
% exact rank-sum distribution by dynamic programming without ties (15 runs),
% normal approximation with tie correction otherwise
x = x(:);  y = y(:);  nx = numel(x);  nt = nx + numel(y);
r = tied_ranks([x; y]);  W = sum(r(1:nx));
if nchoosek(nt, nx) <= 2e5
    Wall = sum(reshape(r(nchoosek(1:nt, nx)), [], nx), 2);
    p = min(1, 2*min(mean(Wall <= W + 1e-9), mean(Wall >= W - 1e-9)));
elseif all(r == round(r))
    smax = sum(nt - nx + 1:nt);  Cn = zeros(nx + 1, smax + 1);  Cn(1, 1) = 1;
    for k = 1:nt
        for j = min(k, nx):-1:1
            Cn(j + 1, k + 1:end) = Cn(j + 1, k + 1:end) + Cn(j, 1:end - k);
        end
    end
    pr = Cn(nx + 1, :)/sum(Cn(nx + 1, :));  s = 0:smax;  w = round(W);
    p = min(1, 2*min(sum(pr(s <= w)), sum(pr(s >= w))));
else
    ny = nt - nx;  mu = nx*(nt + 1)/2;  u = unique(r);  t = arrayfun(@(v) sum(r == v), u);
    s2 = nx*ny/12*((nt + 1) - sum(t.^3 - t)/(nt*(nt - 1)));
    z = (W - mu - 0.5*sign(W - mu))/sqrt(s2);
    p = min(1, erfc(abs(z)/sqrt(2)));
end
end

function a = a12(x, y)
a = (sum(sum(x(:) < y(:)')) + 0.5*sum(sum(x(:) == y(:)')))/(numel(x)*numel(y));
end

function s = a12_label(a)
d = abs(a - 0.5);
if d < 0.06, s = 'negligible'; elseif d < 0.14, s = 'small'; elseif d < 0.21, s = 'medium'; else, s = 'large'; end
end

function q = holm(p)
p = p(:);  m = numel(p);  [ps, o] = sort(p);  q = zeros(m, 1);  run = 0;
for i = 1:m, run = max(run, (m - i + 1)*ps(i)); q(o(i)) = min(1, run); end
end

function [W, p] = brown_forsythe(G)
k = numel(G);  Z = cell(1, k);
for i = 1:k, Z{i} = abs(G{i}(:) - median(G{i})); end
z = cell2mat(Z(:));  N = numel(z);  zm = mean(z);
num = 0;  den = 0;
for i = 1:k
    ni = numel(Z{i});  num = num + ni*(mean(Z{i}) - zm)^2;  den = den + sum((Z{i} - mean(Z{i})).^2);
end
d1 = k - 1;  d2 = N - k;
if den == 0, W = 0; p = 1; return, end
W = (num/d1)/(den/d2);
p = betainc(d2/(d2 + d1*W), d2/2, d1/2);
end

function [X2, p] = fligner_killeen(G)
k = numel(G);  Z = cell(1, k);
for i = 1:k, Z{i} = abs(G{i}(:) - median(G{i})); end
z = cell2mat(Z(:));  N = numel(z);  r = tied_ranks(z);
a = -sqrt(2)*erfcinv(2*(0.5 + r/(2*(N + 1))));
abar = mean(a);  v = var(a);  X2 = 0;  s = 0;
if v == 0, X2 = 0; p = 1; return, end
for i = 1:k
    ni = numel(Z{i});  X2 = X2 + ni*(mean(a(s + 1:s + ni)) - abar)^2;  s = s + ni;
end
X2 = X2/v;  p = chi2_upper(X2, k - 1);
end

function p = fisher_exact(t)
a = t(1, 1);  r1 = sum(t(1, :));  c1 = sum(t(:, 1));  N = sum(t(:));
lo = max(0, r1 + c1 - N);  hi = min(r1, c1);
lp = @(x) gammaln(r1 + 1) + gammaln(N - r1 + 1) + gammaln(c1 + 1) + gammaln(N - c1 + 1) ...
     - gammaln(N + 1) - gammaln(x + 1) - gammaln(r1 - x + 1) - gammaln(c1 - x + 1) - gammaln(N - r1 - c1 + x + 1);
P = exp(arrayfun(lp, lo:hi));  p0 = exp(lp(a));
p = min(1, sum(P(P <= p0*(1 + 1e-7))));
end
