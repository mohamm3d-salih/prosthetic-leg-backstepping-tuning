function p1_step5_revision()
%P1_STEP5_REVISION  Additional analyses (alternatives to metaheuristics and
%practical limits). Needs the results of steps 3 and 4.
%   (A) designs found without a metaheuristic, same budget (104 evaluations):
%       analytical rule (delay and noise limits), systematic manual tuning
%       (equal gains, doubling + golden section) and Nelder-Mead from c3 = c4 = 9
%   (B) share of steady-state samples in which the torque saturates
%   (C) optimum natural frequency against the torque weight, a torque-rate
%       weight and the delay of scenario J2
%   (D) actuator lag and belt compliance (with the reflected motor inertia)
%   (E) velocity noise alone (gyroscope level up to a differentiated encoder)
%   (F) gait-like references: knee-flexion error and passive hip damping
%   (G) every final design re-evaluated on five fresh noise realisations
%   results/p1_revision.json, p1_scan.csv, p1_baseline_hist.csv
S = p1_settings();  S0 = b2_settings();  M = b2_model();  [root, res] = p1_root();  deg = pi/180;
st = jsondecode(fileread(fullfile(res, 'p1_stats.json')));
E4 = jsondecode(fileread(fullfile(res, 'p1_evaluate.json')));
pbest = E4.gains(E4.ibest, :);
Jbest = st.J_best;
fprintf('Best tuned design (step 4): c3 = %.4f, c4 = %.4f, J = %.6f\n', pbest, Jbest);
RV = struct('Jbest', Jbest, 'pbest', pbest);

% =====================================================================================
% (A) baselines with the budget of the optimisers
% =====================================================================================
fprintf('\n(A) baselines with the same budget (%d evaluations)\n', S.budget);
NP = M.NP;  beta_min = 1/NP(3);  tau = 3e-3;  phi = pi/6;
wn_delay = phi/(tau + S0.Ts/2);                         % delay adds at most 30 deg of lag at wn
tq_noise = @(w) sqrt((w.^2*S0.sig_th).^2 + (2*w*S0.sig_w).^2)/beta_min;   % torque noise SD, zeta = 1
wn_noise = fzero(@(w) tq_noise(w) - 0.2*S.umax, [10 1000]); % noise torque at most 20 % of the limit
wn_an = min(wn_delay, wn_noise);  c_an = sqrt(wn_an^2 - 1);
fprintf('   analytical: delay rule wn = %.1f rad/s, noise rule wn = %.1f rad/s -> wn = %.1f, c3 = c4 = %.2f\n', ...
        wn_delay, wn_noise, wn_an, c_an);
RV.analytical = struct('wn_delay', wn_delay, 'wn_noise', wn_noise, 'wn', wn_an, 'c', c_an, ...
                       'phi_deg', phi/deg, 'tau_ms', 1e3*tau, 'noise_frac', 0.2, 'beta_min', beta_min);

% systematic manual tuning: equal gains, doubling then golden section on log(c)
[pm, Jm, hm] = manual_tuning(S);
fprintf('   manual (equal gains): c3 = c4 = %.2f, J = %.6f, %d evaluations\n', pm(1), Jm, numel(hm));
% Nelder-Mead in log coordinates from the trial-and-error gains
[pn, Jn, hn] = nelder_mead(S);
fprintf('   Nelder-Mead: c3 = %.2f, c4 = %.2f, J = %.6f, %d evaluations\n', pn, Jn, numel(hn));
fid = fopen(fullfile(res, 'p1_baseline_hist.csv'), 'w');
fprintf(fid, 'method,eval,J,Jbest_so_far\n');
for k = 1:numel(hm), fprintf(fid, 'manual,%d,%.10g,%.10g\n', k, hm(k), min(hm(1:k))); end
for k = 1:numel(hn), fprintf(fid, 'nelder,%d,%.10g,%.10g\n', k, hn(k), min(hn(1:k))); end
fclose(fid);
e1 = @(h) first_below(h, 1.01*Jbest);  e01 = @(h) first_below(h, 1.001*Jbest);
RV.manual = struct('p', pm, 'J', Jm, 'evals', numel(hm), 'e1', e1(hm), 'e01', e01(hm));
RV.nelder = struct('p', pn, 'J', Jn, 'evals', numel(hn), 'e1', e1(hn), 'e01', e01(hn));

% designs compared in the rest of this step
D = struct('name', {'hand', 'best', 'analytical', 'manual', 'nelder'}, ...
           'p', {S.hand, pbest, [c_an c_an], pm, pn});
nd = numel(D);
for i = 1:nd
    D(i).C = b2_params_to_controller('BSR', D(i).p);
    [D(i).J, D(i).Js] = p1_cost(D(i).p, S.scen);
    D(i).wn = sqrt(1 + prod(D(i).p));  D(i).zeta = sum(D(i).p)/(2*D(i).wn);
    D(i).gap = 100*(D(i).J/Jbest - 1);
end
fprintf('\n   %-11s %9s %9s %8s %7s %11s %9s\n', 'design', 'c3', 'c4', 'wn', 'zeta', 'J', 'gap %');
for i = 1:nd
    fprintf('   %-11s %9.2f %9.2f %8.1f %7.3f %11.6f %9.2f\n', D(i).name, D(i).p, D(i).wn, D(i).zeta, D(i).J, D(i).gap);
end

% nominal test leg, 20 noise realisations, for every design
mets = {'RMSE_ss', 'ess', 'urms', 'umax', 'du_rms', 'hip_pp', 'RMSE_meas'};
Cs = {};  op = {};  sd = {};
for i = 1:nd
    for s = S.eval_seeds, Cs{end+1} = D(i).C; op{end+1} = {}; sd{end+1} = s; end %#ok<AGROW>
end
Rn = n2_run_jobs(Cs, op, sd);
ns = numel(S.eval_seeds);
fprintf('\n   nominal test leg, mean over %d realisations\n   %-11s', ns, 'design');  fprintf('%13s', mets{:});  fprintf('\n');
for i = 1:nd
    V = NaN(ns, numel(mets));
    for s = 1:ns
        r = Rn{(i - 1)*ns + s};
        if r.ok, for m = 1:numel(mets), V(s, m) = r.(mets{m}); end, end
    end
    fprintf('   %-11s', D(i).name);  fprintf('%13.4g', mean(V, 1));  fprintf('\n');
    for m = 1:numel(mets), D(i).nominal.(mets{m}) = struct('mean', mean(V(:, m)), 'sd', std(V(:, m))); end
end

% Monte Carlo (300 random legs) for the analytical and manual designs
smp = dlmread(fullfile(root, 'data', 'mc_samples.csv'), ',', 1, 0);
N = size(smp, 1);  mcd = [3 4];
Cs = {};  op = {};  sd = {};
for n = 1:N
    m = smp(n, :);
    Pn = M.plant('mass_scale', m(1), 'b1', 0.5*m(2), 'b2', 0.3*m(2), 'd0', m(3), 'dw', m(4));
    o = {'P', Pn, 'sig_th', 1e-3*m(5), 'sig_w', 1e-2*m(5), 'delay', m(6), 'x0', [m(8)*deg; 0; (10 + m(7))*deg; 0]};
    for i = mcd, Cs{end+1} = D(i).C; op{end+1} = o; sd{end+1} = m(9); end %#ok<AGROW>
end
Rm = n2_run_jobs(Cs, op, sd);
mc_old = dlmread(fullfile(res, 'p1_mc.csv'), ',', 1, 0);      % case, hand, best (step 4)
fprintf('\n   Monte Carlo, %d random legs: median / 95th percentile of the steady RMS error (deg)\n', N);
fprintf('   %-11s %10.4f %10.4f\n', 'hand', median(mc_old(:, 2)), prctile_(mc_old(:, 2), 95));
fprintf('   %-11s %10.4f %10.4f\n', 'best', median(mc_old(:, 3)), prctile_(mc_old(:, 3), 95));
for q = 1:numel(mcd)
    v = NaN(N, 1);
    for n = 1:N
        r = Rm{(n - 1)*numel(mcd) + q};
        if r.ok, v(n) = r.RMSE_ss; end
    end
    vb = mc_old(:, 3);  ok = isfinite(v) & isfinite(vb);
    D(mcd(q)).mc = struct('failed', sum(~isfinite(v)), 'median', median(v(isfinite(v))), ...
        'p95', prctile_(v(isfinite(v)), 95), 'ratio_to_best_median', median(v(ok)./vb(ok)), ...
        'p_vs_best', signrank_p(v(ok) - vb(ok)), 'n_lower_than_best', sum(v(ok) < vb(ok)));
    fprintf('   %-11s %10.4f %10.4f   (failed %d; median ratio to best %.3f; lower than best in %d/%d, p = %.2g)\n', ...
        D(mcd(q)).name, D(mcd(q)).mc.median, D(mcd(q)).mc.p95, D(mcd(q)).mc.failed, ...
        D(mcd(q)).mc.ratio_to_best_median, D(mcd(q)).mc.n_lower_than_best, sum(ok), D(mcd(q)).mc.p_vs_best);
end
RV.designs = rmfield(D, 'C');

% =====================================================================================
% (B) torque saturation in steady tracking (nominal leg, 20 realisations)
% =====================================================================================
fprintf('\n(B) share of samples with |uc| > umax (nominal leg, 20 realisations)\n');
sel = [1 2 3];  jobs = zeros(0, 2);
for i = sel, for s = S.eval_seeds, jobs(end+1, :) = [i s]; end, end %#ok<AGROW>
nj = size(jobs, 1);  sat_ss = zeros(nj, 1);  sat_all = zeros(nj, 1);  CC = {D.C};
parfor j = 1:nj
    out = n2_simulate(CC{jobs(j, 1)}, 'seed', jobs(j, 2), 'out_dt', 1e-3);
    lim = CC{jobs(j, 1)}.umax + 1e-9;
    last = out.t >= out.t(end) - 10;
    sat_ss(j) = 100*mean(abs(out.uc(last)) > lim);
    sat_all(j) = 100*mean(abs(out.uc) > lim);
end
RV.saturation = struct();
for i = sel
    w = jobs(:, 1) == i;
    RV.saturation.(D(i).name) = struct('steady_pct', mean(sat_ss(w)), 'steady_pct_max', max(sat_ss(w)), 'all_pct', mean(sat_all(w)));
    fprintf('   %-11s steady (last 10 s) %.3f %% (max %.3f %%), whole run %.3f %%\n', D(i).name, mean(sat_ss(w)), max(sat_ss(w)), mean(sat_all(w)));
end

% =====================================================================================
% (C) weight and delay sensitivity along the equal-gain line (zeta ~ 1)
% =====================================================================================
fprintf('\n(C) optimum natural frequency vs cost weights and the delay of scenario J2 (equal gains)\n');
wg = logspace(log10(30), log10(600), 25);
sc = struct('mass', {1.2, 1.2, 1.2, 1.2, 1.5}, 'noise', {1, 1, 1, 1, 2}, 'delay', {0, 1, 3, 5, 0});
jobs = zeros(0, 2);
for a = 1:numel(wg), for b = 1:numel(sc), jobs(end+1, :) = [a b]; end, end %#ok<AGROW>
nj = size(jobs, 1);  CMP = NaN(nj, 5);                   % A, B, Cd, rmse_last10, du_rms_last10
Rf = b2_reference('sine');
parfor j = 1:nj
    w = wg(jobs(j, 1));  c = sqrt(w^2 - 1);  q = sc(jobs(j, 2));
    C = b2_params_to_controller('BSR', [c c]);
    out = n2_simulate(C, 'T', S.T_tune, 'P', M.plant('mass_scale', q.mass), 'R', Rf, 'Ts', S0.Ts, 'nsub', S0.nsub, ...
        'sig_th', S0.sig_th*q.noise, 'sig_w', S0.sig_w*q.noise, 'delay', q.delay, 'seed', S0.noise_seed, 'out_dt', 1e-3);
    if out.ok
        CMP(j, :) = cost_parts(out, Rf, S0.Ts);
    end
end
fid = fopen(fullfile(res, 'p1_scan.csv'), 'w');
fprintf(fid, 'wn,scenario,A,B,Cd,rmse_last10,du_rms_last10\n');
for j = 1:nj, fprintf(fid, '%.6g,%d,%.10g,%.10g,%.10g,%.10g,%.10g\n', wg(jobs(j, 1)), jobs(j, 2), CMP(j, :)); end
fclose(fid);
A = reshape(CMP(:, 1), numel(sc), numel(wg))';  B = reshape(CMP(:, 2), numel(sc), numel(wg))';
Cd = reshape(CMP(:, 3), numel(sc), numel(wg))';  rm1 = reshape(CMP(:, 4), numel(sc), numel(wg))';
du1 = reshape(CMP(:, 5), numel(sc), numel(wg))';
wus = [1e-5 1e-4 1e-3 1e-2];  wdus = [0 1e-10 1e-9 1e-8];  dls = [2 3 4];   % columns of sc: delay 1, 3, 5 ms
dlv = [1 3 5];
SC = struct('wu', {}, 'wdu', {}, 'delay_ms', {}, 'wn_opt', {}, 'rmse', {}, 'du_rms', {});
fprintf('   %-8s %-8s %-6s %10s %12s %14s\n', 'w_u', 'w_du', 'delay', 'wn opt', 'RMSE (deg)', 'du_rms (N.m/s)');
for a = 1:numel(wus)
    for b = 1:numel(wdus)
        for d = 1:numel(dls)
            cols = [1 dls(d) 5];
            Jw = mean(A(:, cols) + wus(a)*B(:, cols) + wdus(b)*Cd(:, cols), 2);
            Jw(any(~isfinite(A(:, cols)), 2)) = Inf;
            [wo, k] = parabolic_min(log(wg(:)), Jw);
            SC(end+1) = struct('wu', wus(a), 'wdu', wdus(b), 'delay_ms', dlv(d), 'wn_opt', exp(wo), ...
                               'rmse', rm1(k, 1), 'du_rms', du1(k, 1)); %#ok<AGROW>
            if b == 1 || d == 2
                fprintf('   %-8.0e %-8.0e %3d ms %10.1f %12.4f %14.1f\n', wus(a), wdus(b), dlv(d), exp(wo), rm1(k, 1), du1(k, 1));
            end
        end
    end
end
RV.scan = struct('wn', wg, 'table', SC);

% =====================================================================================
% (D) actuator lag and belt compliance (3 realisations, nominal leg)
% =====================================================================================
fprintf('\n(D) actuator lag and belt compliance: steady RMS error (deg) / RMS torque rate (N.m/s)\n');
taus = [0 1 2 3 5];  kbs = [Inf 5000 2000 1000 500];  sel = [1 2 3];
Cs = {};  op = {};  lab = zeros(0, 3);
for v = taus, for i = sel, Cs{end+1} = D(i).C; op{end+1} = {'tau_act', v*1e-3}; lab(end+1, :) = [1 v i]; end, end %#ok<AGROW>
for v = kbs,  for i = sel, Cs{end+1} = D(i).C; op{end+1} = {'k_belt', v};      lab(end+1, :) = [2 v i]; end, end %#ok<AGROW>
Rd = n2_run_jobs(Cs, op, 1:3);
RV.actuator = struct('tau_ms', taus, 'k_belt', kbs);
fac = {'tau_act', 'k_belt'};  lev = {taus, kbs};
for f = 1:2
    fprintf('   %s\n', fac{f});
    tab = NaN(numel(lev{f}), numel(sel));  tdu = tab;
    for q = 1:numel(lev{f})
        for ii = 1:numel(sel)
            j = find(lab(:, 1) == f & lab(:, 2) == lev{f}(q) & lab(:, 3) == sel(ii));
            r = Rd{j};
            if r.ok, tab(q, ii) = r.RMSE_ss;  tdu(q, ii) = r.du_rms; end
        end
        fprintf('   %10g', lev{f}(q));  fprintf('   %s %9.4f / %8.1f', D(1).name, tab(q, 1), tdu(q, 1));
        fprintf('   %s %9.4f / %8.1f', D(2).name, tab(q, 2), tdu(q, 2));
        fprintf('   %s %9.4f / %8.1f\n', D(3).name, tab(q, 3), tdu(q, 3));
    end
    RV.actuator.(fac{f}) = struct('rmse', tab, 'du_rms', tdu);
end

% =====================================================================================
% (E) velocity noise alone
% =====================================================================================
fprintf('\n(E) velocity noise alone (angle noise 1e-3 rad): steady RMS error (deg) / RMS torque rate / RMS torque\n');
sws = [0.01 0.03 0.1 0.3 1.0 1.41];
Cs = {};  op = {};  lab = zeros(0, 2);
for v = sws, for i = sel, Cs{end+1} = D(i).C; op{end+1} = {'sig_w', v}; lab(end+1, :) = [v i]; end, end %#ok<AGROW>
Re = n2_run_jobs(Cs, op, 1:3);
tab = NaN(numel(sws), numel(sel), 4);
for q = 1:numel(sws)
    fprintf('   %6.2f rad/s', sws(q));
    for ii = 1:numel(sel)
        r = Re{find(lab(:, 1) == sws(q) & lab(:, 2) == sel(ii))};
        if r.ok, tab(q, ii, :) = [r.RMSE_ss r.du_rms r.urms r.umax]; end
        fprintf('   %s %8.4f / %8.1f / %6.2f', D(sel(ii)).name, tab(q, ii, 1:3));
    end
    fprintf('\n');
end
RV.velnoise = struct('sig_w', sws, 'rmse', tab(:, :, 1), 'du_rms', tab(:, :, 2), 'urms', tab(:, :, 3), 'umax', tab(:, :, 4));

% =====================================================================================
% (F) gait-like references: knee flexion and hip damping
% =====================================================================================
fprintf('\n(F) gait-like references (3 realisations): shank RMSE / knee-flexion RMSE / hip swing (deg)\n');
gp = [1.0 1.2 1.6];  b1s = [0.5 2 5];  sel = [1 2 3];
jobs = zeros(0, 4);
for g = 1:numel(gp), for b = 1:numel(b1s), for i = sel, for s = 1:3, jobs(end+1, :) = [g b i s]; end, end, end, end %#ok<AGROW>
nj = size(jobs, 1);  G = NaN(nj, 4);
parfor j = 1:nj
    jr = jobs(j, :);
    R = b2_reference('gait', gp(jr(1)));
    out = n2_simulate(CC{jr(3)}, 'R', R, 'x0', [0; 0; -5*pi/180; 0], 'P', M.plant('b1', b1s(jr(2))), ...
                      'seed', jr(4), 'out_dt', 1e-3);
    if out.ok
        G(j, :) = gait_parts(out, R);
    end
end
RV.gait = struct('periods', gp, 'b1', b1s);
for g = 1:numel(gp)
    for b = 1:numel(b1s)
        fprintf('   cycle %.1f s, hip damping b1 = %.1f N.m.s/rad\n', gp(g), b1s(b));
        for i = sel
            w = jobs(:, 1) == g & jobs(:, 2) == b & jobs(:, 3) == i;
            v = mean(G(w, :), 1);
            fprintf('      %-11s %8.4f %9.3f %8.2f   (knee flexion range %.1f deg)\n', D(i).name, v);
            RV.gait.(sprintf('g%d_b%d', g, b)).(D(i).name) = struct('shank_rmse', v(1), 'flex_rmse', v(2), 'hip_pp', v(3), 'flex_range', v(4));
        end
    end
end

% =====================================================================================
% (G) every final design on five fresh noise realisations
% =====================================================================================
fprintf('\n(G) final designs of all runs re-evaluated on 5 fresh noise realisations (seeds 201-205)\n');
k = numel(S.algs);  fresh_seeds = 201:205;
X = zeros(0, 2);  Jt = zeros(0, 1);  al = zeros(0, 1);
for a = 1:k
    n = size(st.desc(a).x, 1);
    X = [X; st.desc(a).x]; Jt = [Jt; st.desc(a).J(:)]; al = [al; a*ones(n, 1)]; %#ok<AGROW>
end
nx = size(X, 1);  jobs = zeros(0, 2);
for i = 1:nx, for s = fresh_seeds, jobs(end+1, :) = [i s]; end, end %#ok<AGROW>
nj = size(jobs, 1);  Jf = zeros(nj, 1);
parfor j = 1:nj
    Jf(j) = cost_seed(X(jobs(j, 1), :), jobs(j, 2), S, S0, M);
end
Jfm = zeros(nx, 1);
for i = 1:nx, Jfm(i) = mean(Jf(jobs(:, 1) == i)); end
rho = spearman(Jt, Jfm);
Gs = cell(1, k);  for a = 1:k, Gs{a} = Jfm(al == a); end
[~, pkw] = kruskal_wallis_(Gs);
spread = 100*(max(Jfm)/min(Jfm) - 1);
fprintf('   Spearman rank correlation between tuning cost and fresh-noise cost: %.3f\n', rho);
fprintf('   spread of the fresh-noise cost over all %d designs: %.3f %%; Kruskal-Wallis between optimisers p = %.3f\n', nx, spread, pkw);
fprintf('   %-7s %12s %12s\n', 'alg', 'median J', 'best J');
med = zeros(1, k);
for a = 1:k
    med(a) = median(Gs{a});
    fprintf('   %-7s %12.6f %12.6f\n', S.algs{a}, med(a), min(Gs{a}));
end
RV.fresh = struct('seeds', fresh_seeds, 'rho', rho, 'spread_pct', spread, 'p_kw', pkw, 'median', med, ...
                  'J_tune', Jt, 'J_fresh', Jfm, 'alg', al);

p1_save_json('p1_revision', RV);
end

% ======================================================================================
function [p, J, h] = manual_tuning(S)
% emulates a systematic hand tuning with equal gains: double c while the cost
% decreases, then golden-section search on log(c) inside the last bracket
budget = S.budget;  lo = S.lb(1);  hi = min(S.ub);
h = [];
f = @(c) p1_cost([c c], S.scen);
c = S.hand(1);  Jc = f(c);  h(end+1) = Jc;
while numel(h) < budget && 2*c <= hi
    J2 = f(2*c);  h(end+1) = J2; %#ok<AGROW>
    if J2 < Jc, c = 2*c;  Jc = J2; else, break, end
end
a = log(max(c/2, lo));  b = log(min(2*c, hi));  gr = (sqrt(5) - 1)/2;
x1 = b - gr*(b - a);  x2 = a + gr*(b - a);
if numel(h) < budget, f1 = f(exp(x1)); h(end+1) = f1; else, f1 = Inf; end
if numel(h) < budget, f2 = f(exp(x2)); h(end+1) = f2; else, f2 = Inf; end
while numel(h) < budget && (b - a) > log(1.005)
    if f1 < f2
        b = x2;  x2 = x1;  f2 = f1;  x1 = b - gr*(b - a);  f1 = f(exp(x1));  h(end+1) = f1; %#ok<AGROW>
    else
        a = x1;  x1 = x2;  f1 = f2;  x2 = a + gr*(b - a);  f2 = f(exp(x2));  h(end+1) = f2; %#ok<AGROW>
    end
end
cand = [c, exp(x1), exp(x2)];  Jcand = [Jc, f1, f2];
[J, i] = min(Jcand);  p = [cand(i) cand(i)];
h = h(:);
end

function [p, J, h] = nelder_mead(S)
% Nelder-Mead (fminsearch) in log coordinates, started from the trial-and-error gains,
% limited to the budget of the optimisers; points outside the box are clipped
lb = log(S.lb);  ub = log(S.ub);
hh = zeros(0, 1);  P = zeros(0, 2);
    function v = fobj(z)
        z = min(max(z(:)', lb), ub);
        v = p1_cost(exp(z), S.scen);
        hh(end+1, 1) = v;  P(end+1, :) = exp(z);
    end
opt = optimset('MaxFunEvals', S.budget, 'MaxIter', S.budget, 'TolX', 1e-4, 'TolFun', 1e-7, 'Display', 'off');
fminsearch(@fobj, log(S.hand), opt);
h = hh(1:min(end, S.budget));
[J, i] = min(h);  p = P(i, :);
end

function k = first_below(h, thr)
k = find(h(:) <= thr, 1);
if isempty(k), k = numel(h) + 1; end
end

function v = cost_parts(out, R, Ts)
% A = int t|e3| dt, B = int u^2 dt, Cd = int (du/dt)^2 dt, plus the steady RMS
% error (deg) and RMS torque rate (N.m/s) over the last 10 s
xd = arrayfun(@(tt) b2_ref_eval(tt, R), out.t);
e = out.x(:, 3) - xd;
A = trapz(out.t, out.t.*abs(e));  B = trapz(out.t, out.u.^2);
du = diff(out.u)/Ts;  Cd = sum(du.^2)*Ts;
last = out.t >= out.t(end) - 10;  dul = diff(out.u(last))/Ts;
v = [A, B, Cd, sqrt(mean((e(last)*180/pi).^2)), sqrt(mean(dul.^2))];
end

function v = gait_parts(out, R)
% shank RMSE, knee-flexion RMSE (flexion th1 - th2 against kappa = -x3d), hip swing,
% knee-flexion range; last 10 s of the run, degrees
xd = arrayfun(@(tt) b2_ref_eval(tt, R), out.t);
last = out.t >= out.t(end) - 10;
e = (out.x(last, 3) - xd(last))*180/pi;
flex = (out.x(last, 1) - out.x(last, 3))*180/pi;  kap = -xd(last)*180/pi;
hip = out.x(last, 1)*180/pi;
v = [sqrt(mean(e.^2)), sqrt(mean((flex - kap).^2)), max(hip) - min(hip), max(flex) - min(flex)];
end

function J = cost_seed(p, seed, S, S0, M)
% the tuning cost (13) with another noise realisation
Rf = b2_reference('sine');  C = b2_params_to_controller('BSR', p);  Js = zeros(1, numel(S.scen));
for k = 1:numel(S.scen)
    sc = S.scen(k);
    out = n2_simulate(C, 'T', S.T_tune, 'P', M.plant('mass_scale', sc.mass), 'R', Rf, ...
        'Ts', S0.Ts, 'nsub', S0.nsub, 'sig_th', S0.sig_th*sc.noise, 'sig_w', S0.sig_w*sc.noise, ...
        'delay', sc.delay, 'k_belt', sc.k_belt, 'seed', seed, 'out_dt', 0.02);
    Js(k) = b2_cost(out, Rf, S0.w_u);
end
if any(Js >= 1e6), J = 1e6; else, J = mean(Js); end
end

function [xo, k] = parabolic_min(x, y)
% minimum of y(x) on a grid, refined by a parabola through the three nearest points
[~, k] = min(y);
if k == 1 || k == numel(y) || ~all(isfinite(y(k-1:k+1))), xo = x(k); return, end
x1 = x(k-1); x2 = x(k); x3 = x(k+1);  y1 = y(k-1); y2 = y(k); y3 = y(k+1);
den = (x1 - x2)*(x1 - x3)*(x2 - x3);
a = (x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2))/den;
b = (x3^2*(y1 - y2) + x2^2*(y3 - y1) + x1^2*(y2 - y3))/den;
if a > 0, xo = -b/(2*a); else, xo = x2; end
xo = min(max(xo, x1), x3);
end

function r = spearman(a, b)
ra = tied_ranks_(a(:));  rb = tied_ranks_(b(:));
ra = ra - mean(ra);  rb = rb - mean(rb);
r = sum(ra.*rb)/sqrt(sum(ra.^2)*sum(rb.^2));
end

function r = tied_ranks_(v)
[~, o] = sort(v);  r = zeros(size(v));  r(o) = 1:numel(v);
u = unique(v);
for k = 1:numel(u)
    idx = v == u(k);
    if sum(idx) > 1, r(idx) = mean(r(idx)); end
end
end

function [H, p] = kruskal_wallis_(G)
x = vertcat(G{:});  n = numel(x);  r = tied_ranks_(x);  H = 0;  i0 = 0;
for a = 1:numel(G)
    m = numel(G{a});  H = H + sum(r(i0+1:i0+m))^2/m;  i0 = i0 + m;
end
H = 12/(n*(n + 1))*H - 3*(n + 1);
u = unique(x);  t = arrayfun(@(v) sum(x == v), u);
H = H/(1 - sum(t.^3 - t)/(n^3 - n));
p = gammainc(H/2, (numel(G) - 1)/2, 'upper');
end

function p = prctile_(v, q)
v = sort(v(:));  n = numel(v);
if n == 0, p = NaN; return, end
pos = 1 + (n - 1)*q/100;  lo = floor(pos);  hi = min(lo + 1, n);
p = v(lo) + (pos - lo)*(v(hi) - v(lo));
end

function p = signrank_p(d)
% two-sided Wilcoxon signed-rank test, normal approximation with tie correction
d = d(isfinite(d) & d ~= 0);  n = numel(d);
if n == 0, p = 1; return, end
r = tied_ranks_(abs(d));  W = sum(r(d > 0));
u = unique(abs(d));  t = arrayfun(@(v) sum(abs(d) == v), u);  tie = sum(t.^3 - t);
mu = n*(n + 1)/4;  s2 = n*(n + 1)*(2*n + 1)/24 - tie/48;
z = (W - mu - 0.5*sign(W - mu))/sqrt(s2);
p = min(1, erfc(abs(z)/sqrt(2)));
end
