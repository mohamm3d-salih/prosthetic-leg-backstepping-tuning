function p1_step4_evaluate()
%P1_STEP4_EVALUATE  Tuned gains against the trial-and-error gains (c3 = c4 = 9).
%Designs: trial and error and the best run of every optimiser.
%  (a) robust cost and its scenario costs; natural frequency and damping of
%      e3'' + (c3+c4) e3' + (1+c3 c4) e3 = 0
%  (b) test leg, 20 noise realisations (seeds 101-120): steady RMS and peak
%      error, RMS and peak torque, RMS torque rate, hip swing, IAE; paired exact
%      Wilcoxon signed-rank tests of the best tuned design against trial and error
%  (c) transients after large initial errors (3 realisations)
%  (d) one-factor robustness (shank mass, friction, disturbance amplitude and
%      frequency, noise, measurement delay 0-10 ms with delay margin)
%  (e) gait-like references (1.0, 1.2, 1.6 s cycles)
%  (f) Monte Carlo over 300 random legs and conditions (data/mc_samples.csv)
%  (g) time histories for the figures (test 1 and a 30 deg initial error)
%   results/p1_evaluate.json, p1_timehist_*.csv, p1_robustness.csv, p1_mc.csv
S = p1_settings();  [root, res] = p1_root();  M = b2_model();  deg = pi/180;
st = jsondecode(fileread(fullfile(res, 'p1_stats.json')));
names = [{'Hand'}, cellfun(@tag, S.algs, 'UniformOutput', false)];
labels = [{'Trial and error'}, S.algs];
P = zeros(numel(names), 2);  P(1, :) = S.hand;
for a = 1:numel(S.algs)
    d = st.desc(a);  [~, i] = min(d.J);  P(a + 1, :) = d.x(i, :);
end
[~, ib] = min(arrayfun(@(a) min(st.desc(a).J), 1:numel(S.algs)));  ibest = ib + 1;
fprintf('Best tuned design: %s, c3 = %.4f, c4 = %.4f\n', labels{ibest}, P(ibest, :));
nd = numel(names);  CC = cell(1, nd);
for i = 1:nd, CC{i} = b2_params_to_controller('BSR', P(i, :)); end
E = struct('names', {names}, 'labels', {labels}, 'gains', P, 'ibest', ibest);

% (a) robust cost and error dynamics -----------------------------------------
fprintf('\n(a) design            c3        c4       J (robust)  J1 nominal  J2 delay 3 ms  J3 heavy+noise   wn (rad/s)  zeta\n');
E.J = zeros(nd, 1);  E.Js = zeros(nd, 3);  E.wn = zeros(nd, 1);  E.zeta = zeros(nd, 1);  E.poles = zeros(nd, 2);
for i = 1:nd
    [E.J(i), js] = p1_cost(P(i, :), S.scen);  E.Js(i, 1:numel(js)) = js;
    c3 = P(i, 1);  c4 = P(i, 2);
    E.wn(i) = sqrt(1 + c3*c4);  E.zeta(i) = (c3 + c4)/(2*E.wn(i));
    E.poles(i, :) = sort(real(roots([1, c3 + c4, 1 + c3*c4])))';
    fprintf('   %-16s %8.3f %9.3f %10.4f %11.4f %13.4f %15.4f %12.2f %6.3f\n', labels{i}, c3, c4, E.J(i), E.Js(i, :), E.wn(i), E.zeta(i));
end

% (b) nominal, 20 realisations -----------------------------------------------
mets = {'RMSE_ss', 'ess', 'urms', 'umax', 'du_rms', 'hip_pp', 'IAE', 'RMSE_meas'};
Cs = {};  op = {};  sd = {};
for i = 1:nd
    for s = S.eval_seeds, Cs{end+1} = CC{i}; op{end+1} = {}; sd{end+1} = s; end %#ok<AGROW>
end
Rn = n2_run_jobs(Cs, op, sd);
ns = numel(S.eval_seeds);
V = NaN(nd, ns, numel(mets));
for i = 1:nd
    for s = 1:ns
        r = Rn{(i - 1)*ns + s};
        if r.ok, for m = 1:numel(mets), V(i, s, m) = r.(mets{m}); end, end
    end
end
fprintf('\n(b) nominal real leg, mean +- SD over %d noise realisations\n   %-16s', ns, 'design');
fprintf('%22s', mets{:});  fprintf('\n');
E.nominal = struct();
for i = 1:nd
    fprintf('   %-16s', labels{i});
    for m = 1:numel(mets)
        v = squeeze(V(i, :, m));
        fprintf('%22s', sprintf('%.4g +- %.2g', mean(v), std(v)));
        E.nominal.(names{i}).(mets{m}) = struct('mean', mean(v), 'sd', std(v), 'values', v);
    end
    fprintf('\n');
end
fprintf('\n   best tuned (%s) vs trial and error: change of the mean and exact signed-rank p (20 pairs)\n', labels{ibest});
E.nominal_test = struct();
for m = 1:numel(mets)
    vb = squeeze(V(ibest, :, m));  vh = squeeze(V(1, :, m));
    pc = 100*(mean(vb)/mean(vh) - 1);  p = signrank_p(vb - vh);  nb = sum(vb < vh);
    fprintf('   %-10s %+9.2f %%   p = %.2e   tuned lower in %d/%d\n', mets{m}, pc, p, nb, ns);
    E.nominal_test.(mets{m}) = struct('change_pct', pc, 'p', p, 'n_lower', nb);
end

% (c) transients ---------------------------------------------------------------
tc = {
 'shank error 10 deg (nominal)', {}
 'shank error 30 deg',           {'x0', [10*deg; 0; 30*deg; 0]}
 'shank error 45 deg',           {'x0', [10*deg; 0; 45*deg; 0]}
 'hip angle 30 deg',             {'x0', [30*deg; 0; 10*deg; 0]}
 };
jobs = zeros(0, 3);
for c = 1:size(tc, 1), for i = 1:nd, for s = 1:3, jobs(end+1, :) = [c, i, s]; end, end, end %#ok<AGROW>
nj = size(jobs, 1);  TR = cell(nj, 1);  topt = tc(:, 2);
parfor j = 1:nj
    jr = jobs(j, :);
    TR{j} = p1_transient(CC{jr(2)}, jr(3), topt{jr(1)});
end
tf = {'peak05', 'settle', 'IAE5', 'tsat', 'umax'};
fprintf('\n(c) transients (mean of 3 realisations): peak error after 0.5 s (deg) / settling to 0.5 deg (s) / IAE first 5 s (rad.s) / saturated time (s) / peak torque (N.m)\n');
E.transient = struct('conditions', {tc(:, 1)'});
for c = 1:size(tc, 1)
    fprintf('   %s\n', tc{c, 1});
    for i = 1:nd
        sel = find(jobs(:, 1) == c & jobs(:, 2) == i);
        row = zeros(1, numel(tf));
        for f = 1:numel(tf), row(f) = mean(cellfun(@(r) r.(tf{f}), TR(sel))); end
        fprintf('      %-16s %8.3f %8.3f %8.4f %8.3f %8.2f\n', labels{i}, row);
        E.transient.(sprintf('c%d', c)).(names{i}) = cell2struct(num2cell(row), tf, 2);
    end
end

% (d) one-factor robustness ------------------------------------------------------
F = {'mass',     [1.0 1.2 1.4 1.6 1.8 2.0]
     'friction', [0 0.5 1 2 3 5]
     'dist_amp', [0 0.5 1 2 3]
     'dist_freq', [1 3 5 8 10]
     'noise',    [0 0.5 1 2 4 8]
     'delay',    [0 1 2 3 4 5 6 7 8 9 10]};
Cs = {};  op = {};  lab = [];
for f = 1:size(F, 1)
    for v = F{f, 2}
        switch F{f, 1}
            case 'mass',      o = {'P', M.plant('mass_scale', v)};
            case 'friction',  o = {'P', M.plant('b1', 0.5*v, 'b2', 0.3*v)};
            case 'dist_amp',  o = {'P', M.plant('d0', v)};
            case 'dist_freq', o = {'P', M.plant('dw', v)};
            case 'noise',     o = {'sig_th', 1e-3*v, 'sig_w', 1e-2*v};
            case 'delay',     o = {'delay', v};
        end
        for i = 1:nd, Cs{end+1} = CC{i}; op{end+1} = o; lab(end+1, :) = [f v i]; end %#ok<AGROW>
    end
end
Rr = n2_run_jobs(Cs, op, 1:3);
rm = cellfun(@(r) r.RMSE_ss, Rr);  ok = cellfun(@(r) r.ok, Rr);  rm(~ok) = NaN;
ur = cellfun(@(r) r.urms, Rr);  dr = cellfun(@(r) r.du_rms, Rr);  ur(~ok) = NaN;  dr(~ok) = NaN;
fid = fopen(fullfile(res, 'p1_robustness.csv'), 'w');
fprintf(fid, 'factor,level,design,rmse,urms,du_rms\n');
for j = 1:size(lab, 1), fprintf(fid, '%s,%g,%s,%.10g,%.10g,%.10g\n', F{lab(j, 1), 1}, lab(j, 2), names{lab(j, 3)}, rm(j), ur(j), dr(j)); end
fclose(fid);
fprintf('\n(d) one-factor robustness: steady RMS error (deg), 3 realisations; NaN = diverged\n');
E.robust = struct();
for f = 1:size(F, 1)
    fprintf('   %s\n   %10s', F{f, 1}, 'level');  fprintf('%17s', labels{:});  fprintf('\n');
    tab = NaN(numel(F{f, 2}), nd);
    for q = 1:numel(F{f, 2})
        v = F{f, 2}(q);
        sel = lab(:, 1) == f & lab(:, 2) == v;
        tab(q, :) = rm(sel)';
        fprintf('   %10g', v);  fprintf('%17.4f', tab(q, :));  fprintf('\n');
    end
    E.robust.(F{f, 1}) = struct('levels', F{f, 2}, 'rmse', tab);
end
dl = F{6, 2};  E.delay_margin = zeros(nd, 1);
for i = 1:nd
    mg = 0;
    for q = 1:numel(dl)
        if isfinite(E.robust.delay.rmse(q, i)) && E.robust.delay.rmse(q, i) < 0.5, mg = dl(q); else, break, end
    end
    E.delay_margin(i) = mg;
end
fprintf('   delay margin (largest delay with error < 0.5 deg, all smaller delays passing):');
for i = 1:nd, fprintf('  %s %d ms', labels{i}, E.delay_margin(i)); end
fprintf('\n');

% (e) gait-like references -------------------------------------------------------
gp = [1.0 1.2 1.6];  Cs = {};  op = {};
for g = gp
    for i = 1:nd, Cs{end+1} = CC{i}; op{end+1} = {'R', b2_reference('gait', g), 'x0', [0; 0; -5*deg; 0]}; end %#ok<AGROW>
end
Rg = n2_run_jobs(Cs, op, 1:3);
fprintf('\n(e) gait-like references, 3 realisations: steady RMS error (deg) / RMS torque (N.m) / hip swing (deg)\n');
E.gait = struct('periods', gp);
for q = 1:numel(gp)
    fprintf('   cycle %.1f s\n', gp(q));
    for i = 1:nd
        r = Rg{(q - 1)*nd + i};
        if r.ok, row = [r.RMSE_ss r.urms r.hip_pp r.du_rms]; else, row = NaN(1, 4); end
        fprintf('      %-16s %10.4f %8.3f %8.2f\n', labels{i}, row(1:3));
        E.gait.(sprintf('g%d', q)).(names{i}) = struct('RMSE_ss', row(1), 'urms', row(2), 'hip_pp', row(3), 'du_rms', row(4));
    end
end

% (f) Monte Carlo ----------------------------------------------------------------
smp = dlmread(fullfile(root, 'data', 'mc_samples.csv'), ',', 1, 0);
N = size(smp, 1);  mcd = [1 ibest];
Cs = {};  op = {};  sd = {};
for n = 1:N
    m = smp(n, :);  % mass fric d0 dw noise delay e0 h0 seed
    Pn = M.plant('mass_scale', m(1), 'b1', 0.5*m(2), 'b2', 0.3*m(2), 'd0', m(3), 'dw', m(4));
    o = {'P', Pn, 'sig_th', 1e-3*m(5), 'sig_w', 1e-2*m(5), 'delay', m(6), 'x0', [m(8)*deg; 0; (10 + m(7))*deg; 0]};
    for i = mcd, Cs{end+1} = CC{i}; op{end+1} = o; sd{end+1} = m(9); end %#ok<AGROW>
end
Rm = n2_run_jobs(Cs, op, sd);
VM = NaN(N, 2);
for n = 1:N
    for q = 1:2
        r = Rm{(n - 1)*2 + q};
        if r.ok, VM(n, q) = r.RMSE_ss; end
    end
end
fid = fopen(fullfile(res, 'p1_mc.csv'), 'w');
fprintf(fid, 'case,hand,best\n');
for n = 1:N, fprintf(fid, '%d,%.10g,%.10g\n', n, VM(n, :)); end
fclose(fid);
fprintf('\n(f) Monte Carlo, %d random legs and conditions (delays 0-3 ms)\n', N);
E.mc = struct('N', N);
for q = 1:2
    v = VM(:, q);  v = v(isfinite(v));
    ci = boot_ci_q(v, 50);  ci95 = boot_ci_q(v, 95);
    fprintf('   %-16s failed %d, median %.4f [%.4f, %.4f], P95 %.4f [%.4f, %.4f]\n', labels{mcd(q)}, sum(~isfinite(VM(:, q))), ...
            median(v), ci, prctile_(v, 95), ci95);
    E.mc.(names{mcd(q)}) = struct('failed', sum(~isfinite(VM(:, q))), 'median', median(v), 'median_ci', ci, ...
        'p95', prctile_(v, 95), 'p95_ci', ci95, 'mean', mean(v));
end
okb = all(isfinite(VM), 2);
pm = signrank_p(VM(okb, 2) - VM(okb, 1));  nw = sum(VM(okb, 2) < VM(okb, 1));
fprintf('   best tuned lower in %d of %d cases (signed-rank p = %.2e); median ratio hand/tuned = %.1f\n', nw, sum(okb), pm, ...
        median(VM(okb, 1)./VM(okb, 2)));
E.mc.wins = nw;  E.mc.n_paired = sum(okb);  E.mc.p = pm;  E.mc.median_ratio = median(VM(okb, 1)./VM(okb, 2));

% (g) time histories for the figures ---------------------------------------------
cases = {'nominal', {}; 'e30', {'x0', [10*deg; 0; 30*deg; 0]}};
R0 = b2_reference('sine');
for c = 1:size(cases, 1)
    outs = cell(1, 2);
    for q = 1:2, outs{q} = n2_simulate(CC{mcd(q)}, cases{c, 2}{:}, 'seed', S.eval_seeds(1), 'out_dt', 5e-3); end
    t = outs{1}.t;  xd = arrayfun(@(tt) b2_ref_eval(tt, R0), t);
    fid = fopen(fullfile(res, sprintf('p1_timehist_%s.csv', cases{c, 1})), 'w');
    fprintf(fid, 't,xd_deg,e_hand_deg,u_hand,uc_hand,hip_hand_deg,flex_hand_deg,e_best_deg,u_best,uc_best,hip_best_deg,flex_best_deg\n');
    for k = 1:numel(t)
        row = xd(k)*180/pi;
        for q = 1:2
            o = outs{q};
            row = [row, (o.x(k, 3) - xd(k))*180/pi, o.u(k), o.uc(k), o.x(k, 1)*180/pi, (o.x(k, 1) - o.x(k, 3))*180/pi]; %#ok<AGROW>
        end
        fprintf(fid, '%.4f', t(k));  fprintf(fid, ',%.8g', row);  fprintf(fid, '\n');
    end
    fclose(fid);
end
fprintf('\ntime histories written: p1_timehist_nominal.csv, p1_timehist_e30.csv\n');
p1_save_json('p1_evaluate', E);
end

% ======================================================================
function t = tag(name)
t = strrep(strrep(name, '+', '_'), '-', '');
end

function p = prctile_(v, q)
v = sort(v(:));  n = numel(v);
if n == 0, p = NaN; return, end
pos = 1 + (n - 1)*q/100;  lo = floor(pos);  hi = min(lo + 1, n);
p = v(lo) + (pos - lo)*(v(hi) - v(lo));
end

function ci = boot_ci_q(v, q)
rs = RandStream('mt19937ar', 'Seed', 0);  nb = 2000;  m = zeros(nb, 1);  n = numel(v);
for b = 1:nb, m(b) = prctile_(v(ceil(n*rand(rs, n, 1))), q); end
m = sort(m);  ci = [m(round(0.025*nb)) m(round(0.975*nb))];
end

function p = signrank_p(d)
% two-sided Wilcoxon signed-rank test; exact for n <= 25 without ties, else
% normal approximation with tie correction
d = d(isfinite(d) & d ~= 0);  n = numel(d);
if n == 0, p = 1; return, end
[~, ord] = sort(abs(d));  r = zeros(n, 1);  r(ord) = 1:n;
a = abs(d);  ua = unique(a);  tie = 0;
for k = 1:numel(ua)
    idx = a == ua(k);  t = sum(idx);
    if t > 1, r(idx) = mean(r(idx)); tie = tie + t^3 - t; end
end
W = sum(r(d > 0));
if n <= 25 && tie == 0
    smax = n*(n + 1)/2;  cnt = zeros(1, smax + 1);  cnt(1) = 1;
    for k = 1:n
        cnt(k+1:end) = cnt(k+1:end) + cnt(1:end-k);
    end
    pr = cnt/2^n;  w = round(W);
    pl = sum(pr(1:w+1));  pu = sum(pr(w+1:end));
    p = min(1, 2*min(pl, pu));
else
    mu = n*(n + 1)/4;  s2 = n*(n + 1)*(2*n + 1)/24 - tie/48;
    z = (W - mu - 0.5*sign(W - mu))/sqrt(s2);
    p = min(1, erfc(abs(z)/sqrt(2)));
end
end
