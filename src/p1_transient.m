function m = p1_transient(C, seed, opts)
%P1_TRANSIENT  Transient measures of one run (true shank angle, degrees):
%   peak05  largest |error| after 0.5 s,  settle  time after which |error| stays
%   below 0.5 deg (s),  IAE5  integral of |error| over the first 5 s (rad.s),
%   tsat  time with a saturated command |uc| > umax (s),  umax  peak torque (N.m).
out = n2_simulate(C, opts{:}, 'seed', seed, 'out_dt', 1e-3);
R = b2_reference('sine');
for k = 1:2:numel(opts)
    if strcmp(opts{k}, 'R'), R = opts{k+1}; end
end
m = struct('ok', out.ok, 'peak05', NaN, 'settle', NaN, 'IAE5', NaN, 'tsat', NaN, 'umax', NaN);
if ~out.ok, return, end
t = out.t;  xd = arrayfun(@(tt) b2_ref_eval(tt, R), t);
e = (out.x(:, 3) - xd)*180/pi;
m.peak05 = max(abs(e(t >= 0.5)));
idx = find(abs(e) > 0.5, 1, 'last');
if isempty(idx), m.settle = 0; else, m.settle = t(min(idx + 1, numel(t))); end
w = t <= 5;
m.IAE5 = trapz(t(w), abs(e(w)))*pi/180;
m.tsat = sum(abs(out.uc) > C.umax + 1e-9)*(t(2) - t(1));
m.umax = max(abs(out.u));
end
