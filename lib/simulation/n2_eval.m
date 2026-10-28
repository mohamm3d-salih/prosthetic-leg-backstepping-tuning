function r = n2_eval(C, seeds, varargin)
%N2_EVAL  Run one condition for several noise seeds and average the metrics.
%   r = n2_eval(C, 1:5, 'P', P, 'R', R, ...)   (options of n2_simulate)
%   r.ok = false if any seed diverged. Metrics over the last 10 s of a 30 s run:
%   RMSE_ss, ess (deg), urms, umax (N.m), du_rms (RMS torque rate, N.m/s),
%   hip_pp (deg), flex_min (deg, min knee flexion; < 0 = hyperextension),
%   tstop_max (N.m), Wmax (peak weight norm), IAE (rad.s) and RMSE_meas
%   (RMS of the measured shank angle minus the reference, deg).
opts = varargin;
R = b2_reference('sine');
for k = 1:2:numel(opts)
    if strcmp(opts{k}, 'R'), R = opts{k+1}; end
end
names = {'RMSE_ss', 'ess', 'urms', 'umax', 'du_rms', 'hip_pp', 'flex_min', 'tstop_max', 'Wmax', 'IAE', 'RMSE_meas'};
acc = zeros(numel(seeds), numel(names));
r.ok = true;
for i = 1:numel(seeds)
    out = n2_simulate(C, opts{:}, 'seed', seeds(i), 'out_dt', 1e-3);
    if ~out.ok, r.ok = false; break, end
    m = b2_metrics(out, R);
    last = out.t >= out.t(end) - 10;
    du = diff(out.u(last))/(out.t(2) - out.t(1));
    xd = arrayfun(@(tt) b2_ref_eval(tt, R), out.t);
    em = (out.xm3(last) - xd(last))*180/pi;
    acc(i, :) = [m.RMSE_ss, m.ess, m.urms, m.umax, sqrt(mean(du.^2)), m.hip_pp, ...
                 m.flex_min, max(out.tstop), max(out.Wnorm), m.IAE, sqrt(mean(em.^2))];
end
for j = 1:numel(names)
    if r.ok, r.(names{j}) = mean(acc(:, j)); else, r.(names{j}) = NaN; end
end
end
