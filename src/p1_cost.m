function [J, Js] = p1_cost(p, scen, T)
%P1_COST  Robust tuning cost of the backstepping gains p = [c3 c4].
%   Mean over the scenarios of J = int_0^T t|e3| dt + 1e-4 int_0^T u^2 dt, on the
%   1 kHz digital loop with the fixed noise realisation used for tuning.
S = p1_settings();  S0 = b2_settings();  M = b2_model();  Rf = b2_reference('sine');
if nargin < 2 || isempty(scen), scen = S.scen; end
if nargin < 3, T = S.T_tune; end
C = b2_params_to_controller('BSR', p);
Js = zeros(1, numel(scen));
for k = 1:numel(scen)
    sc = scen(k);
    out = n2_simulate(C, 'T', T, 'P', M.plant('mass_scale', sc.mass), 'R', Rf, ...
        'Ts', S0.Ts, 'nsub', S0.nsub, 'sig_th', S0.sig_th*sc.noise, 'sig_w', S0.sig_w*sc.noise, ...
        'delay', sc.delay, 'k_belt', sc.k_belt, 'seed', S0.noise_seed, 'out_dt', 0.02);
    Js(k) = b2_cost(out, Rf, S0.w_u);          % 1e6 for a failed run
    if Js(k) >= 1e6, break, end
end
if any(Js >= 1e6), J = 1e6; else, J = mean(Js); end
end
