function out = n2_simulate(C, varargin)
%N2_SIMULATE  Closed-loop simulation of the belt-driven leg on a digital loop.
%   out = n2_simulate(C, 'name', value, ...)
%
%   1 kHz zero-order hold, noisy and delayed measurements, forward-Euler
%   controller states and an RK4 plant. Options (defaults in brackets):
%     'T' [30 s], 'x0', 'P' (plant, b2_model), 'R' (reference), 'Ts' [1e-3 s],
%     'nsub' [5], 'sig_th' [1e-3 rad], 'sig_w' [1e-2 rad/s], 'delay' [0 samples],
%     'quant' [0], 'seed' [1], 'out_dt' [0.02 s]
%     'tau_act'  [0 s]    first-order actuator lag: du_act/dt = (u - u_act)/tau_act
%     'k_belt'   [Inf]    belt stiffness (N.m/rad); Inf = rigid belt
%     'J_m'      [0.005]  reflected inertia of the motor and pulley A (kg.m^2)
%     'zeta_belt' [0.1]   damping ratio of the belt mode, c_belt = 2*zeta_belt*sqrt(k_belt*J_m)
%     'b_m'      [0.01]   motor viscous friction (N.m.s/rad)
%     'backlash' [0 rad]  belt backlash (dead zone), compliant belt only
%     'bias_w'   [0 rad/s] constant bias of the measured velocities
%     'pelvis'   [none]   struct('Ax', m, 'Ay', m, 'period', s): prescribed pelvis motion
%     'knee_stop' [false], 'k_stop' [200], 'c_stop' [2], 'flex_min' [0 deg]
%
%   Controller extensions used in related work (off by default): ADRC (types 5
%   and 6), initial RBF weights C.W0, a state predictor over the delay C.pred
%   and a gyroscope-bias observer C.bobs.
%
%   out: t, x (n x 4), u (applied command), uc, tau (torque delivered to the
%        shank), tstop (stop torque), Wnorm, xm3 (measured shank angle), ok

M = b2_model();
o = struct('T', 30, 'x0', [10*pi/180; 0; 10*pi/180; 0], ...
           'P', M.plant(), 'R', b2_reference('sine'), 'Ts', 1e-3, 'nsub', 5, ...
           'sig_th', 1e-3, 'sig_w', 1e-2, 'delay', 0, 'quant', 0, ...
           'seed', 1, 'out_dt', 0.02, ...
           'tau_act', 0, 'k_belt', Inf, 'J_m', 0.005, 'zeta_belt', 0.1, 'b_m', 0.01, ...
           'pelvis', [], 'knee_stop', false, 'k_stop', 200, 'c_stop', 2, 'flex_min', 0, ...
           'backlash', 0, 'bias_w', 0);
for k = 1:2:numel(varargin)
    if ~isfield(o, varargin{k}), error('unknown option "%s"', varargin{k}); end
    o.(varargin{k}) = varargin{k+1};
end
NP = M.NP;
o.rigid = isinf(o.k_belt);
if ~o.rigid, o.c_belt = 2*o.zeta_belt*sqrt(o.k_belt*o.J_m); end

nk = round(o.T/o.Ts);
oe = max(1, round(o.out_dt/o.Ts));
nout = floor(nk/oe) + 1;
t = zeros(nout, 1);  X = zeros(nout, 4);  U = t;  UC = t;  TAU = t;  TS = t;  WN = t;  XM = t;
s = [o.x0(:); o.x0(3); 0; 0];          % [th1 th1d th2 th2d psi psid u_act]
W = zeros(r2_nodes(C), 1);  chi = 0;
if isfield(C, 'W0') && ~isempty(C.W0), W = C.W0(:); end   % frozen-weight tests
  I = 0;  Z = zeros(3, 1);
buf = repmat(s(1:4)', o.delay + 1, 1);
dpred = 0;  bobs = 0;                                     % predictor horizon, bias observer
if isfield(C, 'pred') && ~isempty(C.pred)
    if C.pred < 0, dpred = o.delay; else, dpred = round(C.pred); end
end
if isfield(C, 'bobs') && ~isempty(C.bobs), bobs = C.bobs; end
UH = zeros(max(dpred, 1), 1);  BO = zeros(4, 1);
if dpred > 0
    Pn = struct('a1', NP(1), 'a2', NP(2), 'a3', NP(3), 'k1', NP(4), 'k2', NP(5), 'b1', 0, 'b2', 0, ...
                'd0', 0, 'dw', 0, 'd_step', 0, 't_step', 0, 'joint_dist', 0);
    onom = struct('P', Pn, 'knee_stop', false, 'pelvis', []);
end
dt = o.Ts/o.nsub;
rs = make_stream(o.seed);
ok = true;  io = 0;
for k = 0:nk
    tk = k*o.Ts;
    buf = [s(1:4)'; buf(1:end-1, :)];                  % delay line
    xs = buf(end, :)';
    xm = xs + [o.sig_th; o.sig_w; o.sig_th; o.sig_w].*draw(rs, 4);
    if o.quant > 0
        xm([1 3]) = o.quant*round(xm([1 3])/o.quant);
    end
    if o.bias_w ~= 0, xm([2 4]) = xm([2 4]) + o.bias_w; end   % gyroscope bias
    if bobs > 0                                        % remove the estimated gyroscope bias
        if k == 0, BO = [xm(1); xm(3); 0; 0]; end      % [th1_hat; th2_hat; b1_hat; b2_hat]
        xm([2 4]) = xm([2 4]) - BO(3:4);
    end
    xc = xm;                                           % corrected measurement (before prediction)
    if dpred > 0                                       % predict the state over the delay
        xm = predict_state(xm, UH, dpred, W, C, onom, o.Ts);
    end
    if C.type == 5
        if k == 0, Z = [xm(3); xm(4); 0]; end
        [u, uc] = adrc_law(tk, xm, Z, C, NP, o.R);
        s4 = 0;  beta2 = 0;  h = [];
    elseif C.type == 6
        if k == 0, Z = [xm(4); 0; 0]; end
        [u, uc] = adrc_vel_law(tk, xm, Z, C, NP, o.R);
        s4 = 0;  beta2 = 0;  h = [];
    else
        [u, uc, s4, ~, beta2, h] = b2_control(tk, xm, W, chi, I, C, NP, o.R);
    end
    if mod(k, oe) == 0
        io = io + 1;
        [~, tq, tst] = plant_deriv(tk, s, u, o);
        t(io) = tk;  X(io, :) = s(1:4)';  U(io) = u;  UC(io) = uc;
        TAU(io) = tq;  TS(io) = tst;  WN(io) = norm(W);  XM(io) = xc(3);
    end
    if k == nk, break; end
    if dpred > 0, UH = [UH(2:end); u]; end             % torques u(k-d+1) ... u(k)
    if bobs > 0                                        % bias observer, forward Euler
        eo = [xc(1); xc(3)] - BO(1:2);
        BO = BO + o.Ts*[xc([2 4]) + 2*bobs*eo; -bobs^2*eo];
    end
    % controller states, forward Euler (digital implementation)
    if C.type == 1 || C.type == 2
        W = W + o.Ts*C.gamma*(h*s4 - C.sigma*W);
    end
    if C.type == 2
        chi = chi + o.Ts*(-C.kchi*chi + beta2*(u - uc));
    end
    if C.type == 4 && ~((uc > C.umax && s4 < 0) || (uc < -C.umax && s4 > 0))
        I = I + o.Ts*s4;
    end
    if C.type == 5
        [f2, b2n] = b2_knee_terms(xm, NP);
        e = Z(1) - xm(3);
        Z = Z + o.Ts*[Z(2) - 3*C.wo*e; Z(3) + f2 + b2n*u - 3*C.wo^2*e; -C.wo^3*e];
    end
    if C.type == 6                                     % reduced-order ESO on the velocity
        [f2, b2n] = b2_knee_terms(xm, NP);
        e = Z(1) - xm(4);
        Z = Z + o.Ts*[Z(2) + f2 + b2n*u - 2*C.wo*e; -C.wo^2*e; 0];
    end
    % plant, RK4 with u held
    for m = 0:o.nsub-1
        tt = tk + m*dt;
        s = rk4(tt, s, u, dt, o);
    end
    if abs(s(1)) > 2*pi || abs(s(3)) > 2*pi || ~all(isfinite(s))
        ok = false;  break
    end
end
out = struct('t', t(1:io), 'x', X(1:io, :), 'u', U(1:io), 'uc', UC(1:io), ...
             'tau', TAU(1:io), 'tstop', TS(1:io), 'Wnorm', WN(1:io), 'xm3', XM(1:io), 'ok', ok);
end

% ======================================================================
function [u, uc] = adrc_law(t, xm, Z, C, NP, R)
% model-assisted linear ADRC: PD on the observer states, cancel f2 and the
% estimated total disturbance Z(3), divide by the nominal input gain
[xd, dxd, ddxd] = b2_ref_eval(t, R);
[f2, beta2] = b2_knee_terms(xm, NP);
uc = (ddxd + C.wc^2*(xd - Z(1)) + 2*C.wc*(dxd - Z(2)) - f2 - Z(3))/beta2;
u = min(max(uc, -C.umax), C.umax);
end

function x = predict_state(x, UH, d, W, C, onom, Ts)
% propagate the delayed measurement over d samples with the nominal model,
% the RBF estimate of Delta (shank channel) and the stored torques (oldest first)
useW = (C.type == 1 || C.type == 2);
for j = 1:d
    acc = n2_plant_acc(0, x, UH(j), onom);
    if useW, acc(2) = acc(2) + W(:)'*b2_rbf(x, C.b, r2_nodes(C)); end
    x(2) = x(2) + Ts*acc(1);  x(4) = x(4) + Ts*acc(2);
    x(1) = x(1) + Ts*x(2);    x(3) = x(3) + Ts*x(4);
end
end

function [u, uc] = adrc_vel_law(t, xm, Z, C, NP, R)
% ADRC with the measured knee angle, the observer velocity Z(1) and the
% estimated total disturbance Z(2) (reduced-order ESO driven by the gyroscope)
[xd, dxd, ddxd] = b2_ref_eval(t, R);
[f2, beta2] = b2_knee_terms(xm, NP);
uc = (ddxd + C.wc^2*(xd - xm(3)) + 2*C.wc*(dxd - Z(1)) - f2 - Z(2))/beta2;
u = min(max(uc, -C.umax), C.umax);
end

function s = rk4(t, s, u, dt, o)
k1 = plant_deriv(t, s, u, o);
k2 = plant_deriv(t + 0.5*dt, s + 0.5*dt*k1, u, o);
k3 = plant_deriv(t + 0.5*dt, s + 0.5*dt*k2, u, o);
k4 = plant_deriv(t + dt, s + dt*k3, u, o);
s = s + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function [ds, tq, tst] = plant_deriv(t, s, u, o)
x = s(1:4);
if o.tau_act > 0
    ua = s(7);  dua = (u - ua)/o.tau_act;
else
    ua = u;  dua = 0;
end
if o.rigid
    tq = ua;  dpsi = 0;  ddpsi = 0;
else
    if o.backlash > 0                   % belt backlash (dead zone of width o.backlash, rad)
        dl = s(5) - x(3);  gap = o.backlash/2;
        if abs(dl) <= gap
            tq = 0;
        else
            tq = o.k_belt*(dl - sign(dl)*gap) + o.c_belt*(s(6) - x(4));
        end
    else
        tq = o.k_belt*(s(5) - x(3)) + o.c_belt*(s(6) - x(4));
    end
    dpsi = s(6);
    ddpsi = (ua - tq - o.b_m*s(6))/o.J_m;
end
[acc, tst] = n2_plant_acc(t, x, tq, o);
ds = [x(2); acc(1); x(4); acc(2); dpsi; ddpsi; dua];
end

function rs = make_stream(seed)
if exist('RandStream', 'class') || exist('RandStream', 'file')
    rs = RandStream('mt19937ar', 'Seed', seed);
else
    randn('state', seed);  rs = [];
end
end

function r = draw(rs, n)
if isempty(rs), r = randn(n, 1); else, r = randn(rs, n, 1); end
end
