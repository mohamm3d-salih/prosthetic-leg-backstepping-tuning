function S = b2_settings()
%B2_SETTINGS  Shared settings of the simulator and of the tuning problems.
%
%   Digital controller at 1 kHz with zero-order hold, noisy measurements,
%   RK4 plant with 5 sub-steps, torque weight of the cost, fixed noise seed of
%   the tuning, search boxes and the scenarios of the robust cost.

S.T_tune  = 40;          % simulation length inside the tuning (s)
S.w_u     = 1e-4;        % torque weight in the cost
S.Ts      = 1e-3;        % controller sampling period (s)
S.nsub    = 5;           % plant RK4 sub-steps per sample
S.sig_th  = 1e-3;        % angle noise std (rad)
S.sig_w   = 1e-2;        % velocity noise std (rad/s)
S.noise_seed = 2026;     % fixed noise realisation during tuning
S.kchi    = 100;         % anti-windup gain
S.n_pop   = 40;          % optimizer population
S.n_iter  = 50;          % iterations  -> 40*51 = 2040 cost evaluations per run

% search ranges ------------------------------------------------------------
% backstepping + RBF, narrow range:  c3 c4 log10(gamma) log10(sigma) b kr
S.lb.proposed = [1 1 1 -6 0.5 0];     S.ub.proposed = [60 150 5 -1 6 10];
% wide ranges used for controller comparisons (optima are interior)
S.lb.proposedW = [1 1 1 -6 0.5 0];    S.ub.proposedW = [500 1500 6 -1 6 100];
S.lb.BSW  = [1 1];                    S.ub.BSW  = [500 1500];              % c3 c4
S.lb.SMCW = [1 0 0];                  S.ub.SMCW = [500 1500 200];          % lambda ks K
S.lb.PIDW = [0 0 0];                  S.ub.PIDW = [5000 50000 300];        % Kp Ki Kd
% robust (final) tuning: same wide boxes, cost averaged over scenarios
for k = {'proposed', 'BS', 'SMC', 'PID'}
    S.lb.([k{1} 'R']) = S.lb.([k{1} 'W']);  S.ub.([k{1} 'R']) = S.ub.([k{1} 'W']);
end
S.lb.proposedR6 = S.lb.proposedW;  S.ub.proposedR6 = S.ub.proposedW;

% scenarios of the robust cost (kinds ending in R): J = mean over scenarios.
% mass = shank mass/inertia scale, noise = noise multiplier, delay in samples (ms)
S.scenarios_R = struct('name', {'nominal', 'delay 3 ms', 'heavy shank + 2x noise'}, ...
                       'mass', {1.2, 1.2, 1.5}, 'noise', {1, 1, 2}, 'delay', {0, 3, 0});
% delay-tolerant variant (kind proposedR6): the three scenarios plus a 6 ms delay
S.scenarios_R6 = [S.scenarios_R, struct('name', 'delay 6 ms', 'mass', 1.2, 'noise', 1, 'delay', 6)];
end
