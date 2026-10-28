function M = b2_model()
%B2_MODEL  Constants of the nominal model and of the simulated (test) leg.
%
%   M = b2_model() returns
%     M.NP          nominal constants [a1 a2 a3 k1 k2] used by the controller
%     M.plant(...)  function that builds the constants of the simulated leg
%
%   Geometry and masses from Salman & Kadhim (2022), Table 1. Angles are
%   absolute, measured from the downward vertical. Equal pulleys (n = 1), so
%   the motor torque acts on the shank only, B = [0; 1].
%
%   Test leg (default of M.plant):
%     P = M.plant();      shank mass and inertia +20 %, b1 = 0.5, b2 = 0.3 N.m.s/rad,
%                         external torque d(t) = 1.0*sin(3t) N.m on the shank
%   Other conditions, e.g. a 40 % heavier shank and double friction:
%     P = M.plant('mass_scale', 1.4, 'b1', 1.0, 'b2', 0.6);

M.NP = leg_consts(5.28, 2.23, 0.033, 0.033);
M.plant = @plant;
end

function c = leg_consts(m1, m2, I1, I2)
L1 = 0.302; r1 = 0.236; r2 = 0.189; g = 9.81;
c = [I1 + m1*r1^2 + m2*L1^2, ...      % a1
     m2*L1*r2, ...                    % a2
     I2 + m2*r2^2, ...                % a3
     (m1*r1 + m2*L1)*g, ...           % k1
     m2*r2*g];                        % k2
end

function P = plant(varargin)
% name-value options (defaults = test leg)
o = struct('mass_scale', 1.2, 'b1', 0.5, 'b2', 0.3, 'd0', 1.0, 'dw', 3.0, ...
           'd_step', 0.0, 't_step', 10.0, 'joint_dist', 0);
for k = 1:2:numel(varargin)
    o.(varargin{k}) = varargin{k+1};
end
c = leg_consts(5.28, 2.23*o.mass_scale, 0.033, 0.033*o.mass_scale);
P.a1 = c(1); P.a2 = c(2); P.a3 = c(3); P.k1 = c(4); P.k2 = c(5);
P.b1 = o.b1; P.b2 = o.b2;
P.d0 = o.d0; P.dw = o.dw;              % d(t) = d0*sin(dw*t) + d_step*(t >= t_step)
P.d_step = o.d_step; P.t_step = o.t_step;
P.joint_dist = o.joint_dist;           % 0: external torque on the shank (default)
                                       % 1: knee-joint torque (reaction -d on the thigh)
end
