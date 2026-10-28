function [acc, tstop] = n2_plant_acc(t, x, tq, o)
%N2_PLANT_ACC  Joint accelerations of the simulated leg.
%   tq is the torque delivered to the shank (the motor torque for a rigid belt).
%
%   Optional effects (all off unless the fields/options are set):
%   - Coulomb/Stribeck friction in both joints, when P has the fields fc
%     (Coulomb torque, N.m), fs (static torque, N.m) and vs (Stribeck velocity, rad/s).
%   - Pelvis motion: the hip follows a prescribed acceleration (ax, ay), which
%     adds an inertial field -m*a to gravity:
%       G1 = (k1/g)*[(g + ay)*sin(th1) + ax*cos(th1)]
%       G2 = (k2/g)*[(g + ay)*sin(th2) + ax*cos(th2)]
%     (angles from the downward vertical, x forward, y up).
%   - Knee stop: below the knee flexion flex_min a unilateral spring-damper
%     torque tau_s = max(0, k_stop*(phi_min - phi) - c_stop*dphi) acts between
%     thigh and shank (phi = th1 - th2).

P = o.P;
c = P.a2*cos(x(1) - x(3));
s = P.a2*sin(x(1) - x(3));
d = P.d0*sin(P.dw*t) + P.d_step*(t >= P.t_step);
rel = x(4) - x(2);
Q1 = -P.b1*x(2) + P.b2*rel - d*P.joint_dist;
Q2 = tq - P.b2*rel + d;
if isfield(P, 'fc') && P.fc > 0           % Coulomb + Stribeck friction, smoothed
    fr = @(v) (P.fc + (P.fs - P.fc)*exp(-(v/P.vs)^2))*tanh(v/0.01);
    f2 = fr(rel);
    Q1 = Q1 - fr(x(2)) + f2;
    Q2 = Q2 - f2;
end
tstop = 0;
if o.knee_stop
    phi = x(1) - x(3);  dphi = x(2) - x(4);
    pmin = o.flex_min*pi/180;
    if phi < pmin
        tstop = max(0, o.k_stop*(pmin - phi) - o.c_stop*dphi);
        Q1 = Q1 + tstop;
        Q2 = Q2 - tstop;
    end
end
if isempty(o.pelvis)
    G1 = P.k1*sin(x(1));
    G2 = P.k2*sin(x(3));
else
    g = 9.81;
    w2 = 2*(2*pi/o.pelvis.period);                 % two pelvis oscillations per cycle
    ax = -o.pelvis.Ax*w2^2*sin(w2*t);
    ay = -o.pelvis.Ay*w2^2*cos(w2*t);
    G1 = P.k1/g*((g + ay)*sin(x(1)) + ax*cos(x(1)));
    G2 = P.k2/g*((g + ay)*sin(x(3)) + ax*cos(x(3)));
end
r1 = Q1 - s*x(4)^2 - G1;
r2 = Q2 + s*x(2)^2 - G2;
D = P.a1*P.a3 - c^2;
acc = [(P.a3*r1 - c*r2)/D; (P.a1*r2 - c*r1)/D];
end
