function [u, uc, s4, dhat, beta2, h] = b2_control(t, x, W, chi, I, C, NP, R)
%B2_CONTROL  Control laws for the shank angle th2 of the belt-driven leg.
%   [u, uc, s4, dhat, beta2, h] = b2_control(t, x, W, chi, I, C, NP, R)
%
%   C.type = 0  backstepping (the controller of this paper)
%            1  backstepping + online RBF compensation
%            2  backstepping + online RBF + anti-windup auxiliary state
%            3  model-based sliding-mode control with a boundary layer
%            4  PID with conditional-integration anti-windup
%   Types 1-4 are kept because the simulator is shared with related work;
%   this repository only uses type 0.
%
%   Outputs: u = sat(uc) applied torque, uc = unsaturated command,
%            s4 = signal driving the adaptation (e4, or e4 - chi for type 2),
%            dhat = RBF estimate of the lumped disturbance, beta2 = nominal
%            input gain, h = RBF activations (empty if no network is used).
%
%   Backstepping (type 0):
%     e3 = th2 - xd,  alpha = -c3*e3 + dxd,  e4 = dth2 - alpha
%     uc = (ddxd - e3 - c3*de3 - c4*e4 - f2)/beta2,   u = sat(uc, +-umax)

[xd, dxd, ddxd] = b2_ref_eval(t, R);
e3 = x(3) - xd;
de3 = x(4) - dxd;
[f2, beta2] = b2_knee_terms(x, NP);
dhat = 0;  h = [];  s4 = 0;
switch C.type
    case {0, 1, 2}
        alpha = -C.c3*e3 + dxd;
        e4 = x(4) - alpha;
        if C.type == 2, s4 = e4 - chi; else, s4 = e4; end
        if C.type >= 1
            h = b2_rbf(x, C.b, r2_nodes(C));
            dhat = W(:)'*h;
        end
        uc = (ddxd - e3 - C.c3*de3 - C.c4*s4 - f2 - dhat - C.kr*tanh(s4/C.eps_r))/beta2;
        if C.type == 2
            uc = uc - C.kchi*chi/beta2;
        end
    case 3
        s4 = de3 + C.lam*e3;
        uc = (ddxd - C.lam*de3 - f2 - C.ks*s4 - C.K*tanh(s4/C.phi))/beta2;
    case 4
        uc = -(C.Kp*e3 + C.Ki*I + C.Kd*de3);
        s4 = e3;
    otherwise
        error('unknown controller type %d', C.type);
end
u = min(max(uc, -C.umax), C.umax);
end
