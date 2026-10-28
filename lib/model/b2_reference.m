function R = b2_reference(kind, varargin)
%B2_REFERENCE  Shank-angle reference th2_d(t).
%
%   R = b2_reference('sine')               th2_d = sin(t)
%   R = b2_reference('sine', A, w)         th2_d = A*sin(w*t)
%   R = b2_reference('gait', period)       gait-like reference, th2_d = -kappa(t)
%
%   Evaluate with [xd, dxd, ddxd] = b2_ref_eval(t, R).
%
%   kappa(t) is a 6-harmonic Fourier series of the knee flexion through
%   typical landmarks of level walking: 5 deg at initial contact, 18 deg at
%   15 % of the cycle, 3 deg at 40 %, 38 deg at toe-off (60 %), 62 deg peak
%   swing flexion (73 %) and 30 deg at 87 %. The coefficients were fitted once
%   to these landmarks and are stored below. It is a synthetic gait-like
%   trajectory, not a recorded one.

switch lower(kind)
    case 'sine'
        A = 1; w = 1;
        if numel(varargin) >= 1, A = varargin{1}; end
        if numel(varargin) >= 2, w = varargin{2}; end
        R = struct('type', 0, 'A', A, 'w', w, 'off', 0);
    case 'gait'
        period = 1.2;
        if numel(varargin) >= 1, period = varargin{1}; end
        c = gait_coeffs();
        R = struct('type', 1, 'period', period, 'a0', -c.a0, 'a', -c.a, 'b', -c.b);
    otherwise
        error('unknown reference "%s"', kind);
end
end

function c = gait_coeffs()
% knee-flexion Fourier coefficients (rad); kappa = a0 + sum a_k cos(k W t) + b_k sin(k W t)
% fitted coefficients of the gait profile (do not edit by hand)
c.a0 = 0.40717290441868975;
c.a = [-0.05019766634599601 -0.24217393535364365 -0.023574300222480905 0.0030889615200351164 -0.0055122222655246906 -0.000492791832635231];
c.b = [-0.36969645524847145 0.127534177929969 0.031048842911684547 -0.006897278216224858 -0.005042917119077543 8.049055722781581e-05];
end
