function [f2, beta2] = b2_knee_terms(x, NP)
%B2_KNEE_TERMS  Shank dynamics of the nominal model: ddth2 = f2 + beta2*u.
%
%   With equal pulleys (n = 1):
%     f2    = [c*(N1 + G1) - a1*(N2 + G2)] / D
%     beta2 = a1 / D,          D = a1*a3 - c^2,   c = a2*cos(th1 - th2)
%   beta2 stays in [1/a3, a1/(a1*a3 - a2^2)] = [8.876, 12.18] in every posture.

a1 = NP(1); a2 = NP(2); a3 = NP(3); k1 = NP(4); k2 = NP(5);
c = a2*cos(x(1) - x(3));
s = a2*sin(x(1) - x(3));
N1 = s*x(4)^2;  N2 = -s*x(2)^2;
G1 = k1*sin(x(1));  G2 = k2*sin(x(3));
D = a1*a3 - c^2;
f2 = (c*(N1 + G1) - a1*(N2 + G2))/D;
beta2 = a1/D;
end
