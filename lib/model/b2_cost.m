function J = b2_cost(out, R, w_u)
%B2_COST  Tuning cost of one run: J = int t|e3| dt + w_u int u^2 dt (trapezoid rule).
%   A run that failed (leg fell over or the simulation diverged) returns J = 1e6.
if nargin < 3, w_u = 1e-4; end
if ~out.ok, J = 1e6; return, end
xd = arrayfun(@(tt) b2_ref_eval(tt, R), out.t);
e = out.x(:, 3) - xd;
J = trapz(out.t, out.t.*abs(e)) + w_u*trapz(out.t, out.u.^2);
end
