function C = b2_params_to_controller(kind, p, umax)
%B2_PARAMS_TO_CONTROLLER  Map a search vector p to a controller struct.
%   kind: 'BS', 'BSW', 'BSR'           backstepping, p = [c3 c4]  (used here)
%         'proposed', 'proposedW', ...  backstepping + RBF, p = [c3 c4 log10(gamma) log10(sigma) b kr]
%         'SMC', 'SMCW', 'SMCR'         sliding mode, p = [lambda ks K]
%         'PID', 'PIDW', 'PIDR'         PID, p = [Kp Ki Kd]
%   The suffix only names the tuning problem; the mapping is the same.
if nargin < 3, umax = 10; end
S = b2_settings();
switch kind
    case {'proposed', 'proposedW', 'proposedR', 'proposedR6'}
        C = b2_controller_params(2, 'c3', p(1), 'c4', p(2), 'gamma', 10^p(3), ...
            'sigma', 10^p(4), 'b', p(5), 'kr', p(6), 'kchi', S.kchi, 'umax', umax);
    case {'BS', 'BSW', 'BSR'}
        C = b2_controller_params(0, 'c3', p(1), 'c4', p(2), 'umax', umax);
    case {'SMC', 'SMCW', 'SMCR'}
        C = b2_controller_params(3, 'lam', p(1), 'ks', p(2), 'K', p(3), 'phi', 0.05, 'umax', umax);
    case {'PID', 'PIDW', 'PIDR'}
        C = b2_controller_params(4, 'Kp', p(1), 'Ki', p(2), 'Kd', p(3), 'umax', umax);
    otherwise
        error('unknown controller kind "%s"', kind);
end
end
