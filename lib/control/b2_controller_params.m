function C = b2_controller_params(type, varargin)
%B2_CONTROLLER_PARAMS  Controller parameter struct with default values.
%
%   C = b2_controller_params(0, 'c3', 214, 'c4', 97)
%   Fields that a control law does not use are simply ignored.
C = struct('type', type, 'c3', 9, 'c4', 9, 'gamma', 0, 'sigma', 0, 'b', 1, ...
           'kr', 0, 'eps_r', 0.05, 'kchi', 100, 'umax', 10, ...
           'lam', 0, 'ks', 0, 'K', 0, 'phi', 0.05, 'Kp', 0, 'Ki', 0, 'Kd', 0);
for k = 1:2:numel(varargin)
    if ~isfield(C, varargin{k})
        error('unknown controller parameter "%s"', varargin{k});
    end
    C.(varargin{k}) = varargin{k+1};
end
end
