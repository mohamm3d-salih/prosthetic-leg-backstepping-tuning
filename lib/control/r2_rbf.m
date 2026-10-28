function h = r2_rbf(x, b, nodes)
%R2_RBF  Gaussian RBF activations h_j = exp(-||xi - c_j||^2 / b^2).
%   NODES selects the network:
%     8, 27, 64, 125    inputs xi = [th1; th2; dth2/3], 2-5 centres per input on [-1, 1]
%     1                 constant regressor (adaptive bias), h = 1
%     16, 81, 256, 625  inputs xi = [th1; dth1/3; th2; dth2/3], 2-5 centres per input
%   The first input varies fastest in h.
if nargin < 3 || isempty(nodes), nodes = 8; end
if nodes == 1, h = 1; return, end     % adaptive bias (one constant basis function)
x = x(:);
switch nodes
    case {8, 27, 64, 125}
        z = [x(1); x(3); x(4)/3];
        m = round(nodes^(1/3));
    case {16, 81, 256, 625}
        z = x ./ [1; 3; 1; 3];
        m = round(nodes^(1/4));
    otherwise
        error('r2_rbf: unsupported number of nodes %d', nodes);
end
if m == 3, cc = [-1; 0; 1]; else, cc = linspace(-1, 1, m)'; end
h = 1;
for i = 1:numel(z)
    h = kron(exp(-(z(i) - cc).^2/b^2), h);    % earlier inputs vary fastest
end
end
