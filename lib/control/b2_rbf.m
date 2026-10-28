function h = b2_rbf(x, b, nodes)
%B2_RBF  Gaussian RBF activations used by the adaptive controllers.
%   Thin wrapper around r2_rbf; the default is the 8-node network.
%   Not used by the backstepping controller of this paper.
if nargin < 3 || isempty(nodes), nodes = 8; end
h = r2_rbf(x, b, nodes);
end
