function n = r2_nodes(C)
%R2_NODES  Number of RBF nodes of controller C (C.nodes if set, otherwise 8).
%   Only relevant for the adaptive controllers; the simulator calls it to size
%   the (unused) weight vector.
if isstruct(C) && isfield(C, 'nodes') && ~isempty(C.nodes)
    n = C.nodes;
else
    n = 8;
end
end
