function J = p1_popcost(X, scen)
%P1_POPCOST  Cost of every row of X, evaluated in parallel. Used by all optimisers.
J = zeros(size(X, 1), 1);
parfor i = 1:size(X, 1)
    J(i) = p1_cost(X(i, :), scen);
end
end
