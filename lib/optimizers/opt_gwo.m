function R = opt_gwo(f, lb, ub, n_pop, n_iter, seed)
%OPT_GWO  Grey wolf optimiser, a decreasing linearly from 2 to 0.
%   R = opt_gwo(f, lb, ub, n_pop, n_iter, seed)
rng(seed);  nv = numel(lb);
X = lb + rand(n_pop, nv).*(ub - lb);  J = f(X);
[~, idx] = sort(J);  L = X(idx(1:3), :);  JL = J(idx(1:3));  Jhist = zeros(n_iter, 1);
for it = 1:n_iter
    a = 2 - 2*(it - 1)/(n_iter - 1);  Xnew = zeros(n_pop, nv);
    for j = 1:3
        Aw = 2*a*rand(n_pop, nv) - a;  C = 2*rand(n_pop, nv);
        D = abs(C.*L(j, :) - X);  Xnew = Xnew + (L(j, :) - Aw.*D);
    end
    X = min(max(Xnew/3, lb), ub);  J = f(X);
    allX = [X; L];  allJ = [J; JL];  [~, idx] = sort(allJ);
    L = allX(idx(1:3), :);  JL = allJ(idx(1:3));  Jhist(it) = JL(1);
end
R = struct('x', L(1, :), 'J', JL(1), 'Jhist', Jhist, 'evals', n_pop*(2:n_iter + 1)');
end
