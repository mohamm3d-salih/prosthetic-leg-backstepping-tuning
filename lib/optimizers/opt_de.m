function R = opt_de(f, lb, ub, n_pop, n_iter, seed, F, CR)
%OPT_DE  DE/rand/1/bin with greedy selection.
%   R = opt_de(f, lb, ub, n_pop, n_iter, seed, F, CR)
%   f maps the rows of X to a column of costs (it may use parfor inside).
if nargin < 7, F = 0.5; end
if nargin < 8, CR = 0.9; end
rng(seed);  nv = numel(lb);
X = lb + rand(n_pop, nv).*(ub - lb);  J = f(X);
LB = repmat(lb, n_pop, 1);  UB = repmat(ub, n_pop, 1);  Jhist = zeros(n_iter, 1);
for it = 1:n_iter
    U = X;
    for i = 1:n_pop
        r = randperm(n_pop - 1, 3);  r(r >= i) = r(r >= i) + 1;
        v = X(r(1), :) + F*(X(r(2), :) - X(r(3), :));
        cross = rand(1, nv) < CR;  cross(randi(nv)) = true;
        U(i, cross) = v(cross);
    end
    low = U < LB;  U(low) = (X(low) + LB(low))/2;
    high = U > UB;  U(high) = (X(high) + UB(high))/2;
    Ju = f(U);  better = Ju <= J;  X(better, :) = U(better, :);  J(better) = Ju(better);
    Jhist(it) = min(J);
end
[Jb, k] = min(J);
R = struct('x', X(k, :), 'J', Jb, 'Jhist', Jhist, 'evals', n_pop*(2:n_iter + 1)');
end
