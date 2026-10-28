function R = opt_ba(f, lb, ub, n_pop, n_iter, seed)
%OPT_BA  Bat algorithm (bounded, minimisation).
%   R = opt_ba(f, lb, ub, n_pop, n_iter, seed); f maps the rows of X to a column of costs.
rng(seed);  nv = numel(lb);  f_min = 0;  f_max = 2;  A0 = 1;  r0 = 0.5;  alpha = 0.9;  gam = 0.9;
X = lb + rand(n_pop, nv).*(ub - lb);  V = zeros(n_pop, nv);  vmax = 0.5*(ub - lb);
A = A0*ones(n_pop, 1);  r = zeros(n_pop, 1);
J = f(X);  [Jg, k] = min(J);  xbest = X(k, :);  Jhist = zeros(n_iter, 1);
for it = 1:n_iter
    fr = f_min + (f_max - f_min)*rand(n_pop, 1);
    V = V + (X - xbest).*fr;  V = min(max(V, -vmax), vmax);  Xnew = X + V;
    walk = rand(n_pop, 1) > r;
    Xnew(walk, :) = xbest + 0.1*mean(A)*(2*rand(nnz(walk), nv) - 1).*(ub - lb);
    Xnew = min(max(Xnew, lb), ub);
    Jnew = f(Xnew);
    accept = (Jnew <= J) & (rand(n_pop, 1) < A);
    X(accept, :) = Xnew(accept, :);  J(accept) = Jnew(accept);
    A(accept) = alpha*A(accept);  r(accept) = r0*(1 - exp(-gam*it));
    [Jmin, k] = min(Jnew);  if Jmin < Jg, Jg = Jmin;  xbest = Xnew(k, :); end
    Jhist(it) = Jg;
end
R = struct('x', xbest, 'J', Jg, 'Jhist', Jhist, 'evals', n_pop*(2:n_iter + 1)');
end
