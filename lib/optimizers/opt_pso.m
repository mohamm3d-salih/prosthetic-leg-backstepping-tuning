function R = opt_pso(f, lb, ub, n_pop, n_iter, seed)
%OPT_PSO  Particle swarm optimisation, inertia 0.9 -> 0.4, c1 = c2 = 2.
%   R = opt_pso(f, lb, ub, n_pop, n_iter, seed)
rng(seed);  nv = numel(lb);  w_max = 0.9;  w_min = 0.4;  c1 = 2;  c2 = 2;
X = lb + rand(n_pop, nv).*(ub - lb);  V = zeros(n_pop, nv);  vmax = 0.5*(ub - lb);
J = f(X);  pbest = X;  Jp = J;  [Jg, k] = min(Jp);  gbest = pbest(k, :);  Jhist = zeros(n_iter, 1);
for it = 1:n_iter
    w = w_max - (w_max - w_min)*(it - 1)/(n_iter - 1);
    V = w*V + c1*rand(n_pop, nv).*(pbest - X) + c2*rand(n_pop, nv).*(gbest - X);
    V = min(max(V, -vmax), vmax);  X = min(max(X + V, lb), ub);
    J = f(X);  better = J < Jp;  pbest(better, :) = X(better, :);  Jp(better) = J(better);
    [Jmin, k] = min(Jp);  if Jmin < Jg, Jg = Jmin;  gbest = pbest(k, :); end
    Jhist(it) = Jg;
end
R = struct('x', gbest, 'J', Jg, 'Jhist', Jhist, 'evals', n_pop*(2:n_iter + 1)');
end
