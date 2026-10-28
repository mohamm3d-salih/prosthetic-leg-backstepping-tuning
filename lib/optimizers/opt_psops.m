function R = opt_psops(f, lb, ub, n_pop, n_iter, seed)
%OPT_PSOPS  PSO for 80 % of the budget, then a compass (pattern) search.
%   R = opt_psops(f, lb, ub, n_pop, n_iter, seed)
rng(seed);  nv = numel(lb);  w_max = 0.9;  w_min = 0.4;  c1 = 2;  c2 = 2;
n_iter_pso = round(0.8*n_iter);  n_poll = floor(n_pop*(n_iter - n_iter_pso)/(2*nv));  step0 = 0.05;  step_min = 1e-4;
X = lb + rand(n_pop, nv).*(ub - lb);  V = zeros(n_pop, nv);  vmax = 0.5*(ub - lb);
J = f(X);  pbest = X;  Jp = J;  [Jg, k] = min(Jp);  gbest = pbest(k, :);  Jhist = zeros(n_iter_pso + n_poll, 1);
for it = 1:n_iter_pso
    w = w_max - (w_max - w_min)*(it - 1)/(n_iter_pso - 1);
    V = w*V + c1*rand(n_pop, nv).*(pbest - X) + c2*rand(n_pop, nv).*(gbest - X);
    V = min(max(V, -vmax), vmax);  X = min(max(X + V, lb), ub);
    J = f(X);  better = J < Jp;  pbest(better, :) = X(better, :);  Jp(better) = J(better);
    [Jmin, k] = min(Jp);  if Jmin < Jg, Jg = Jmin;  gbest = pbest(k, :); end
    Jhist(it) = Jg;
end
y = (gbest - lb)./(ub - lb);  step = step0;
for p = 1:n_poll
    if step >= step_min
        Yp = [repmat(y, nv, 1) + step*eye(nv); repmat(y, nv, 1) - step*eye(nv)];
        Yp = min(max(Yp, 0), 1);  kk = find(any(Yp ~= y, 2));
        Jt = inf(2*nv, 1);
        if ~isempty(kk), Jt(kk) = f(lb + Yp(kk, :).*(ub - lb)); end
        [Jmin, k] = min(Jt);
        if Jmin < Jg, Jg = Jmin;  y = Yp(k, :);  step = min(2*step, step0); else, step = step/2; end
    end
    Jhist(n_iter_pso + p) = Jg;
end
R = struct('x', lb + y.*(ub - lb), 'J', Jg, 'Jhist', Jhist, ...
           'evals', [n_pop*(2:n_iter_pso + 1), n_pop*(n_iter_pso + 1) + 2*nv*(1:n_poll)]');
end
