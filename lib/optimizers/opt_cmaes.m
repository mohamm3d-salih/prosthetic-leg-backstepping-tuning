function R = opt_cmaes(f, lb, ub, n_pop, n_iter, seed)
%OPT_CMAES  (mu/mu_w, lambda) CMA-ES in the scaled box, with clipping repair.
%   R = opt_cmaes(f, lb, ub, n_pop, n_iter, seed)
rng(seed);  N = numel(lb);  lambda = n_pop;  mu = floor(lambda/2);  n_gen = n_iter + 1;  s = 0.3;
w = log(mu + 1/2) - log(1:mu)';  w = w/sum(w);  mueff = sum(w)^2/sum(w.^2);
cc = (4 + mueff/N)/(N + 4 + 2*mueff/N);  cs = (mueff + 2)/(N + mueff + 5);
c1 = 2/((N + 1.3)^2 + mueff);  cmu = min(1 - c1, 2*(mueff - 2 + 1/mueff)/((N + 2)^2 + mueff));
damps = 1 + 2*max(0, sqrt((mueff - 1)/(N + 1)) - 1) + cs;  chiN = sqrt(N)*(1 - 1/(4*N) + 1/(21*N^2));
m = rand(N, 1);  pc = zeros(N, 1);  ps = zeros(N, 1);  B = eye(N);  D = ones(N, 1);  C = eye(N);  invsqrtC = eye(N);
Jg = inf;  xg = lb;  Jhist = zeros(n_gen, 1);
for gen = 1:n_gen
    Z = randn(N, lambda);  Y = m + s*(B*(D.*Z));  Y = min(max(Y, 0), 1);
    Xc = lb + Y'.*(ub - lb);  J = f(Xc);
    [Jmin, k] = min(J);  if Jmin < Jg, Jg = Jmin;  xg = Xc(k, :); end
    Jhist(gen) = Jg;
    [~, idx] = sort(J);  m_old = m;  m = Y(:, idx(1:mu))*w;
    ps = (1 - cs)*ps + sqrt(cs*(2 - cs)*mueff)*invsqrtC*(m - m_old)/s;
    hsig = norm(ps)/sqrt(1 - (1 - cs)^(2*gen))/chiN < 1.4 + 2/(N + 1);
    pc = (1 - cc)*pc + hsig*sqrt(cc*(2 - cc)*mueff)*(m - m_old)/s;
    Yd = (Y(:, idx(1:mu)) - m_old)/s;
    C = (1 - c1 - cmu)*C + c1*(pc*pc' + (1 - hsig)*cc*(2 - cc)*C) + cmu*Yd*diag(w)*Yd';
    s = min(s*exp((cs/damps)*(norm(ps)/chiN - 1)), 1);
    C = triu(C) + triu(C, 1)';  [B, D2] = eig(C);  D = sqrt(max(diag(D2), 1e-20));  invsqrtC = B*diag(1./D)*B';
end
R = struct('x', xg, 'J', Jg, 'Jhist', Jhist, 'evals', lambda*(1:n_gen)');
end
