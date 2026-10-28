function S = p1_settings()
%P1_SETTINGS  All settings of the study in one place.
%   Controller : backstepping (type 0), parameters p = [c3 c4]
%   Cost       : mean over three scenarios of J = int t|e3| dt + 1e-4 int u^2 dt,
%                T = 20 s, 1 kHz digital controller, noisy measurements (fixed
%                noise seed 2026):
%                  1) test leg, 2) test leg + 3 ms measurement delay,
%                  3) shank 50 % heavier than nominal + doubled noise
%   Optimisers : DE, PSO, BA, GWO, CMA-ES, PSO+PS with the same population,
%                budget, box, cost and seeds
%   Budget     : population 8, 12 generations -> 8 x 13 = 104 cost evaluations
%                per run, 15 runs per optimiser
B = b2_settings();
S.lb = [1 1];   S.ub = [500 1500];       % search box of c3, c4
S.hand = [9 9];                          % trial-and-error gains of Salman & Kadhim (2022)
S.n_pop = 8;    S.n_iter = 12;           % 104 evaluations per run
S.T_tune = 20;                           % simulated time of one tuning run (s)
S.budget = S.n_pop*(S.n_iter + 1);
S.seeds = 1:15;                          % 15 runs per optimiser
S.algs  = {'DE', 'PSO', 'BA', 'GWO', 'CMA-ES', 'PSO+PS'};
S.funs  = {@opt_de, @opt_pso, @opt_ba, @opt_gwo, @opt_cmaes, @opt_psops};
base = B.scenarios_R;
S.scen = struct('mass', {base.mass}, 'noise', {base.noise}, 'delay', {base.delay}, 'k_belt', Inf);
S.scen_names = {base.name};
S.eval_seeds = 101:120;                  % 20 noise realisations, never used in the tuning
S.grid_n = 20;                           % cost map: 20 x 20 log-spaced points
S.umax = 10;
end
