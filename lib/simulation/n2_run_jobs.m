function R = n2_run_jobs(Cs, opts, seeds)
%N2_RUN_JOBS  Evaluate many (controller, condition) pairs in parallel.
%   Cs{j}: controller struct, opts{j}: cell of n2_simulate options,
%   seeds{j} (or one vector for all jobs): noise seeds. R{j} = n2_eval result.
n = numel(Cs);
if ~iscell(seeds), seeds = repmat({seeds}, 1, n); end
R = cell(1, n);
parfor j = 1:n
    o = opts{j};
    R{j} = n2_eval(Cs{j}, seeds{j}, o{:});
end
end
