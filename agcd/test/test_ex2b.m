%% Clean up and set paths
addpath ..;
addpath fastgcd;
clear all;

%% Define methods
add_field = @(s,name,value) cell2struct([struct2cell(s); value], ...
                                        [fieldnames(s); name],1);

gcd_cofe_1 = @(p, w, d, opt) gcd_cofe(p, w, d, add_field(opt,'maxiter', 1));
gcd_cofe_5 = @(p, w, d, opt) gcd_cofe(p, w, d, add_field(opt,'maxiter', 5));
methods = {@gcd_nothing, gcd_cofe_5, @gcd_cofe, @gcd_fastgcd, @gcd_fastgcd_onestep};
names = {'LRA', 'VP$_{g}$(5)', 'VP$_{g}$', 'FASTGCD', 'FASTGCD(1)'};
methods2 = {gcd_cofe_1, @gcd_fastgcd_onestep};

normalize_vec = @(v) (v / norm(v,2));

%% Run a modified example from Section 4.2 of Bini, Boito
g10 = [1; 1; 1; 1];
g20 = [1; -1; 1; -1];
ds = [50; 100; 200; 500; 1000];
rng(54321);
testpoly = cell(1);
for k=1:length(ds)
  h0 = floor(11 * rand(ds(k)+1,1))-5;
  testpoly{k, 1} = normalize_vec(conv(g10, h0)); 
  testpoly{k, 2} = normalize_vec(conv(g20, h0));
end

nreps = 20;
sigma = 1e-4;
opt.epsgrad = 1e-14;
%opt.disp = 'iter';

[res, iters, times] = run_time_tests(testpoly, methods, nreps, ...
                         ds, opt, sigma, methods2);
idx = cellfun(@length, testpoly);
ntest = size(testpoly, 1);
idx = [(1:ntest)' idx(:,1) ds];

save_text_table('res_ex2b_matlab.txt', [{'k','n', 'd'} names], [idx res]);
save_text_table('iters_ex2b_matlab.txt', [{'k','n', 'd'} names], [idx iters]);
save_text_table('times_ex2b_matlab.txt', [], [idx times]);



