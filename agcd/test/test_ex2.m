%% Clean up and set paths
addpath ..;
addpath fastgcd;
clear all;

%% Define methods
add_field = @(s,name,value) cell2struct([struct2cell(s); value], ...
                                        [fieldnames(s); name],1);

gcd_nls_complex_1 = @(p, w, d, opt) gcd_nls_complex(p, w, d, add_field(opt,'maxiter', 1));
gcd_nls_complex_2 = @(p, w, d, opt) gcd_nls_complex(p, w, d, add_field(opt,'maxiter', 2));
methods = {@gcd_nothing, gcd_nls_complex_2, @gcd_nls_complex, @gcd_fastgcd, @gcd_fastgcd_onestep};
names = {'LRA', 'VP$_{h}$(2)', 'VP$_{h}$', 'FASTGCD', 'FASTGCD(1)'};
methods2 = {gcd_nls_complex_1, @gcd_fastgcd_onestep};

normalize_vec = @(v) (v / norm(v,2));

%% Run a modified example from Section 4.6 of Bini, Boito
h0 = [-1; 1; 0; 10; 1];
d = length(h0) - 1;
testpoly = cell(1);
for k=1:8
  n1 = 25*k; n2 = 15*k; n3 = 10*k;
  g1 = conv(conv([-1; zeros(n1 - 1, 1); 1], ...
      [-2; zeros(n2 - 1, 1) ; 1]), [-3; zeros(n3 - 1, 1) ; 1]);
  g2 = conv(conv([1i; zeros(n1 - 1, 1); 1], ...
      [5; zeros(n2 - 1,1) ; 1]), [2; zeros(n3 - 1, 1) ; 1]);
  testpoly{k, 1} = normalize_vec(conv(g1, h0)); 
  testpoly{k, 2} = normalize_vec(conv(g2, h0));
end

nreps = 20;
sigma = 1e-4;
opt.epsgrad = 1e-10;

[res, iters, time] = run_time_tests(testpoly, methods, nreps, ...
                         d * ones(size(testpoly,1),1), opt, sigma, methods2);
idx = cellfun(@length, testpoly);
ntest = size(testpoly, 1);
idx = [(1:ntest)' idx(:,1) ds];

save_text_table('res_ex2_matlab.txt', [{'k','n', 'd'} names], [idx res]);
save_text_table('iters_ex2_matlab.txt', [{'k','n', 'd'} names], [idx iters]);
save_text_table('times_ex2_matlab.txt', [], [idx times]);



