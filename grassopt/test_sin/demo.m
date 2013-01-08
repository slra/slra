addpath '..';
addpath '../../doc';
addpath '../../test_c';
gensin

opt.maxiter = 200;
opt.epsabs = 0;
opt.gradtol = 1e-5;
opt.disp = 'off';
opt.Display = 'off';

testnos = (1:10)';
methods = {'slra_mex'  'slra_mex_chp'  'slra_grass'  'slra_fmincon' 'slra_reg'};
[fields, res] = run_tests(testnos, methods, opt);
methnames = {'slra-lm'  'slra-perm' 'slra-grass' 'slra-mat' 'slra-reg'};

fits = [testnos res{5}];
fid = fopen('sinfits.txt', 'wt');
fprintf(fid, '   {test}');
fprintf(fid, '   {%s}', methnames{:});
fprintf(fid, '\n');
fclose(fid);
save('sinfits.txt', 'fits', '-ascii', '-append');

iters = [testnos res{1}];
fid = fopen('siniters.txt', 'wt');
fprintf(fid, '   {test}');
fprintf(fid, '   {%s}', methnames{:});
fprintf(fid, '\n');
fclose(fid);
save('siniters.txt', 'iters', '-ascii', '-append');

