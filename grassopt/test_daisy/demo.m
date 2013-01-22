addpath '..';
addpath '../..';
addpath '../../test_c';


ells = [1 1 2 5 5];
test_examples = { 'erie_n20' 'destill_n30' 'heating_system'  'dryer' 'flutter'};

gendata(ells, test_examples);

opt.maxiter = 200;
opt.epsrel = 1e-5;
opt.epsabs = 1e-5;
opt.epsgrad = 1e-5;
opt.gradtol = 1e-5;
opt.disp = 'off';
opt.Display = 'off';

testnos = (1:length(ells))';
methods = {'slra'  'slra_mex_chp'  'slra_grass'  'slra_fmincon' 'slra_reg'};
[fields, res] = run_tests(testnos, methods, opt);
methnames = {'slra-lm'  'slra-perm' 'slra-grass' 'slra-mat' 'slra-reg'};

fits = [res{5}];
fid = fopen('fits.txt', 'wt');
fprintf(fid, '   {%s}', methnames{:});
fprintf(fid, '\n');
fclose(fid);
save('fits.txt', 'fits', '-ascii', '-append');

switches = res{4};
iters = [res{1} switches(:,2)];
fid = fopen('iters.txt', 'wt');
fprintf(fid, '   {%s}', methnames{:});
fprintf(fid, '   {switches}');
fprintf(fid, '\n');
fclose(fid);
save('iters.txt', 'iters', '-ascii', '-append');

