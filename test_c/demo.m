addpath '..';
addpath '../doc';
opt.maxiter = 500;
opt.epsabs = 0;
opt.gradtol = 1e-5;
opt.disp = 'off';
opt.Display = 'off';

methods = {'slra'};
[fields, res] = run_tests(1:7, methods, opt);
fields
cell2mat(res)
%{:}

