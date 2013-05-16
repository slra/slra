addpath '..';
addpath '../..';
addpath '../../test_c';

testno = 1;
opt.maxiter = 200;
opt.checkgradient = 'yes';
opt.epsrel = 1e-5;
opt.epsabs = 1e-5;
opt.epsgrad = 1e-5;
opt.gradtol = 1e-5;
opt.disp = 'off';


for testno=1:5
  figure;
  eval(['info = run_test(@slra_grass, testno, opt);']);
  info1 = run_test(@slra, testno, opt);
 
  figure;
  loglog(info.iterinfo(1,:), info.iterinfo(2,:), 'b.-', info1.iterinfo(1,:), info1.iterinfo(2,:), 'rx-');
  title(['Test #' num2str(testno)]);
  xlabel('time');
  ylabel('fmin');
  decades_equal(gca);
end
