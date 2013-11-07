addpath '..';
addpath '../..';
addpath '../../test_c';

testno = 1;
opt.maxiter = 200;
%opt.checkgradient = 'no';
opt.epsrel = 1e-3;
opt.epsabs = 1e-5;
opt.epsgrad = 1e-5;
opt.gradtol = 1e-5;
opt.disp = 'off';


for testno=1:5
  eval(['info = run_test(@slra_grass, testno, opt);']);
  info1 = run_test(@slra_mex_chp, testno, opt);
  info2 = run_test(@slra, testno, opt);
  info3 = run_test(@slra_reg, testno, opt);
  info4 = run_test(@slra_fmincon, testno, opt);
  
  i_gr = info.iterinfo;
  i_perm = info1.iterinfo;
  i_perm0 = info2.iterinfo;
  i_reg = info3.iterinfo;
  i_fmin = info4.iterinfo;
  
  yLimits = [min([i_gr(2,end), i_perm(2,end), i_perm0(2,end), i_reg(2,end), i_fmin(2,end)]) max([i_gr(2,1), i_perm(2,1), i_perm0(2,1), i_reg(2,1), i_fmin(2,1)])];
  xLimits = [min([i_gr(1,1), i_perm(1,1), i_perm0(1,1), i_reg(1,1), i_fmin(1,1)]) max([i_gr(1,end), i_perm(1,end), i_perm0(1,end), i_reg(1,end), i_fmin(1,end)]) ];
  yLimits = (10.^(log10(yLimits) + [-0.125 0.125]))';
  xLimits = (10.^(log10(xLimits) + [-0.125 0.125]))';
   
 
 
  hFig = figure;
  graph1 = loglog(i_gr(1,:), i_gr(2,:), 'm.-', i_perm(1,:), i_perm(2,:), 'rs-', i_perm0(1,:), i_perm0(2,:), 'bo-', i_reg(1,:), i_reg(2,:), 'gs-',i_fmin(1,:), i_fmin(2,:), 'k*-');
  title(['Test #' num2str(testno)]);
  xlabel('time, s.');
  ylabel('fmin');
  decades_equal(gca, xLimits, yLimits);
  set(gcf, 'Position', get(gcf, 'Position') - [0 0 0 140]);
  
  i_gr = info.iterinfo';
  i_perm = info1.iterinfo';
  i_perm0 = info2.iterinfo';
  i_reg = info3.iterinfo';
  i_fmin = info4.iterinfo';
   
  save(['i_gr'  num2str(testno) '.txt'], 'i_gr', '-ascii');
  save(['i_perm'  num2str(testno) '.txt'], 'i_perm', '-ascii');
  save(['i_perm0'  num2str(testno) '.txt'], 'i_perm0', '-ascii');
  save(['i_reg'  num2str(testno) '.txt'], 'i_reg', '-ascii');
  save(['i_fmin'  num2str(testno) '.txt'], 'i_fmin', '-ascii');
%  pba = get(gca, 'PlotBoxAspectRatio');
%  pf = get(gcf, 'Position');
%  print(hFig, '-depsc', '-loose', ['test' num2str(testno) '.eps']);
end
