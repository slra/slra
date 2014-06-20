clear all;
addpath '..';
addpath '../..';
addpath '../../test_c';

opt.maxiter = 200;
%opt.checkgradient = 'no';
opt.epsrel = 1e-4;
opt.epsabs = 1e-5;
opt.epsgrad = 1e-7;
opt.gradtol = 1e-5;
opt.disp = 'off';
%opt1 = opt;
%opt1.method = 'ps';
%opt1.avoid_xi = 1;
%opt2 = opt;
opt.method = 'ps';

optini = opt;
optini.maxiter = 0;

%tests = { 'manopt'  'run_test(@slra_grass, testno, opt)'	'm.-';
%          'chp gsl' 'run_test(@slra_mex_chp, testno, opt)'  'rs-';
%          'fxp gsl' 'run_test(@slra, testno, opt)'          'bs-';
%          'chp svd' 'run_test(@slra_mex_chp, testno, opt2)' 'r*-'; 
%          'fxp svd' 'run_test(@slra, testno, opt2)'         'b*-';
%          'lm pinv' 'run_test(@slra, testno, opt1)'         'cv-'    
%          };
      
tests = { 'gr'      'run_test(@slra_grass, testno, opt)'	'm.-';
          'perm'    'run_test(@slra_mex_chp, testno, opt)'  'rs-';
          'perm0'   'run_test(@slra, testno, opt)'          'bo';
          'reg'     'run_test(@slra_reg, testno, opt)'      'gs-'; 
          'ddlc'    'run_test(@slra_ddlc, testno, opt)'     'cv-' };      
      
for testno=1:7
  infoini = run_test(@slra, testno, optini);
  [Q,~] = qr(infoini.Rh',0);
  opt.Rini = Q';
    
  info = cell(size(tests, 1), 1);
  xMax = 0; xMin = 1e15; yMax = 0; yMin = 1e15;
  for i=1:length(info)
    eval(['info{i} = ' tests{i,2} ';']);
    yMin = min(yMin,  min(info{i}.iterinfo(2,:)));
    yMax = max(yMax,  max(info{i}.iterinfo(2,:)));
    xMin = min(xMin,  min(info{i}.iterinfo(1,:)));
    xMax = max(xMax,  max(info{i}.iterinfo(1,:)));
  end    
  
  % equalize all starting points
  for i=1:length(info)
    info{i}.iterinfo(1,:) = info{i}.iterinfo(1,:) - min(info{i}.iterinfo(1,:))+xMin;
  end    
  
  
  yLimits = (10.^(log10([yMin yMax]) + [-0.125 0.125]))';
  xLimits = (10.^(log10([xMin xMax]) + [-0.125 0.125]))';
   
  hFig = figure;
  loglogarg = {};
  for i=1:length(info)
    loglogarg = [ loglogarg {info{i}.iterinfo(1,:), info{i}.iterinfo(2,:), tests{i,3}} ];
  end

  loglog(loglogarg{:});
  title(['Test #' num2str(testno)]);
  xlabel('time, s.');
  ylabel('fmin');
  decades_equal(gca, xLimits, yLimits);
  set(gcf, 'Position', get(gcf, 'Position') - [0 0 0 140]);
  legend(tests{:,1}, 'Location', 'EastOutside');
  
  for i=1:length(info)
    iterinfo = info{i}.iterinfo';
    save(['i_' tests{i,1}  num2str(testno) '.txt'], 'iterinfo', '-ascii');
  end
%  pba = get(gca, 'PlotBoxAspectRatio');
%  pf = get(gcf, 'Position');
%  print(hFig, '-depsc', '-loose', ['test' num2str(testno) '.eps']);
end
