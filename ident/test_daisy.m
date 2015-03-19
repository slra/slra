%% Test system identification methods on examples from DAISY
clear all, close all, warning off, 
addpath ../doc/ident/data, 
% run('~/mfiles/FrequencyDomainToolbox/setup')

f  = fopen('daisy-results.txt', 'w');
ff = fopen('daisy-exec-time.txt', 'w');

%% control parameters
cmp = 1; % compare criterion (1) or normalized percentage misfit (0)
bar_plots = 1;
save_data = 0;
plots = 0;  % plots of the outputs and printed figures
stop  = 0;  % stopping after each method
eiv   = 1;  % EIV or not
order_selection = 0; ell = [];
idt   = .6;  % w(1:(idt * T)) - data for identification
val   = 1 - idt; % w(end - (val * T):end) - data for validation
first_ident = 0; % first part of the data is used for identification
detrended   = 1;
opt_l.maxiter =  500; opt_l.disp = 'off'; opt_l.method = 'll'; % for slra with LM method
opt_q.maxiter =  500; opt_q.disp = 'off'; opt_q.method = 'qb'; % for slra with QN method
opt_n.maxiter = 1000; opt_n.disp = 'off'; opt_n.method = 'nn'; % for slra with NM method
opt_p.maxiter =  500; opt_p.disp = 'off'; opt_n.method = 'ps'; % for slra with permutations method

opt_ls.maxiter =  500; opt_ls.disp = 'off'; opt_ls.method = 'ls'; % for slra with LM method
opt_q2.maxiter =  500; opt_q2.disp = 'off'; opt_q2.method = 'q2'; % for slra with QN method
opt_qp.maxiter =  500; opt_qp.disp = 'off'; opt_qp.method = 'qp'; % for slra with QN method
opt_qf.maxiter =  500; opt_qf.disp = 'off'; opt_qf.method = 'qf'; % for slra with QN method
opt_n2.maxiter = 1000; opt_n2.disp = 'off'; opt_n2.method = 'n2'; % for slra with NM method
opt_nr.maxiter = 1000; opt_nr.disp = 'off'; opt_nr.method = 'nr'; % for slra with NM method
opt_pu.maxiter =  500; opt_pu.disp = 'off'; opt_nu.method = 'pu'; % for slra with permutations method


%% data sets
test_examples = {'#' 'name'                                                               'file name'          'l'
                 1   'Data of a simulation of the western basin of Lake Erie'             'erie'                1 
                 2   'Data of a simulation of the western basin of Lake Erie'             'erie_n10'            1 
                 3   'Data of a simulation of the western basin of Lake Erie'             'erie_n20'            1 
                 4   'Data of a simulation of the western basin of Lake Erie'             'erie_n30'            1 
                 5   'Data of Ethane-ethylene destillation column'                        'destill'             1 
                 6   'Data of Ethane-ethylene destillation column'                        'destill_n10'         1 
                 7   'Data of Ethane-ethylene destillation column'                        'destill_n20'         1 
                 8   'Data of Ethane-ethylene destillation column'                        'destill_n30'         1       
                 9   'Data of a 120 MW power plant'                                       'powerplant'          2 
                10   'Heating system'                                                     'heating_system'      2 
                11   'Data from an industrial dryer (supplied by Cambridge Control Ltd)'  'dryer2'              1 
                12   'Data of a laboratory setup acting like a hair dryer'                'dryer'               5 
                13   'Data of the ball-and-beam setup in SISTA'                           'ballbeam'            2 
                14   'Wing flutter data'                                                  'flutter'             5 
                15   'Data from a flexible robot arm'                                     'robot_arm'           4 
                16   'Data of a glass furnace (Philips)'                                  'glassfurnace'        1 
                17   'Heat flow density through a two layer wall'                         'thermic_res_wall'    2 
                18   'Simulation data of a pH neutralization process in a stirring tank'  'pHdata'              2 
                19   'Data of a CD-player arm'                                            'CD_player_arm'       1 
                20   'Data from a test setup of an industrial winding process'            'winding'             2  
                21   'Liquid-saturated steam heat exchanger'                              'exchanger'           2 
                22   'Data from an industrial evaporator'                                 'evaporator'          1 
                23   'Continuous stirred tank reactor'                                    'cstr'                1 
                24   'Model of a steam generator at Abbott Power Plant in Champaign IL'   'steamgen'            1                  
};
fprintf('\nTest examples:\n')
for i = 2:size(test_examples, 1), fprintf('%3d - %s\n', test_examples{i, 1}, test_examples{i, 2}), end 
%I = input('\nWhich examples to test? (Matlab array or Enter for all.) ... '); 
%I = [1:4,10:13,15:17,19:23]; % siso = [10 12 13 14 15 21 23]; I = intersect(I, siso);
%I = [1, 5, 9:15 17:24];
I = [];
if isempty(I) 
  I = 2:size(test_examples, 1);
else
  I = I + 1;
end

%% methods
methods = {'#' 'name'      'command'                                                        'line style'
           1  'ident'    'opt_l.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_l);' '--b'
           2  'ident-qn'    'opt_q.exct = exct; sys = ident([ui, yi], nu, l, opt_q);' '--b'
           3  'ident-nm'    'opt_n.exct = exct; sys = ident([ui, yi], nu, l, opt_n);' '--b'
           4  'ident-p'    'opt_n.exct = exct; sys = ident([ui, yi], nu, l, opt_p);' '--b'
           5  'pem'          'sys = pem(iddata(yi, ui), ny * l, ''nk'', zeros(1, nu), ''dist'', ''none'', ''ss'', ''can'', ''InitialState'', ''Estimate'', ''LimitError'', 0); x0 = sys.x0; sys = ss(sys.a, sys.b, sys.c, sys.d, -1);' '--r'
           6  'moesp'        'sys = n4sid(iddata(yi, ui), l * ny, ''nk'', zeros(1, nu), ''dist'', ''none'', ''InitialState'', ''Estimate'', ''Focus'', ''Simulation'', ''N4Weight'', ''MOESP''); sys = ss(sys.a, sys.b, sys.c, sys.d, -1);' '--g'
           7  'cva'        'sys = n4sid(iddata(yi, ui), l * ny, ''nk'', zeros(1, nu), ''dist'', ''none'', ''InitialState'', ''Estimate'', ''Focus'', ''Simulation'', ''N4Weight'', ''CVA''); sys = ss(sys.a, sys.b, sys.c, sys.d, -1);' '--g'
           8 'stls' '[sys, info_stls] = stlsident([ui, yi], nu, l, opt_stls);' '--r'           
           9 'ident-m'  'opt.exct = exct; opt.solver = ''m''; [sys, info] = ident([ui, yi], nu, l, opt);'  '--r'
          10 'frq-dom'  'sys = sysid_rik([ui yi], nu, l);' '--b'
          11  'ident-ls'    'opt_ls.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_ls);' '--b'           
          12  'ident-q2'    'opt_q2.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_q2);' '--b'           
          13  'ident-qp'    'opt_qp.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_qp);' '--b'           
          14  'ident-qf'    'opt_qf.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_qf);' '--b'
          15  'ident-n2'    'opt_n2.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_n2);' '--b'
          16  'ident-nr'    'opt_nr.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_nr);' '--b'           
          17  'ident-pu'    'opt_pu.exct = exct; [sys, info_slra] = ident([ui, yi], nu, l, opt_pu);' '--b'           
};
    
fprintf('\nMethods:\n')
disp(methods(:, 1:2));
    %J = input('\nWhich methods to use? (Matlab array or Enter for all.) ... '); 
    J = [2 12:14];    %J = [1 7];
if isempty(J) 
  J = 2:size(methods, 1);
else
  J = J + 1;
end

%% evaluation criteria
eval_criteria = {'#' 'name'    'command'
                 1   'idt fit'     '[M, wh] = misfit([ui yi], sys, opt_l); c = m2c(M, ui, yi, eiv, cmp); uh = wh(:, 1:nu); yh = wh(:, nu + 1:end);'
                 2   'val fit'     '[M, wh] = misfit([uv yv], sys, opt_l); c = m2c(M, uv, yv, eiv, cmp); uh = wh(:, 1:nu); yh = wh(:, nu + 1:end);'
                 3   'compare'     '[yh, c] = compare(iddata(yv, uv), idss(sys)); c = mean(c);'
                };
fprintf('\nEvaluation criteria:\n')
disp(eval_criteria(:, 1:2))
%K = input('\nWhich evaluation criteria to show? (Matlab array or Enter for all.) ... '); 
K = [1 2]; %K = [1];
if isempty(K) 
  K = 2:size(eval_criteria, 1);
else
  K = K + 1;
end

%% run all selected methods on all selected data sets
exct = []; % if eiv == 0, exct will be 1:set to nu
res = zeros(length(I), length(J) * (1 + length(K)));
for i = 1:length(I)
  exm_name = test_examples{I(i), 2};
  dat_name = test_examples{I(i), 3};
  l        = test_examples{I(i), 4};
  col      = 1; % pointer for the current column in res
  fprintf('\n------------------------------------------------------------\n')
  fprintf('Test example #%d: %s', I(i) - 1, exm_name)
  fprintf('\n------------------------------------------------------------\n')
  fprintf(f, '%d & ', i);
  eval(dat_name)
  if detrended
    zd = detrend(iddata(y, u)); u = zd.u; y = zd.y;
  end
    
  % split the data for identification and validation
  T = size(u, 1); nu = size(u, 2); ny = size(y, 2); 
  Tidt = ceil(idt * T); Tval = floor(val * T); if idt == 1, Tval = T; end
  TT = 1:T;
  if first_ident
        TTidt = TT(1:Tidt);
        TTval = TT(end - Tval + 1:end);
  else
        TTval = TT(1:Tval);
        TTidt = TT(end - Tidt + 1:end);
  end
  ui = u(TTidt, :); yi = y(TTidt, :); 
  uv = u(TTval, :); yv = y(TTval, :); 
    
  % plots
  if plots
    for j = 1:ny
      figure(j)
      plot(TT, y(:, j), '-k'); hold on
      xlabel('x'), ylabel('y'), title('t')
      set(gca, 'fontsize', 15)
      ax = axis; axis([1 T ax(3:4)])
      plot([TTidt(end) TTidt(end)], ax(3:4), ':')
    end
  end

  if order_selection
      sys_ell = pem(iddata(yi, ui)); ell = [ell size(sys_ell, 'order') / size(y, 2)]; 
  end
  
  % run methods
  if ~eiv, exct = 1:nu; end
  for j = J
    ls     = methods{j, 4};
    m_name = methods{j, 2};
    fprintf('\nMethod %8s : ', m_name)
    try
       tic, eval(methods{j, 3}), t = toc;
    catch 
       sys = NaN; t = NaN; info_slra.time = NaN; 
    end
    fprintf('exec. time %6.2f, ', t), fprintf(f, '%5.2f & ', t); 
    if strfind(methods{j, 3}, 'slra')
            fprintf(ff, '%6.2f & ', info_slra.time);
    elseif strfind(methods{j, 3}, 'stls')
            fprintf(ff, '%6.2f & ', info_stls.time);
    end
    res(i, col) = t; col = col + 1;
    % evaluation
    for k = K 
      try 
          eval(eval_criteria{k, 3}), 
      catch
          c = NaN; yh = NaN;
      end
      fprintf('%s %9.5f, ', eval_criteria{k, 2}, c), fprintf(f, '%9.5f & ', c);
      res(i, col) = c; col = col + 1;
      % add plots
      if plots
          for jj = 1:ny
              figure(jj)
              if k == 2
                  plot(TTidt, yh(:, jj), ls);
              else
                  plot(TTval, yh(:, jj), ls);
              end
          end
      end
    end
  end
  if stop
      fprintf('\nPress any key to continue.'), pause, 
  end
  fprintf('\n'), fprintf(f, '\\\\\n'); fprintf(ff, '\\\\\n'); 
  close all
  
  % save data
  if save_data
      save(exm_name, 'info')
  end
end

idt_misfit_ind = 2:length(K) + 1:size(res, 2);
idt_misfits = res(:, idt_misfit_ind);

val_misfit_ind = 3:length(K) + 1:size(res, 2);
val_misfits = res(:, val_misfit_ind);

time_ind = 1:length(K) + 1:size(res, 2);
times = res(:, time_ind);

if bar_plots
    %% bar plots
    file_str = 'all-';
    if detrended
        file_str = [file_str 'dtr-'];
    end
    if first_ident
        file_str = [file_str 'end-'];
    end

    figure
    b = bar(idt_misfits);
    title('t'), xlabel('x'), ylabel('y')
    axis(gca, [0 size(idt_misfits, 1) + 1 0 100])
    set(gca, 'xtick', [1:size(idt_misfits, 1)])
    set(gca, 'fontsize', 15)
    %    legend(methods(J, 2)), legend boxoff
    eval(['print -depsc ' file_str 'idt-' int2str(idt*100) '.eps'])
    m_idt_misfits = mean(idt_misfits);

    figure
    b = bar(val_misfits);
    title('t'), xlabel('x'), ylabel('y')
    axis(gca, [0 size(val_misfits, 1) + 1 0 100])
    set(gca, 'xtick', [1:size(val_misfits, 1)])
    set(gca, 'fontsize', 15)
    %    legend(methods(J, 2)), legend boxoff
    eval(['print -depsc ' file_str 'val-' int2str(idt * 100) '.eps'])
    m_val_misfits = mean(val_misfits);

    figure
    b = bar(times);
    title('t'), xlabel('x'), ylabel('y')
    ax = axis;
    axis(gca, [0 size(times, 1) + 1 0 5])
    set(gca, 'xtick', [1:size(times, 1)])
    set(gca, 'fontsize', 15)
    legend(methods(J, 2), 2), legend boxoff
    eval(['print -depsc ' file_str 'time-' int2str(idt * 100) '.eps'])
        
    mean_times = mean(times)
end

fclose(f);