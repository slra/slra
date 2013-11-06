% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_grass(p, s, r, opts)
  addpath ..;
  
  obj = slra_mex_obj('new', p,s,r);

  % Currently without psi
  if isempty(opts) | ~isfield(opts, 'Rini') | isempty(opts.Rini)
    opts.Rini = slra_mex_obj('getRini', obj); 
  end
  m = slra_mex_obj('getM', obj);

  if ~exist('opts', 'var') || isempty(opts)
    opts = struct();
  end
  if isfield(opts, 'maxiter'),      
    params.maxiter = opts.maxiter; 
  end
  if isfield(opts, 'epsgrad'),      
    params.tolgradnorm = opts.epsgrad;
  end
  if ~isfield(opts, 'epsabs'),      
    params.epsabs = 0; 
  end
  if ~isfield(opts, 'epsrel'),      
    params.epsrel = 1e-5; 
  end

  if isfield(opts, 'disp') 
    if (strcmp(opts.disp, 'iter'))
      params.verbosity = 2;
    else if (strcmp(opts.disp, 'notify'))  
      params.verbosity = 1;
    else  
      params.verbosity = 0;
    end  
  end

  manifold = grassmannfactory(sum(s.m), sum(s.m)-r,1);
  problem.M = manifold;
  problem.cost = @(x) slra_mex_obj('func', obj, x');
  problem.grad = @(x) problem.M.egrad2rgrad(x, slra_mex_obj('grad', obj, x')');
  problem.hess = @(x, eta) problem.M.ehess2rhess(x, slra_mex_obj('grad', obj, x')', slra_mex_obj('mhess', obj, x', eta')', eta);
  
  [x0, ignore] = qr(opts.Rini',0);

  problem.stopnow = @(problem,x,info,last) ((last >=2) && prod(double(abs(info(last-1).x - x) < params.epsabs + params.epsrel * abs(x))) ) ;


  if isfield(opts, 'checkgradient') 
    checkgradient(problem);
    pause
    checkhessian(problem);
  end

  [x xcost stats] = trustregions(problem, x0, params);
%  params.linesearch = @manopt.solvers.linesearch.linesearch;
%  params.beta_type = 'D';
%  params.orth_value = Inf;
%   [x xcost stats] = conjugategradient(problem, x0, params);
%   [x xcost stats] = steepestdescent(problem, x0, params);

%  figure;  
 % plotprofile(problem, x0, problem.grad(x0)/norm(problem.grad(x0),2), linspace(-2,2, 1000));
 % figure;
 % plotprofile(problem, x, problem.grad(x)/norm(problem.grad(x0),2), linspace(-2,2, 1000));
  
  info.Rh = x';
  ph = slra_mex_obj('getPh', obj, info.Rh);
  info.iter = length(stats)-1;
  info.time = sum(stats(end).time);
  info.fmin = xcost;
  info.iterinfo = zeros(3, length(stats));
  info.iterinfo(1,:) = [stats.time];
  info.iterinfo(2,:) = [stats.cost];
  info.iterinfo(3,:) = [stats.gradnorm];
  
  slra_mex_obj('delete', obj);
  [varargout{1:nargout}] = deal(ph, info);
end
