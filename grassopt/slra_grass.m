% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_grass(p, s, r, opts)
  addpath ..;
  addpath '~/src/manopt/';
  import manopt.solvers.trustregions.*;
  import manopt.manifolds.grassmann.*;
  import manopt.tools.*;
  
  importmanopt;

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
  problem.grad = @(x) manifold.proj(x, slra_mex_obj('grad', obj, x')');
  x0 = opts.Rini';

  problem.stopnow = @(problem,x,info,last) ((last >=2) && prod(double(abs(info(last-1).x - x) < params.epsabs + params.epsrel * abs(x))) ) ;


  if isfield(opts, 'checkgradient') 
    checkgradient(problem);
  end

  [x xcost stats] = trustregions(problem, x0, params);
  
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
