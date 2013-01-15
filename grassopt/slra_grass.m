% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_grass(p, s, r, opts)
  obj = slra_mex_obj('new', p,s,r);

  % Currently without psi
  if isempty(opts) | ~isfield(opts, 'Rini') | isempty(opts.Rini)
    opts.Rini = slra_mex_obj('getRini', obj); 
  end
  m = slra_mex_obj('getM', obj);

  if ~exist('opts', 'var') || isempty(opts)
    opts = struct();
  end
  if ~isfield(opts, 'maxiter'),      opts.maxiter = 100; end
  if ~isfield(opts, 'maxinner'),      opts.maxinner = 40; end
  if ~isfield(opts, 'gradtol'),       opts.gradtol = 1e-5; end
  if ~isfield(opts, 'verbosity'),     opts.verbosity = 0; end
  if ~isfield(opts, 'order'),         opts.order = 2; end
  if ~isfield(opts, 'computeRMSE'),   opts.computeRMSE = false; end
  if ~isfield(opts, 'kappa'),         opts.kappa = 0.1; end
  if ~isfield(opts, 'theta'),         opts.theta = 1.0; end

  if isfield(opts, 'disp') & strcmp(opts.disp, 'iter'),     opts.verbosity = 1; end

  function y = grRetr(x, v, t)
    if nargin < 3 || isempty(t)  t = 1;  end
    [y R] = qr(x'+t*(v'), 0);         
    y = y';
  end

  fns.R = @grRetr;
  fns.g = @(U, H1, H2) sum(H1(:) .* H2(:));
  fns.f = @(U) slra_mex_obj('func', obj, U);
  fns.proj =  @(U, H) H - H *(U.')*(U);
  fns.fgrad = @(U) fns.proj(U, slra_mex_obj('grad', obj, U));
  
  
  function Hh = hessM(U, H)
    fdeps = 1e-12;
    G0 = slra_mex_obj('grad', obj, U);
    G1 = slra_mex_obj('grad', obj, (U+fdeps*H));
    Hh = fns.proj(U, (G1 - G0)/fdeps);
  end
  
  fns.fhess = @hessM; %@(U, H) fns.proj(U, slra_mex_obj('hessMult', obj, U', H')');

  params.x0 = opts.Rini;
  params.Delta_bar = sqrt((m-r))*pi/2;
  params.Delta0   = params.Delta_bar/8;
  params.max_inner = opts.maxinner;
  params.max_outer = opts.maxiter;
  params.useRand = 0;
  params.testgh = 1;
  params.debug = 0;
  params.verbosity = opts.verbosity;
  params.epsilon = opts.gradtol;
  params.kappa = opts.kappa;
  params.theta = opts.theta;

  addpath '~/software/GenRTR/solvers';
  [U stats] = rtr(fns, params);
  
  
  info.Rh = U;
  ph = slra_mex_obj('getPh', obj, info.Rh);
  info.iter = length(stats)-1;
  info.time = sum([stats.time]);
  stats = stats(end);
  info.fmin = stats.fx;
%  info = struct('Rh',Rh, 'Vh', [], 'fmin', stats.fx, 'iter', stats.k, 'time', stats.time);
  [varargout{1:nargout}] = deal(ph, info);
end
