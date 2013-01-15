% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_fmincon(p, s, r, opt)
  obj = slra_mex_obj('new', p,s,r);

  % Currently without psi


  if isempty(opt) | ~isfield(opt, 'Rini') | isempty(opt.Rini)
    opt.Rini = slra_mex_obj('getRini', obj); 
  end
  m = slra_mex_obj('getM', obj);

  if isempty(opt) | ~isfield(opt, 'psi') | isempty(opt.psi)
    opt.psi = eye((m-r) * m); 
  end

  if ~isempty(opt) & isfield(opt, 'maxiter') opt.MaxIter = opt.maxiter; end
  if ~isempty(opt) & isfield(opt, 'disp') opt.Display = opt.disp; end

  t2R = @(th) reshape(th * opt.psi, m - r, m); 
  R2t = @(R)  (opt.psi' \  R(:))';
    
  prob.options = opt; 
  prob.x0 = R2t(opt.Rini); 
  prob.solver = 'fmincon'; 
  prob.objective = @(th) slra_mex_obj('func', obj, t2R(th));
  prob.nonlcon = @(th) deal([], t2R(th) * t2R(th)' - eye(m - r));
  tic, [x, fval, flag, info] = fmincon(prob); t_slra = toc; 
  Rh = t2R(x);
  ph = slra_mex_obj('getPh', obj, Rh);
  [varargout{1:nargout}] = deal(ph, struct('Rh', Rh, 'Vh', [], 'fmin', fval, 'iter', info.iterations, 'time', t_slra));
end
