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
    
  t_slra_start = 0;
  iterinfo = zeros(3, opt.MaxIter+1);
  function  stop = outfun(x, optimValues, state)
    iterinfo(1, optimValues.iteration + 1) = toc(t_slra_start);
    iterinfo(2, optimValues.iteration + 1) = optimValues.fval;
    stop = false;
  end

  function [f, g, Hinfo] = myobjective(th)
    f = slra_mex_obj('func', obj, t2R(th));
    if nargout > 1
      g = R2t(slra_mex_obj('grad', obj, t2R(th)));
    end
  end
    
  prob.options = optimset(opt, 'OutputFcn',@outfun); 
%  prob.options = opt; 
  prob.x0 = R2t(opt.Rini); 
  prob.solver = 'fmincon'; 
  prob.objective = @myobjective;
  prob.nonlcon = @(th) deal([], t2R(th) * t2R(th)' - eye(m - r));
  t_slra_start = tic; [x, fval, flag, info] = fmincon(prob); t_slra = toc(t_slra_start); 
  Rh = t2R(x);
  ph = slra_mex_obj('getPh', obj, Rh);
  [varargout{1:nargout}] = deal(ph, struct('Rh', Rh, 'Vh', [], 'fmin', fval, 'iter', info.iterations, 'time', t_slra,'iterinfo', iterinfo(:,1:info.iterations+1)));
end
