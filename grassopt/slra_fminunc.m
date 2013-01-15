% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_fminunc(p, s, r, opt)
  obj = slra_mex_obj('new', p,s,r);

  % Currently without psi


  if isempty(opt) | ~isfield(opt, 'Rini') | isempty(opt.Rini)
    opt.Rini = slra_mex_obj('getRini', obj); 
  end
  m = slra_mex_obj('getM', obj);

  if isempty(opt) | ~isfield(opt, 'psi') | isempty(opt.psi)
    opt.psi = eye(m); 
  end
  d = m-r;

  t2R = @(th) [reshape(th, d, size(opt.psi,1)-d)  -eye(d)] * opt.psi; 
  PQ2X = @(PQ) PQ(1:end, (size(PQ,2)-size(PQ,1)+1):end) \ PQ(1:end, 1:(size(PQ,2) - size(PQ,1)));
  R2t = @(R) (PQ2X(reshape(kron(opt.psi', eye(d)) \  R(:), d, size(opt.psi,1))))(:);
    
  prob.options = opt; 
  prob.x0 = R2t(opt.Rini); 
  prob.solver = 'fminunc'; 
  prob.objective = @(th) slra_mex_obj('func', obj, t2R(th));
  

  tic, [x, fval, flag, info] = fminunc(prob); t_slra = toc, 
  Rh = t2R(x);
  ph = slra_mex_obj('getPh', obj, Rh);
  [varargout{1:nargout}] = deal(ph, struct('Rh', Rh, 'Vh', [], 'fmin', fval, 'iter', info.iterations, 'time', t_slra));
end
function X = PQ2x(PQ) 
  X = 
end


