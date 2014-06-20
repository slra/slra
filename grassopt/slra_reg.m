% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_reg(p, s, r, opt)
  obj = slra_mex_obj('new', p,s,r);

  % Currently without psi
  if isempty(opt) | ~isfield(opt, 'Rini') | isempty(opt.Rini)
    opt.Rini = slra_mex_obj('getRini', obj); 
  end
  m = slra_mex_obj('getM', obj);

  if  ~isfield(opt, 'psi') || isempty(opt.psi)
    opt.psi = eye(m); 
  end
  if  ~isfield(opt, 'method') || isempty(opt.method)
    opt.method = 'off'; 
  end
  d = m-r;

  if ~isempty(opt) & isfield(opt, 'maxiter') opt.MaxIter = opt.maxiter; end
  if ~isempty(opt) & isfield(opt, 'disp') opt.Display = opt.disp; end

  if isempty('opt') || ~isfield(opt, 'g') || isempty(opt.g), opt.g = norm(p) ^ 2; end

  function R = t2R(th)
    R = reshape(th, m - r, m); 
  end;

  opt.g = opt.g;
  R2t = @(R) (R(:));
  C = @(th) t2R(th) * t2R(th)' - eye(m - r);

  function [f, g, Hinfo] = myobjective(th)
    f = slra_mex_obj('func', obj, t2R(th)) + opt.g * norm(C(th), 'fro') ^ 2;
    if nargout > 1
      g = R2t(slra_mex_obj('grad', obj, t2R(th))+ opt.g * 4 * C(th) * t2R(th));
    end
    if nargout > 2
      Hinfo = th;
    end  
  end
  
  function res = HM(Hinfo, th)  
    res = [];
    for i=1:size(th,2)  
      R =  t2R(Hinfo);
      E = t2R(th(:,i));
      addreg = opt.g * 4 * (C(R2t(R)) * E + (E * R'  + R * E') * R);
      res = [res R2t(slra_mex_obj('mhess', obj, R, E) + addreg)];
    end  
  end;  

  % Measuring time
  t_slra_start = 0;
  iterinfo = zeros(3, opt.MaxIter+1);
  function  stop = outfun(x, optimValues, state)
    iterinfo(1, optimValues.iteration + 1) = toc(t_slra_start);
    iterinfo(2, optimValues.iteration + 1) = optimValues.fval;
    stop = false;
  end

prob.options = optimset(opt, 'GradObj', 'on',  'Hessian', 'on', 'LargeScale','on', 'HessMult', @HM, 'OutputFcn', @outfun);
  prob.x0 = R2t(opt.Rini); 
  prob.solver = 'fminunc'; 
  prob.objective = @myobjective;
  
  t_slra_start = tic; [x, fval, flag, info] = fminunc(prob); t_slra = toc(t_slra_start);
%  info, pause
  Rh = t2R(x);
  ph = slra_mex_obj('getPh', obj, Rh);
  [varargout{1:nargout}] = deal(ph, struct('Rh', Rh, 'Vh', [], 'fmin', fval, 'iter', info.iterations, 'time', t_slra, 'iterinfo', iterinfo(:,1:info.iterations+1)));
end


