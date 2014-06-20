% SLRA - solves the structured low-rank approximation problem 
function varargout = slra_ddlc(p, s, r, opt)
  if isempty(opt) 
    opt = struct(); 
  end
  opt.avoid_xi = 1;
  [ph, info] = slra(p,s,r,opt);
  [varargout{1:nargout}] = deal(ph, info);
end
