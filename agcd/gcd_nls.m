function [ph, info] = gcd_nls(p, w, d, opt)
  if ~exist('opt', 'var')
    opt = struct();
  end
  if ~prod(cellfun(@isreal,p))
    error('Only real polynomials are allowed. Use gcd_nls_complex instead');      
  end    
  
  p = reshape(p,length(p), 1);
  bfp = cell2mat(p);
  bfn = cellfun(@length, p) - 1;
  bfell = bfn - d;
  if (~isfield(opt, 'hini') || isempty(opt.hini) )
  if (~isfield(opt, 'gini') || isempty(opt.gini) )
    opt.gini = g_ini(p, d);
  end    
    opt.hini = lsdivmult(p, d, opt.gini);
  end    
  opt.Rini = opt.hini';
  s = struct('m', d+1, 'n', bfn-d+1);
  s.gcd = 1;
    if isempty(w) 
      [ph, info] = slra(bfp, s, sum(s.m)-1, opt);
    else
      bfw = cell2mat(w);
      s.w = 1./bfw;
      [ph, info] = slra(bfw.*bfp, s, sum(s.m)-1, opt);
    end 
  ph = mat2cell(ph, bfn+1, [1]);
end