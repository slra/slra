function [ph, info] = gcd_nls_complex(p, w, d, opt)
  if ~exist('opt', 'var')
    opt = struct();
  end
 
  p = reshape(p, 1, length(p));
  pext = [cellfun(@real, p, 'UniformOutput', false); ... 
          cellfun(@imag, p, 'UniformOutput', false)];
  p = reshape(p, length(p), 1);
  pext = reshape(pext, 2 * length(p), 1);     
      
  bfp = cell2mat(pext);
  bfn = cellfun(@length, p) - 1;
  bfell = bfn - d;
  if (~isfield(opt, 'hini') || isempty(opt.hini) )
    if (~isfield(opt, 'gini') || isempty(opt.gini) )
      opt.gini = g_ini(p, d);
    end    
    opt.hini = lsdivmult(p, d, opt.gini);
  end    
  opt.psi = [eye(2 * d + 2) kron(eye(d+1), [0 1; -1 0])];
  opt.Rini = [real(opt.hini(:)') imag(opt.hini(:)')] * opt.psi;
  opt.Rini = reshape(opt.Rini, 2, 2*(d+1));
  
  disp(pext)
  disp(bfp)
  disp(opt.Rini)
  
  
  s = struct('m', [d+1 d+1], 'n', bfn-d+1);
  s.gcd = 1;
    if isempty(w) 
      [ph, info] = slra(bfp, s, sum(s.m)-2, opt);
    else
      bfw = cell2mat(w);
      s.w = 1./bfw;
      [ph, info] = slra(bfw.*bfp, s, sum(s.m)-2, opt);
    end 
  ph = mat2cell(ph, 2*(bfn+1), [1]);
end