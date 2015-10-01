function [h, resid] = lsdivmult(p, d, g)
if ~exist('opt', 'var')
  opt = struct();
end
p = reshape(p,length(p), 1);
bfp = cell2mat(p);
bfn = cellfun(@length, p) - 1;
bfell = bfn - d;
  st2 = cumsum([0;bfell]+1);
  mm = [];
  for k=1:length(bfn)
    mm = [mm; multmat(g(st2(k):st2(k+1)-1), d)];
  end
  h = mm \ bfp;
  resid = mm * h - bfp;
end