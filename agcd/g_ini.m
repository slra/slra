function [g] = g_ini(p, d)
if ~exist('opt', 'var')
  opt = struct();
end
p = reshape(p,length(p), 1);
bfp = cell2mat(p);
bfn = cellfun(@length, p) - 1;
bfell = bfn - d;
  Sd = [];
  SdBdiag = cell(1, length(p)-1);
  for i=2:length(p)
    Sd = [Sd; multmat(p{i}, bfell(1))];
    SdBdiag{i-1} = -multmat(p{1}, bfell(i));
  end
  Sd = [Sd blkdiag(SdBdiag{:})];
  [U,S,V] = svd(Sd);
  g = V(:,end);
end