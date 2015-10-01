function [hini, sv] = h_ini_mp(p, d, kadd)
  if ~exist('kadd', 'var')
    opt = 0;
  end
  p = reshape(p,length(p), 1);
  bfp = cell2mat(p);
  bfn = cellfun(@length, p) - 1;
  maxn = max(bfn);
  Sp = [];
  
  for i=1:length(p)
    Sp = [Sp multmat(p{i}, 2 * maxn - length(p{i})+kadd)];
  end
  [U,S,V] = svd(Sp, 'econ');
  sv = diag(S);
  P = U(:, (end-d+1):end);
  hini = fliplr(poly(eig(P(1:end-1,:) \ P(2:end,:))));
  hini = -hini(:);  
end