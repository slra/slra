function s = s_old2s(s_old, bfn)
if isfield(s_old, 'K')
    N = s_old.K; s_old = s_old.a;
else
    N = 1;
end
s.m = []; s.n = ones(N, 1) * bfn; 
add_phi = any(s_old(:, 1) == 1); if add_phi, s.phi = []; end
add_w   = any(s_old(:, 1) == 4); if add_w, s.w = [];     end
while ~isempty(s_old)
    switch s_old(1)
      case 1,
        m = (s_old(1, 2) / s_old(1, 3)) * ones(s_old(1, 3), 1); 
        phi = fliplr(eye(sum(m))); w = 1;
        np = s_old(1, 3) * N * (s_old(1, 2) / s_old(1, 3) + bfn - 1);
      case 2,
        m = (s_old(1, 2) / s_old(1, 3)) * ones(s_old(1, 3), 1); 
        phi = eye(sum(m)); w = 1;
        np = s_old(1, 3) * N * (s_old(1, 2) / s_old(1, 3) + bfn - 1);
      case 3,
        m = ones(s_old(1, 2), 1); 
        phi = eye(sum(m)); np = sum(m) * sum(s.n); w = 1;
      case 4,
        m = ones(s_old(1, 2), 1);
        phi = eye(sum(m)); np = sum(m) * sum(s.n); w = inf;
      otherwise, 
        error('Incorrect structure S_OLD.')
    end
    if add_phi, s.phi = blkdiag(s.phi, phi);, end
    if add_w, s.w = [s.w; w * ones(np, 1)], end
    s.m = [s.m; m]; s_old(1, :) = [];
end
