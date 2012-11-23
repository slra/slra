function [sys, hh] = h2ss(h, n, tol, i ,j)
if length(size(h)) == 2
  [p, T] = size(h); if p > T, h = h'; [p, T] = size(h); end
  h = reshape(h, p, 1, T);
end
[p, m, T] = size(h);
if ~exist('i', 'var') | ~isreal(i) | isempty(i) 
  i = ceil(T * m / (m + p)); 
end
if ~exist('j', 'var') | ~isreal(j) | isempty(j) 
  j = T - i; 
elseif j > T - i
  error('Not enough data.')
end  
[U, S, V] = svd(blkhank(h(:, :, 2:end), i, j), 0); s = diag(S);  
if ~exist('n', 'var') | isempty(n)
  if ~exist('tol') || isempty(tol), tol = 1e-12; end , n = sum(s > tol);
else
  if n > min(i * p - 1, j * m), error('Not enough data'), end
end
sqrt_s = sqrt(s(1:n))';
O =  sqrt_s(ones(size(U, 1), 1), :) .* U(:, 1:n);
C = (sqrt_s(ones(size(V, 1), 1), :) .* V(:, 1:n))'; 
b = C(:, 1:m); c = O(1:p, :); 
a = O(1:end - p, :) \ O((p + 1):end, :);, sys = ss(a, b, c, h(:, :, 1), -1);
if nargout > 1, hh = shiftdim(impulse(sys, T), 1); end 
