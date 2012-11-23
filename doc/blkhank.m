function H = blkhank(w, i, j)
if length(size(w)) == 3 
  [q, N, T]  = size(w);
  if nargin < 3 | isempty(j), j = T - i + 1; end
  if j <= 0, error('Not enough data.'), end
  H = zeros(i * q, j * N);
  for ii = 1:i
    for jj = 1:j
      H(((ii - 1) * q + 1):(ii * q), ...
        ((jj - 1) * N + 1):(jj * N)) = w(: ,:, ii + jj - 1);
    end
  end
else
  [q, T] = size(w); if T < q, w = w'; [q, T] = size(w); end
  if nargin < 3 | isempty(j), j = T - i + 1; end
  if j <= 0, error('Not enough data.'), end
  H = zeros(i * q, j); 
  for ii = 1:i
    H(((ii - 1) * q + 1):(ii * q), :) = w(:, ii:(ii + j - 1));
  end
end
