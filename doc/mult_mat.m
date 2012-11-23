function S = mult_mat(R, q, n)
[g, nc] = size(R); S = zeros(n * g, nc + (n - 1) * q);
for i = 1:n
  S((1:g) + (i - 1) * g, (1:nc) + (i - 1) * q) = R;
end
