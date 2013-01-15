function S = s2s(s, np)
if nargin > 1, p = 1:np; end
if ~isfield(s, 'phi'), s.phi = eye(sum(s.m)); end
if ~isfield(s, 'n'), s.n = (length(p) - sum(s.m)) + 1; end
[m, mp] = size(s.phi); q = length(s.m); N = length(s.n); n = sum(s.n); 
s2np = @(s) sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                   - length(s.m) * length(s.n);, np = s2np(s); % = N * mp + q * n - q * N; 
p = 1:np; 
tmp = cumsum([1; s.m(:)]); Imb = tmp(1:end - 1); Ime = tmp(2:end) - 1;
tmp = cumsum([1; s.n(:)]); Inb = tmp(1:end - 1); Ine = tmp(2:end) - 1;
S = zeros(mp, n); ind = 1;
for j = 1:N
  for i = 1:q
    npij = s.m(i) + s.n(j) - 1;
    pij = p(ind:(ind + npij - 1)); ind = ind + npij;
    Hij = hankel(pij(1:s.m(i)), pij(s.m(i):end));
    S(Imb(i):Ime(i), Inb(j):Ine(j)) = Hij;
  end
end
