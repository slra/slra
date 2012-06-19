function [M, ph] = misfit_ext(R, tts, p, w, bfs, phi)
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end
if ~exist('w', 'var') | isempty(w), w = ones(np, 1); end
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('bfs') | isempty(bfs), vec_tts = tts(:); NP = 1:np;
                                 bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);, end
Im = unique([find(isnan(p)); find(w == 0)]); 
If = find(isinf(w)); 
Ibm = setdiff(1:np, Im); Imf = [Im; If(:)]; Ibmf = setdiff(1:np, Imf);
g = reshape(R * phi * reshape(bfs, mp, n * np), size(R, 1) * n, np);  
gm = g(:, Im); perp = @(a) null(a')'; gmp = perp(gm); gr = gmp * g(:, Ibmf); 
h  = gmp * (g(:, Ibmf) * p(Ibmf) + g(:, If) * p(If, :)); 
wr = w(Ibmf); inv_Wr = diag(1 ./ wr);
dpr  = inv_Wr * gr' * (pinv(gr * inv_Wr * gr') * h); 
ph(Im) = - pinv(gm) * (g(:, Ibm) * p(Ibm, :) - g(:, Ibmf) * dpr );
ph(Ibmf) = p(Ibmf) - dpr; ph(If) = p(If); ph = ph(:); 
M = dpr' * diag(wr) * dpr; % = h' * (pinv(gr * inv_Wr * gr') * h); 
