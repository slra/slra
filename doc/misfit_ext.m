function [M, ph] = misfit_ext(R, tts, p, w, bfs, phi, s0)
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('s0') || isempty(s0), s0 = zeros(m, n); end
if ~exist('bfs') | isempty(bfs), vec_tts = tts(:); NP = 1:np;
                                 bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);, end
Im = find(isnan(p)); Ig = setdiff(1:np, Im); 
if exist('w') & ~isempty(w)
  if any(size(w) == 1), w = diag(w); end
  if size(w, 1) == np, w = w(Ig, Ig); end  
  If = find(isinf(diag(w)));
  if ~isempty(If)
    pf = p(Ig(If));
    s0 = s0 + reshape(bfs(:, Ig(If)) * pf, m, n);
    w(If, :) = []; w(:, If) = []; p(Ig(If)) = []; 
    bfs(:, Ig(If)) = []; 
    Ig_ = Ig; np_ = np; np = length(p); 
    tts = reshape(bfs * vec(1:np), m, n);
    Im = find(isnan(p)); Ig = setdiff(1:np, Im); 
  end
  sqrt_w = sqrtm(w); inv_sqrt_w = pinv(sqrt_w);
  bfs = double(bfs); p(Ig) = sqrt_w * p(Ig); 
  bfs(:, Ig) = bfs(:, Ig) * inv_sqrt_w;
end
g = reshape(R * phi * reshape(bfs, mp, n * np), size(R, 1) * n, np);  
perp = @(a) null(a')';, perp_gm = perp(g(:, Im)); bg = perp_gm * g(:, Ig); 
dpg = bg' * pinv(bg * bg') * (bg * p(Ig) + perp_gm * vec(R * s0)); 
M = dpg' * dpg; ph(Ig) = p(Ig) - dpg; ph = ph(:);
ph(Im) = - pinv(g(:, Im)) * (vec(R * s0) + g(:, Ig) * ph(Ig));
if exist('w') & ~isempty(w) 
  ph(Ig) = inv_sqrt_w * ph(Ig); 
  if ~isempty(If)
    ph_ = Inf * ones(np_, 1);
    ph_(Ig_(If)) = pf; ph_(isinf(ph_)) = ph;
    ph = ph_;
  end  
end
