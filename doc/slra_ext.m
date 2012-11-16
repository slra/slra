function [ph, info] = slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C, s0)
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end 
vec_tts = tts(:); NP = 1:np;
bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('s0') || isempty(s0), s0 = zeros(m, n); end
if ~exist('psi', 'var') | isempty(psi), psi = eye(m * (m - r)); end
if ~exist('th2R') || isempty(th2R), th2R = @(th) reshape(th * psi, m - r, m); end 
if ~exist('C') || isempty(C), C = @(th) th2R(th) * th2R(th)' - eye(m - r); end  
if ~exist('Rini') | isempty(Rini)
  pext = [0; p];
  Rini = lra(phi * (s0 + pext(tts + 1)), r); 
end
prob = optimset(); 
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
if reg
  prob.solver = 'fminunc';
else
  prob.solver = 'fmincon'; 
end
prob.options = optimset('disp', 'iter'); 
pext = [0; p];
prob.x0 = R2th(Rini, phi * (s0 + pext(tts + 1)), psi); 
Im = find(isnan(p)); Ig = setdiff(1:np, Im); 
if exist('w') & ~isempty(w)
  if any(size(w) == 1), w = diag(w); end
  if size(w, 1) == np, w = w(Ig, Ig); end  
  If = isinf(diag(w)); 
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
if reg
  if ~exist('opt') || ~isfield(opt, 'g') || isempty(opt.g), opt.g = norm(p(Ig)) ^ 2; end
  prob.objective = @(th) misfit_ext(th2R(th), tts, p, [], bfs, phi, s0) ...
                         + opt.g * norm(C(th), 'fro') ^ 2;
else
  prob.objective = @(th) misfit_ext(th2R(th), tts, p, [], bfs, phi, s0);
  prob.nonlcon = @(th) deal([], C(th));
end  
if reg
  [x, fval, flag, info] = fminunc(prob);
else
  [x, fval, flag, info] = fmincon(prob); 
end
info.fmin = fval;, info.Rh = th2R(x);
[M, ph] = misfit_ext(info.Rh, tts, p, [], bfs, phi, s0);
if exist('w') & ~isempty(w) 
  ph(Ig) = inv_sqrt_w * ph(Ig); 
  if ~isempty(If)
    ph_ = Inf * ones(np_, 1);
    ph_(Ig_(If)) = pf; ph_(isinf(ph_)) = ph;
    ph = ph_;
  end  
end
function [R, P, dh] = lra(d, r)
d(find(isnan(d))) = 0;
[u, s, v] = svd(d); R = u(:, (r + 1):end)'; P = u(:, 1:r);
if nargout > 2, dh = u(:, 1:r) * s(1:r, 1:r) * v(:, 1:r)'; end
function th = R2th(R, d, psi)
if size(psi, 1) == size(psi, 2)
  th = R(:)' / psi;
else
  P = null(R); dh = P * (P \ d);
  th = lra(psi * kron(dh, eye(size(R, 1))), size(psi, 1) - 1);
end
