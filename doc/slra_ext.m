function [ph, info] = slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C)
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end 
vec_tts = tts(:); NP = 1:np;
bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);
if ~exist('w', 'var') | isempty(w), w = ones(np, 1); end 
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('psi', 'var') | isempty(psi), psi = eye(m * (m - r)); end
if ~exist('th2R'), th2R = @(th) reshape(th * psi, m - r, m); end 
if ~exist('C'), C = @(th) th2R(th) * th2R(th)' - eye(m - r); end  
if ~exist('Rini') | isempty(Rini), Rini = lra(phi * p(tts), r); end
prob = optimset(); 
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
if reg
  prob.solver = 'fminunc'; 
else
  prob.solver = 'fmincon'; 
end
prob.options = optimset('disp', 'off'); 
prob.x0 = R2th(Rini, phi * p(tts), psi); 
if reg
  g = norm(p(find(~isnan(p))));
  prob.objective = @(th) misfit_ext(th2R(th), tts, p, w, bfs, phi) + g * norm(C(th));
else
  prob.objective = @(th) misfit_ext(th2R(th), tts, p, w, bfs, phi);
  prob.nonlcon = @(th) deal([], C(th));
end  
if reg
  [x, fval, flag, info] = fminunc(prob); 
else
  [x, fval, flag, info] = fmincon(prob); 
end
info.fmin = fval;, info.Rh = th2R(x);
[M, ph] = misfit_ext(info.Rh, tts, p, w, bfs, phi);
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
