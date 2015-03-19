%% SLRA_EXT - MATLAB solver for structured low-rank approximation
%             (requires Optimization Toolbox)
%
%  minimize over ph |p - ph|^2_w subject to rank(S(ph)) <= r
%
%  where  S(ph) is the (m x n) structured matrix 
%         p is a given vector (of length np)
%         r is bound on the rank
%         |p|^2_w = sum(w .* (p.^2)) - the weighted semi-norm
%             defined by a vector of weights w (of length np)
%             w(i) = Inf  <=>  ph(i) = p(i) (constraint on the approximation)
%
%% Syntax
%  [ph, info] = slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C, s0)
%
%% Input 
%     tts  - structure specification S(p) = phi * (s0 + p(tts)) 
%     p    - data vector
%     r    - bound on the rank
%     w    - vector of weights 
%     Rini - initial approximation (default unstructured LRA)
%     phi  - full row rank matrix in the structure specification (default I)
%     psi  - 
%     opt  - options [optional]
%     - opt.method - optimization method:
%          'fmincon' -  uses fmincon
%          'fminunc' -  uses the exact penalty method and fminunc function
%     - other options of fmincon or fminunc can be passed in this structure    
%     th2R - structure on the kernel R = th2R(th), R S(ph) = 0 (default unstr.)
%     C    - constraint on R, C(th) = 0, (default RR' = I)
%     s0   - the constant term in the structure specification (default 0)
%
%% Reference
% I. Markovsky and K. Usevich. Structured low-rank approximation 
% with missing data. SIAM J. Matrix Anal. Appl., pages 814-830, 2013.
%
%% Output
%  ph        - approximation of p, corresponding to low-rank matrix
%  info      - information on optimization returned by MATLAB
%
%% See also
%   slra
function [ph, info] = slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C, s0)
p = p(:); if ~exist('opt'), opt = []; end
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end 
vec_tts = tts(:); NP = 1:np;
bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('s0') || isempty(s0), s0 = zeros(mp, n); end
if ~exist('psi', 'var') | isempty(psi), psi = eye(m * (m - r)); end
if ~isfield(opt, 'R0'), opt.R0 = 0; end
if ~exist('th2R') || isempty(th2R) 
  th2R = @(th) reshape(th * psi + opt.R0(:)', m - r, m); 
end 
if ~exist('C') || isempty(C), C = @(th) th2R(th) * th2R(th)' - eye(m - r); end  
if ~exist('Rini') | isempty(Rini)
  pext = [0; p];
  Rini = lra(phi * (s0 + pext(tts + 1)), r);
end
prob.options = optimset(opt); 
if isempty(prob.options.Display)
   prob.options = optimset(prob.options, 'disp', 'off'); 
end
global Th F; Th = []; F = []; 
prob.options = optimset(prob.options, 'OutputFcn', @OutputFcn);
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
if reg
  prob.solver = 'fminunc';
else
  prob.solver = 'fmincon'; 
end
pext = [0; p]; 
%Rini = - Rini(:, end - size(Rini, 1) + 1:end) \ Rini; % quick and dirty fix needed for slra_armax
prob.x0 = R2th(Rini, phi * (s0 + pext(tts + 1)), psi, opt.R0); 
Im = find(isnan(p)); 
if exist('w') && length(w(:)) == length(p), 
  Im = unique([Im(:); find(w(:) == 0)]); 
  p(Im) = NaN;
end
Ig = setdiff(1:np, Im); 
if exist('w') & ~isempty(w)
  if any(size(w) == 1), w = diag(w); end
  if size(w, 1) == np, w = w(Ig, Ig); end  
  If = find(isinf(diag(w))); 
  if ~isempty(If)
    pf = p(Ig(If));
    s0 = s0 + reshape(bfs(:, Ig(If)) * pf, mp, n);
    w(If, :) = []; w(:, If) = []; p(Ig(If)) = []; 
    bfs(:, Ig(If)) = []; 
    Ig_ = Ig; np_ = np; np = length(p); 
    tts = reshape(bfs * vec(1:np), mp, n);
    Im = find(isnan(p)); 
    if exist('w') && length(w(:)) == length(p), 
      Im = unique([Im(:); find(w(:) == 0)]); 
      p(Im) = NaN;
    end
    Ig = setdiff(1:np, Im); 
  end
  sqrt_w = sqrtm(w); inv_sqrt_w = pinv(sqrt_w);
  bfs = double(bfs); p(Ig) = sqrt_w * p(Ig); 
  bfs(:, Ig) = bfs(:, Ig) * inv_sqrt_w;
end
if reg
  if ~exist('opt') || ~isfield(opt, 'g') || isempty(opt.g) 
    opt.g = norm(p(Ig)) ^ 2; 
  end
  prob.objective = @(th) Mslra_ext(th2R(th), tts, p, [], bfs, phi, s0) ...
                         + opt.g * norm(C(th), 'fro') ^ 2;
else
  prob.objective = @(th) Mslra_ext(th2R(th), tts, p, [], bfs, phi, s0);
  prob.nonlcon = @(th) deal([], C(th));
end  
if reg
  [x, fval, flag, info] = fminunc(prob);
else
  prob.options = optimset(prob.options, 'alg', 'sqp');
  [x, fval, flag, info] = fmincon(prob); 
end
info.fmin = fval;, info.Rh = th2R(x); 
info.Th = Th; info.F = F;
[M, ph] = Mslra_ext(info.Rh, tts, p, [], bfs, phi, s0);
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

function th = R2th(R, d, psi, R0)
if (size(psi, 1) ~= size(psi, 2)) && (all(all(R0 == 0)))
  P = null(R); dh = P * (P \ d); size(dh);
  th = lra(psi * kron(dh, eye(size(R, 1))), size(psi, 1) - 1);
else
  th = vec(R - R0)' / psi;
end
function stop = OutputFcn(x, optimValues, state) 
global Th F; Th = [Th x']; F = [F optimValues.fval']; stop = 0;
