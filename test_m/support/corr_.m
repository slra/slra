% CORR_ - finds the STLS correction without calling 
% DECODE_STRUCT and 

function [ dp, ch, ph ] = corr( x, c, struct )

[S,V] = decode_struct( struct );

[m,n_d] = size(c);
[n,d]   = size(x);
if isstruct(struct)
  K = struct.k;
  struct = struct.a;
else
  K = 1;
end
t = K;
s = size(V,3);
q = size(struct,1);

% Find yr
if t*d == 1
  f = zeros(1,s);
  for k = 1:s
    f(k) = [x;-1]' * V(:,:,k)' * [x;-1];
  end
  r  = c(:,1:n)*x-c(:,n+1);
  yr = mb02gd(f,r);
else
  x_ext = [x;-eye(d)];
  if t > 1
    x_ext = kron(eye(t),x_ext);
  end
  f = zeros(t*d,t*d,s);
  for k = 1:s
    f(:,:,k)  = x_ext' * V(:,:,k)' * x_ext;
  end
  F  = reshape(f,size(f,1),size(f,2)*size(f,3));
  r  = c(:,1:n)*x-c(:,n+1:end);
  rt = r';
  yr = mb02gd(F,rt(:));
end

% Compute Vp G^T yr
dp = [];
Yr = reshape(yr,d,m)';
col_k_b = 0; % first column of X_ext,k in X_ext 
col_k_e = 0; %  last column of X_ext,k in X_ext  
x_ext   = [x;-eye(d)];
% first compute G^T yr
for k = 1:q
  col_k_b = col_k_e + 1;
  col_k_e = col_k_e + struct(k,2);
  x_ext_k = x_ext(col_k_b:col_k_e,:);
  switch struct(k,1)
   case {1,2} % Toeplitz/Hankel
    nk = struct(k,2);
    bm = K;
    bn = struct(k,3);
    mb = m/bm;
    nb = nk/bn;
    % Form the matrix gk    
    if bm == 1 & bn == 1 % scalar case
      gk = zeros((m+nk-1),d*nk);
      for l = 1:nk
        gk(l:l-1+m,(l-1)*d+1:l*d) = Yr;
      end
    else % block case
      Yr_ = zeros(bn*size(Yr));
      for l = 1:d
        Yr_(:,l:d:end) = kron(Yr(:,l),eye(bn));
      end
      gk = zeros((mb+nb-1)*bm*bn,d*nk);
      for l = 1:nb
        gk((l-1)*bn*bm+1:(l-1)*bn*bm+m*bn, (l-1)*d*bn+1:l*d*bn) = Yr_; 
      end
    end
    % find dp_k
    if struct(k,1) == 1
      x_ext_k = fliplr(reshape(x_ext_k',bn*d,nb));      
    else
      x_ext_k = x_ext_k';
    end
    dp_k = gk * x_ext_k(:);
   case 3     % Unstructured
    dp_k = x_ext_k * Yr';
   case 4
    dp_k = [];
  end
  dp = [dp; dp_k(:)];
end
% then multiply by Vp

if nargout > 1
  % Find CH
  dp_ = [0; dp];
  S   = decode_struct( struct, [], m );
  ch  = c - reshape(dp_(S(:)+1),m,n+d);
  if nargout > 2
    % Extract p from c and check consistency of C and STRUCT
    np = max(max(S));
    p  = zeros(np,1);
    for i = 1:np
      ci   = c(S==i);
      ci11 = ci(1,1);
      if ~all(all(ci == ci11(ones(size(ci)))))
        error('The structure of [A B] does not match the specification.')
      end
      p(i) = ci11;
    end
    ph  = p - dp;
  end
end