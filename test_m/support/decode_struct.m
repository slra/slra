% DECODE_STRUCT - Decodes the structure specification and forms the matrices Vk

function [ S, V ] = decode_struct( struct, vp, m )

if isstruct(struct)
  K = struct.k;
  struct = struct.a;
else
  K = 1;
end

n_d = sum(struct(:,2)); % = n+d
q   = size(struct,1);   % number of blocks

% Determine the s-dependence parameter indices of the Toeplitz/Hankel structured blocks
TH  = (struct(:,1) == 1) | (struct(:,1) == 2);
if (TH == 0)
  s = 1;
else
  temp = struct(find(TH),:);
  [max_n,max_i] = max(temp(:,2));
  s    = K * max_n / temp(max_i,3); % s-dependence parameter
end

if nargin < 2
  vp = [];
end

for i = 1:q
  if rem(struct(i,2),struct(i,3))
    error('Incorrect structure specification.')
  end
end

% default number of rows is s
if nargin < 3
  m = s;
end

% Form S
S   = zeros(m,n_d);   % the explicit structure description
np  = 0;              % parameter counter
col = 1;              % column counter 
for i = 1:q
  ni = struct(i,2);   % # of column of Si
  Si = zeros(s,ni);   % i-th block of S 
  switch struct(i,1)
   case {1,2} % Toeplitz/Hankel block
    bm  = K;             % # of rows for the repeated block
    bn  = struct(i,3);   % # of columns for the repeated block
    mb  = m/bm;          % # of block rows
    nb  = ni/bn;         % # of block columns
    rb  = np + reshape(1:bm*bn,bn,bm)'; % repeated block
    if struct(i,1) == 1 % Toeplitz block
      for i = 1:mb
        for j = 1:nb
          Si((i-1)*bm+1:i*bm,(j-1)*bn+1:j*bn) = rb + (i - 1 + nb - j)*bm*bn;
        end
      end
    else % Hankel block
      for i = 1:mb
        for j = 1:nb
          Si((i-1)*bm+1:i*bm,(j-1)*bn+1:j*bn) = rb + (i + j - 2)*bm*bn;
        end
      end
    end
    np = np + (nb+mb-1)*bm*bn;
   case 3 % Unstructured block
    Si = reshape((np+1):(np+m*ni),ni,m)';
    np = np + m*ni;
   case 4 % Noise-free block
    Si = zeros(m,ni);
   otherwise
    error(sprintf('Unknown structure type ''%s''.',s(i,1)))
  end
  S(:,col:col+ni-1) = Si;
  col = col + ni; 
end

if nargout == 2
  t = K;
  % Construct {V1,...Vs}  
  if t == 1
    V  = zeros(n_d,n_d,s);
    s1 = S(1,:);
    if isempty(vp)
      s1(find(s1 == 0)) = -1;
      for k = 1:s
        sk  = S(k,:);
        V(:,:,k) = (sk(ones(n_d,1),:)' == s1(ones(n_d,1),:));
      end
    else
      vp = kron(eye(s),vp);
      nv = size(vp,1);
      vp = [zeros(1,nv+1); [zeros(nv,1) vp]];
      np = size(vp,1);
      for k = 1:s
        sk = S(k,:);
        V(:,:,k) = vp( s1(ones(n_d,1),:) + 1 + np * (sk(ones(n_d,1),:)') );
      end
    end
  else
    tn_d = t * n_d;
    st   = s/t;
    V  = zeros(tn_d,tn_d,st);
    s1 = S(1:t,:)'; 
    s1 = s1(:)';
    if isempty(vp)
      s1(find(s1 == 0)) = -1;
      for k = 1:st 
        sk = S(t*(k-1)+1:t*k,:)'; 
        sk = sk(:)';
        V(:,:,k) = (sk(ones(tn_d,1),:)' == s1(ones(tn_d,1),:));
      end    
    else
      vp = kron(eye(s),vp);
      nv = size(vp,1);
      vp = [zeros(1,nv+1); [zeros(nv,1) vp]];
      np = size(vp,1);
      for k = 1:s
        sk = S(k*(t-1)+1:k*t,:)';
        sk = sk(:)';
        V(:,:,k) = vp( s1(ones(n_d,1),:) + 1 + np * (sk(ones(n_d,1),:)') );
      end
    end
  end
end

