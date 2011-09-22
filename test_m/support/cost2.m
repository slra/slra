% ---------------------------------------------------------------------------
% COST2 - STLS cost function
function cost = cost2( x, a, b, struct  )

[S,V] = decode_struct(struct);

[m,n] = size(a);
d     = size(b,2);
s     = size(V,3);
t     = size(V,1) / (n+d);
x_ext = [x;-eye(d)];
if t > 1
  x_ext = kron(eye(t),x_ext);
end

% Construct the first row of Vr
f = zeros(t*d,t*d,s);
for k = 1:s
  f(:,:,k)  = x_ext' * V(:,:,k)' * x_ext;
end
F  = reshape(f,size(f,1),size(f,2)*size(f,3));
r  = a*x-b;
rt = r';

% Solve the system Vr*yr = R(:)
yr = mb02gd(F,rt(:));
cost = rt(:)'*yr;
