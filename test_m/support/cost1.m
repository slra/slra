% ---------------------------------------------------------------------------
% COST1 - STLS cost function; t = 1, d = 1
function cost = cost1( x, a, b, struct )

[S,V] = decode_struct(struct);

% Construct the first row of Vr
s = size(V,3);
f = zeros(1,s);
for k = 1:s
  f(k) = [x;-1]' * V(:,:,k)' * [x;-1];
end
r = a*x-b;

% Solve the system Vr*yr = r
yr   = mb02gd(f,r); % Slicot solver
cost = r'*yr;
