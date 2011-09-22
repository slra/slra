load test-blkhank.mat

opt.maxiter = 100;
opt.epsrel  = 1e-7; 
opt.epsabs  = 0; 
opt.epsgrad = 1e-5; 
opt.disp    = 'disp';
opt.maxiter = 400;

function H = blkhankel(arr, L)
  [p,q,N] = size(arr);
	K = N - L + 1;
	
	if L <= 0 || K <= 0 
	  error('Incorrect params') 
	end
	
	arr2 = reshape(arr, p, q*N);
	H = zeros(p*L, q*K);
	
	for i=1:L
	  H(((i-1)*p+1):(i*p), :) = arr2(:,((i-1)*q+1):((i+K-1)*q));
	endfor

endfunction

[p, m, T] = size(y);
y = y0 + 0.02 * randn(p, m, T);
h = blkhankel(y, 3);

a = h(1:6,:)';
b = h(7:9,:)';
a0 = h0(1:6,:)';
b0 = h0(7:9,:)';


s.k = 2;
s.q = 1;
s.a = [2 9 3];
[xh, info, v] =  stls(a,b,s,[],opt);
[xh0, info, v] =  stls(a0,b0,s,opt);

xh
xh0


rk = 6;
L = 4;

g = blkhankel(y, L);
g0 = blkhankel(y0, L);

ag = g(1:rk,:)';
bg = g((rk+1):(L*3),:)';
ag0 = g0(1:rk,:)';
bg0 = g0((rk+1):(L*3),:)';

sg.k = 2;
sg.q = 1;
sg.a = [2 (L*3) 3];

%[zh, info, v] =  stls(ag,bg,sg);
%[zh0, info, v] =  stls(ag0,bg0,sg);

zh0
zh


