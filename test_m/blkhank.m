load test-blkhank.mat

addpath support

opt.maxiter = 100;
opt.epsrel  = 1e-7; 
opt.epsabs  = 0; 
opt.epsgrad = 1e-5; 
opt.disp    = 'disp';
opt.maxiter = 50;



function arr = hankelize(H, p, q) 
  [pl, qk] = size(H);
  L = int8(pl / p)
  K = int8(qk / q)
  N = L+K -1;
  
  arr2 = zeros(p, q * N)

  w = zeros(1,N);

	if L <= 0 || K <= 0 
	  error('Incorrect params') 
	end
	
	for i=1:L
	  arr2(:,((i-1)*q+1):((i+K-1)*q)) =  arr2(:,((i-1)*q+1):((i+K-1)*q)) + H(((i-1)*p+1):(i*p), :);
	  w(1,i:(i+K-1)) = w(1,i:(i+K-1)) + ones(1,K);
	endfor

	for i=1:N
	  arr2(:,((i-1)*q+1):(i*q)) =  arr2(:,((i-1)*q+1):(i*q)) / w(1, i);
	endfor

	arr = reshape(arr2, p, q, N);

endfunction

function arrout = getssa(arrin, L, rk)  
  [p,q,N] = size(arrin);
	
  g = blkhankel(arrin, L);
  [Ug,Sg,Vg] = svd(g);
  grec = Ug(:,1:rk) * Sg(1:rk,1:rk) * Vg(:,1:rk)';

  arrout = hankelize(grec,p,q);  

endfunction




[p, m, T] = size(y);
y = y0 + 0.01 * randn(p, m, T);
h = blkhankel(y, 3);


a0 = h0(1:6,:)';
b0 = h0(7:9,:)';
a = h(1:6,:)';
b = h(7:9,:)';

%yssa = getssa(y, 20, 6);

%y(:,:, 1:3)
%yssa(:,:, 1:3)



%hssa = blkhankel(yssa, 3);
%assa = hssa(1:6,:)';
%bssa = hssa(7:9,:)';
%x0ssa = tls(assa,bssa);


s.k = 2;
s.a = [2 9 3];
[xh0, info, v] =  stls(a0,b0,s,[],opt);
[xh, infoh, v] =  stls(a,b,s,[],opt);
%[xhssa, infohssa, v] =  stls(a,b,s,x0ssa,opt);

xh0
xh
%xhssa

info
infoh
%infohssa



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

[zh, info, v] =  stls(ag,bg,sg);
[zh0, info, v] =  stls(ag0,bg0,sg);

zh0
zh


