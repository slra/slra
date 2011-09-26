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
	end

end
