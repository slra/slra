function [f] = gcd_nls_evaluate(p, h_array)
  d = size(h_array, 1) - 1;
  p = reshape(p,length(p), 1);
  bfp = cell2mat(p);
  bfn = cellfun(@length, p) - 1;
  s = struct('m', d+1, 'n', bfn-d+1);
  s.gcd = 1;
  
  f = zeros(1, size(h_array,2));
  
  obj =  slra_mex_obj('new', bfp, s, d);

  for i=1:length(f) 
    f(i) = slra_mex_obj('func', obj, h_array(:,i)');
  end    
  slra_mex_obj('delete', obj);
end

