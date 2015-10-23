function [ ph, info ] = gcd_uvgcd(p, w, d, opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
  bfn = cellfun(@length, p) - 1
  bfell = bfn - d;
  bfg = mat2cell(opt.gini, bfell + 1);
  
  if (length(bfn) ~= 2)
    error('Only 2 polynomials are allowed');  
  end    
  
  [hh, g1h, g2h,  uvres, uvcond] = uvGCDfixedDegree(p{1},p{2},d)  
  ph = {conv(hh,g1h) ; conv(hh,g2h)};
  info.iter = -1;
end



