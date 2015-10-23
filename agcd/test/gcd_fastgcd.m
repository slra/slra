function [ ph, info ] = gcd_fastgcd(p, w, d, opt)
  bfn = cellfun(@length, p) - 1
  bfell = bfn - d;
  bfg = mat2cell(opt.gini, bfell + 1);
  
  if (length(bfn) ~= 2)
    error('Only 2 polynomials are allowed');  
  end    

  hh =    opt.hini.'
  g1h = bfg{2}.'
  g2h = bfg{1}.'
  
  [hh,res1, g1h,g2h,num_iter] = c_f_newton_iter(p{1}.',p{2}.', ...
           hh, g1h, g2h, 100, 1);
  hh
  g1h
  g2h
       
  ph = {conv(hh,g1h).' ; conv(hh,g2h).'};
  info.iter = num_iter;
end

