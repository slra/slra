function [ ph, info  ] = gcd_fastgcd_onestep(p, w, d, opt)
  bfn = cellfun(@length, p) - 1;
  bfell = bfn - d;
  bfg = mat2cell(opt.gini, bfell + 1);
  
  if (length(bfn) ~= 2)
    error('Only 2 polynomials are allowed');  
  end    

  hh =    opt.hini.';
  g1h = bfg{1}.';
  g2h = bfg{2}.';
  
  N = length(p{1}) - 1;
  M = length(p{2}) - 1;
  altro_Z=c_iterfast(p{1}.',p{2}.', g1h, g2h, hh, 100, 0);
  hh=fliplr((altro_Z(1:d+1)).'); % new gcd
  g1h=fliplr((altro_Z(d+2:2+N)).');
  g2h=fliplr((altro_Z(N+3:N+3+M-d)).'); % new cofactors
  ph = {conv(hh,g1h).' ; conv(hh,g2h).'};
  info.iter = 1;
end

