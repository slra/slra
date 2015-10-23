function [ph, info] = gcd_nothing(p, w, d, opt)
  bfn = cellfun(@length, p) - 1;
  bfell = bfn - d;
  bfg = mat2cell(opt.gini, bfell + 1);
  
  ph = cell(length(bfn),1);
  for i=1:length(bfn)
    ph{i} = multmat(bfg{i}, d) * opt.hini;
  end    
  info.iter = -1;
end

