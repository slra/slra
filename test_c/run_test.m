function myinf = run_test(slra_fun, num, opt) 
  num = num2str(num);
  s_raw = load_raw([pwd '/s' num '.txt']);
  [m, r] = deal(s_raw(3), s_raw(4));
  s = struct('n', s_raw(4 + (1:s_raw(1))), 'm', s_raw(4 + s_raw(1) + (1:s_raw(2))));
  w_raw = load_raw([pwd '/w' num '.txt']);
  if (~isempty(w_raw))
    s.w = w_raw;
  end
  p = load_raw([pwd '/p' num '.txt']);
  p = p(1:(sum(s.n) * length(s.m) + sum(s.m) * length(s.n) - length(s.n) * length(s.m)));
  phi_raw = load_raw([pwd '/phi' num '.txt']);
  if (~isempty(phi_raw))
    s.phi = reshape(phi_raw, m, sum(s.m));
  end
  r_raw = load_raw([pwd '/r' num '.txt']);
  if (~isempty(r_raw))
    opt.Rini = reshape(r_raw, m-r, m);
  end
  
  st = tic;
  [ph, myinf] = slra_fun(p, s, r, opt);
  myinf.time = toc(st);
  myinf.m = m;
  myinf.d = m-r;
  myinf.q = s_raw(2);
  
  pos = 1;
  mnrm = 0;
  p = p(:);
  for l=1:length(s.n)
    for k=1:length(s.m)
      T = s.m(k) + s.n(l)-1;
      tp = p(pos:(pos+T-1));
      mnrm = mnrm + norm(tp - mean(tp))^2;
      pos = pos+T;
    end
  end    
    
  myinf.fit = (1 - sqrt(myinf.fmin) / sqrt(mnrm)) * 100;
end
