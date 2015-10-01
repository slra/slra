p = {[-1i + 0.05; (1-1i); 1], [1i;(-1-1i);1]}
[ph, info] = gcd_nls_complex(p,[],1, struct('disp', 'iter'))
info.Rh
ph{:}