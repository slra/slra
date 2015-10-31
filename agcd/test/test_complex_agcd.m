u = conv([-1i;1], [1i;1]);
v = conv([-1i;1], [-1;1]);

%p = {u + [0.01;0;0], v};
p = {u + [0;0;0], v};
[ph, info] = gcd_nls_complex(p,[],1, struct('disp', 'iter'));
info.hh
ph{:}



