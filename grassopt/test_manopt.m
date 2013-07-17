     
% Generate the problem data.
T = 100;
p = 1:T;
p = p + 0.05 * randn(1, T);
r = 2;
s = struct('m', r+1, 'n', T-r);

[ph, info] = slra_grass(p, s, r, struct('maxiter', 100, 'disp', 'iter', 'checkgradient', 'yes'));
[ph1, info1] = slra(p, s, r, struct('maxiter', 100, 'disp', 'iter'));
     
% Display some statistics.
figure;
loglog(info.iterinfo(1,:), info.iterinfo(2,:), 'b.-', info1.iterinfo(1,:), info1.iterinfo(2,:), 'rx-');
xlabel('time');
ylabel('fmin');



