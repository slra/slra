function [res, iters, times] = run_time_tests(testpoly, methods, nreps, ds, opt, sigma, methods2)
  %UNTITLED2 Summary of this function goes here
  %   Detailed explanation goes here

  %% setup arrays
  ntest = size(testpoly,1);
  times = zeros(ntest, length(methods2));
  iters = zeros(ntest, length(methods));
  res = zeros(ntest, length(methods));


  %% Run the loop
  for k=1:ntest
    k
    p0 = testpoly{k, 1};
    q0 = testpoly{k, 2};
    for i=1:nreps
      i
      p = p0 + sigma * randn(size(p0));% sqrt(12) * (rand(size(p)) -0.5);
      q = q0 + sigma * randn(size(q0));% sqrt(12) * (rand(size(q)) -0.5);
      d = ds(k);
      opt.gini = g_ini({p,q}, d);
      opt.hini = lsdivmult({p,q}, d, opt.gini);

      for j=1:length(methods)
        [ph, info] = methods{j}({p,q}, [], d, opt);
        res(k, j) = res(k, j) + norm([p;q] -cell2mat(ph),2)^2;
        iters(k, j) = iters(k, j) + info.iter;
      end
      for j=1:length(methods2)
        j
        times(k,j) = times(k,j) + meas_func(methods2{j}, {{p,q}, [], d, opt});
      end
    end
  end

  res = sqrt(res / nreps);
  iters = iters / nreps;
  times = times / nreps;
end

