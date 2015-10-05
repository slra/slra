function t = meas_func(func, pars)
  tic; eltm = toc;  n = 0;
  while (eltm < 1)
    func(pars{:});
    n = n+ 1;
    eltm = toc;
  end  
  t = eltm / n;
end