function M = multmat(h, m)
  M = toeplitz([h(:);zeros(m,1)], [h(1); zeros(m,1)]);
end