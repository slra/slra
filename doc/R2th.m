function th = R2th(R, d, psi)
if size(psi, 1) == size(psi, 2)
  th = R(:)' / psi;
else
  P = null(R); dh = P * (P \ d);
  th = lra(psi * kron(dh, eye(size(R, 1))), size(psi, 1) - 1);
end
