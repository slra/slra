
syms u0 u1 v0 v1 w0 w1;
mmu = [u0 0 0; u1 u0 0; 0 u1 u0; 0 0 u1]
mmv = [v0 0 0; v1 v0 0; 0 v1 v0; 0 0 v1]
mmw = [w0 0 0; w1 w0 0; 0 w1 w0; 0 0 w1]
sylv = [-mmv mmu zeros(4,3); -mmw zeros(4,3) mmu]
det(sylv * (sylv'))
