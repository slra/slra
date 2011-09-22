%Use:             yey(m)
%This function returns a square matrix
%of size mxm, with ones on the main
%anti-diagonal.
function [M]=yey(m)

M = zeros(m);
M([1:m:m^2]+[m-1:-1:0]) = 1;