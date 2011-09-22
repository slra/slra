%function [Apost]=colhouse(A,v)
%
%This function performs a Householder Post-Multiplication on the matrix A,
%using the Householder transform built from the Householder vector v, so
%Apost=A(I-2*v*v.'/(v.'*v))
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999

function [Apost]=colhouse(A,v)
beta=-2/(v.'*v);
w=beta*A*v;
Apost=A+w*v.';
