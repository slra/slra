%function [Apre]=rowhouse(A,v)
%
%This function performs a Householder Pre-Multiplication on the matrix A,
%using the Householder transform built from the Householder vector v, so
%Apre=(I-2*v*v.'/(v.'*v))A
%
%Author: Philippe Lemmerling
%Last modified: 04/09/1999

function [Apre]=rowhouse(A,v)
beta=-2/(v.'*v);
w=beta*A.'*v;
Apre=A+v*w.';
