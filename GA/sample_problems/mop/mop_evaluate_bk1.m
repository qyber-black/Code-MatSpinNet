function [ fx ] = mop_evaluate_bk1( x )
% Fun��o teste BK1: 2 entradas e 2 sa�das
% Dom�nio: -5 a 10

fx(2) = (x(1)-5)^2+(x(2)-5)^2;
fx(1) = x(1)^2+x(2)^2;

end

