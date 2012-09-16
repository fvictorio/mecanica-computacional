function [ phi ] = exacta( x )
%Solucion exacta a la ODE T'' + T = 4(x-1/2)^2 -1
%Con T'(0) = -10 y T'(1) = 0

k1 = -6;
k2 = -(6*cos(1)-4)/(sin(1));

phi = k1*sin(x) + k2*cos(x) + 4*x.^2 - 4*x -8;

end

