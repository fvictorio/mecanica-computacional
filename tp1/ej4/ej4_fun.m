function [ phi ] = ej4_fun(x)
%ej4_fun diferencias finitas en una dimension para mallas uniformes
%Problema: phi'' + phi  - 10 x (x+1) = 0

N = length(x);

K = zeros(N);
b = zeros(N, 1);

K(1, 1) = 1;
b(1)    = 1;

for i = 2 : N-1
   hmenos = x(i) - x(i-1);
   hmas = x(i+1) - x(i);
   
   K(i, i-1) = 2/(hmenos*(hmenos+hmas));
   K(i, i)   = -2/(hmas*hmenos) + 1;
   K(i, i+1) = 2/(hmas*(hmenos+hmas));
   
   b(i) = 10*x(i)*(x(i)+1);
end

K(N, N) = 1;
b(N)    = 0;

phi = K\b;

end

