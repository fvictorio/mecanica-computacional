function [ T, x, K, f ] = ej4dn( N, phi1, q2)
% Solucion analitica
% y1 = -x^2/2 + x/2 + 1
% y2 = 9/8

h = 1/N;
x = 0 : h : 1;
N_nodos = N+1;

K = zeros(N_nodos);
f = zeros(N_nodos,1);

Kelem = [1/h -1/h; -1/h 1/h];
Nelem = length(Kelem);

for i = 1 : (N+1) - Nelem + 1
    K(i:i+Nelem-1,i:i+Nelem-1) = K(i:i+Nelem-1,i:i+Nelem-1) + Kelem;
end

K(1,:) =       [1 zeros(1, N_nodos - 1)];
f(1) = phi1;

f(N_nodos) = -q2;

for i = 2 : N_nodos - 1
    if (abs(x(i) - 0.5) < 1e-3)
        f(i) = h/2;
        break;
    end
    f(i) = h;
end

T = K\f;

end