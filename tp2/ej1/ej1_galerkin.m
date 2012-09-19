function [ error, a, K, x_error, y_ex, y_ap] = ej1_galerkin( M )
%EJ1 [ error, a, K, x_error, y_ex, y_ap] = ej1_galerkin( M )
% Aproximacion de 1 + sin(pi/2*x) usando residuos ponderados,
% con Galerkin.

syms x;

phi_ex = inline('1 + sin(pi/2*x)', 'x');

psi = inline('1 + x', 'x');
N_m = inline('sin(m*pi*x)', 'm', 'x');

K = zeros(M);
f = zeros(M, 1);

for l = 1 : M
    f(l) = double(int(N_m(l, x)*(phi_ex(x)-psi(x)), x, 0, 1));
    K(l,l) = double(int(N_m(l, x)*N_m(l, x), x, 0, 1));
end

a = K\f;

x_error = linspace(0, 1, 1000)';
y_ex = phi_ex(x_error);
y_ap = psi(x_error);
for i = 1 : M
    y_ap = y_ap + a(i) * N_m(i, x_error);
end

error = norm(y_ex - y_ap)/norm(y_ex);

end

