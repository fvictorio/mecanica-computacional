function [ error, a, K, x_error, y_ex, y_ap] = ej2( M )
%EJ1 [ error, a, K, x_error, y_ex, y_ap] = ej2( M )
% Resolucion de phi'' + phi + 1 = 0 en una dimension usando residuos% ponderados con Galerkin. (Con debilitacion.)

syms x;

% Solucion analitica de la ecuacion diferencial
phi_ex = inline('((1 + sin(1) - cos(1))/(cos(1) + sin(1))) * sin(x) + cos(x) - 1', 'x');

psi = inline('0', 'x');
N_m = inline('x.^(m)', 'm', 'x');
DN_m = inline('(m)*x.^(m-1)', 'm', 'x');

K = zeros(M);
f = zeros(M, 1);

for l = 1 : M
    f(l) = -double(int(N_m(l, x), x, 0, 1));
    for m = 1 : M
        K(l,m) = -double(int(DN_m(l,x) * DN_m(m,x), x, 0, 1)) + double(int(N_m(l,x)*N_m(m,x),x,0,1)) - N_m(l, 1)*N_m(m,1);
    end
end

a = K\f;

x_error = linspace(0, 1, 1000);
y_ex = phi_ex(x_error);
y_ap = psi(x_error);
for i = 1 : M
    y_ap = y_ap + a(i) * N_m(i, x_error);
end

error = norm(y_ex - y_ap)/norm(y_ex);

end

