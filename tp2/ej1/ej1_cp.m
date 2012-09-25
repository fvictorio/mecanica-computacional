function [ error, a, K, x_error, y_ex, y_ap ] = ej1_cp( M )
%EJ1 [ error, a, K, x_error, y_ex, y_ap ] = ej1_cp( M )
% Aproximacion de 1 + sin(pi/2*x) usando residuos ponderados,
% con colocacion puntual.

% Funcion a aproximar
phi_ex = inline('1 + sin(pi/2*x)', 'x');

% Psi que satisface las condiciones de borde
psi = inline('1 + x', 'x');
N_m = inline('sin(m*pi*x)', 'm', 'x');

K = zeros(M);
f = zeros(M, 1);

% Los puntos donde va a evaluarse el error
% (O sea, los x_l de los w_l = d(x - x_l))
x = linspace(0, 1, M+2)';
x = x(2 : M+1);

for l = 1 : M
    f(l) = phi_ex(x(l)) - psi(x(l));
    for m = 1 : M
        K(l,m) = N_m(m, x(l));
    end
end

% Calculo los coeficientes
a = K\f;

% Calculo el error
x_error = linspace(0, 1, 1000)';
y_ex = phi_ex(x_error);
y_ap = psi(x_error);
for i = 1 : M
    y_ap = y_ap + a(i) * N_m(i, x_error);
end

error = norm(y_ex - y_ap);

end

