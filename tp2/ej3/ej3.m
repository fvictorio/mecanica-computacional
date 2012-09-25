function [ a, K, f ] = ej3( M )
%EJ3 Summary of this function goes here
%   Detailed explanation goes here

syms x y;

psi = inline('2 - x.^2 - y.^2', 'x', 'y');
N   = inline('sin(m*pi*(x+1)/2) .* sin(m*pi*(y+1)/2)', 'm', 'x', 'y');
N2 = inline('-m^2 * pi^2 * sin(m*pi*(x+1)/2) .* sin(m*pi*(y+1)/2) / 4 ', 'm', 'x', 'y');
%N   = inline('(1-x^(2*m)) * (1-y^(2*m))', 'm', 'x', 'y');
%N2x = inline('', 'm', 'x', 'y');

K = zeros(M);
f = zeros(M, 1);

for l = 1 : M
    f(l) = 4 * double(int(int(N(l, x, y), x, -1, 1), y, -1, 1));
    for m = 1 : M
        K(l, m) = double(int(int(N(l,x,y)*(N2(m,x,y) + N2(m,x,y)), x, -1, 1), y, -1, 1));
    end
end

a = K \ f;
x = [-1 : 0.05 : 1];
y = [-1 : 0.05 : 1];

plot(x, N(1,x,y)); return;

[X, Y] = meshgrid(x, y);
phi = psi(X, Y);
for m = 1 : M
    phi = phi + a(m) * N(m, X, Y);
end

surf(X,Y,phi);


end

