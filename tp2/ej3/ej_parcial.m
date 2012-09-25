function [ a, K, f ] = ej_parcial( M )
%EJ3 Summary of this function goes here
%   Detailed explanation goes here

syms r t;

%psi = inline('2 - x.^2 - y.^2', 'x', 'y');
N   = inline('sin(m*pi*(r-1)) .* sin(m*(t-pi/6)/(2-1/6))', 'm', 't', 'r');

K = zeros(M);
f = zeros(M, 1);

for l = 1 : M
    f(l) = -10*double( int(N(l,t,2), t, pi/6, 2*pi) );
    for m = 1 : M
        K(l, m) = 0;
        K(l, m) = K(l, m) - double( int( int(N(l, t, r) * diff(N(l, t, r), r) * diff(N(m, t, r), r),r, 1, 2),t,pi/6,2*pi ));
        K(l, m) = K(l, m) - double( int( int(N(l,t,r) * (1/r) * diff(N(m, t, r), r),r, 1, 2),t,pi/6,2*pi ));
        K(l, m) = K(l, m) - double( int( int(N(l,t,r) * (1/r^2) * diff(N(m, t, r), t, 2),r, 1, 2),t,pi/6,2*pi ));
    end
end

a = K \ f;

end

