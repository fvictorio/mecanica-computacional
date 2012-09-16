function [phi, phi_ex, error] = ej4(x)
    if nargin == 0
        x = [0 0.05 0.15 0.35 0.5 0.75 1]';
    end

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

    phi_ex = exacta(x);

    error = norm(phi-phi_ex, 2)/norm(phi_ex,2);