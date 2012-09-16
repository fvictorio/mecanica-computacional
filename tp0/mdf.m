function [x, T] = mdf(N, T0, TL)
    dx = 1/N;
    T = zeros(N-1, 1);
    M = zeros(N-1, N-1);
    b = zeros(N-1, 1);
    
    M(1, 1) = -2 - dx^2;
    M(1, 2) = 1;
    b(1) = -dx^2 - T0;
    for i = 2 : N-2
        M(i, i-1) = 1;
        M(i, i)   = -2 - dx^2;
        M(i, i+1) = 1;
        b(i) = -dx^2;
    end
    M(N-1, N-2) = 1;
    M(N-1, N-1) = -2 - dx^2;
    b(N-1) = -dx^2 - TL;
    
    x = [1:N-1]'*dx;
    T = M^(-1) * b;