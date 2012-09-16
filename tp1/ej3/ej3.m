function [phi, K, X, Y] = ej3(dx, dy, Lx, Ly, Tinf, Tsup, Tizq, Tder)
    x = 0 : dx : Lx;
    y = 0 : dy : Ly;
    
    Nx = length(x);
    Ny = length(y);
    N  = Nx * Ny;
    
    if nargin == 4
        Tinf = 0;
        Tsup = 0;
        Tizq = 0;
        Tder = 100;
    end
    
    Q = inline('2*(x^2 + y^2)');
    
    if Nx <= Ny
        orden = 'x';
        Nord = Nx;
    else
        orden = 'y';
        Nord = Ny;
    end
    
    K = zeros(N);
    b = zeros(N, 1);
    
    for i = 1 : Nx
        for j = 1 : Ny
            p = pos(i, j, Nord, orden);
            if i == 1
                K(p, p) = 1;
                b(p) = Tizq;
            elseif i == Nx
                K(p, p) = 1;
                b(p) = Tder;
            elseif j == 1
                K(p, p) = 1;
                b(p) = Tinf;
            elseif j == Ny
                K(p, p) = 1;
                b(p) = Tsup;
            else
                K(p,p) = -2/(dx*dx) - 2/(dy*dy);
                p2 = pos(i, j+1, Nord, orden);
                K(p,p2) = 1/(dy*dy);
                p2 = pos(i, j-1, Nord, orden);
                K(p,p2) = 1/(dy*dy);
                p2 = pos(i+1, j, Nord, orden);
                K(p,p2) = 1/(dx*dx);
                p2 = pos(i-1, j, Nord, orden);
                K(p,p2) = 1/(dx*dx);
                b(p) = Q((i-1)*dx, (j-1)*dy);
            end
        end
    end
    
    phi_vec = K \ b;
    [m, n] = size(phi_vec);
    
    phi = reshape(phi_vec, Ny, Nx);
    
    [X, Y] = meshgrid(x, y);