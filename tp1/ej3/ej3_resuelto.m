function [phi, K, Nx, Ny] = ej3_resuelto(dx, dy, Lx, Ly)
    x = 0:dx:Lx;
    y = 0:dy:Ly;

    Q = inline('2*(x^2 + y^2)');

    T = 0;

    Nx = length(x) - 1; % floor(Lx/dx)
    Ny = length(y) - 1;

    K   = zeros((Nx+1)*(Ny+1)); % matriz global
    phi = zeros((Nx+1)*(Ny+1), 1);
    b   = zeros((Nx+1)*(Ny+1), 1);

    p = pos((Nx+1), (Ny+1), Ny);

    for i = 1:Nx+1
        for j = 1:Ny+1
            p = pos(i, j, Ny);
            if ((i==1)||(i==Nx+1)||(j==1)||(j==Ny+1))
                K(p,p) = 1;
                b(p) = T;
            else
                K(p,p) = -2/(dx*dx) - 2/(dy*dy);
                p2 = pos(i, j+1, Ny);
                K(p,p2) = 1/(dy*dy);
                p2 = pos(i, j-1, Ny);
                K(p,p2) = 1/(dy*dy);
                p2 = pos(i-1, j, Ny);
                K(p,p2) = 1/(dx*dx);
                p2 = pos(i+1, j, Ny);
                K(p,p2) = 1/(dx*dx);
                b(p) = -Q(i*dx, j*dy);
                end
        end
    end

    K
    b
    phi = K\b;
    [m n] = size(phi);

    phir = reshape(phi, Ny+1, Nx+1);

    [X,Y] = meshgrid(x,y);


