function [phi, X, Y, K, b] = mdf2D(dx, dy, Lx, Ly, cond_inf, cond_sup, cond_izq, cond_der, promediar_esquinas)
% [phi, X, Y] = mdf2D(dx, dy, Lx, Ly, cinf, csup, cizq, cder)
% Resolucion de la ecuacion diferencial 
% a^2phi/ax^2 + a^2phi/ay^2 = 2(x^2 + y^2)
% en un recinto rectangular de dimenxiones Lx y Ly, utilizando una malla
% con espaciamiento homogeneo de dx y dy.
% Las condiciones de borde son constantes y se aplican a cada borde
% en su totalidad, pero pueden ser Dirichlet o Neumann.
% Las condiciones se especifican con un string en el que el primer caracter
% representa el tipo de condicion y el resto el valor de esta. Por ejemplo,
% una condicion Dirichlet de valor 100 se especifica como "d100" o "D100",
% mientras que una condicion Neumann nula es "n0" o "N0".
    x = 0 : dx : Lx;
    y = 0 : dy : Ly;
    
    Nx = length(x);
    Ny = length(y);
    N  = Nx * Ny;
    
    tipo_cond_inf = lower(cond_inf(1));
    tipo_cond_sup = lower(cond_sup(1));
    tipo_cond_izq = lower(cond_izq(1));
    tipo_cond_der = lower(cond_der(1));
    
    val_cond_inf = str2num(cond_inf(2:length(cond_inf)));
    val_cond_sup = str2num(cond_sup(2:length(cond_sup)));
    val_cond_izq = str2num(cond_izq(2:length(cond_izq)));
    val_cond_der = str2num(cond_der(2:length(cond_der)));
    
    if nargin < 9 
        promediar_esquinas = 0;
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
            %AD HOC, ARREGLAR XXX
            if i == 2 && j == 2
                K(p, p) = 1;
                b(p) = 100;
                continue
            end % XXX            
            if i == 1 % Izquierda
                if tipo_cond_izq == 'd'
                    K(p, p) = 1;
                elseif tipo_cond_izq == 'n'
                    K(p, p) = -3/(2*dy);
                    p2 = pos(i+1, j, Nord, orden);
                    K(p, p2) = 2/dy;
                    p2 = pos(i+2, j, Nord, orden);
                    K(p, p2) = -1/(2*dy);
                end
                if promediar_esquinas
                    if (j==1) && (tipo_cond_izq == tipo_cond_inf)
                        b(p) = (val_cond_izq + val_cond_inf)/2;
                    elseif (j==Ny) && (tipo_cond_izq == tipo_cond_sup)
                        b(p) = (val_cond_izq + val_cond_sup)/2;
                    else
                        b(p) = val_cond_izq;
                    end
                else
                    b(p) = val_cond_izq;
                end
            elseif i == Nx % Derecha
                if tipo_cond_der == 'd'
                    K(p, p) = 1;
                elseif tipo_cond_der == 'n'
                    K(p, p) = 3/(2*dy);
                    p2 = pos(i-1, j, Nord, orden);
                    K(p, p2) = -2/dy;
                    p2 = pos(i-2, j, Nord, orden);
                    K(p, p2) = 1/(2*dy);
                end
                if promediar_esquinas
                    if (j==1) && (tipo_cond_der == tipo_cond_inf)
                        b(p) = (val_cond_der + val_cond_inf)/2;
                    elseif (j==Ny) && (tipo_cond_der == tipo_cond_sup)
                        b(p) = (val_cond_der + val_cond_sup)/2;
                    else
                        b(p) = val_cond_der;
                    end
                else
                    b(p) = val_cond_der;
                end
            elseif j == 1 % Abajo
                if tipo_cond_inf == 'd'
                    K(p, p) = 1;
                elseif tipo_cond_inf == 'n'
                    K(p, p) = -3/(2*dy);
                    p2 = pos(i, j+1, Nord, orden);
                    K(p, p2) = 2/dy;
                    p2 = pos(i, j+2, Nord, orden);
                    K(p, p2) = -1/(2*dy);
                end
                b(p) = val_cond_inf;
            elseif j == Ny % Arriba
                if tipo_cond_sup == 'd'
                    K(p, p) = 1;
                elseif tipo_cond_sup == 'n'
                    K(p, p) = 3/(2*dy);
                    p2 = pos(i, j-1, Nord, orden);
                    K(p, p2) = -2/dy;
                    p2 = pos(i, j-2, Nord, orden);
                    K(p, p2) = 1/(2*dy);
                end
                b(p) = val_cond_sup;
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
