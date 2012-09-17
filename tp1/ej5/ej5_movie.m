function [F] = ej5_movie(metodo, Lx, Ly, Lt, dx, dy, dt)
%EJ5_MOVIE Summary of this function goes here
%   Detailed explanation goes here



if lower(metodo) == 'e'
    [u, Nx, Ny, Nt, X, Y] = ej5_expl(Lx, Ly, Lt, dx, dy, dt);
elseif lower(metodo) == 'i'
    [u, Nx, Ny, Nt, X, Y] = ej5_impl(Lx, Ly, Lt, dx, dy, dt);
elseif lower(metodo) == 'c'
    [u, Nx, Ny, Nt, X, Y] = ej5_cn(Lx, Ly, Lt, dx, dy, dt);
end

F(Nt) = struct('cdata', [], 'colormap', []);
figure;

for t = 1 : Nt
    surf(X, Y, u([1:Ny]+Ny*(t-1),1:Nx));

    axis([0,Lx,0,Ly,0,50]);
    grid on;
    set(gca, 'xtick', [0:dx:Lx])
    set(gca, 'ytick', [0:dy:Ly])
    set(gca, 'ztick', [0:10:50])
    xlabel('x')
    ylabel('y')
    zlabel('u')
    if lower(metodo) == 'e'
        title('Ejercicio 5 - Explicito')
    elseif lower(metodo) == 'i'
        title('Ejercicio 5 - Implicito')
    elseif lower(metodo) == 'c'
        title('Ejercicio 5 - Crank-Nicholson')
    end

    F(t) = getframe(gcf);
end

end

