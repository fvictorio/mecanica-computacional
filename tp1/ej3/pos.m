function [p] = pos(i, j, N, orden)
% Dados los indices de una malla bidimensional, devuelve el indice
% correspondiente al vector asociado si se enumeran por el lado con menos
% nodos. Por defecto se asume que se enumeran empezando por las x.
    if nargin == 4
       if (orden == 'y') || (orden == 'Y')
           temp = i;
           i = j;
           j = temp;
       end
    end
    p = (j-1) * N + i;
    return;