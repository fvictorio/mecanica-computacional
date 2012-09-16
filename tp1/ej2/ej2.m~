dsafcvxfunction [error, x, phi_mdf, phi_ex] = ej2(N, q0, q1)
%El error al resolver T'' + T = 4(x-1/2)^2 -1
%por el metodo de diferencias finitas al compararlo
%con la solucion exacta utilizando norma euclidea

if nargin == 1
    q0 = 10;
    q1 = 0;
end

[phi_mdf, x] = mdf(N, q0, q1);
phi_ex  = exacta(x);
error = norm(phi_mdf - phi_ex)/norm(phi_ex);


end

