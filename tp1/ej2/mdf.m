function [phi, x] = mdf(N, q0, q1)
%Ecuacion del calor en una dimension
%T'' + T = 4(x-1/2)^2 -1
%q0 y q1 son condiciones Neumann en x=0 y x=1

x = linspace(0, 1, N+1)';

h = x(2) - x(1);

N = length(x) + 2;
K = zeros(N);
f = zeros(N, 1);

K(1, 1) = -1/(2*h);
K(1, 3) =  1/(2*h);
f(1) = -q0;

for i = 2 : N-1
   K(i, i) = 1 - 2/h^2; 
   K(i, i-1) = 1/h^2;
   K(i, i+1) = 1/h^2;
   f(i) = 4*(x(i-1)-1/2)^2 - 1;
end

K(N, N-2) = -1/(2*h);
K(N, N)   =  1/(2*h);
f(N) = -q1;

phi_con_ficticios = K\f;
phi = phi_con_ficticios(2:N-1);

end

