L = 1;
Ne = 10;
bx = 1;
t = 1;
D = 1;

x = linspace(0, L, Ne+1);
Nn = length(x);
h = L/Ne;

K = zeros(Nn);
f = zeros(Nn, 1);

Kelem = D*[1/h -1/h; -1/h 1/h];
for i = 1 : Nn - 1
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + Kelem;
end

for l = 1 : Nn
    f(l) = f(l) + (h/2) * bx;
end

K(1,:) = zeros(1, Nn);
K(1,1) = 1;
f(1) = 0;

f(Nn) = f(Nn) + t;

u = K \ f;

plot(x,u);
xlabel('x');
ylabel('y');
title('Desplazamientos');

figure;
t = zeros(Nn,1);
t(1) = dot(TT([0,1,2], 1), u(1:3))/h;
coef = TT([-1, 0, 1], 1);
for i = 2 : Nn - 1
    t(i) = dot(coef, u(i-1:i+1))/h;
end
t(Nn) = dot(TT([-2,-1,0], 1), u(Nn-2:Nn))/h;

plot(x,t);
xlabel('x');
ylabel('y');
title('Tensiones');












