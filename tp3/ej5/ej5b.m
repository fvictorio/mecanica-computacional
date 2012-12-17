L = 1;
Ne = 80;
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
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + Kelem * (-exp(-x(i+1)) + exp(-x(i)));
end

for l = 2 : Nn-1
    %fun = eval(strcat('@(x,y) (',al,'+',bl,'*x+',cl,'*y).*(',am,'+',bm,'*x+',cm,'*y)'));
    Nj = eval(strcat('@(x) (x-', num2str(x(l-1)), ')/(', num2str(h), ') .* exp(-x)'));
    Ni = eval(strcat('@(x) (', num2str(x(l+1)), '-x)/(', num2str(h), ') .* exp(-x)'));
    f(l) = f(l) + bx * quad(Nj, x(l-1), x(l));
    f(l) = f(l) + bx * quad(Ni, x(l), x(l+1));
end

Nj = eval(strcat('@(x) (x-', num2str(x(Nn-1)), ')/(', num2str(h), ') .* exp(-x)'));
f(Nn) = f(Nn) + bx * quad(Nj, x(Nn-1), x(Nn));

K(1,:) = zeros(1, Nn);
K(1,1) = 1;
f(1) = 0;

f(Nn) = f(Nn) + t * quad(Nj, x(Nn-1), x(Nn));

u = K \ f;

plot(x,u);
xlabel('x');
ylabel('y');
title('Desplazamientos');

figure;
t = zeros(Nn,1);
t(1) = dot(TT([0,1,2,3,4], 1), u(1:5))/h;
t(2) = dot(TT([0,1,2,3,4], 1), u(2:6))/h;
coef = TT([-2,-1, 0, 1,2], 1);
for i = 3 : Nn - 2
    t(i) = dot(coef, u(i-2:i+2))/h;
end
t(Nn-1) = dot(TT([-4, -3, -2,-1,0], 1), u(Nn-5:Nn-1))/h;
t(Nn) = dot(TT([-4, -3, -2,-1,0], 1), u(Nn-4:Nn))/h;

plot(x,t);
xlabel('x');
ylabel('y');
title('Tensiones');












