nodos = [0 0; 2.5 0; 5 0; 0 2.5; 2.5 2.5; 5 2.5; 0 5; 2.5 5; 5 5];
elems = [1 2 4; 2 5 4; 2 3 5; 3 6 5; 4 5 7; 5 8 7; 5 6 8; 6 9 8];

Nnodos = length(nodos(:,1));
Nelems = length(elems(:,1));

nodos_d = [2 0; 3 50; 6 100];
elems_q = [5 1 3 2; 1 1 3 2];
elems_h = [6 2 3; 8 2 3];

Qsup = 1.2;
elems_fuente = [4 6 7 8];

Qpuntual = 5;
pos_fuente_puntual = [1 1];
elems_fuente_puntual = [1];

k_conduct = 2;
T_amb = 30;
h_value = 1.2;

K = zeros(Nnodos);
f = zeros(Nnodos, 1);

% K correspondiente al laplaciano
% integrar en el elemento
% diff(Nl, x)*diff(Nm,x) + diff(Nl, y)*diff(Nm,y)
for i_elem = 1 : Nelems
    Kelem = zeros(3);
    nodos_elem = elems(i_elem,:);
    x = nodos(nodos_elem, 1)';
    y = nodos(nodos_elem, 2)';
    b = [];
    c = [];
    AA = det([1 1 1; x; y]);
    b(1) = y(2) - y(3);
    b(2) = y(3) - y(1);
    b(3) = y(1) - y(2);
    c(1) = x(3) - x(2);
    c(2) = x(1) - x(3);
    c(3) = x(2) - x(1);
    b = b/AA;
    c = c/AA;
    
    for l = 1 : 3
        for m = 1 : 3
            Kelem(l, m) = k_conduct * (b(l)*b(m) + c(l)*c(m));
        end
    end
    
    K(nodos_elem,nodos_elem) = K(nodos_elem,nodos_elem) + Kelem;
end

% K correspondiente a la condicion Robin
for i = 1 : size(elems_h, 1)
    Kelem = zeros(3);
    i_elem = elems_h(i, 1);
    nodos_elem = elems(i_elem,:);
    nodo_i = elems_h(i, 2);
    nodo_j = elems_h(i, 3);
    x = nodos(nodos_elem, 1)';
    y = nodos(nodos_elem, 2)';
    a = [];
    b = [];
    c = [];
    AA = det([1 1 1; x; y]);
    a(1) = x(2)*y(3) - x(3)*y(2);
    a(2) = x(3)*y(1) - x(1)*y(3);
    a(3) = x(1)*y(2) - x(2)*y(1);
    b(1) = y(2) - y(3);
    b(2) = y(3) - y(1);
    b(3) = y(1) - y(2);
    c(1) = x(3) - x(2);
    c(2) = x(1) - x(3);
    c(3) = x(2) - x(1);
    a = a/AA;
    b = b/AA;
    c = c/AA;
    al = num2str(a(nodo_i));
    bl = num2str(b(nodo_i));
    cl = num2str(c(nodo_i));
    am = num2str(a(nodo_j));
    bm = num2str(b(nodo_j));
    cm = num2str(c(nodo_j));
    fun = eval(strcat('@(x,y) (',al,'+',bl,'*x+',cl,'*y).*(',am,'+',bm,'*x+',cm,'*y)'));
    Kelem(nodo_i, nodo_j) = -(h_value / k_conduct) * int_triangle(x, y, fun);
    fun = eval(strcat('@(x,y) (',al,'+',bl,'*x+',cl,'*y).*(',al,'+',bl,'*x+',cl,'*y)'));
    Kelem(nodo_i, nodo_i) = -(h_value / k_conduct) * int_triangle(x, y, fun);
    fun = eval(strcat('@(x,y) (',am,'+',bm,'*x+',cm,'*y).*(',am,'+',bm,'*x+',cm,'*y)'));
    Kelem(nodo_j, nodo_j) = -(h_value / k_conduct) * int_triangle(x, y, fun);
    K(nodos_elem, nodos_elem) = K(nodos_elem, nodos_elem) + Kelem;
end

%====== F ========

% Q en triangulo superior

for i = 1 : length(elems_fuente)
    i_elem = elems_fuente(i);
    nodos_elem = elems(i_elem, :);
    x = nodos(nodos_elem, 1)';
    y = nodos(nodos_elem, 2)';
    A = det([1 1 1; x; y]) / 2;
    for l = nodos_elem
        f(l) = f(l) + Qsup * A / 3;
    end
end

% Q puntual

for i = 1 : length(elems_fuente_puntual)
    i_elem = elems_fuente_puntual(i);
    nodos_elem = elems(i_elem, :);
    x = nodos(nodos_elem, 1)';
    y = nodos(nodos_elem, 2)';
    a = [];
    b = [];
    c = [];
    AA = det([1 1 1; x; y]);
    a(1) = x(2)*y(3) - x(3)*y(2);
    a(2) = x(3)*y(1) - x(1)*y(3);
    a(3) = x(1)*y(2) - x(2)*y(1);
    b(1) = y(2) - y(3);
    b(2) = y(3) - y(1);
    b(3) = y(1) - y(2);
    c(1) = x(3) - x(2);
    c(2) = x(1) - x(3);
    c(3) = x(2) - x(1);
    a = a/AA;
    b = b/AA;
    c = c/AA;
    f(nodos_elem(1)) = (a(1) + b(1) * pos_fuente_puntual(1) + c(1) * pos_fuente_puntual(2)) * Qpuntual;
    f(nodos_elem(2)) = (a(2) + b(2) * pos_fuente_puntual(1) + c(2) * pos_fuente_puntual(2)) * Qpuntual;
    f(nodos_elem(3)) = (a(3) + b(3) * pos_fuente_puntual(1) + c(3) * pos_fuente_puntual(2)) * Qpuntual;
end

% Condicion Neumann

for i = 1 : size(elems_q, 1)
    i_elem = elems_q(i, 1);
    nodos_elem = elems(i_elem, :);
    nodo_i = nodos_elem(elems_q(i, 2));
    nodo_j = nodos_elem(elems_q(i, 3));
    val_q = elems_q(i, 4);
    f(nodo_i) = f(nodo_i) + (val_q / k_conduct) / 2;
    f(nodo_j) = f(nodo_j) + (val_q / k_conduct) / 2;
end

% Condicion Robin

for i = 1 : size(elems_h, 1)
    i_elem = elems_h(i, 1);
    nodos_elem = elems(i_elem, :);
    nodo_i = nodos_elem(elems_h(i, 2));
    nodo_j = nodos_elem(elems_h(i, 3));
    f(nodo_i) = f(nodo_i) - (h_value * T_amb / k_conduct) / 2;
    f(nodo_j) = f(nodo_j) - (h_value * T_amb / k_conduct) / 2;
end

% Imponer temperaturas dirichlet
for i = 1 : size(nodos_d, 1)
    ind_nodo = nodos_d(i,1);
    K(ind_nodo, :) = zeros(1, Nnodos);
    K(ind_nodo, ind_nodo) = 1;
    f(ind_nodo) = nodos_d(i, 2);
end

phi = K \ f;

if 1
    x = nodos(:, 1);
    y = nodos(:, 2);
    tri = delaunay(x, y);
    trisurf(tri, x, y, zeros(Nnodos, 1), phi, 'FaceColor', 'interp');
    xlabel('x');
    ylabel('y');
    view(2);
    colorbar;
end












