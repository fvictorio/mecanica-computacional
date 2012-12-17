cosote;

tini = 0;
tend = 800;
tdelta = 0.01;


Nnodes = size(coordinates, 1);
Nelems = size(elements, 1);

T_ev = zeros(Nnodes, length(tini:tdelta:tend));
T_ini = zeros(Nnodes, 1);

% Matriz de masa
C = zeros(Nnodes);
for i_elem = 1 : Nelems
    disp(strcat('Elemento ', num2str(i_elem), ' de ', num2str(Nelems)));
    Celem = zeros(3);
    nodos_elem = elements(i_elem,:);
    x = coordinates(nodos_elem, 1)';
    y = coordinates(nodos_elem, 2)';
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

    for l = 1 : 3
        al = num2str(a(l));
        bl = num2str(b(l));
        cl = num2str(c(l));
        for m = 1 : 3
            am = num2str(a(m));
            bm = num2str(b(m));
            cm = num2str(c(m));
            fun = eval(strcat('@(x,y) (',al,'+',bl,'*x+',cl,'*y).*(',am,'+',bm,'*x+',cm,'*y)'));
            Celem(l, m) = int_triangle(x, y, fun);
        end
    end

    C(nodos_elem,nodos_elem) = C(nodos_elem,nodos_elem) + Celem;
end

paso_tiempo = 0;
for t = tini : tdelta : tend
    paso_tiempo = paso_tiempo + 1;
    disp(strcat('Paso ', num2str(paso_tiempo), ' de ', num2str(size(T_ev,2))));
    K = zeros(Nnodes);
    f = zeros(Nnodes, 1);

    for i_elem = 1 : Nelems
        Kelem = zeros(3);
        nodos_elem = elements(i_elem,:);
        x = coordinates(nodos_elem, 1)';
        y = coordinates(nodos_elem, 2)';
        b = [];
        c = [];
        AA = det([1 1 1; x; y]);
        A = AA/2;
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
                Kelem(l, m) = b(l)*b(m) + c(l)*c(m);
            end
        end

        K(nodos_elem,nodos_elem) = K(nodos_elem,nodos_elem) + Kelem*A;
    end
    
    % Neumann
    for i = 1 : size(pointload, 1)
        n_load = pointload(i,1);
        val_load = pointload(i,2);
        
        f(n_load) = val_load;
    end
    
    % Dirichlet
    % (Creo que siempre deberia ir al final de lo de la K)
    for i = 1 : size(fixnodes, 1)
        n_fix = fixnodes(i,1);
        val_fix = fixnodes(i,2);
        
        f(n_fix) = val_fix;
        K(n_fix,:) = [zeros(1,n_fix-1) 1 zeros(1, Nnodes-n_fix)];
    end
    

    
    T = C\(tdelta*(f - K*T_ini) + C*T_ini); % explicito
    T_ini = T;
    T_ev(:,paso_tiempo) = T;
end


































