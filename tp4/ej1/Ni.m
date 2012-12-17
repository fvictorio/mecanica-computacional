function [ z ] = Ni( i, x, y, xs, ys )

j = mod(i, 3) + 1;
k = mod(i+1, 3) + 1;

xi = xs(i);
xj = xs(j);
xk = xs(k);
yi = ys(i);
yj = ys(j);
yk = ys(k);

ai = xj*yk - xk*yj;
bi = yj - yk;
ci = xk - xj;

A = det([1 1 1; xi xj xk; yi yj yk]);

z = (ai + bi*x + ci*y)/A;

end

