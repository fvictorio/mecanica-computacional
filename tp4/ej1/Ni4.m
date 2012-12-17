function [ z ] = Ni4( i, x, y, xs, ys )

    xsCell = num2cell(xs);
    ysCell = num2cell(ys);
    
    [xi xk] = xsCell{:};
    [yi yk] = ysCell{:};
    
    dx = xk - xi;
    dy = yk - yi;
    A = dx * dy;
    
    z = 1;
    
    if i == 1
        z = z * (xk - x) * (yk -  y);
    elseif i == 2
        z = z * (x - xi) * (yk -  y);
    elseif i == 3
        z = z * (x - xi) * (y - yi);
    elseif i == 4
        z = z * (xk - x) * (y - yi);
    end
    
    z = z / A;
end

