function [ y ] = N_alpha( alpha, h, x )
%N_ALPHA Summary of this function goes here
%   Detailed explanation goes here
if alpha == 1
    y = (2/h^2)*x.^2 - (1/h)*x;
elseif alpha == 2
    y = -(4/h^2)*x.^2 + 1;
else
    y = (2/h^2)*x.^2 + (1/h)*x;
end

end

