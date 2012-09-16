function [y] = ej2_temp_exacta(x)
   A = -(1+e^-1)/(e + e^-1);
   B = (1-e)/(e+e^-1);
   y = A*exp(x) + B*exp(-x) + 1;
