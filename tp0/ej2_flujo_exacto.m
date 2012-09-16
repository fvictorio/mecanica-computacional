function [y] = ej2_flujo_exacto(x)
   A = -(1+e^-1)/(e + e^-1);
   B = (1-e)/(e+e^-1);
   y = B*exp(-x) - A*exp(x);
