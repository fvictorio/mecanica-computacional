function [y] = ej1_flujo_exacto(x)
   A = (exp(-1))/(exp(1) - exp(-1));
   B = (exp(1))/(exp(-1) - exp(1));
   y = A*exp(x) - B*exp(-x);
