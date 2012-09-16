function [y] = ej1_temp_exacta(x)
   A = (exp(-1))/(exp(1) - exp(-1));
   B = (exp(1))/(exp(-1) - exp(1));
   y = A*exp(x) + B*exp(-x) + 1;
