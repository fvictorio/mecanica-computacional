%=======================================================================
% MAT-femcCal 1.0  - MAT-femCal is a learning tool for undestanding 
%                    the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  kx =              1.00000 ;
  ky =              1.00000 ;
 heat=              0.00000 ;
%
% Coordinates
%
global coordinates
coordinates = [
        -5.83333   ,         3.29897  ;
        -5.83333   ,         2.29897  ;
        -4.83333   ,         3.29897  ;
        -4.94705   ,         1.89897  ;
        -4.11371   ,         2.39897  ;
        -5.83333   ,         1.29897  ;
        -3.83333   ,         3.29897  ;
        -4.96730   ,         0.79897  ;
        -4.10128   ,         1.29897  ;
        -5.83333   ,         0.29897  ;
        -2.83333   ,         3.29897  ;
        -2.83333   ,         2.29897  ;
        -4.10128   ,         0.29897  ;
        -4.96730   ,        -0.20103  ;
        -2.83333   ,         1.29897  ;
        -5.83333   ,        -0.70103  ;
        -2.83333   ,         0.29897  ;
        -3.83333   ,        -0.70103  ;
        -4.96730   ,        -1.20103  ;
        -5.83333   ,        -1.70103  ;
        -2.83333   ,        -0.70103  ;
        -1.83333   ,         0.29897  ;
        -3.83333   ,        -1.70103  ;
        -4.96730   ,        -2.20103  ;
        -1.83333   ,        -0.70103  ;
        -0.83333   ,         0.29897  ;
        -5.83333   ,        -2.70103  ;
        -3.83333   ,        -2.70103  ;
        -1.83333   ,        -1.70103  ;
        -4.96730   ,        -3.20103  ;
        -0.85566   ,        -1.11770  ;
         0.16667   ,         0.29897  ;
        -2.83333   ,        -2.70103  ;
        -5.83333   ,        -3.70103  ;
        -1.83333   ,        -2.70103  ;
         0.16667   ,        -0.70103  ;
        -3.90066   ,        -3.89825  ;
        -0.83705   ,        -2.27047  ;
        -4.96730   ,        -4.20103  ;
        -2.83333   ,        -3.70103  ;
         0.16667   ,        -1.70103  ;
        -5.83333   ,        -4.70103  ;
        -1.83333   ,        -3.70103  ;
        -3.96938   ,        -4.88436  ;
         0.16667   ,        -2.70103  ;
        -2.83333   ,        -4.70103  ;
        -0.83333   ,        -3.70103  ;
        -4.92506   ,        -5.33159  ;
        -5.83333   ,        -5.70103  ;
         0.16667   ,        -3.70103  ;
        -4.11371   ,        -5.80103  ;
        -2.83333   ,        -5.70103  ;
        -5.83333   ,        -6.70103  ;
        -4.83333   ,        -6.70103  ;
        -3.83333   ,        -6.70103  ;
        -2.83333   ,        -6.70103  ] ; 
%
% Elements
%
global elements
elements = [
      1   ,      2   ,      3   ; 
     40   ,     43   ,     35   ; 
     25   ,     29   ,     31   ; 
     35   ,     43   ,     47   ; 
     31   ,     29   ,     38   ; 
     25   ,     31   ,     26   ; 
     38   ,     29   ,     35   ; 
     38   ,     35   ,     47   ; 
     42   ,     49   ,     48   ; 
     48   ,     49   ,     54   ; 
     42   ,     48   ,     39   ; 
     39   ,     48   ,     44   ; 
     42   ,     39   ,     34   ; 
     44   ,     48   ,     51   ; 
     34   ,     39   ,     30   ; 
     44   ,     51   ,     52   ; 
     30   ,     39   ,     37   ; 
     37   ,     39   ,     44   ; 
     30   ,     37   ,     28   ; 
     37   ,     44   ,     46   ; 
     26   ,     22   ,     25   ; 
     25   ,     22   ,     21   ; 
     16   ,     20   ,     19   ; 
     19   ,     20   ,     24   ; 
     16   ,     19   ,     14   ; 
     24   ,     20   ,     27   ; 
     19   ,     24   ,     23   ; 
     14   ,     19   ,     18   ; 
     16   ,     14   ,     10   ; 
     24   ,     27   ,     30   ; 
     10   ,     14   ,      8   ; 
     30   ,     27   ,     34   ; 
     24   ,     30   ,     28   ; 
      8   ,     14   ,     13   ; 
     13   ,     14   ,     18   ; 
      8   ,     13   ,      9   ; 
      9   ,     13   ,     15   ; 
      8   ,      9   ,      4   ; 
      4   ,      9   ,      5   ; 
      8   ,      4   ,      6   ; 
      5   ,      9   ,     12   ; 
      4   ,      5   ,      3   ; 
      6   ,      4   ,      2   ; 
      8   ,      6   ,     10   ; 
      2   ,      4   ,      3   ; 
     45   ,     41   ,     38   ; 
     28   ,     23   ,     24   ; 
     55   ,     56   ,     52   ; 
     12   ,     11   ,      7   ; 
     47   ,     50   ,     45   ; 
     35   ,     33   ,     40   ; 
     40   ,     33   ,     28   ; 
     53   ,     54   ,     49   ; 
     17   ,     15   ,     13   ; 
     36   ,     32   ,     26   ; 
     18   ,     21   ,     17   ; 
     17   ,     21   ,     22   ; 
     52   ,     46   ,     44   ; 
      7   ,      3   ,      5   ; 
     31   ,     38   ,     41   ; 
     41   ,     36   ,     31   ; 
     23   ,     18   ,     19   ; 
     54   ,     55   ,     51   ; 
     15   ,     12   ,      9   ; 
     46   ,     40   ,     37   ; 
     51   ,     48   ,     54   ; 
     28   ,     37   ,     40   ; 
     13   ,     18   ,     17   ; 
      7   ,      5   ,     12   ; 
     51   ,     55   ,     52   ; 
     45   ,     38   ,     47   ; 
     31   ,     36   ,     26   ] ; 
%
% Fixed Nodes
%
fixnodes = [
     18  ,   100.00000  ;
     21  ,   100.00000  ;
     23  ,   100.00000  ;
     25  ,   100.00000  ;
     28  ,   100.00000  ;
     29  ,   100.00000  ;
     32  ,     0.00000  ;
     33  ,   100.00000  ;
     35  ,   100.00000  ;
     36  ,     0.00000  ;
     41  ,     0.00000  ;
     45  ,     0.00000  ;
     50  ,     0.00000  ] ;
%
% Punctual Fluxes
%
pointload = [
      1  ,    0.00000  ;
      2  ,    0.00000  ;
      3  ,    0.00000  ;
      6  ,    0.00000  ;
      7  ,    0.00000  ;
     10  ,    0.00000  ;
     11  ,    0.00000  ;
     12  ,    0.00000  ;
     15  ,    0.00000  ;
     16  ,    0.00000  ;
     17  ,    0.00000  ;
     20  ,    0.00000  ;
     22  ,    0.00000  ;
     26  ,    0.00000  ;
     27  ,    0.00000  ;
     34  ,    0.00000  ;
     40  ,    0.00000  ;
     42  ,    0.00000  ;
     43  ,    0.00000  ;
     46  ,    0.00000  ;
     47  ,    0.00000  ;
     49  ,    0.00000  ;
     52  ,    0.00000  ;
     53  ,    0.00000  ;
     54  ,    0.00000  ;
     55  ,    0.00000  ;
     56  ,    0.00000  ] ;
%
% Side loads
%
sideload = [
      3  ,      1  ,    0.00000   ;
      1  ,      2  ,    0.00000   ;
     11  ,      7  ,    0.00000   ;
     12  ,     11  ,    0.00000   ;
      7  ,      3  ,    0.00000  ];

