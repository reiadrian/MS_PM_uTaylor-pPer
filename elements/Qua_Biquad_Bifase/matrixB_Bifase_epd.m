function [Bd,DetJe_d,DerCae_p,DetJe_p]= matrixB_Bifase_epd(coord_n,xg,e_DatElemSet,e_VG)

%*********************************************************************************
%*  ARGUMENTOS DE ENTRADA:                                                       *
%*  ���������������������                                                        *
%*  coord_n : Coordenadas de los nodos que componen el elemento.                 *
%*            coord_n = [x1  x2  x3 x4 x5  x6  x7 x8                             *
%*                       y1  y2  y3 y4 y5  y6  y7 y8]                            *
%*  numeraci�n de nodos: 4------7------3                                         *
%*                       |             |                                         *
%*                       8             6                                         *
%*                       |             |                                         *
%*                       1------5------2                                         *
%*  xg      : Vector de posicion del punto de Gauss en estudio.                  *
%*            xg = [E n]                                                         *
%*  npe     : Numero de nodos por elemento.                                      *
%*                                                                               *
%*  LISTA DE VARIABLES:                                                          *
%*  ������������������                                                           *
%*  E       : coordenada natural                                                 *
%*  n       : coordenada natural                                                 *
%*  Elementos de 8 nodos                                                         *
%*  dN_E8   : Derivadas de las funciones de forma respecto a "E"                 *
%*            dN_E = [dN1_E  dN2_E  dN3_E dN4_E dN5_E  dN6_E  dN7_E dN8_E]       *
%*  dN_n8   : Derivadas de las funciones de forma respecto a "n"                 *
%*            dN_n = [dN1_n  dN2_n  dN3_n dN4_n N5_n  dN6_n  dN7_n dN8_n]        *
%*  dN_En8   : Derivadas de las funciones de forma respecto a "E" y a "n"        *
%*  J8       : matriz jacobiana de la transformacion                             *
%*            J = [dx_E8 dy_E8 ; dx_n8 dy_n8]                                    *
%*  dN_xy8   : Derivadas de las funciones de forma respecto a "x" y a "y"        *
%*                                                                               *
%*  Elementos de 4 nodos                                                         *
%*  dN_E4    : Derivadas de las funciones de forma respecto a "E"                *
%*            dN_E = [dN1_E  dN2_E  dN3_E dN4_E]                                 *
%*  dN_n4    : Derivadas de las funciones de forma respecto a "n"                *
%*            dN_n = [dN1_n  dN2_n  dN3_n dN4_n]                                 *
%*  dN_En4   : Derivadas de las funciones de forma respecto a "E" y a "n"        *
%*  J4       : matriz jacobiana de la transformacion                             *
%*            J = [dx_E dy_E ; dx_n dy_n]                                        *
%*                                                                               *
%*  ARGUMENTOS DE SALIDA:                                                        *
%*  ��������������������                                                         *
%*  Bd      : Matriz B para elemento cuadrangular de 8 nodos en estado plano.    *
%*  DetJe_d : determinante de la matriz jacobiana para elemento de 8 nodos.      *
%*  DerCae_p: Derivadas de las funciones de forma respecto a "x" y a "y" (Elem 4)*
%*  DetJe_p : determinante de la matriz jacobiana para elemento de 8 nodos.      *
%*********************************************************************************

%ndn = e_VG.ndn;
ndn_d = e_DatElemSet.ndn_d;
ndn_p = e_DatElemSet.ndn_p;
ntens = e_VG.ntens;
%dofpe = e_DatElemSet.dofpe;
dofpe_d = e_DatElemSet.dofpe_d;
dofpe_p=  e_DatElemSet.dofpe_p;

E = xg(1);
n = xg(2);

% Desplazamientos
% Derivada de funciones de forma: Elemento cuadrilatero de 8 nodos
% Derivadas de Ni respecto "E" y "n":
% ����������������������������������
   dN_E8 = [-1/4*((1-n).*(-E-n-1)+(1-E).*(1-n))  1/4*((1-n).*(E-n-1)+(1+E).*(1-n)) ...
            1/4*((1+n).*(E+n-1)+(1+E).*(1+n))  -1/4*((1+n).*(-E+n-1)+(1-E).*(1+n))...
            -E+E.*n                             1/2*(1-n.^2)                      ...
            -E-E.*n                            -1/2*(1-n.^2)                     ];
        
   dN_n8 = [-1/4*((1-E).*(-E-n-1)+(1-E).*(1-n))  -1/4*((1+E).*(E-n-1)+(1+E).*(1-n))...
            1/4*((1+E).*(E+n-1)+(1+E).*(1+n))   1/4*((1-E).*(-E+n-1)+(1-E).*(1+n))...
            -1/2*(1-E.^2)                       -n-E.*n                           ...
             1/2*(1-E.^2)                       -n+E.*n                          ];
         
   dN_En8 = [dN_E8;dN_n8];

% Matriz Jacobiana:
% ����������������
   J11_8 = coord_n(1,:)*dN_E8';
   J12_8 = coord_n(2,:)*dN_E8';
   J21_8 = coord_n(1,:)*dN_n8';
   J22_8 = coord_n(2,:)*dN_n8';
   J8 = [J11_8 J12_8 ; J21_8 J22_8];
   DetJe_d = det(J8);

% Matriz B:
% ��������
   Bd = zeros(ntens,dofpe_d);
   dN_xy8 = J8\dN_En8;
   Bd(1,1:ndn_d:dofpe_d) = dN_xy8(1,:);
   Bd(2,2:ndn_d:dofpe_d) = dN_xy8(2,:);
   Bd(4,1:ndn_d:dofpe_d) = dN_xy8(2,:);
   Bd(4,2:ndn_d:dofpe_d) = dN_xy8(1,:);

% Poropresiones
% Derivada de funciones de forma: Elemento cuadrilatero de 4 nodos
% Derivadas de Ni respecto "E" y "n":
% ����������������������������������
   dN_E4 = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
   dN_n4 = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
   dN_En4 = [dN_E4;dN_n4];

% Matriz Jacobiana:
% ����������������
%VERIFICAR QUE LAS COORDENADAS SEAN SOLO LAS DE ESQUINA
   J11_4 = coord_n(1,1:ndn_p:dofpe_p)*dN_E4'; 
   J12_4 = coord_n(2,1:ndn_p:dofpe_p)*dN_E4';
   J21_4 = coord_n(1,1:ndn_p:dofpe_p)*dN_n4';
   J22_4 = coord_n(2,1:ndn_p:dofpe_p)*dN_n4';
   J4 = [J11_4 J12_4 ; J21_4 J22_4];
   DetJe_p = det(J4);
% Derivadas cartesianas
   DerCae_p = J4\dN_En4;   
  