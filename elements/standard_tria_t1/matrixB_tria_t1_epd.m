function [B,detJ,dN_xy] = matrixB_tria_t1_epd(coord_n,e_DatElemSet,e_VG)

%*********************************************************************************
%*  Obtiene la matriz B del tri�ngulo lineal de 3 nodos para EPD � EPT           *
%*  evaluada en xg.                                                              *
%*                                                                               *
%*  ARGUMENTOS DE ENTRADA:                                                       *
%*  ���������������������                                                        *
%*  coord_n : Coordenadas de los nodos que componen el elemento.                 *
%*            coord_n = [x1  x2  x3                                              *
%*                       y1  y2  y3]                                             *
%*  numeraci�n de nodos: 3                                                       *
%*                       |\                                                      *
%*                       | \                                                     *
%*                       |  \                                                    *
%*                       |   \                                                   *
%*                       1----2                                                  *
%*  xg      : Vector de posici�n del punto de Gauss en estudio.                  *
%*            xg = [E n]                                                         *
%*  npe     : N�mero de nodos por elemento.                                      *
%*                                                                               *
%*  LISTA DE VARIABLES:                                                          *
%*  ������������������                                                           *
%*  E       : coordenada natural                                                 *
%*  n       : coordenada natural                                                 *
%*  dN_E    : Derivadas de las funciones de forma respecto a "E"                 *
%*            dN_E = [dN1_E  dN2_E  dN3_E]                                       *
%*  dN_n    : Derivadas de las funciones de forma respecto a "n"                 *
%*            dN_n = [dN1_n  dN2_n  dN3_n]                                       *
%*  dN_En   : Derivadas de las funciones de forma respecto a "E" y a "n"         *
%*  J       : matriz jacobiana de la transformaci�n                              *
%*            J = [dx_E dy_E ; dx_n dy_n]                                        *
%*  dN_xy   : Derivadas de las funciones de forma respecto a "x" y a "y"         *
%*                                                                               *
%*  ARGUMENTOS DE SALIDA:                                                        *
%*  ��������������������                                                         *
%*  B       : Matriz B para elemento triangular de 3 nodos en estado plano.      *
%*  detJ    : determinante de la matriz jacobiana.                               *
%*********************************************************************************

%global ndn dofpe ntens
ndn = e_VG.ndn;
ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;
%xg = e_DatElemSet.xg;

%La derivada de las funci�n de forma es uniforme, por lo que no depende de la coordenada.
% E = xg(1);
% n = xg(2);

% Derivadas de Ni respecto "E" y "n":
% ����������������������������������
dN_E = [-1,1,0];
dN_n = [-1,0,1];
dN_En = [dN_E;dN_n];

% Matriz Jacobiana:
% ����������������
J11 = coord_n(1,:)*dN_E';
J12 = coord_n(2,:)*dN_E';
J21 = coord_n(1,:)*dN_n';
J22 = coord_n(2,:)*dN_n';
J = [J11,J12;J21,J22];
detJ = det(J);

% Matriz B:
% ��������
B = zeros(ntens,dofpe);
dN_xy = J\dN_En;
B(1,1:ndn:dofpe) = dN_xy(1,:);
B(2,2:ndn:dofpe) = dN_xy(2,:);
B(4,1:ndn:dofpe) = dN_xy(2,:);
B(4,2:ndn:dofpe) = dN_xy(1,:);