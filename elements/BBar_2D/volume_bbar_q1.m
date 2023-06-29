function [vol] = volume_bbar_q1(coord_n,Area,e_VG)

%global npg wg xg
npg = e_VG.npg;
wg = e_VG.wg;
xg = e_VG.xg;

vol = 0;

for i = 1:npg;
   for j = 1:npg;
          
      E = xg(1);
      n = xg(2);

      % Derivadas de Ni respecto "E" y "n":
      % ����������������������������������
      dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
      dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
      dN_En = [dN_E ; dN_n];

      % Matriz Jacobiana:
      % ����������������
      J11 = coord_n(1,:)*dN_E';
      J12 = coord_n(2,:)*dN_E';
      J21 = coord_n(1,:)*dN_n';
      J22 = coord_n(2,:)*dN_n';
      J = [J11 J12 ; J21 J22];
      detJ = det(J);
          
      vol = vol + wg(i)*wg(j)*detJ;
     
   end
end

vol = vol*Area;