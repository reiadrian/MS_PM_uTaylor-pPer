function [kt] = elem_stiff_bbar_hexahedral(coord_n,u,hvar_old,aux_var,Eprop,ce,e_VG)

%global_var
struhyp = e_VG.struhyp;
conshyp = e_VG.conshyp;
dofpe = e_VG.dofpe;
npg = e_VG.npg;
xg = e_VG.xg;
wg = e_VG.wg;
ntens = e_VG.ntens;
sihvarpg = e_VG.sihvarpg;
siavarpg = e_VG.siavarpg;
sitvare = e_VG.sitvare;

% Inicializaciones
% ����������������
kt      = zeros(dofpe,dofpe);
eps_new = zeros(sitvare,1);

% C�lculo de la matriz Bh
% �����������������������
BH = matrixBh_bbar_hexahedral(coord_n,e_VG);

for i = 1:npg;
   for j = 1:npg;
      for k = 1:npg;
      
      % C�lculo de �ndices
      % ������������������
      ini =    ntens*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf =    ntens*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; itv = ini:inf;
      ini = sihvarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf = sihvarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; ihv = ini:inf;
      ini = siavarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf = siavarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; iav = ini:inf;

      % Deformaci�n
      % �����������
      if struhyp == 3
         [B,detJ] = matrixB_bbar_hexahedral(coord_n,[xg(i),xg(j),xg(k)],BH,e_VG);
      end
      eps_new(itv) = B*u;
      
      % Modelo constitutivo
      % �������������������
      switch conshyp
          case 2
              ct = rmapct_plasJ2(eps_new(itv),hvar_old(ihv),aux_var(iav),Eprop,ce,e_VG); 
          case 7
              ct = rmapct_danio_plasdp(eps_new(itv),hvar_old(ihv),aux_var(iav),Eprop,ce,e_VG);
      end
      
      % C�lculo de Kep
      % ��������������
      kt  = kt + B.'*ct*B*wg(i)*wg(j)*wg(k)*detJ;
      
      end
   end
end