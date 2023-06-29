function m_FFe = f_FF_tria_t1(e_DatElemSet,e_VG)
   
   ndn = e_VG.ndn;
   struhyp = e_VG.struhyp;
   dofpe = e_DatElemSet.dofpe;
   nPG = e_DatElemSet.npg;
   xg = e_DatElemSet.xg;
   
   % Funciones de forma
   E = xg(:,1);
   n = xg(:,2);
   N1 = 1-n-E;
   N2 = E;
   N3 = n;
   
   m_FFe = zeros(ndn,dofpe,nPG);
   switch struhyp
      case {1,2,20} %Deformaci�n y tensi�n plana, y Large deformations with plane deformation.
         %m_FFpg = [N1,0,N2,0,N3,0;0,N1,0,N2,0,N3];
         m_FFe(1,1,:) = N1;
         m_FFe(1,3,:) = N2;
         m_FFe(1,5,:) = N3;
         m_FFe(2,2,:) = N1;
         m_FFe(2,4,:) = N2;
         m_FFe(2,6,:) = N3;


      otherwise
         error('Elemento tria_t1: Matriz de forma: Hip�tesis estructural no definido.');
   end

end